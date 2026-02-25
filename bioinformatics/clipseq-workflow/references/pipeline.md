# CLIP/eCLIP Pipeline Cheatsheet

Minimal, runnable templates. Replace `$` vars. Keep samples in an array to loop.

## 0) Project variables
```bash
samples=(Sample1 Sample2)                      # base name without _R1/_R2
threads=${SLURM_CPUS_PER_TASK:-16}
genome_fa=/path/to/genome.fa                   # e.g., mm10.fa
genome_index=/path/to/index_prefix             # STAR or Bowtie index prefix
gtf=/path/to/annotation.gtf                    # for STAR/RCAS
chrom_sizes=/path/to/chrom.sizes               # for bigWig if needed
out=work
mkdir -p $out/trimmed $out/aln $out/logs $out/qc
```

## 1) QC raw reads
```bash
fastqc -t $threads -o $out/qc raw/*fastq.gz
```

## 2) Adapter + UMI handling (eCLIP-style example)
- Common eCLIP adapters (example; replace if protocol differs):  
  3' AACTTGTAGATCGGA, 3' AGGACCAAGATCGGA; 5' CTTCCGATCTACAAGTT, 5' CTTCCGATCTTGGTCCT  
- eCLIP often has a 5 bp UMI on R2 5' end; bleed-through may appear at R1 tail.

```bash
# Trim adapters; drop last 5 bp of R1 if UMI bleed-through
cutadapt -u -5 \
  -a AACTTGTAGATCGGA -a AGGACCAAGATCGGA \
  -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT \
  -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA \
  -G CTTCCGATCTACAAGTT -G CTTCCGATCTTGGTCCT \
  -m 10 \
  -o $out/trimmed/${s}_R1.trim.fq.gz \
  -p $out/trimmed/${s}_R2.trim.fq.gz \
  raw/${s}_R1.fastq.gz raw/${s}_R2.fastq.gz

# Extract 5 bp UMI from R2 5' (adjust pattern/side to protocol)
umi_tools extract \
  --bc-pattern=NNNNN --stdin=$out/trimmed/${s}_R2.trim.fq.gz --stdout=$out/trimmed/${s}_R2.umi.fq.gz \
  --read2-in=$out/trimmed/${s}_R1.trim.fq.gz --read2-out=$out/trimmed/${s}_R1.umi.fq.gz \
  --log=$out/logs/${s}_extract.log
```
Skip UMI steps if library lacks UMIs.

## 3) Alignment
### STAR (recommended for CLIP)
```bash
STAR --runThreadN $threads \
  --genomeDir ${genome_index} \
  --readFilesIn $out/trimmed/${s}_R2.umi.fq.gz $out/trimmed/${s}_R1.umi.fq.gz \  # eCLIP: R2 forward
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --alignEndsType EndToEnd \
  --outFileNamePrefix $out/aln/${s}. \
  --sjdbGTFfile $gtf
samtools sort -@ $threads -o $out/aln/${s}.sorted.bam $out/aln/${s}.Aligned.out.bam
samtools index -@ $threads $out/aln/${s}.sorted.bam
```

### Bowtie2 local (short inserts / heavy adapter contamination)
```bash
bowtie2 --very-sensitive-local -p $threads -x ${genome_index} \
  -1 $out/trimmed/${s}_R1.trim.fq.gz -2 $out/trimmed/${s}_R2.trim.fq.gz \
  | samtools sort -@ $threads -o $out/aln/${s}.sorted.bam
samtools index -@ $threads $out/aln/${s}.sorted.bam
```

### Bowtie1 legacy (as in some older CLIP pipelines)
```bash
bowtie --sam --chunkmbs 512 -p $threads ${genome_index} \
  -1 <(zcat $out/trimmed/${s}_R1.trim.fq.gz) -2 <(zcat $out/trimmed/${s}_R2.trim.fq.gz) \
  | samtools view -bS -T ${genome_fa} - \
  | samtools sort -@ $threads -o $out/aln/${s}.sorted.bam
samtools index -@ $threads $out/aln/${s}.sorted.bam
```

## 4) UMI-based deduplication
```bash
umi_tools dedup -I $out/aln/${s}.sorted.bam -S $out/aln/${s}.dedup.bam \
  --paired --method=adjacency --edit-distance-threshold=1 \
  --output-stats=$out/logs/${s}_dedup
samtools index -@ $threads $out/aln/${s}.dedup.bam
```

## 5) QC after mapping
```bash
samtools flagstat $out/aln/${s}.dedup.bam > $out/logs/${s}_flagstat.txt
plotFingerprint -b $out/aln/*.dedup.bam -p $threads --binSize 100 -o $out/qc/fingerprint.png
multiBamSummary bins --binSize 1000 -b $out/aln/*.dedup.bam -out $out/qc/mbs.npz
plotCorrelation -in $out/qc/mbs.npz --corMethod spearman --skipZeros --plotTitle "CLIP vs Ctrl" \
  --whatToPlot heatmap -o $out/qc/corr_heatmap.png
```

## 6) Coverage tracks (strand-specific CPM)
```bash
bamCoverage -b $out/aln/${s}.dedup.bam -o $out/aln/${s}.plus.bw \
  --filterRNAstrand forward --normalizeUsing CPM -bs 25 -p $threads
bamCoverage -b $out/aln/${s}.dedup.bam -o $out/aln/${s}.minus.bw \
  --filterRNAstrand reverse --normalizeUsing CPM -bs 25 -p $threads
```

## 7) Peak calling
### PEAKachu (with control)
```bash
peakachu callpeaks \
  --exp $out/aln/CLIP1.dedup.bam,$out/aln/CLIP2.dedup.bam \
  --ctrl $out/aln/CTRL1.dedup.bam,$out/aln/CTRL2.dedup.bam \
  --format bam --paired --max-insert 200 \
  --norm DESeq2 --mode adaptive --fc 2 --padj 0.05 \
  --out-prefix $out/peaks/clip
```

### PureCLIP (no control, nucleotide resolution)
```bash
PureCLIP -i $out/aln/${s}.dedup.bam -bai $out/aln/${s}.dedup.bam.bai \
  -g ${genome_fa} -o $out/peaks/${s}.pureclip.bed
```

## 8) Motif discovery (MEME-ChIP) and annotation (RCAS)
```bash
mkdir -p $out/motif
# extend peaks ±20 bp and get fasta
bedtools slop -i $out/peaks/clip_peaks.bed -g $chrom_sizes -b 20 \
  | bedtools getfasta -fi $genome_fa -bed - -s > $out/motif/peaks.fa
meme-chip -oc $out/motif/meme \
  -dna -meme-nmotifs 20 -meme-minw 5 -meme-maxw 20 -meme-p 4 \
  -dreme-e 0.05 -meme-maxsize 1000000 \
  $out/motif/peaks.fa

# Functional annotation
RCAS --bed $out/peaks/clip_peaks.bed --gtf $gtf --genome ${genome_fa} --out $out/rcas
```

## 9) Troubleshooting low mapping rate
- Confirm species/index matches; try mapping 10k reads to contamination DB (Kraken2).  
- Re-run cutadapt with correct adapters; inspect FastQC adapter content and read length.  
- Try local alignment (`--very-sensitive-local` Bowtie2) or relaxed clipping in STAR.  
- Inspect unaligned reads for low complexity; filter with `prinseq`/`fastp`.  
- If warnings like “Exhausted best-first chunk memory”, raise Bowtie1 `--chunkmbs` (e.g., 512–1024).

---

# nf-core/clipseq (Nextflow)
**Status:** DSL1; requires Nextflow ≤22.10.6.

## Samplesheet
Single-end only. CSV with headers:
```
sample,fastq
exp1_rep1,/data/exp1_R1.fastq.gz
exp1_rep2,/data/exp1_R2.fastq.gz
ctrl_rep1,/data/ctrl_R1.fastq.gz
```
One row per replicate. Controls are provided separately via params (see below).

## Quick run
```bash
nextflow run nf-core/clipseq -r 1.0.0 \
  -profile docker,test \          # choose docker/singularity/conda or institute profile
  --input design.csv \
  --fasta /ref/genome.fa \
  --gtf /ref/genes.gtf \          # optional but recommended
  --clip_samples exp1_rep1,exp1_rep2 \
  --control_samples ctrl_rep1 \
  --outdir results
```

## Major steps (defaults)
- Cutadapt trim adapters
- Optional rRNA/tRNA depletion (Bowtie2 pre-map)
- STAR alignment
- UMI-tools deduplication
- Crosslink site calling and bedGraph generation
- Peak calling: iCount, Paraclu, PureCLIP, Piranha (selectable)
- Motif: DREME
- QC: FastQC, Preseq, RSeQC, MultiQC

## Useful params
- `--aligner {star}` (default)  
- `--skip_rrna_removal` to disable pre-mapping  
- `--umitools_bc_pattern` if UMIs differ from default  
- `--outdir` where results go  
- `--max_cpus / --max_memory` to fit cluster limits  
- `-resume` to continue failed runs

## Tips
- Keep Nextflow version at/below 22.10.6 (DSL1).  
- For Singularity on HPC: `-profile singularity` and set `NXF_SINGULARITY_CACHEDIR`.  
- Use `-profile test` first to validate environment.  
- Ensure fasta/gtf/index versions match; mismatches drive low mapping and spurious peaks.
