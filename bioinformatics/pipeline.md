# LACE / PARN-LACE-seq Processing (single-end)

Concise, reproducible pipeline to trim adapters/poly(A), remove rRNA, map with Bowtie1, and deliver QC/bigWig outputs. Assumes single-end ~50 bp reads (typical for LACE-seq). Adapt assemblies/paths as needed.

## 0) Inputs & environment
- Reads: `*.fq.gz`, single-end.
- Assembly: use the assembly stated by the study (e.g., mm10 or hg38). Mixing mm9/mm10 or hg19/hg38 is a common source of off-chrom "isolated gene" false positives.
- Conda env needs `cutadapt`, `bowtie`, `samtools`, `deeptools`.

```bash
conda activate lace_seq   # rename to your env
```

## 1) Identify adapters (sanity check)
Run on 1–2k reads to confirm ligated adapters and random bases:
```bash
gzip -cd sample_1.fq.gz | head -n 2000 | grep -o "ATCTCGTATGCCGTCTTCTGCTT" | head
```
Typical layout (Illumina P7 + LACE 3' adapter + polyA):
```
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTAAAAAAAA
^^^^ Illumina          ^^^^^^^^^^^^^^^^^^^^^ LACE 3' adapter
```
Some libraries add 4 random nt at read 5' to reduce ligation bias.

## 2) Adapter & poly(A) trimming (cutadapt)
Two variants: remove or keep the leading 4 nt. Keep both to compare mapping.
```bash
# A) Drop 5' random 4 nt (recommended if present)
cutadapt -j 16 \
  -u 4 \                             # trim first 4 nt (omit if none)
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTT \
  -q 30,0 --max-n 0.25 --trim-n -m 18 \
  -o sample.noAdapt.Cut4bp.fq sample_1.fq.gz

cutadapt -j 16 -a "A{15}" -n 2 -m 18 \
  -o sample.clean.Cut4bp.fq sample.noAdapt.Cut4bp.fq

# B) Keep 5' bases (control)
cutadapt -j 16 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTT \
  -q 30,0 --max-n 0.25 --trim-n -m 18 \
  -o sample.noAdapt.Keep4bp.fq sample_1.fq.gz

cutadapt -j 16 -a "A{15}" -n 2 -m 18 \
  -o sample.clean.Keep4bp.fq sample.noAdapt.Keep4bp.fq
```
Notes
- Remove `-u 4` if your library has no random bases.
- Use `--report=full` if you need trimming stats.

## 3) rRNA depletion with Bowtie1
Prepare species-matched rRNA FASTA (e.g., pre-rRNA from NCBI or Ensembl) and build index once:
```bash
bowtie-build rRNA.fa rRNA
```
Filter rRNA:
```bash
bowtie -v 2 -p 12 --un sample.no_rRNA.fq rRNA sample.clean.Cut4bp.fq > /dev/null
```
`sample.no_rRNA.fq` carries reads that did **not** align to rRNA.

## 4) Genome alignment (Bowtie1, single-end)
```bash
GENOME=mm10   # or hg38/hg19/mm9; must match your rRNA source and annotations
bowtie -v 2 -m 10 --best --strata -p 12 -S \
  $GENOME sample.no_rRNA.fq > sample.sam

samtools view -@ 12 -bS sample.sam | samtools sort -@ 12 -o sample.sorted.bam
samtools index sample.sorted.bam
rm sample.sam
```
- `-v 2`: allow ≤2 mismatches; adjust if read length differs.
- `-m 10`: discard reads with >10 genomic hits (repeats).

## 5) QC & deliverables
- Mapping: `samtools flagstat sample.sorted.bam`, `samtools idxstats`.
- Replicate concordance: `multiBamSummary bins` + `plotCorrelation` (deepTools).
- Tracks: strand-agnostic coverage (LACE is not strand-specific by default):
```bash
bamCoverage -b sample.sorted.bam -o sample.CPM.bw --normalizeUsing CPM --binSize 10
```
- Peak / site calling: follow study-specific tool (e.g., PureCLIP/CLIPper) if needed.

## 6) Repeat element enrichment (optional)
1) Download RepeatMasker BED for your assembly (UCSC table browser -> rmsk -> BED6).
2) Intersect peaks or aligned reads with repeats:
```bash
bedtools intersect -wa -u -a peaks.bed -b rmsk.bed \
  | cut -f4 \
  | sort | uniq -c | sort -k1,1nr > repeat_enrichment.txt
```
3) For read-level enrichment use `bedtools coverage -a rmsk.bed -b sample.sorted.bam`.

## 7) Metaplot around features (deepTools)
Use peaks or TSS as reference regions:
```bash
# build matrix over TSS (example Gencode)
computeMatrix reference-point \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  -S sample.CPM.bw \
  -R genes.tss.bed \
  -out tss.mat.gz

plotProfile -m tss.mat.gz -out tss_metaplot.pdf --perGroup
```
For gene-body plots, switch to `scale-regions` with `-m 3000 -b 2000 -a 2000`.

## 8) Motif analysis around binding sites
1) Generate summit-centered FASTA (e.g., +/-50 bp):
```bash
bedtools slop -i peaks.bed -g chrom.sizes -b 50 > peaks.100bp.bed
bedtools getfasta -fi genome.fa -bed peaks.100bp.bed -fo peaks.100bp.fa
```
2) Run MEME-ChIP:
```bash
meme-chip -oc meme_out -meme-nmotifs 10 -db JASPAR.meme peaks.100bp.fa
```
   or HOMER:
```bash
findMotifsGenome.pl peaks.100bp.bed genome.fa homer_out -len 6,8,10
```
3) Use central enrichment (centrimo) or E-value to pick candidates; cross-check with known RNA-binding motifs if applicable.

## 6) SLURM example (trimming → rRNA → genome)
```bash
#!/bin/bash
#SBATCH --job-name=LACEseq
#SBATCH --partition=Fnode2
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --output=%j.out
#SBATCH --error=%j.err
set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate lace_seq

READ=sample_1.fq.gz
GENOME=/path/to/mm10    # bowtie prefix
RRNA=/path/to/mouse_pre_rRNA

cutadapt -j 20 -u 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTT \
  -q 30,0 --max-n 0.25 --trim-n -m 18 -o sample.noAdapt.fq $READ
cutadapt -j 20 -a "A{15}" -n 2 -m 18 -o sample.clean.fq sample.noAdapt.fq

bowtie -v 2 -p 20 --un sample.no_rRNA.fq $RRNA sample.clean.fq > /dev/null
bowtie -v 2 -m 10 --best --strata -p 20 -S $GENOME sample.no_rRNA.fq > sample.sam
samtools view -@ 20 -bS sample.sam | samtools sort -@ 20 -o sample.sorted.bam
samtools index sample.sorted.bam
rm sample.sam
```

## 9) Common pitfalls
- Adapter mismatch (e.g., mm9 vs mm10 index) → off-chrom/"isolated gene" hits; always match assemblies across genome, rRNA, and annotation.
- Forgetting polyA trimming leaves long tails, causing alignment failures.
- Missing rRNA filter inflates duplicate/low-complexity reads and lowers unique mapping.
- Wrong strand assumption: LACE-seq often treated as unstranded; check library notes before making stranded bigWigs.
