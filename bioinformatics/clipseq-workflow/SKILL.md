---
name: clipseq-workflow
description: "End-to-end CLIP/eCLIP (and related CLIP/iCLIP/PAR-CLIP) processing from raw FASTQ to peaks/motifs: QC, adapter/UMI handling, alignment (STAR/Bowtie2/Bowtie1), UMI-based deduplication, strand-aware coverage tracks, peak calling (e.g., PEAKachu/PureCLIP), and motif/functional analysis. Trigger when analyzing CLIP-seq data, troubleshooting low mapping rates, or generating bigWig/peak outputs."
---

# CLIP-seq Workflow

Concise operating guide for CLIP/eCLIP-style datasets. Use this when you need to trim adapters/UMIs, map reads, deduplicate, call peaks, make bigWigs, and run motif/annotation steps.

## Quick Start (happy path)
1. **Set project vars**: samples, genome fasta + index, annotation, chrom sizes.  
2. **QC + trim**: FastQC → cutadapt (remove 5'/3' adapters; trim UMI if present) → optional low-complexity filter.  
3. **Align**: STAR (spliced) or Bowtie2/Bowtie1 (short inserts); remember eCLIP orientation: `--readFilesIn R2 R1` (R2 = forward).  
4. **Deduplicate**: `umi_tools dedup` (paired, adjacency, edit-distance 1).  
5. **QC**: `samtools flagstat`, `plotFingerprint`, correlation heatmap.  
6. **Peaks**: PEAKachu (with control) or PureCLIP (no control).  
7. **Motifs/annotation**: extend peaks ±20 bp → fasta → MEME-ChIP; annotate with RNA-Centric Annotation System (RCAS).  
8. **Tracks**: `bamCoverage --filterRNAstrand {forward,reverse} --normalizeUsing CPM` to bigWig.

Full command templates and parameters: see `references/pipeline.md`.

## nf-core/clipseq one-liner
- Pipeline: Nextflow DSL1 (requires Nextflow ≤22.10.6).  
- Samplesheet (single-end only): CSV with headers `sample,fastq` (one row per replicate). Example:  
  ```
  sample,fastq
  exp1_rep1,clip0001_01.fastq.gz
  exp1_rep2,clip0001_02.fastq.gz
  ```
- Run (uses containers):  
  ```bash
  nextflow run nf-core/clipseq -r 1.0.0 \
    -profile docker,test        # replace test with docker/singularity/... or institute profile
    --input design.csv \
    --fasta /path/genome.fa \
    --outdir results
  ```
- Defaults: Cutadapt → rRNA/tRNA pre-map (Bowtie2) → STAR align → UMI-tools dedup → crosslink/bedgraph → peak calling (iCount/Paraclu/PureCLIP/Piranha) → motif (DREME) → QC (FastQC/Preseq/RSeQC/MultiQC).

## Decision Notes
- **Adapters/UMIs**: Required for CLIP. If UMIs are in R2 5' (common eCLIP), extract with `umi_tools extract --bc-pattern=NNNNN --read2-in`. Drop first 5 bp of R1 to avoid UMI bleed-through if protocol matches ENCODE eCLIP.
- **Aligner choice**:  
  - STAR: best for transcriptome-aware spliced alignment; use for standard eCLIP/iCLIP.  
  - Bowtie2 local: good for very short inserts/adapter bleed.  
  - Bowtie1 legacy: only if pipeline depends on it; increase `--chunkmbs` if “Exhausted best-first chunk memory”.
- **Orientation**: eCLIP libraries are often “R2 forward / R1 reverse”. Swap mates accordingly in STAR/Bowtie inputs to keep strand interpretation correct.
- **Low mapping rate triage**:  
  1) Confirm species/index match; 2) rerun FastQC for adapter/length; 3) map a small subset locally with relaxed params (`--very-sensitive-local` in Bowtie2/STAR `--clip3pNbases`); 4) check contamination via Kraken2 on unmapped; 5) ensure UMI/adapter trimming was correct.
- **Peak caller selection**:  
  - PEAKachu: uses control, DESeq2 normalization; paired-end supported.  
  - PureCLIP: single-sample, single-nucleotide resolution; good when no control.  
  - CLIPper / Piranha are alternatives but not bundled here.

## Outputs to deliver
- Sorted/indexed BAM + QC (flagstat, duplication, plotFingerprint).  
- Strand-specific CPM bigWigs (plus/minus).  
- Peak file (BED/GFF) with statistics; optionally motif logos and RCAS report.

## References
Load `references/pipeline.md` for step-by-step commands, parameter defaults, and example snippets (cutadapt, STAR, umi_tools, PEAKachu, MEME-ChIP, bigWig generation).
