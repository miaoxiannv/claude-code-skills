---
name: clipseq-workflow
description: "End-to-end CLIP/eCLIP (and related CLIP/iCLIP/PAR-CLIP) processing from raw FASTQ to peaks/motifs: QC, adapter/UMI handling, alignment (STAR/Bowtie2/Bowtie1), UMI-based deduplication, strand-aware coverage tracks, peak calling (e.g., PEAKachu/PureCLIP), and motif/functional analysis. Trigger when analyzing CLIP-seq data, troubleshooting low mapping rates, or generating bigWig/peak outputs."
---

# CLIP-seq Workflow

Concise operating guide for CLIP/eCLIP-style datasets **without Nextflow**. Run the steps yourself with the listed tools: trim adapters/UMIs, map reads, deduplicate, call peaks, make bigWigs, and run motif/annotation.

## Quick Start (manual flow)
1. **Set project vars**: sample names, genome fasta/index, GTF, chrom sizes, output folders.  
2. **QC + trim**: FastQC → cutadapt (5'/3' adapters + UMI trim) → optional low-complexity filter (fastp/prinseq).  
3. **Align**: STAR (spliced) or Bowtie2 local (short inserts). eCLIP orientation: `R2 forward / R1 reverse`; pass to STAR/Bowtie as `--readFilesIn R2 R1`.  
4. **Dedup UMIs**: `umi_tools extract` (if present) → `umi_tools dedup --paired --method adjacency`.  
5. **QC after map**: `samtools flagstat`, deepTools `plotFingerprint` and `plotCorrelation`.  
6. **Peaks**: PEAKachu (needs control) or PureCLIP (no control, nucleotide resolution).  
7. **Motifs/annotation**: extend peaks ±20 bp → `bedtools getfasta` → MEME-ChIP; annotate with RCAS.  
8. **Tracks**: deepTools `bamCoverage --filterRNAstrand {forward,reverse} --normalizeUsing CPM` to bigWig.

Command templates for each step live in `references/pipeline.md` (bash-friendly loops over samples).

## Tool picks (swap as needed)
- **Adapters/UMIs**: cutadapt for trimming; `umi_tools extract` when UMIs exist. For heavy adapter content/low-complexity, add fastp/prinseq.  
- **Aligner**: STAR for spliced/transcriptome-aware; Bowtie2 local for short inserts or heavy adapter bleed; Bowtie1 only if legacy pipeline required.  
- **Peak callers**: PEAKachu (requires control; DESeq2-based), PureCLIP (single-sample, nucleotide-level), CLIPper/Piranha as alternates.  
- **Motif/annotation**: MEME-ChIP for motifs; RCAS for functional annotation.  
- **Tracks**: deepTools `bamCoverage` strand-specific CPM.

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
