---
name: laceseq-workflow
description: "End-to-end LACE/PARN-LACE-seq preprocessing: adapter + poly(A) trimming (cutadapt), optional 5' random-base clipping, Bowtie rRNA filtering, genome mapping, and QC (flagstat, bigWig, replicate correlation). Use when handling LACE-seq single-end reads or troubleshooting low mapping/off-chrom \"isolated gene\" hits caused by assembly mismatch or incomplete trimming."
---

# LACE-seq Workflow

Concise steps to clean LACE/PARN-LACE-seq single-end reads (~50 bp), remove rRNA, map with Bowtie1, and produce QC/coverage outputs, plus repeat/metaplot/motif downstreams. Optimized for single-end LACE libraries; swap parameters explicitly if your library differs.

## Inputs & prerequisites
- Reads: single-end FASTQ (gz ok), ~50 bp typical; if paired, adjust Bowtie flags.
- Indexes: Bowtie1 genome + rRNA built from the **same assembly** (mm10↔mm10, hg38↔hg38).
- Tools: cutadapt, bowtie (v1), samtools, deeptools, bedtools, MEME Suite or HOMER.

## Quick Start (per sample)
- Inspect ~2k reads to verify adapters and optional 5' 4 nt.
- Trim adapters + poly(A) with cutadapt; if unsure about 5' random bases, run both “Cut4bp” and “Keep4bp” branches.
- rRNA filter with Bowtie1 (`--un` keeps desired reads) → genome align (`-v 2 -m 10 --best --strata`), then sort/index BAM.
- QC: `flagstat`, `idxstats`, deepTools correlation; generate CPM bigWig via `bamCoverage`.
- Downstream (optional): repeat enrichment, metaplot with outlier filtering, motif discovery (MEME-ChIP/HOMER).

## Core Commands (per sample)
```bash
# trim (drop 5' 4 nt if present)
cutadapt -j 16 -u 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTT \
  -q 30,0 --max-n 0.25 --trim-n -m 18 -o sample.noAdapt.fq sample_1.fq.gz
cutadapt -j 16 -a "A{15}" -n 2 -m 18 -o sample.clean.fq sample.noAdapt.fq

# rRNA filter (Bowtie1)
bowtie -v 2 -p 12 --un sample.no_rRNA.fq /path/to/rRNA sample.clean.fq > /dev/null

# genome align (Bowtie1, single-end)
bowtie -v 2 -m 10 --best --strata -p 12 -S /path/to/genome sample.no_rRNA.fq > sample.sam
samtools view -@12 -bS sample.sam | samtools sort -@12 -o sample.sorted.bam
samtools index sample.sorted.bam
```

## Outputs & QC
- Mapping stats: `samtools flagstat`, `samtools idxstats`.
- Replicate concordance: deepTools `multiBamSummary` + `plotCorrelation`.
- Coverage: `bamCoverage -b sample.sorted.bam -o sample.CPM.bw --normalizeUsing CPM --binSize 10` (unstranded unless library says otherwise).
- Repeat enrich: intersect peaks with RepeatMasker BED; see reference.
- Metaplot: deepTools `computeMatrix` + `plotProfile`; optionally `computeMatrixOperations filterValues --max 0.2` to drop outlier bins.
- Motif: extract +/-50 bp around peaks, run MEME-ChIP or HOMER.

## Troubleshooting
- Off-chrom spikes: usually assembly mismatch (mm9/mm10, hg19/hg38) or untrimmed poly(A); realign with consistent indexes and confirm trimming.
- Low unique mapping: check adapter correctness, ensure rRNA FASTA matches species, and confirm `-m 10` isn’t over-filtering very repetitive regions.
- Strand assumptions: LACE-seq often treated as unstranded; only make stranded tracks if protocol notes specify.

## Reference
- Detailed pipeline (trim variants, rRNA/genome map), repeat enrichment, metaplot with outlier filtering, motif steps, SLURM script, and pitfalls: [references/pipeline.md](references/pipeline.md)
