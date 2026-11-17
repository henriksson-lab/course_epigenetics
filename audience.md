# Audience

A **final-year MSc course** (second cycle) at Umeå University. Students take it
near the end of their master's, so we pitch it as a specialist course: they
can already analyse sequencing data, and the course teaches them the biology of
chromatin and how to analyse the assays that measure it.

## What we can expect them to know

### Biology
- Solid molecular biology, genetics and cell biology: DNA, transcription,
  transcription factors, the central dogma, basic gene regulation.
- What a promoter is; probably the word "enhancer"; probably not much about
  chromatin beyond "DNA wrapped around histones".
- *Recap, don't assume:* nucleosome structure, histone modifications, DNA
  methylation chemistry, 3D genome organisation, imprinting, X inactivation.

### Sequencing and analysis
- **Have done an RNA-seq analysis in R**: count matrix → DESeq2/edgeR (or
  limma-voom) → differential expression → volcano plot / heatmap → gene-set
  enrichment. We can reuse design formulas, size factors, dispersion,
  multiple-testing correction and log fold changes without re-teaching them.
- Know what a FASTQ file and an aligned BAM are, in outline; have probably run
  or seen an aligner (STAR, HISAT2, or a pseudo-aligner like salmon).
- *Recap, don't assume:* SAM flags, MAPQ, duplicates, paired-end fragment
  logic, coordinate systems (0- vs 1-based), interval files (BED), signal
  tracks (bigWig), genome browsers.

### Computing
- **Comfortable on the Linux command line**: navigating, pipes, `grep`, `awk`
  basics, running tools with flags, writing a small shell script, working on a
  remote server / HPC over ssh.
- Comfortable in R and Bioconductor (installing packages, data frames,
  `SummarizedExperiment` at least as a user via DESeq2).
- *Recap, don't assume:* conda/mamba environments, containers, workflow managers
  (nf-core / Snakemake), job schedulers (SLURM), `GenomicRanges`.

### Statistics
- Hypothesis testing, multiple testing (FDR), linear models with design
  matrices, the negative binomial as "the count model DESeq2 uses".
- **No machine learning assumed.** Supervised learning, train/test splits,
  overfitting and neural networks are all taught from scratch, conceptually.
- *Recap, don't assume:* the Poisson model for background in peak calling,
  binomial / beta-binomial models for methylation proportions, hidden Markov
  models (chromatin states), distance-decay and matrix balancing (Hi-C).

## Pedagogical implications
- **Biology and analysis together.** Each assay gets the biology it measures,
  the wet-lab principle (because artefacts come from the chemistry), the
  processing on the command line, and the downstream analysis in R.
- **Build on RNA-seq.** Every time a step is "just like RNA-seq" (differential
  binding with DESeq2, size factors) say so, and spend the time on what is
  *different* (no gene annotation to count on, peaks must be called, background
  and input controls, global changes that break normalisation).
- **One shared foundation** (file formats, intervals, alignment for enrichment
  assays, browsers, QC metrics) taught once, then reused by every assay.
- **Small data, real data.** Exercises use real public data (ENCODE, GEO,
  4DN, Roadmap) subset to one chromosome or region, so a laptop or a short
  HPC job is enough. The full-size pipeline is shown and, where possible, run
  once on HPC.
- **Judgement.** Final-year students should leave able to read a methods
  section critically: is there an input/IgG control, are there replicates, is
  the normalisation valid under a global change, is "epigenetic inheritance"
  really shown?
