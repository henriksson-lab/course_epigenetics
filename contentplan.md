# Content plan: epigenetics

Working list of the course's teaching units. It is a planning document, not a
contract. **We do not plan days or weeks here.** We write all the material we
need and decide the timetable later. Each unit is one folder (`<slug>/` with
`lec.qmd`, `exercise.qmd`, optional `quiz.qmd`). Order is set in `_quarto.yml`.

Status legend:
- **core**: decided, part of the course (all 31 units below, as of 2026-11-17)
- **suggested**: proposed, waiting for a decision
- **further**: candidate for the "Further reading" section (not examined)
- **written**: lecture, exercise and quiz exist (update as units are built)

Open questions are marked **[?]**.

---

## Block A: Foundations

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `chromatin-intro` | What is epigenetics? Chromatin, nucleosomes, histone marks, DNA methylation | written | Definitions (mitotic heritability vs "everything on top of the sequence"); writers/readers/erasers; euchromatin vs heterochromatin; the assays of the course in one map |
| `gene-regulation` | Promoters, enhancers, insulators, transcription factors | written | Pioneer factors; TF binding motifs; enhancer–promoter communication; sets up ATAC, ChIP and Hi-C |
| `genomic-intervals` | Genomic intervals on the command line and in R | written | BED/narrowPeak/bigWig/BAM; 0- vs 1-based; `bedtools`; `GenomicRanges`; genome builds and liftOver |
| `enrichment-alignment` | From reads to coverage for enrichment assays | written | Bowtie2; MAPQ and multi-mappers; duplicates; blacklist; paired-end fragments; coverage and normalisation (CPM, RPGC) with deepTools; IGV |
| `epigenomics-design` | Experimental design and QC for epigenomic assays | written | Replicates, input/IgG controls, spike-ins; QC metrics shared across assays (FRiP, fingerprint, library complexity); ENCODE standards |

## Block B: Chromatin accessibility (ATAC-seq)

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `atacseq-intro` | ATAC-seq: Tn5, fragment sizes, nucleosome ladder | written | Also DNase-seq and MNase-seq in brief; TSS enrichment; mitochondrial reads |
| `atacseq-analysis` | ATAC-seq processing and peak calling | written | Tn5 shift; MACS3 for ATAC; consensus peaks; nf-core/atacseq |
| `differential-accessibility` | Differential accessibility | written | Counts in peaks → DESeq2 (what they know); normalisation pitfalls (background vs peak-based); DiffBind/csaw |
| `motifs-footprinting` | Motifs, motif enrichment and TF footprinting | written | PWMs; HOMER/MEME; chromVAR deviations; footprinting (TOBIAS) and its caveats |

## Block C: Protein–DNA binding and histone marks (ChIP-seq and successors)

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `chipseq-intro` | ChIP-seq: the assay, TFs vs histone marks | written | Crosslinking, sonication, antibodies; narrow vs broad signals; input control |
| `chipseq-peaks` | Peak calling | written | The Poisson background model (local lambda); MACS3; broad-domain callers; IDR for replicates |
| `cutrun-cuttag` | CUT&RUN and CUT&Tag | written | Why they replaced ChIP for many uses; low background; SEACR; different QC |
| `chipseq-downstream` | Annotation and differential binding | written | Peak-to-gene annotation (ChIPseeker); profiles/heatmaps (deepTools); DiffBind; spike-in normalisation for global changes (e.g. EZH2 inhibitor) |
| `chromatin-states` | Combining marks into chromatin states | written | ChromHMM (an HMM over mark combinations); ENCODE/Roadmap annotations; using public reference epigenomes |

## Block D: DNA methylation

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `methylation-biology` | DNA methylation biology | written | DNMTs, TETs, 5mC/5hmC; CpG islands; imprinting; X inactivation; methylation in development and cancer; plants (CHG/CHH, RdDM) |
| `methylation-bisulfite` | Bisulfite and enzymatic sequencing (WGBS, RRBS, EM-seq) | written | Bisulfite chemistry; Bismark alignment; M-bias; coverage; calling per-CpG methylation |
| `methylation-differential` | Differential methylation | written | Beta-binomial models; DMLs vs DMRs; smoothing (BSmooth/DSS/methylKit) |
| `methylation-arrays` | Methylation arrays (Illumina EPIC) | written | Most common in clinical and population studies; beta vs M values; normalisation; cell-type composition as a confounder (link to their linear-model knowledge); minfi/sesame |
| `epigenetic-clocks` | Epigenetic clocks and biomarkers | written | Horvath/GrimAge as penalised regression; liquid biopsy / cfDNA methylation; critical view |
| `longread-epigenomics` | Methylation and chromatin from long reads | written | Nanopore/PacBio 5mC calling (modkit); haplotype-resolved and allele-specific methylation; Fiber-seq |

## Block E: 3D genome

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `genome-3d-intro` | 3D genome organisation | written | Territories, A/B compartments, TADs, loops; cohesin/CTCF loop extrusion; imaging (FISH) vs sequencing views |
| `hic-processing` | Hi-C: from reads to contact matrices | written | Proximity ligation; pairs; resolution vs depth; `.cool`/`.hic`; balancing (ICE); distance decay |
| `hic-features` | Calling compartments, TADs and loops | written | Eigenvector of the correlation matrix (compartments); insulation score; loop callers; differential Hi-C; Micro-C, HiChIP, Capture-C in brief |

## Block F: Integration and interpretation

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `integration-rnaseq` | Linking regulatory elements to genes | written | Combining peaks with RNA-seq; peak-to-gene links; ABC model; contacts from Hi-C |
| `noncoding-variants` | Epigenomics and genetic variation | written | Why GWAS hits sit in enhancers; caQTL/meQTL; allele-specific accessibility; enrichment of variants in cell-type annotations |
| `epigenetic-inheritance` | Epigenetic inheritance: what is really shown? | core | Transgenerational inheritance in plants vs mammals; reprogramming waves; critical reading of claims. A good "judgement" unit |
| `cancer-epigenetics` | Epigenetics in disease and therapy | core | IDH/TET, CIMP, histone mutations (H3K27M), EZH2/DNMT inhibitors; methylation-based tumour classification |

## Block G: Machine learning for prediction

Conceptual, not technical: **what are people trying to predict from DNA
sequence and epigenomic data, how does a model like a CNN do it, and how do we
know a prediction is real?** No specific deep-learning packages (Enformer,
chromBPNet, …) are taught as tools; they appear as examples of what is being
attempted. Students train small CNNs themselves on simple tasks.

| Slug | Unit | Status | Notes |
|---|---|---|---|
| `ml-basics` | Machine learning from scratch | written | Students have **no prior ML**. Supervised learning: features, labels, a model with parameters fitted to data; a first classifier (logistic regression, which extends the linear models they know); training vs test data; overfitting and why we hold data out; accuracy, ROC and PR curves. Exercise: predict open vs closed regions from GC content and a few motif counts in R |
| `ml-prediction-tasks` | What people try to predict, and how to tell if it works | written | A map of the tasks: TF binding and accessibility from sequence; expression from sequence (long-range models); effects of non-coding variants; enhancer–gene links; imputing missing marks; 3D contacts from sequence; age and tumour type from methylation. Genomics-specific evaluation traps: hold out whole chromosomes, leakage between overlapping regions, always compare with a simple baseline |
| `ml-cnn-dna` | How a CNN reads DNA | written | One-hot encoding; a filter is a motif scanner (link to PWMs in `motifs-footprinting`); pooling; from logistic regression to a neural network (layers of the same idea); training by gradient descent in pictures. Exercise: train a small CNN to tell bound from unbound sequences (e.g. CTCF peaks vs GC-matched background), look at what the filters learned |
| `ml-interpretation-variants` | What did the model learn? Variant effects | written | In-silico mutagenesis (doable by hand on a toy model); predicted effect = prediction(alt) − prediction(ref); why the big models exist and what they still get wrong (e.g. differences between individuals). Exercise: mutate sequences through the CNN from `ml-cnn-dna` |

Notes for this block:
- **Language.** R with the [`torch`](https://torch.mlverse.org/) package
  (decided 2026-11-17): trains small CNNs on a CPU, no Python needed.
- **Compute.** A small CNN on one chromosome trains on a CPU in minutes; no GPU
  needed.
- **Exam.** Concept questions: what a filter detects, which split leaks, read a
  ROC/PR curve, compute an in-silico mutagenesis effect from given scores.

---

## Cross-cutting

- **Key papers** per unit as optional further reading (designprinciples §9);
  seed list in `plan/key-papers.md`.
- **Glossary** (`dictionary.qmd`): grow it as units are written.
- **Datasets** (decided 2026-11-17): the shared system is the **ENCODE cell
  lines GM12878 and K562** (ATAC, histone marks, TFs incl. CTCF, RNA-seq, WGBS,
  Hi-C), subset to one chromosome for exercises. Record each file's ENCODE
  accession, genome build (GRCh38) and download command in the unit that first
  uses it. **[?]** For differential analyses with biological replicates
  (differential accessibility/binding/methylation) we still need a second,
  treatment-or-time-course dataset; choose it when writing
  `differential-accessibility`.
- **Compute**: students have access to an HPC cluster with SLURM (decided
  2026-11-17). Exercises are written as plain shell commands that run the
  same on a laptop, a server or inside a SLURM job; show once how to wrap
  them in an `sbatch` script, and don't make units depend on a particular cluster.
- **Examination**: an individual **written exam only** (decided 2026-11-17).
  Every unit therefore needs pen-and-paper practice: reading figures, small
  calculations, and judging designs (designprinciples §6). Computer exercises
  build understanding but are not examined directly.
