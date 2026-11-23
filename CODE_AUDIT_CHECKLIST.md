# Code correctness audit checklist

Scope: every source document with an executable code fence (R or shell), the
figure scripts (`<unit>/make-*.R`), the demos, and the build configuration.
Generated `docs/` output is excluded; it is audited by rendering the source.

Mark an item complete only after a reviewer has checked every executable chunk
for syntax, dependencies (is the tool/package in `intro/install.qmd`?),
data-shape assumptions, genome build and coordinate conventions (0- vs 1-based),
numerical correctness, and consistency with its surrounding teaching material.
Record fixes in the Git history.

## Review on 2026-11-23

A parallel review corrected scientific descriptions, worked answers and code
across the course. This is a partial audit: unchecked items below still need
the complete independent review specified above.

- Corrected interval/summit and Tn5 insertion coordinates, local-background
  Poisson calculations, an impossible peak-selection task, inconsistent QC
  questions, and empirical permutation-p-value reporting.
- Distinguished standard PBAT from four-strand non-directional libraries;
  corrected dmrseq's model and support for two replicates per group, methylation
  averaging, array coverage, and germline methylation-reprogramming exceptions.
- Updated the cancer unit's tazemetostat status using the
  [June 2026 FDA withdrawal notice](https://www.govinfo.gov/content/pkg/FR-2026-06-22/pdf/2026-12367.pdf)
  and qualified the claim that tumour MLH1 methylation excludes inherited risk.
- Corrected tied-score PR calculations and regenerated the integration PR plot;
  corrected CNN arithmetic, sequence-shuffling and attribution explanations,
  and claims about allele-specific binding and mapping-bias correction.
- Corrected glossary definitions, a draft-syllabus label, missing preparation
  dependencies, reference/index creation order, and topic-specific dispatch of
  the methylation and Hi-C preparation scripts.

Validation includes syntax checks of 205 R, 163 shell and 18 Python examples
and 28 R figure scripts; TypeScript type checking; full runs of the chromatin-intro, methylation-biology
and epigenetic-inheritance R exercises; genomic-intervals R1–R4; and the full
integration figure script. The downloader passed isolated fresh-directory
checks for ten affected topic selections and the whole-course selection,
including missing-index recovery, using a tiny real bgzipped FASTA and mocked
downloads/preparation jobs. Full public-data downloads were not repeated.
Several optional R/Python packages are absent locally; syntax checks do not
establish that every exercise runs end to end.
All 96 site pages render; the generated site has no missing local file links,
assets or cross-page anchors. All source SVGs parse as XML.
The footprint/bias plotting blocks were also rerun after realigning the
supplied ATAC reads, using the official JASPAR MA0139.2 matrix; the corrected
figures and their quoted summary statistics were regenerated together.

Sources for method corrections include the [Bismark library options](https://felixkrueger.github.io/Bismark/options/alignment/),
[dmrseq vignette](https://bioconductor.org/packages/release/bioc/vignettes/dmrseq/inst/doc/dmrseq.html),
[chromVAR vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/chromVAR/inst/doc/Introduction.html),
[Infinium assay documentation](https://support-docs.illumina.com/ARR/Inf_HD_Methylation/Content/ARR/Methylation/Overview_mMeth.htm),
[average precision definition](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html),
[DeepLIFT paper](https://arxiv.org/abs/1704.02685), and
[WASP paper](https://doi.org/10.1038/nmeth.3582).

The allele-specific-binding exercise retains its illustrative ChIP 42:18 and
input 31:29 counts. At the teacher's request, solution B5 now explains both
tests: the binomial test against 0.5 gives p = 0.00267, whereas a direct
two-sided Fisher comparison with input gives p = 0.06093. The solution explains
the different null hypotheses and the uncertainty in both sets of counts.

Also unresolved: the differential-methylation hand rules are described as DSS
defaults, but the [current upstream implementation](https://raw.githubusercontent.com/haowulab/DSS/master/R/DMR.R)
uses different run/gap logic and strict thresholds. The DSS version behind the
reported 5,061 DMRs is unspecified and is not installed here. Pin and rerun that
version before changing the real-data answers; discuss whether the hand rules
should instead be explicitly labelled as a teaching simplification.

The motif workflow is now standardised at the teacher's request: it reads the
saved differential_peaks.tsv from Differential accessibility instead of
refitting DESeq2 with different shrinkage/filtering settings. Both target groups
use the same `ns` background; GC-matched partners are sampled without replacement
from that pool. `ns` means neither differential class under the chosen thresholds,
not proof of equal accessibility. The differential figure script also exports
this table for the motif figure script, preserving `make figs` dependency order.

The canonical run reproduces 211,824 tested peaks: 19,404 B_up, 16,897 T_up and
175,523 ns. Its figure-script export matches the exercise R11 table exactly.
Enrichment against ns yields 511 enriched motifs in B_up and 507 in T_up at BH
5%; the affected figures, worked answers and SVG table are regenerated from
those results. The stricter liftOver filter removes one interchromosomal peak,
but the ordered top 100,000 peaks used by chromVAR are identical, so its existing
figures retain the same inputs.

## Executable teaching documents

<!-- - [ ] `<unit>/lec.qmd` — added YYYY-MM-DD; needs audit -->
- [ ] `chromatin-intro/lec.qmd` — added 2026-11-17; needs audit
- [ ] `chromatin-intro/exercise.qmd` — added 2026-11-17; all R chunks run as written against data/ (2026-11-17); needs independent audit
- [ ] `chromatin-intro/make-chromatin-figs.R` — added 2026-11-17; needs audit
- [ ] `genomic-intervals/lec.qmd` — added 2026-11-17; needs audit (bash code slide)
- [ ] `genomic-intervals/exercise.qmd` — added 2026-11-17; bash and R chunks run against data/, numbers in solutions from those runs (2026-11-17); needs independent audit
- [ ] `genomic-intervals/make-intervals-figs.R` — added 2026-11-17; needs audit (calls bedtools)
- [x] `data/download.sh` — reviewed 2026-11-23; dependency order and topic dispatch fixed, isolated preparation/index-recovery checks passed
- [ ] `gene-regulation/lec.qmd` — added 2026-11-18; needs audit
- [ ] `gene-regulation/exercise.qmd` — added 2026-11-18; chunks run as written by its author; needs independent audit
- [ ] `gene-regulation/make-generegulation-figs.R` — added 2026-11-18; needs audit
- [ ] `enrichment-alignment/lec.qmd` — added 2026-11-18; needs audit
- [ ] `enrichment-alignment/exercise.qmd` — added 2026-11-18; chunks run as written by its author; needs independent audit
- [ ] `enrichment-alignment/make-enrichalign-figs.R` — added 2026-11-18; needs audit
- [ ] `epigenomics-design/lec.qmd` — added 2026-11-18; needs audit
- [ ] `epigenomics-design/exercise.qmd` — added 2026-11-18; chunks run as written by its author; command-line and R parts re-run independently 2026-11-19, numbers agree; needs independent audit
- [ ] `epigenomics-design/make-design-figs.R` — added 2026-11-18; needs audit
- [ ] `atacseq-intro/lec.qmd` — added 2026-11-18; needs audit
- [ ] `atacseq-intro/exercise.qmd` — added 2026-11-18; chunks run as written by its author; needs independent audit
- [ ] `atacseq-intro/make-atacintro-figs.R` — added 2026-11-18; needs audit
- [ ] `atacseq-analysis/lec.qmd` — added 2026-11-18; needs audit
- [ ] `atacseq-analysis/exercise.qmd` — added 2026-11-18; chunks run as written by its author; needs independent audit
- [ ] `atacseq-analysis/make-atacanalysis-figs.R` — added 2026-11-18; needs audit
- [ ] `differential-accessibility/lec.qmd` — added 2026-11-18; needs audit
- [ ] `differential-accessibility/exercise.qmd` — added 2026-11-18; chunks run as written by its author; needs independent audit
- [ ] `differential-accessibility/make-diffacc-figs.R` — added 2026-11-18; needs audit
- [ ] `motifs-footprinting/lec.qmd` — added 2026-11-18; needs audit
- [ ] `motifs-footprinting/exercise.qmd` — added 2026-11-18; chunks run as written by its author; needs independent audit
- [ ] `motifs-footprinting/make-motifs-figs.R` — added 2026-11-18; needs audit
- [ ] `data/prepare-reads.sh`, `data/filter_fastq_by_names.py` — added 2026-11-18 (read-name `/1` fix 2026-11-18)
- [ ] `chipseq-intro/lec.qmd` — added 2026-11-19; needs audit
- [ ] `chipseq-intro/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `chipseq-intro/make-chipintro-figs.R` — added 2026-11-19; needs audit
- [ ] `chipseq-peaks/lec.qmd` — added 2026-11-19; needs audit
- [ ] `chipseq-peaks/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `chipseq-peaks/make-chippeaks-figs.R` — added 2026-11-19; needs audit
- [ ] `cutrun-cuttag/lec.qmd` — added 2026-11-19; needs audit
- [ ] `cutrun-cuttag/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `cutrun-cuttag/make-cuttag-figs.R` — added 2026-11-19; needs audit
- [ ] `chipseq-downstream/lec.qmd` — added 2026-11-19; needs audit
- [ ] `chipseq-downstream/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `chipseq-downstream/make-chipdown-figs.R` — added 2026-11-19; needs audit
- [ ] `chromatin-states/lec.qmd` — added 2026-11-19; needs audit
- [ ] `chromatin-states/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `chromatin-states/make-chromstates-figs.R` — added 2026-11-19; needs audit
- [ ] `methylation-biology/lec.qmd` — added 2026-11-19; needs audit
- [ ] `methylation-biology/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `methylation-biology/make-methbio-figs.R` — added 2026-11-19; needs audit
- [ ] `methylation-bisulfite/lec.qmd` — added 2026-11-19; needs audit
- [ ] `methylation-bisulfite/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `methylation-bisulfite/make-bisulfite-figs.R` — added 2026-11-19; needs audit
- [ ] `methylation-differential/lec.qmd` — added 2026-11-19; needs audit
- [ ] `methylation-differential/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `methylation-differential/make-diffmeth-figs.R` — added 2026-11-19; needs audit
- [ ] `methylation-arrays/lec.qmd` — added 2026-11-19; needs audit
- [ ] `methylation-arrays/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `methylation-arrays/make-arrays-figs.R` — added 2026-11-19; needs audit
- [ ] `epigenetic-clocks/lec.qmd` — added 2026-11-19; needs audit
- [ ] `epigenetic-clocks/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `epigenetic-clocks/make-clocks-figs.R` — added 2026-11-19; needs audit
- [ ] `longread-epigenomics/lec.qmd` — added 2026-11-19; needs audit
- [ ] `longread-epigenomics/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `longread-epigenomics/make-longread-figs.R` — added 2026-11-19; needs audit
- [ ] `genome-3d-intro/lec.qmd` — added 2026-11-19; needs audit
- [ ] `genome-3d-intro/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `genome-3d-intro/make-genome3d-figs.py` — added 2026-11-19; needs audit
- [ ] `hic-processing/lec.qmd` — added 2026-11-19; needs audit
- [ ] `hic-processing/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `hic-processing/make-hicproc-figs.py` — added 2026-11-19; needs audit
- [ ] `hic-features/lec.qmd` — added 2026-11-19; needs audit
- [ ] `hic-features/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `hic-features/make-hicfeatures-figs.py` — added 2026-11-19; needs audit
- [ ] `integration-rnaseq/lec.qmd` — added 2026-11-19; needs audit
- [ ] `integration-rnaseq/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `integration-rnaseq/make-integ-figs.R` — added 2026-11-19; needs audit
- [ ] `noncoding-variants/lec.qmd` — added 2026-11-19; needs audit
- [ ] `noncoding-variants/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `noncoding-variants/make-variants-figs.R` — added 2026-11-19; needs audit
- [ ] `ml-basics/lec.qmd` — added 2026-11-19; needs audit
- [ ] `ml-basics/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `ml-basics/make-mlbasics-figs.R` — added 2026-11-19; needs audit
- [ ] `ml-prediction-tasks/lec.qmd` — added 2026-11-19; needs audit
- [ ] `ml-prediction-tasks/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `ml-prediction-tasks/make-mltasks-figs.R` — added 2026-11-19; needs audit
- [ ] `ml-cnn-dna/lec.qmd` — added 2026-11-19; needs audit
- [ ] `ml-cnn-dna/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `ml-cnn-dna/make-mlcnn-figs.R` — added 2026-11-19; needs audit
- [ ] `ml-interpretation-variants/lec.qmd` — added 2026-11-19; needs audit
- [ ] `ml-interpretation-variants/exercise.qmd` — added 2026-11-19; chunks run as written by its author; needs independent audit
- [ ] `ml-interpretation-variants/make-mlinterp-figs.R` — added 2026-11-19; needs audit
- [ ] `data/prepare-chip.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/prepare-cuttag.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/prepare-hic.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/prepare-integ.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/prepare-meth.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/prepare-meth2.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/prepare-ml2.sh` — added 2026-11-19; run end to end by its author
- [ ] `data/bam_slice_remote.py` — added 2026-11-19
- [ ] `plan/merge-manifests.py`, `plan/build-key-papers.py` — added 2026-11-19

## Build and presentation code

- [ ] `Makefile`
- [ ] `_quarto.yml`
- [ ] `canvas-quiz/build.py` — adapted from quantbio-3MB00A to one quiz per topic slug (2026-11-17)
- [ ] `demos/build.sh`

## Interactive demos (TypeScript)
