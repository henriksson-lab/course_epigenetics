# Code correctness audit checklist

Scope: every source document with an executable code fence (R or shell), the
figure scripts (`<unit>/make-*.R`), the demos, and the build configuration.
Generated `docs/` output is excluded; it is audited by rendering the source.

Mark an item complete only after a reviewer has checked every executable chunk
for syntax, dependencies (is the tool/package in `intro/install.qmd`?),
data-shape assumptions, genome build and coordinate conventions (0- vs 1-based),
numerical correctness, and consistency with its surrounding teaching material.
Record fixes in the Git history.

## Executable teaching documents

<!-- - [ ] `<unit>/lec.qmd` — added YYYY-MM-DD; needs audit -->
- [ ] `chromatin-intro/lec.qmd` — added 2026-11-17; needs audit
- [ ] `chromatin-intro/exercise.qmd` — added 2026-11-17; all R chunks run as written against data/ (2026-11-17); needs independent audit
- [ ] `chromatin-intro/make-chromatin-figs.R` — added 2026-11-17; needs audit
- [ ] `genomic-intervals/lec.qmd` — added 2026-11-17; needs audit (bash code slide)
- [ ] `genomic-intervals/exercise.qmd` — added 2026-11-17; bash and R chunks run against data/, numbers in solutions from those runs (2026-11-17); needs independent audit
- [ ] `genomic-intervals/make-intervals-figs.R` — added 2026-11-17; needs audit (calls bedtools)
- [ ] `data/download.sh` — added 2026-11-17
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
