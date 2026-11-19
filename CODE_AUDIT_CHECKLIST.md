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

## Build and presentation code

- [ ] `Makefile`
- [ ] `_quarto.yml`
- [ ] `canvas-quiz/build.py` — adapted from quantbio-3MB00A to one quiz per topic slug (2026-11-17)
- [ ] `demos/build.sh`

## Interactive demos (TypeScript)
