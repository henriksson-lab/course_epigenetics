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

## Build and presentation code

- [ ] `Makefile`
- [ ] `_quarto.yml`
- [ ] `canvas-quiz/build.py` — adapted from quantbio-3MB00A to one quiz per topic slug (2026-11-17)
- [ ] `demos/build.sh`

## Interactive demos (TypeScript)
