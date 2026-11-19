# Parallel authoring: environment and rules for unit agents

Read together with `plan/unit-writing-brief.md`. Several agents write units at
the same time, so these rules keep them out of each other's way.

## Machine (shared, 40 cores, ~170 GB RAM, **disk is tight**)
- `/big` has ~150 GB free in total for everyone. Budget per agent: **≤ 10 GB of
  final data in `data/`** and ≤ 20 GB temporary, deleted when done. Prefer
  processed files, chr21 (or one region) subsets, and remote/indexed access.
- Use at most 8 threads per job.
- Tools: bowtie2, bwa, fastp (`~/bin/fastp`), bedtools, `/usr/bin/samtools`
  (1.21; the conda `samtools` first on PATH is 1.3.1 — too old), bgzip/tabix,
  deepTools (needs `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`), java 19,
  quarto, google-chrome, cairosvg, montage.
- MACS3: `$OLD/venv/bin/macs3`. Extra R packages (JASPAR2024, TFBSTools,
  motifmatchr, chromVAR, universalmotif): `R_LIBS=$OLD/Rlib`, where
  `OLD=/big/henriksson/temp/claude-tmp/claude-1001/-data-henriksson-github-course-epigenetics/b14f07af-1a89-496e-98e1-f6268cc854b7/scratchpad`.
- Install anything else you need **into your own scratch dir** (a python venv,
  an R library via `R_LIBS`, a downloaded jar/binary) — never system-wide or
  into the user's libraries. For R packages needing GSL, build with miniconda
  removed from PATH. Tell students the conda/Bioconductor install line.

## Data
- Existing course data: `data/manifest.tsv` and `data/README.md` (GM12878
  ENCODE peaks/bigWigs/GTF/genome, prepared chr21+chrM reads for CTCF ChIP,
  input and ATAC, Corces 2016 counts, hg38.fa.bgz). Reuse it where it fits.
- New downloads: put them in `data/` with a clear name, and list each one in
  **your own fragment** `plan/manifest-<block>.tsv` (same 4 columns as the
  manifest: id, url | `derived` | `prepare-<block>.sh`, description,
  used_by units). The maintainer merges fragments into `data/manifest.tsv`.
- If files must be built (subsetting reads, converting formats), write
  `data/prepare-<block>.sh` (new file, yours alone), runnable by a teacher.
- Never check large files into git (`data/*` is ignored, except scripts).

## Rendering and checking
Never render in the repo. Private copy: `rsync -a --exclude docs --exclude data
--exclude .git <repo>/ $SCR/site/ && ln -s <repo>/data $SCR/site/data`, re-sync
your unit folders after edits, `cd $SCR/site && quarto render <unit>/lec.qmd`
(and exercise, quiz). Screenshot every slide:
`timeout 90 google-chrome --headless=new --no-sandbox --user-data-dir=$(mktemp -d -p $SCR)
--disable-gpu --hide-scrollbars --window-size=1280,720 --screenshot=$SCR/s.png
"file://$SCR/site/docs/<unit>/lec.html#/N"`; tile with montage; look; fix.

## Shared files you must not edit
`_quarto.yml`, `contentplan.md`, `plan/key-papers.md`, `data/manifest.tsv`,
`data/download.sh`, `intro/install.qmd`, `CODE_AUDIT_CHECKLIST.md`,
`index.qmd`, `dictionary.qmd`, other agents' units. Report what they need.

## Final report (your last message)
Per unit: files; slide count and number with figures; key numbers from runs;
verified key papers with DOIs; tools/packages students need (exact install
lines); glossary terms (with one-line definitions); data added (fragment
file); runtimes; uncertainties.
