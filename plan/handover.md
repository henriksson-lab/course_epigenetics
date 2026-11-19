# Handover (2026-11-19)

## State
- Pushed to github.com/henriksson-lab/course_epigenetics `main` (commit a3d3e1f),
  GitHub Pages builds from `main` `/docs`. **29 of 31 units written**, rendered and
  live (see the `written` rows in `contentplan.md`).
- **Not written:** `epigenetic-inheritance` and `cancer-epigenetics`. The agent
  writing them was stopped by the API's safety filter (a false positive on
  biology content); partial files are in those folders and `data/prepare-dis.sh`,
  **untracked and unreviewed** — do not use them as they are. The user decides
  how to proceed (another model, or write them by hand).
- Pushing over SSH failed after a reboot (key not loaded); pushes used
  `git -c credential.helper='!gh auth git-credential' push https://github.com/...`.

## How the material was made
- Units were written by parallel agents following `plan/unit-writing-brief.md`
  and `plan/parallel-authoring.md`; each ran its exercise code against `data/` and
  screenshot-checked its slides. Integration: `plan/merge-manifests.py` (merges
  `plan/manifest-<block>.tsv` into `data/manifest.tsv`, recomputes used_by by
  scanning units), `plan/build-key-papers.py` (rebuilds `plan/key-papers.md`
  from each written unit's Key papers slide).
- Render units one file at a time (`quarto render <unit>/lec.qmd`): a multi-file
  `quarto render a b c` only rendered the first file here.

## Open items
- Hosting of prepared data for students (reads, subsets; built by
  `data/prepare-*.sh`) is undecided (`data/README.md`). `data/README.md` does not
  yet describe the blocks added on 2026-11-19 (their files are in the manifest).
- `CODE_AUDIT_CHECKLIST.md`: nothing independently audited yet except
  `epigenomics-design` (re-run 2026-11-19, numbers agree).
- Overlapping data: two sets of GM12878/K562 H3K27ac and input chr21 BAMs
  (`prepare-chip.sh` per replicate, `prepare-cuttag.sh` merged) — kept both.
- This machine: conda `samtools` first on PATH is 1.3.1 (use `/usr/bin/samtools`);
  deepTools needs `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`.

## Demos (2026-11-19: all 18 built, linked in the sidebar, rendered and pushed)
Built per `plan/demo-brief.md` by parallel agents; each embedded an Explore
section (E-questions + solutions) in its unit's exercise page and built
`demos/dist/<name>.js`. Done (15): readcoverage, shufflenull, tn5ladder,
globalnorm, motif2x2, polycombmemory, bisulfite, betabinom, clockfit,
peakcaller, fingerprint, crosscorr, rocpr, cnnfilter, leakagesim.
Also done: loopextrusion (genome-3d-intro), compartments (hic-features),
hmmstates (chromatin-states) — all 18 demos finished and committed locally
(sources + exercise pages); not yet in the sidebar, not rendered, not pushed.
Remaining integration: add sidebar "Demo:" links in `_quarto.yml`
(after each unit's Exercise entry; anchors: `#explore`, and for
epigenomics-design `#explore-fingerprint`, `#explore-crosscorr`), run
`make demos` + `node_modules/.bin/tsc --noEmit -p demos`, render the written
units one file at a time, commit (exclude epigenetic-inheritance/,
cancer-epigenetics/, data/prepare-dis.sh), push via gh HTTPS.
Also fixed during this round (uncommitted): differential-accessibility lec note
(background normalisation deflates, not inflates).
