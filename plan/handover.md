# Handover (2026-11-18, before reboot)

## State
- Committed and pushed (origin = github.com/henriksson-lab/course_epigenetics, branch `main`,
  GitHub Pages on from `main` `/docs`): scaffold + `chromatin-intro` + `genomic-intervals`
  (commit 057d976).
- **Not yet committed:** seven new units and integration edits (see `git status`):
  `gene-regulation`, `enrichment-alignment`, `epigenomics-design`, `atacseq-intro`,
  `atacseq-analysis`, `differential-accessibility`, `motifs-footprinting`; plus
  `data/prepare-reads.sh`, `data/filter_fastq_by_names.py`, new `data/download.sh`
  (handles `prepare-reads.sh` and `derived` sources), `data/manifest.tsv` (used_by
  lists + prepared/derived rows), `_quarto.yml` sidebar (Foundations + ATAC block),
  `intro/install.qmd`, `index.qmd`, `.gitignore`, `plan/unit-writing-brief.md`.
- The units were written by parallel agents; each ran its code against `data/` and
  checked its slides by screenshot. Six reported complete. **`epigenomics-design`
  was still in its final checks at reboot**: files are all present (lec, exercise,
  quiz, 24 figs, make-design-figs.R) but it has not reported. Review it first:
  render, screenshot the slides, check that the numbers in the solutions match a run,
  and that the referenced figures exist.

## Remaining integration TODO (then `make`, commit, push — user asked to push when done)
1. Review `epigenomics-design` (above). Check its cross-references by title.
2. `contentplan.md`: mark the seven units `written`.
3. `plan/key-papers.md`: tick and add the verified DOIs reported by the agents
   (listed in each unit's Key papers section; all verified via doi.org).
4. `CODE_AUDIT_CHECKLIST.md`: add the seven units' lec/exercise/make-*-figs.R entries.
5. `dictionary.qmd`: add glossary terms (agents suggested lists; at least the core:
   PWM, motif, pioneer factor, enhancer, MAPQ, SAM FLAG, duplicate, blacklist, CPM/RPGC,
   FRiP, fingerprint, NSC/RSC, NRF/PBC, IDR, Tn5 shift, TSS enrichment, local lambda,
   consensus peaks, size factor/TMM, chromVAR, footprint).
6. `make render`, then screenshot-check a few slides of each new deck, commit
   (no LLM attribution), push, and confirm the Pages build succeeds.

## Facts worth knowing
- Datasets: ENCODE GM12878 (GRCh38) for processing units; reads are chr21+chrM subsets
  built by `data/prepare-reads.sh` (fixed: strips `/1` from BAM read names).
  CTCF ChIP/input FASTQs are **Phred+64** (the unit converts them with perl).
  Corces et al. 2016 (GSE74912 ATAC counts, GSE74246 RNA) for differential and motif
  units: B cells vs CD4 T cells, donors 1022/4983/5483, `~ donor + celltype`,
  technical reps summed; lifted hg19→hg38.
- `differential-accessibility` uses the apeglm-shrunk LFC (19 404 B-up / 16 897 T-up);
  `motifs-footprinting` uses the MLE LFC (19 439 / 16 972) and says so in its R1 solution.
- Host for the prepared reads (~1 GB) for students: **undecided** (data/README.md).
- This machine: the conda `samtools` first on PATH is 1.3.1 (too old for markdup); use
  `/usr/bin/samtools`. System deepTools needs `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`.
- Tools installed only in the session scratchpad (may be gone after reboot): MACS3
  venv and an R library with JASPAR2024/TFBSTools/motifmatchr/chromVAR/universalmotif.
  Reinstall if needed (system GSL: build with conda removed from PATH).
- Next units in `contentplan.md` after this: the ChIP-seq block.
