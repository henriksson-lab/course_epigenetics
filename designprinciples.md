# Design principles

The principles below govern how we build course material in this repository.
Check every new lecture, exercise and quiz against them.

See [audience.md](audience.md) for who we design for, [elo.md](elo.md) for what
we teach, and [contentplan.md](contentplan.md) for the topic list.

Adapted from the quantbio-3MB00A course, whose scaffold this repo reuses.

## 1. Biology ↔ assay ↔ analysis

Every assay is taught as one chain: **the biological question → what the
assay physically measures (and what the chemistry gets wrong) → processing on
the command line → analysis and interpretation in R**. We never teach a tool in
the abstract, and never teach an assay without running its data.

- The unit of design is the **topic**, not the day.
- Most assays get two units: one on the biology and the assay, one on the
  analysis. Split further if a unit grows too large (principle 7).

## 2. Show the reason for each step, and the code that does it

Every processing step is tied to the artefact or question that motivates it
(why remove duplicates? why a blacklist? why an input control?), and every
claim about the data is shown with code a student can run.

- First encounter: transparent tools and small steps (`samtools`, `bedtools`,
  `GenomicRanges`, a few lines of R). Then point to the production pipeline
  (nf-core, ArchR, …) that does the same thing at scale.
- Build intuition with pictures (genome-browser views, a toy fragment-length
  distribution) before formulas.

## 3. Reuse what they know from RNA-seq

Students have done RNA-seq in R. Say explicitly when a step is the same
(DESeq2 on a count matrix of peaks) and spend the teaching on what differs
(peak calling, controls, normalisation under global changes, interval logic).

## 4. Exercises over stretch goals

Prefer a good number of core exercises that every student completes over a few
exercises with long optional tails. A short neutral **"Where this leads"** note
at the end is fine.

## 5. Linked hints over named hints

Every tool, function, file format or concept named in a hint links to its
documentation (manual, vignette, format spec, or canonical paper), so a student
never has to guess search terms to finish a core exercise.

## 6. Pen and paper, because the exam is written

The only examination is an **individual written exam**. Computer exercises
build understanding, but what is examined must be doable without a computer.
Every unit therefore includes pen-and-paper questions (in its quiz and in a
"By hand" section of its exercise) of three kinds:

- **Read the figure:** a genome-browser view, a fingerprint / TSS-enrichment /
  fragment-size plot, a Hi-C contact map, a methylation profile, an M-bias plot,
  a volcano or MA plot of peaks, a ROC/PR curve.
- **Small calculations:** a Poisson peak p-value from a local lambda, FRiP,
  a methylation level from C/T counts, a Tn5 shift, an observed/expected Hi-C
  value, precision and recall from a confusion matrix.
- **Judge the design or the claim:** missing controls or replicates,
  normalisation under a global change, leakage between training and test
  data, causal vs correlational claims.

Never examine tool flags or syntax from memory; examine what a step does and
why.

## 7. Slide decks: one topic, up to ~30 slides

Each lecture deck is recorded as a **video**, so a deck is a self-contained
snippet on one topic.

- **Soft limit: about 30 slides** (content slides; title slide not counted).
  Aim lower when the topic allows.
- **At least half of the slides carry an explanatory figure** (more is fine):
  a diagram that explains a mechanism or a step, or a plot of real data. Not
  decoration: the figure should carry the point of the slide. Follow
  `plan/figure-qa-checklist.md`.
- When a deck would grow past that, first ask whether the new material is
  **its own topic** and should become a separate unit (its own folder with
  `lec.qmd` and `exercise.qmd`, placed next to the original in the sidebar).
  Prefer splitting over cutting good material.

## 8. Names are topics, never days or weeks

- One folder per unit, named by topic in kebab-case (`chipseq-peaks`,
  `methylation-arrays`). **No day, week or sequence numbers** in folder or
  file names, and no "Day N" / "Week N" in slide subtitles or cross-references:
  refer to other units by their title.
- The teaching order lives **only** in `_quarto.yml` (sidebar) and `index.qmd`.
  Reordering means editing those, never renaming a folder.
- Quizzes for Canvas are likewise named by topic (`canvas-quiz/private/questions-<slug>.yml`).

## 9. Key papers as optional further reading

This is a final-year course: students should know the primary literature
behind each method. Every unit ends with a **"Key papers"** list, on the last
slide of `lec.qmd` and as a section at the end of `exercise.qmd`.

- 3–6 papers per unit: the paper that introduced the assay or method, a
  landmark biological finding made with it, and where useful one good review
  or a critical/benchmark paper.
- One line per paper saying **why** it is worth reading.
- Clearly marked **optional**: not required reading, not examined.
- Every citation is **verified** (authors, year, journal, DOI) and linked by
  DOI. Never cite from memory without checking. Seed list: `plan/key-papers.md`.

## 10. Honest about simplifications

When an exercise uses a shortcut (one chromosome, a pre-aligned BAM, a
downsampled Hi-C matrix, a toy peak caller), say so plainly and name the unit
where the real thing is shown.
