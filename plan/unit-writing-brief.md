# Brief: writing a teaching unit

How to write one unit of this course (a folder `<slug>/` with `lec.qmd`,
`exercise.qmd`, `quiz.qmd`, figures and figure scripts). Written for whoever
writes the unit, human or LLM.

## Read first
- `designprinciples.md` (the rules), `audience.md` (who the students are),
  `contentplan.md` (the unit's row: its scope).
- `plan/figure-style.md` and `plan/figure-qa-checklist.md` (figures).
- The finished units `chromatin-intro/` and `genomic-intervals/`: copy their
  **format exactly** (YAML headers, speaker notes, exercise sections, callouts,
  quiz layout, key-papers slide and section).

## The lecture (`lec.qmd`)
- reveal.js header as in the finished units; `subtitle` describes the topic,
  never a day or week.
- About **20–28 content slides**, soft maximum 30. **At least half carry an
  explanatory figure** (more is fine). Every slide has `::: notes` speaker notes
  (the lecture is recorded as a video): what to say, numbers, caveats.
- First slide "About this lecture"; then the content; then "Summary" (table)
  and "Key papers (optional further reading)".
- Figures: hand-drawn SVGs (`figNN-name.svg`, style guide, QA'd by rendering
  with `cairosvg file.svg -o x.png --output-width 1100` and looking at the
  PNG), and real-data PNGs made by `make-<unit>-figs.R` (run from the unit
  folder, reading `../data/`). Reproducible, seeded.
- Use `{.smaller}` on dense slides, but never on slides that are mostly code.

## The exercise (`exercise.qmd`)
Sections in this order: Goal · Background · Data · By hand · On the command
line (when relevant) · In R · Hints (collapsed callouts) · Solutions (collapsed
callouts: by hand / command line / R) · Where this leads · Key papers.
- **By hand** matters most: the exam is written. Small calculations, reading a
  figure, judging a design.
- Every command and R chunk is **run exactly as written** against `data/`;
  numbers in the solutions come from those runs. Commands are plain shell that
  works on a laptop or on an HPC login node / SLURM job.
- Data paths: `D=/path/to/course/data` (shell) and `D <- "/path/to/course/data"` (R).
- Hints link every tool/function to its documentation.

## The quiz (`quiz.qmd`)
8–10 short questions with collapsed hint/answer callouts, answerable on paper.

## Key papers
3–6 per unit, on the last slide and at the end of the exercise, marked optional
and not examined. Verify every citation (authors, year, journal, DOI), e.g.
`curl -sL -H "Accept: text/x-bibliography; style=apa" https://doi.org/<doi>`.
Never cite from memory without checking.

## Cross-references and names
Refer to other units by title (*Genomic intervals*), never by day/week or slug.

## Do not edit shared files
`_quarto.yml`, `contentplan.md`, `plan/key-papers.md`, `data/manifest.tsv`,
`intro/install.qmd`, `CODE_AUDIT_CHECKLIST.md`, `index.qmd`, `dictionary.qmd`,
and other units are integrated by the course maintainer. Report what they need
(packages, tools, glossary terms, verified papers, data files).
