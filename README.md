# Epigenetics (MSc course, Umeå University): course site sources

Teacher information. This file is not part of the published site.

The scaffold (build system, theme, slide editor, demo framework, Canvas quiz
builder) is reused from the quantbio-3MB00A course.

## Where things are

- `contentplan.md`: the list of teaching units and their status. Start here.
- `designprinciples.md`: rules for writing material (read before adding a unit).
- `audience.md`: who the students are and what they already know.
- `elo.md`: expected learning outcomes (draft).
- `<slug>/`: one folder per teaching unit: `lec.qmd` (reveal.js slides),
  `exercise.qmd`, optional `quiz.qmd`, figures `figNN-*.svg`, and optional
  `make-*-figs.R` scripts for R-made figures.
- `_quarto.yml`: site config; the sidebar defines the teaching order.
- `plan/`: internal working documents (not published).
- `demos/`: interactive TypeScript demos (see `demos/README.md`).
- `tools/slide-editor/`: edit slide text in the browser (`make edit-slides`).
- `canvas-quiz/`: Canvas bonus-quiz builder; answers live in the git-ignored
  `canvas-quiz/private/`.

## Building this site

```bash
make           # render the whole site to docs/, then git add
make preview   # start a live-reload preview server
make figs      # re-run every <unit>/make-*.R figure script
make quiz      # build Canvas quizzes from canvas-quiz/private/questions-<slug>.yml
make clean     # delete docs/ and the Quarto cache
```

You will need [Quarto](https://quarto.org/docs/get-started/), a recent
[R](https://cran.r-project.org/), and for the demos `npm i` once.
