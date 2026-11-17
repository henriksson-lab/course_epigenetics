# Figure style guide (SVG lecture diagrams)

All hand-drawn lecture figures follow this guide so the course looks like one
system. QA every figure with `plan/figure-qa-checklist.md`.

## Canvas
- Plain SVG, `xmlns="http://www.w3.org/2000/svg"`, a `viewBox` **1120 wide**
  (height to fit the content, typically 300–560), and `width`/`height`
  attributes equal to the viewBox. White background rect first.
- Slides are 1280×720 and figures are shown at 80–95 % width, often with
  bullets below: keep figures **wide and short** (aspect ≥ 2:1 unless the
  slide puts the figure in a column, then ≈ 4:3 and viewBox 700–800 wide).
- `font-family="Helvetica, Arial, sans-serif"`. Body labels **≥ 15 px**,
  panel titles 17–18 px bold, small annotations ≥ 13 px. Never below 13.
- No external resources, no embedded raster images, no scripts, no CSS
  classes needed (inline attributes are fine).

## Colours
Text and outlines `#1b2a3a`; secondary text `#5b6b7a`; light fills
`#eaf1f7`, `#f4f6f8`; DNA strand `#1b2a3a` (or `#5b6b7a` when in the background).

Histones (used whenever octamer components are shown):
H2A `#f4a261`, H2B `#e9c46a`, H3 `#2a9d8f`, H4 `#8ab17d`, H1 `#9c6644`.

Chromatin marks / tracks (same in R-made figures):

| mark | colour |
|---|---|
| H3K4me3 | `#c0392b` |
| H3K27ac | `#e67e22` |
| H3K4me1 | `#d4a017` |
| H3K36me3 | `#2e8b57` |
| H3K27me3 | `#2a6f97` |
| H3K9me3 | `#6b4fa0` |
| CTCF | `#8c564b` |

Activity: active/open `#2e8b57` (green), repressed/closed `#a8324a` (crimson).
Highlights/warnings: `#c0392b`. Interval sets A/B in interval diagrams:
A `#2a6f97` (blue), B `#e67e22` (orange), result `#2e8b57` (green).

DNA methylation: methylated CpG = **filled black circle** "lollipop",
unmethylated = **open circle** (the field's standard).

## Conventions
- Genes: exons as boxes, introns as thin lines, an arrow at the TSS bent in the
  direction of transcription. + strand genes point right.
- Genome coordinates increase left to right. Label chromosome positions
  with a ruler when positions matter.
- Arrowheads via a `<marker>`; stroke 2–2.5 px.
- Prefer labels next to objects over legends. Each panel has a short title.
- Do not repeat the slide's bullet text inside the figure; the figure shows
  the mechanism, the bullets say it.
