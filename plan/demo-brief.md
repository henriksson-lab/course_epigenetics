# Brief: interactive demos

Small interactive TypeScript demos that build intuition, in the style of the
quantbio-3MB00A course (`../quantbio-3MB00A/demos/`, e.g. `countratio`,
`pca2d`, `shuffle`-like demos; read `demos/README.md` and 2–3 demos there,
and their "Explore" sections in the quantbio exercise pages, before writing).

## Technical rules
- One demo per folder: `demos/src/<name>/main.ts` (+ helper `.ts` files),
  bundled by `demos/build.sh <name>` to `demos/dist/<name>.js` (IIFE).
  Build only your own demos (`demos/build.sh name1 name2`); `npm i` is done.
- Use the shared helpers in `demos/src/lib/ui.ts` (svgEl, htmlEl, slider,
  button, presets, checkbox, scale, rng, autoMount, draggable). You may add
  helpers to your own demo folder; do **not** edit `lib/ui.ts` (shared).
- Auto-mount on `<div data-demo="<name>">`; options via `data-*` attributes.
- No network, no CDN libraries, no fetch: plain DOM + SVG (or canvas).
- Fixed random seeds (`rng(seed)`), so a recorded video is reproducible.
- Prefix CSS with `[data-demo="<name>"]`; inject styles once (see countratio).
  Use the course colours from `plan/figure-style.md` (mark colours etc.).
- Must work in a slide (`.reveal` sizes) and on an exercise page, at 1280 px
  and on a narrow screen. Text ≥ 15 px.
- Type-check: `node_modules/.bin/tsc --noEmit -p demos` must pass for your files.
- Models must be honest: label simulated data as simulated; the numbers a demo
  shows must follow from a correct model of the assay or statistic.

## Embedding (in the unit's exercise page)
Add a section before "## In R" (or before "## Hints" if there is no R part):

````md
## Explore: <title> {#explore}

<one paragraph: what the demo shows and what each control does>

```{=html}
<link rel="stylesheet" href="../demos/dist/demo.css">
<div data-demo="<name>"></div>
<script src="../demos/dist/<name>.js"></script>
```

**E1.** … **E2.** … **E3.** …   (3–4 guided questions that make the student
find the intuition by moving the controls; concrete, checkable)

::: {.callout-note collapse="true" title="Solutions"}
… answers with the numbers the demo actually shows (check them) …
:::
````

If a unit gets two demos use ids `#explore` and `#explore-<short>`. Optionally
add a "Try it" slide in the lecture that embeds the same demo (only if it fits
without overflow; check by screenshot).

## Checking
- Render in a private copy (never in the repo): `rsync -a --exclude docs
  --exclude data --exclude .git <repo>/ $SCR/site/ && ln -s <repo>/data
  $SCR/site/data`; copy your built `demos/dist/<name>.js` too; then
  `cd $SCR/site && quarto render <unit>/exercise.qmd`.
- Screenshot the rendered page section with headless Chrome
  (`timeout 90 google-chrome --headless=new --no-sandbox --user-data-dir=$(mktemp -d -p $SCR)
  --disable-gpu --hide-scrollbars --window-size=1280,1400 --screenshot=… "file://…/exercise.html#explore"`),
  and test interaction logic with a small node/jsdom-free unit test of the model
  code (pure functions in a separate `.ts` file are easiest to test with
  `node_modules/.bin/esbuild --bundle --platform=node` + node).
- Never kill processes by pattern (no `pkill -f chrome`): other agents run in parallel.

## Do not edit shared files
`_quarto.yml` (sidebar), `demos/src/lib/*`, `demos/build.sh`, other units.
Report the sidebar link text you want, e.g. `Demo: Poisson peak caller` →
`chipseq-peaks/exercise.qmd#explore`.

## Final report
Per demo: name, unit and anchor, what it shows, controls, the checks you ran,
the solution numbers, sidebar link text, anything uncertain.
