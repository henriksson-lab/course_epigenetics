// Fit an epigenetic clock: an elastic net on simulated methylation data.
// People aged 20–80; P candidate CpGs, of which K change linearly with age, the rest are noise.
// The elastic net (coordinate descent, model.ts) is fitted on N training people along a
// lambda grid; a fixed set of 200 test people is predicted.
// (a) predicted vs true age at the chosen lambda (test filled, training faint);
// (b) median absolute error against lambda, training and test, with the number of CpGs used.
// Pre-selection: none, top 50 CpGs by |r| with age in the training people, or in ALL people
// (leakage: the test people's ages help choose the CpGs).
// Auto-mounts on <div data-demo="clockfit">. Optional data-*:
//   data-n="100" data-p="1000" data-k="20" data-noise="0.05" data-alpha="0.5"
//   data-lambda="1" data-select="none|train|all" data-seed="1"
import { svgEl, htmlEl, slider, button, presets, scale, autoMount, Slider } from "../lib/ui";
import { simulate, buildData, topCorrelated, path, lambdaGrid, predict, median, ClockSim, Data, PathPoint } from "./model";

const NAME = "clockfit";
const COL = { test: "#2a6f97", train: "#5b6b7a", bad: "#c0392b", ok: "#2e8b57", fg: "#1b2a3a", muted: "#5b6b7a", base: "#e67e22" };
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .cf-ctrl { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: .25em 1.6em; width: 100%; }
[data-demo="${NAME}"] .cf-row { display: flex; flex-wrap: wrap; gap: .4em; align-items: center; }
[data-demo="${NAME}"] .cf-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .cf-plots svg { flex: 1 1 340px; min-width: 0; max-width: min(560px, calc(50% - 6px)); height: auto; }
[data-demo="${NAME}"] svg text { font-size: 16px; fill: ${COL.fg}; }
[data-demo="${NAME}"] svg .ttl { font-size: 17px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15px; fill: ${COL.muted}; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 12em 1fr 3.6em; }
.reveal [data-demo="${NAME}"] { font-size: 17px; gap: .3em; }
.reveal [data-demo="${NAME}"] .cf-plots { justify-content: center; }
.reveal [data-demo="${NAME}"] .cf-plots svg { max-width: 380px; }
.reveal [data-demo="${NAME}"] .qb-note { font-size: 16px; }
@media (max-width: 760px) {
  [data-demo="${NAME}"] .cf-plots svg { max-width: 100%; flex-basis: 100%; }
  [data-demo="${NAME}"] .qb-slider { grid-template-columns: auto 1fr 3.6em; }
}
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const NMAX = 400, NTEST = 200, PS = [20, 50, 100, 200, 500, 1000, 2000], KMAX = 40, NSEL = 50;
const LOG_LO = -2, LOG_HI = 2, NGRID = 25;
const GRID = lambdaGrid(10 ** LOG_LO, 10 ** LOG_HI, NGRID);       // largest first
type Sel = "none" | "train" | "all";

const fmtL = (v: number): string => (v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v >= 0.1 ? v.toFixed(2) : v.toFixed(3));

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (v: string | undefined, d: number) => (v !== undefined && Number.isFinite(Number(v)) ? Number(v) : d);
  const nearestP = (p: number) => PS.reduce((b, v, i) => (Math.abs(Math.log(v / p)) < Math.abs(Math.log(PS[b] / p)) ? i : b), 0);
  const st = {
    n: Math.min(NMAX, Math.max(20, Math.round(num(ds.n, 100)))),
    pi: nearestP(num(ds.p, 1000)),
    k: Math.min(KMAX, Math.max(0, Math.round(num(ds.k, 20)))),
    noise: num(ds.noise, 0.05),
    alpha: Math.min(1, Math.max(0, num(ds.alpha, 0.5))),
    li: 0,
    sel: (["none", "train", "all"].includes(ds.select ?? "") ? ds.select : "none") as Sel,
    seed: num(ds.seed, 1),
  };
  const lamIndex = (lam: number) => Math.round((LOG_HI - Math.log10(lam)) / (LOG_HI - LOG_LO) * (NGRID - 1));
  st.li = Math.min(NGRID - 1, Math.max(0, lamIndex(num(ds.lambda, 1))));

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "cf-ctrl", "", host);
  const c1 = htmlEl("div", "", "", ctr), c2 = htmlEl("div", "", "", ctr);
  const sl: Record<string, Slider> = {};
  sl.n = slider(c1, "training people N", 20, NMAX, 10, st.n, v => { st.n = v; refit(); });
  sl.p = slider(c1, "candidate CpGs P", 0, PS.length - 1, 1, st.pi, v => { st.pi = v; refit(); }, v => String(PS[v]));
  sl.k = slider(c1, "of which age CpGs", 0, KMAX, 1, st.k, v => { st.k = v; refit(); });
  sl.noise = slider(c1, "noise SD (β)", 0.01, 0.15, 0.005, st.noise, v => { st.noise = v; refit(); }, v => v.toFixed(3));
  // slider runs from weak (left) to strong (right) penalty
  sl.lam = slider(c2, "penalty λ", 0, NGRID - 1, 1, NGRID - 1 - st.li, v => { st.li = NGRID - 1 - v; draw(); }, v => fmtL(GRID[NGRID - 1 - v]));
  sl.alpha = slider(c2, "α (0 ridge, 1 lasso)", 0, 1, 0.05, st.alpha, v => { st.alpha = v; refit(); });
  const r1 = htmlEl("div", "cf-row", "", c2);
  htmlEl("span", "", `pre-select the ${NSEL} CpGs with the largest |r| with age:`, r1);
  const r2 = htmlEl("div", "cf-row", "", c2);
  const selBtns = presets(r2, [
    { label: "no selection", apply: () => { st.sel = "none"; refit(); } },
    { label: "in training people", apply: () => { st.sel = "train"; refit(); } },
    { label: "in all people (leaky)", apply: () => { st.sel = "all"; refit(); } },
  ]);
  selBtns[["none", "train", "all"].indexOf(st.sel)].classList.add("active");
  button(r2, "new data", () => { st.seed++; sim = simulate(st.seed, NMAX, NTEST, PS[PS.length - 1], KMAX); refit(); });
  const note = htmlEl("p", "qb-note", "", host);

  const plots = htmlEl("div", "cf-plots", "", host);
  const W = 400, H = 410;
  const svgA = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Predicted against true age" }, plots);
  const svgB = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Error against penalty" }, plots);

  let sim: ClockSim = simulate(st.seed, NMAX, NTEST, PS[PS.length - 1], KMAX);
  let data: Data | null = null, cols: number[] = [], pa: PathPoint[] = [], baseline = 0;

  function compute(): void {
    const P = PS[st.pi];
    const d = buildData(sim, st.n, P, Math.min(st.k, P), st.noise);
    const all = Array.from({ length: P }, (_, j) => j);
    cols = st.sel === "none" || P <= NSEL ? all
      : topCorrelated(d.X, d.y, st.sel === "all" ? [...d.train, ...d.test] : d.train, NSEL);
    pa = path(d, cols, st.alpha, GRID);
    const med = median(d.train.map(i => d.y[i]));
    baseline = median(d.test.map(i => Math.abs(d.y[i] - med)));
    data = d;
  }

  let pending = false;
  function refit(): void {
    if (pending) return;
    pending = true;
    note.textContent = "fitting …";
    setTimeout(() => requestAnimationFrame(() => { pending = false; compute(); draw(); }), 0);
  }

  function draw(): void {
    if (!data) return;
    const P = PS[st.pi], k = Math.min(st.k, P);
    const pt = pa[st.li];
    const selTxt = st.sel === "none" || P <= NSEL ? `all ${P} CpGs` : `the top ${NSEL} of ${P} CpGs, chosen in ${st.sel === "all" ? "all people (training + test)" : "the training people"}`;
    note.textContent = `Simulated: ${st.n} training and ${NTEST} test people aged 20–80; ${k} of ${P} CpGs change with age, the rest are noise. ` +
      `Elastic net on ${selTxt}.`;
    drawA(pt); drawB(pt, k);
  }

  // ---------------- (a) predicted vs true ----------------
  function drawA(pt: PathPoint): void {
    const svg = svgA, d = data!;
    svg.replaceChildren();
    const C = cols.map(j => d.X[j]);
    const ptr = predict(pt.fit, C, d.train), pte = predict(pt.fit, C, d.test);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent = `(a) Predicted age at λ = ${fmtL(pt.lambda)}`;
    const t1 = svgEl("text", { x: 8, y: 46 }, svg);
    const s1 = svgEl("tspan", {}, t1); s1.textContent = `test error ${pt.testMae.toFixed(1)} y`; s1.style.fill = COL.test; s1.style.fontWeight = "700";
    const s2 = svgEl("tspan", {}, t1); s2.textContent = `  training ${pt.trainMae.toFixed(1)} y`; s2.style.fill = COL.train;
    svgEl("text", { x: 8, y: 68, class: "sm" }, svg).textContent =
      `median absolute error · test r = ${pt.testR.toFixed(2)}`;
    const L = 56, R = W - 18, T = 86, B = H - 50;
    const lo = 0, hi = 100;
    const sx = scale(lo, hi, L, R), sy = scale(lo, hi, B, T);
    for (const v of [0, 20, 40, 60, 80, 100]) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: T, y2: B, class: "grid" }, svg);
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: sx(v), y: B + 20, "text-anchor": "middle", class: "sm" }, svg).textContent = String(v);
      svgEl("text", { x: L - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    svgEl("line", { x1: sx(lo), y1: sy(lo), x2: sx(hi), y2: sy(hi), stroke: COL.fg, "stroke-width": 1.5, "stroke-dasharray": "6 4" }, svg);
    const clip = (v: number) => Math.max(lo - 2, Math.min(hi + 2, v));
    const gtr = svgEl("g", { fill: "none", stroke: COL.train, "stroke-opacity": 0.45, "stroke-width": 1.3 }, svg);
    d.train.forEach((i, ii) => svgEl("circle", { cx: sx(d.y[i]).toFixed(1), cy: sy(clip(ptr[ii])).toFixed(1), r: 3.5 }, gtr));
    const gte = svgEl("g", { fill: COL.test, "fill-opacity": 0.75 }, svg);
    d.test.forEach((i, ii) => svgEl("circle", { cx: sx(d.y[i]).toFixed(1), cy: sy(clip(pte[ii])).toFixed(1), r: 3.8 }, gte));
    svgEl("text", { x: (L + R) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = "true age";
    const yt = svgEl("text", { x: 14, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${(T + B) / 2})` }, svg);
    yt.textContent = "predicted age";
    const lg = svgEl("text", { x: R - 4, y: B - 10, "text-anchor": "end", class: "sm" }, svg);
    const l1 = svgEl("tspan", {}, lg); l1.textContent = "● test  "; l1.style.fill = COL.test;
    const l2 = svgEl("tspan", {}, lg); l2.textContent = "○ training"; l2.style.fill = COL.train;
  }

  // ---------------- (b) error vs lambda ----------------
  function drawB(pt: PathPoint, k: number): void {
    const svg = svgB;
    svg.replaceChildren();
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent = "(b) Error against the penalty";
    const t1 = svgEl("text", { x: 8, y: 46 }, svg);
    t1.textContent = `CpGs used: ${pt.nonzero} (${pt.nAge} of ${k} age, ${pt.nonzero - pt.nAge} noise)`;
    const best = pa.reduce((a, b) => (b.testMae < a.testMae ? b : a));
    const bl = svgEl("text", { x: 8, y: 68, class: "sm" }, svg);
    bl.textContent = `best test error ${best.testMae.toFixed(1)} y · median guess ${baseline.toFixed(1)} y`;
    const L = 52, R = W - 20, T = 122, B = H - 50;
    const ymax = Math.max(20, baseline * 1.15, ...pa.map(q => q.testMae));
    const top = Math.ceil(ymax / 5) * 5;
    const sx = scale(LOG_LO, LOG_HI, L, R), sy = scale(0, top, B, T);
    for (let v = 0; v <= top; v += 5) {
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    for (let e = LOG_LO; e <= LOG_HI; e++) {
      svgEl("line", { x1: sx(e), x2: sx(e), y1: T, y2: B, class: "grid" }, svg);
      svgEl("text", { x: sx(e), y: B + 20, "text-anchor": "middle", class: "sm" }, svg).textContent = fmtL(10 ** e);
    }
    // number of CpGs used, along the top
    pa.forEach((q, i) => {
      if (i % 4 !== 0) return;
      svgEl("text", { x: sx(Math.log10(q.lambda)), y: T - 8, "text-anchor": "middle", class: "sm" }, svg).textContent = String(q.nonzero);
    });
    svgEl("text", { x: 8, y: T - 30, class: "sm" }, svg).textContent = "CpGs used at each λ:";
    svgEl("line", { x1: L, x2: R, y1: sy(baseline), y2: sy(baseline), stroke: COL.base, "stroke-width": 2, "stroke-dasharray": "7 5" }, svg);
    const line = (key: "trainMae" | "testMae", col: string, w: number) => {
      const dstr = pa.map((q, i) => `${i ? "L" : "M"}${sx(Math.log10(q.lambda)).toFixed(1)},${sy(Math.min(top, q[key])).toFixed(1)}`).join("");
      svgEl("path", { d: dstr, fill: "none", stroke: col, "stroke-width": w }, svg);
    };
    line("trainMae", COL.train, 2.5);
    line("testMae", COL.test, 3);
    const xc = sx(Math.log10(pt.lambda));
    svgEl("line", { x1: xc, x2: xc, y1: T, y2: B, stroke: COL.fg, "stroke-width": 1.5 }, svg);
    svgEl("circle", { cx: xc, cy: sy(Math.min(top, pt.testMae)), r: 6, fill: COL.test }, svg);
    svgEl("circle", { cx: xc, cy: sy(Math.min(top, pt.trainMae)), r: 6, fill: COL.train }, svg);
    const lg = svgEl("text", { x: L + 6, y: T + 20, class: "sm" }, svg);
    const a1 = svgEl("tspan", {}, lg); a1.textContent = "— test  "; a1.style.fill = COL.test; a1.style.fontWeight = "700";
    const a2 = svgEl("tspan", {}, lg); a2.textContent = "— training  "; a2.style.fill = COL.train; a2.style.fontWeight = "700";
    const a3 = svgEl("tspan", {}, lg); a3.textContent = "- - median guess"; a3.style.fill = COL.base; a3.style.fontWeight = "700";
    svgEl("text", { x: (L + R) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = "λ (log scale): weak ← penalty → strong";
    const yt = svgEl("text", { x: 14, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${(T + B) / 2})` }, svg);
    yt.textContent = "median abs. error (y)";
  }

  compute();
  draw();
}

autoMount(NAME, mount);
