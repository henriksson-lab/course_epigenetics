// Differential methylation with replicates: pooled-read Fisher test vs tests that use the
// replicates. Two groups of replicates, many CpGs; each replicate has its own true level
// (beta-binomial, overdispersion phi). A chosen fraction of CpGs differs by delta beta.
// (a) histogram of p-values for one test, null CpGs and truly differential CpGs stacked;
// (b) calls at BH-FDR 5% for the three tests, true and false.
// Auto-mounts on <div data-demo="betabinom">. Optional data-*:
//   data-cpgs="5000" data-coverage="20" data-reps="3" data-phi="0.1" data-delta="0.2"
//   data-frac="0" (fraction differential, 0–0.3) data-seed="1"
import { svgEl, htmlEl, slider, button, presets, scale, autoMount, Slider } from "../lib/ui";
import { simulate, testFisher, testBB, testT, estimatePhi, summarize, SimParams, Summary } from "./model";

const NAME = "betabinom";
const COL = { nul: "#2a6f97", diff: "#e67e22", ok: "#2e8b57", bad: "#c0392b", fg: "#1b2a3a", muted: "#5b6b7a" };
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .bb-ctrl { display: grid; grid-template-columns: repeat(auto-fit, minmax(290px, 1fr)); gap: .25em 1.6em; width: 100%; }
[data-demo="${NAME}"] .bb-row { display: flex; flex-wrap: wrap; gap: .4em; align-items: center; }
[data-demo="${NAME}"] .bb-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .bb-plots svg { flex: 1 1 340px; min-width: 0; max-width: min(560px, calc(50% - 6px)); height: auto; }
[data-demo="${NAME}"] svg text { font-size: 16px; fill: ${COL.fg}; }
[data-demo="${NAME}"] svg .ttl { font-size: 17px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15px; fill: ${COL.muted}; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 11.5em 1fr 3.4em; }
.reveal [data-demo="${NAME}"] { font-size: 17px; gap: .3em; }
.reveal [data-demo="${NAME}"] .bb-plots { justify-content: center; }
.reveal [data-demo="${NAME}"] .bb-plots svg { max-width: 380px; }
.reveal [data-demo="${NAME}"] .qb-note { font-size: 16px; }
@media (max-width: 760px) {
  [data-demo="${NAME}"] .bb-plots svg { max-width: 100%; flex-basis: 100%; }
  [data-demo="${NAME}"] .qb-slider { grid-template-columns: auto 1fr 3.4em; }
}
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

type Method = "fisher" | "bb" | "t";
const METHODS: { key: Method; label: string; short: string; sub: string }[] = [
  { key: "fisher", label: "Fisher, pooled reads", short: "Fisher", sub: "pooled reads" },
  { key: "bb", label: "beta-binomial", short: "beta-binom.", sub: "replicates" },
  { key: "t", label: "t-test on β", short: "t-test", sub: "replicates" },
];

const fmt = (v: number): string => (v >= 1000 ? Math.round(v).toLocaleString("en-US").replace(/,/g, " ") : String(v));

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (v: string | undefined, d: number) => (v !== undefined && Number.isFinite(Number(v)) ? Number(v) : d);
  const P: SimParams = {
    cpgs: Math.max(500, Math.round(num(ds.cpgs, 5000))),
    coverage: num(ds.coverage, 20), reps: Math.round(num(ds.reps, 3)), phi: num(ds.phi, 0.1),
    delta: num(ds.delta, 0.2), frac: num(ds.frac, 0), seed: num(ds.seed, 1),
  };
  let method: Method = "fisher";

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "bb-ctrl", "", host);
  const c1 = htmlEl("div", "", "", ctr), c2 = htmlEl("div", "", "", ctr);
  const sl: Record<string, Slider> = {};
  sl.cov = slider(c1, "reads per CpG", 5, 100, 5, P.coverage, v => { P.coverage = v; schedule(); });
  sl.reps = slider(c1, "replicates per group", 1, 8, 1, P.reps, v => { P.reps = v; schedule(); });
  sl.phi = slider(c1, "overdispersion φ", 0, 0.2, 0.005, P.phi, v => { P.phi = v; schedule(); }, v => v.toFixed(3));
  sl.frac = slider(c2, "differential CpGs", 0, 0.3, 0.01, P.frac, v => { P.frac = v; schedule(); }, v => `${Math.round(100 * v)}%`);
  sl.delta = slider(c2, "their true Δβ", 0, 0.5, 0.05, P.delta, v => { P.delta = v; schedule(); });
  const r1 = htmlEl("div", "bb-row", "", c2);
  presets(r1, [
    { label: "nothing differs", apply: () => { P.frac = 0; sl.frac.set(0); schedule(); } },
    { label: "10% differ", apply: () => { P.frac = 0.1; sl.frac.set(0.1); if (P.delta === 0) { P.delta = 0.2; sl.delta.set(0.2); } schedule(); } },
  ]);
  button(r1, "new data", () => { P.seed++; schedule(); });
  const r2 = htmlEl("div", "bb-row", "", host);
  htmlEl("span", "", "histogram (a) shows:", r2);
  const mBtns = presets(r2, METHODS.map(m => ({ label: m.label, apply: () => { method = m.key; draw(); } })));
  mBtns[0].classList.add("active");
  const note = htmlEl("p", "qb-note", "", host);

  const plots = htmlEl("div", "bb-plots", "", host);
  const W = 400, H = 380;
  const svgA = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Histogram of p-values" }, plots);
  const svgB = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Calls at FDR 5% per test" }, plots);

  let res: { p: Record<Method, Float64Array>; s: Record<Method, Summary | null>; isDiff: Uint8Array; phiHat: number } | null = null;

  function compute(): void {
    const d = simulate(P);
    const phiHat = estimatePhi(d);
    const p = { fisher: testFisher(d), bb: testBB(d, phiHat), t: testT(d) };
    const ok = P.reps >= 2;
    res = { p, isDiff: d.isDiff, phiHat,
      s: { fisher: summarize(p.fisher, d.isDiff), bb: ok ? summarize(p.bb, d.isDiff) : null, t: ok ? summarize(p.t, d.isDiff) : null } };
  }

  let pending = false;
  function schedule(): void {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => { pending = false; compute(); draw(); });
  }

  function draw(): void {
    if (!res) return;
    const nd = res.s.fisher!.nDiff, nn = res.s.fisher!.nNull;
    note.textContent = `Simulated: ${fmt(P.cpgs)} CpGs, ${P.reps} vs ${P.reps} samples, about ${P.coverage} reads per CpG and sample. ` +
      `${nd ? `${fmt(nd)} CpGs truly differ (Δβ = ${P.delta.toFixed(2)}), ${fmt(nn)} do not.` : "No CpG differs between the groups."} ` +
      (P.reps >= 2 ? `φ estimated from the replicates: ${res.phiHat.toFixed(3)}.` : "With one sample per group, φ cannot be estimated.");
    drawA(); drawB();
  }

  // ---------------- (a) p-value histogram ----------------
  function drawA(): void {
    const svg = svgA;
    svg.replaceChildren();
    const m = METHODS.find(x => x.key === method)!;
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent = `(a) p-values: ${m.label}`;
    const p = res!.p[method], isDiff = res!.isDiff;
    if (Number.isNaN(p[0])) {
      svgEl("text", { x: W / 2, y: H / 2, "text-anchor": "middle" }, svg).textContent = "needs ≥ 2 replicates per group";
      return;
    }
    const NB = 20, h0 = new Array(NB).fill(0), h1 = new Array(NB).fill(0);
    let nNull = 0, nBelow = 0;
    p.forEach((v, i) => {
      const b = Math.min(NB - 1, Math.floor(v * NB));
      if (isDiff[i]) h1[b]++; else { h0[b]++; nNull++; if (v < 0.05) nBelow++; }
    });
    const L = 54, R = W - 10, T = 108, B = H - 50;
    const ymax = Math.max(...h0.map((v, i) => v + h1[i]), nNull / NB * 1.3, 10);
    const step = niceStep(ymax / 5), top = Math.ceil(ymax / step) * step;
    const sx = scale(0, 1, L, R), sy = scale(0, top, B, T);
    const t2 = svgEl("text", { x: 8, y: 46, class: "sm" }, svg);
    t2.textContent = `no-difference CpGs with p < 0.05: ${(100 * nBelow / Math.max(1, nNull)).toFixed(1)}%`;
    t2.style.fill = nBelow / Math.max(1, nNull) > 0.075 ? COL.bad : COL.muted;
    const lg = svgEl("text", { x: 8, y: 68, class: "sm" }, svg);
    const a1 = svgEl("tspan", {}, lg); a1.textContent = "■ no difference  "; a1.style.fill = COL.nul;
    const a2 = svgEl("tspan", {}, lg); a2.textContent = "■ truly differential"; a2.style.fill = COL.diff;
    for (let v = 0; v <= top + 1e-9; v += step) {
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svg).textContent = fmt(v);
    }
    const bw = sx(1 / NB) - sx(0);
    for (let b = 0; b < NB; b++) {
      const x = sx(b / NB) + 0.5;
      if (h0[b]) svgEl("rect", { x, width: bw - 1, y: sy(h0[b]), height: B - sy(h0[b]), fill: COL.nul, "fill-opacity": 0.8 }, svg);
      if (h1[b]) svgEl("rect", { x, width: bw - 1, y: sy(h0[b] + h1[b]), height: sy(h0[b]) - sy(h0[b] + h1[b]), fill: COL.diff, "fill-opacity": 0.85 }, svg);
    }
    const yl = sy(nNull / NB);
    svgEl("line", { x1: L, x2: R, y1: yl, y2: yl, stroke: COL.fg, "stroke-width": 2, "stroke-dasharray": "7 5" }, svg);
    svgEl("text", { x: (L + R) / 2, y: yl - 6, "text-anchor": "middle", class: "sm" }, svg).textContent = "valid test: flat at this line";
    for (const v of [0, 0.25, 0.5, 0.75, 1]) svgEl("text", { x: sx(v), y: B + 20, "text-anchor": "middle", class: "sm" }, svg).textContent = String(v);
    svgEl("text", { x: (L + R) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = "p-value";
    svgEl("text", { x: L - 6, y: T - 16, "text-anchor": "end", class: "sm" }, svg).textContent = "CpGs";
  }

  // ---------------- (b) calls at FDR 5% ----------------
  function drawB(): void {
    const svg = svgB;
    svg.replaceChildren();
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent = "(b) CpGs called at FDR 5% (BH)";
    const lg = svgEl("text", { x: 8, y: 46, class: "sm" }, svg);
    const a1 = svgEl("tspan", {}, lg); a1.textContent = "■ true calls  "; a1.style.fill = COL.ok;
    const a2 = svgEl("tspan", {}, lg); a2.textContent = "■ false calls  "; a2.style.fill = COL.bad;
    const a3 = svgEl("tspan", {}, lg); a3.textContent = "- - all true"; a3.style.fill = COL.fg;
    const nd = res!.s.fisher!.nDiff;
    const L = 54, R = W - 10, T = 100, B = H - 66;
    const ymax = Math.max(nd * 1.1, ...METHODS.map(m => res!.s[m.key]?.calls ?? 0), 10);
    const step = niceStep(ymax / 5), top = Math.ceil(ymax / step) * step;
    const sy = scale(0, top, B, T);
    for (let v = 0; v <= top + 1e-9; v += step) {
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svg).textContent = fmt(v);
    }
    const slot = (R - L) / 3, bw = slot * 0.5;
    METHODS.forEach((m, k) => {
      const x = L + slot * (k + 0.5);
      const s = res!.s[m.key];
      const lab = svgEl("text", { x, y: B + 20, "text-anchor": "middle", "font-weight": m.key === method ? 700 : 400 }, svg);
      lab.textContent = m.short;
      svgEl("text", { x, y: B + 38, "text-anchor": "middle", class: "sm" }, svg).textContent = m.sub;
      if (!s) {
        svgEl("text", { x, y: B - 12, "text-anchor": "middle", class: "sm" }, svg).textContent = "needs ≥ 2";
        return;
      }
      if (s.trueCalls) svgEl("rect", { x: x - bw / 2, width: bw, y: sy(s.trueCalls), height: B - sy(s.trueCalls), fill: COL.ok, "fill-opacity": 0.85 }, svg);
      if (s.falseCalls) svgEl("rect", { x: x - bw / 2, width: bw, y: sy(s.calls), height: sy(s.trueCalls) - sy(s.calls), fill: COL.bad, "fill-opacity": 0.85 }, svg);
      const t1 = svgEl("text", { x, y: sy(s.calls) - 25, "text-anchor": "middle" }, svg);
      t1.textContent = `${fmt(s.calls)} calls`;
      const t2 = svgEl("text", { x, y: sy(s.calls) - 7, "text-anchor": "middle", class: "sm" }, svg);
      t2.textContent = s.calls ? `${fmt(s.falseCalls)} false` : "";
      t2.style.fill = s.calls && s.falseCalls / s.calls > 0.1 ? COL.bad : COL.muted;
      if (nd) svgEl("text", { x, y: B + 58, "text-anchor": "middle", class: "sm" }, svg).textContent = `found ${Math.round(100 * s.trueCalls / nd)}%`;
    });
    if (nd) {
      svgEl("line", { x1: L, x2: R, y1: sy(nd), y2: sy(nd), stroke: COL.fg, "stroke-width": 2, "stroke-dasharray": "7 5" }, svg);
    }
  }

  compute();
  draw();
}

function niceStep(x: number): number {
  const p = 10 ** Math.floor(Math.log10(Math.max(x, 1e-9)));
  const m = x / p;
  return Math.max(1, (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p);
}

autoMount(NAME, mount);
