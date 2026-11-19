// globalnorm: what each normalisation assumes, under a global change and under FRiP differences.
// Simulated libraries (control vs treated, 3 replicates each) over 2000 peaks, 2000 background
// bins and a spike-in. (a) MA plot of log2 fold changes after the chosen normalisation,
// coloured by the simulated truth; (b) size factors of each method against the true ones,
// FRiP per library, and the differential calls split into true and false.
//
// Auto-mounts on <div data-demo="globalnorm">. data-* options (all optional):
//   data-scenario="loss"   start preset: none | loss | lossfrip | frip | gain
//   data-method="mor"      start normalisation: total | mor | tmm | bg | spike
//   data-seed="1"
import { svgEl, htmlEl, slider, button, presets, scale, autoMount, Slider } from "../lib/ui";
import { simulate, sizeFactors, test, tally, METHODS, METHOD_LABEL, UP_LFC, DEFAULTS, Method, Params, Sim } from "./model";

const NAME = "globalnorm";
const C = { unch: "#8a8a8a", glob: "#2a6f97", up: "#e67e22", bad: "#c0392b", fg: "#1b2a3a", ctl: "#5b6b7a", trt: "#a8324a" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { display: block; }
[data-demo="${NAME}"] .gn-controls { display: flex; flex-direction: column; gap: .35em; margin-bottom: .4em; }
[data-demo="${NAME}"] .gn-row { display: flex; flex-wrap: wrap; gap: .35em .6em; align-items: center; }
[data-demo="${NAME}"] .gn-lab { color: var(--qb-muted); min-width: 7.5em; }
[data-demo="${NAME}"] .gn-sliders { display: grid; grid-template-columns: repeat(auto-fill, minmax(440px, 1fr)); gap: .25em 1.6em; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 16em minmax(5em, 1fr) 4.6em; }
[data-demo="${NAME}"] .qb-slider-name { white-space: normal; }
[data-demo="${NAME}"] .gn-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .gn-plots svg { flex: 1 1 460px; min-width: 0; max-width: 640px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 17px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #555; }
[data-demo="${NAME}"] svg .num { font-size: 22px; font-weight: 700; }
.reveal [data-demo="${NAME}"] { font-size: 18px; }
.reveal [data-demo="${NAME}"] .qb-slider { grid-template-columns: 13em minmax(5em, 1fr) 4.2em; }
.reveal [data-demo="${NAME}"] .gn-sliders { grid-template-columns: repeat(auto-fill, minmax(400px, 1fr)); }
@media (max-width: 760px) {
  [data-demo="${NAME}"] .gn-sliders { grid-template-columns: 1fr; }
  [data-demo="${NAME}"] .gn-plots svg { max-width: 100%; flex-basis: 100%; }
  [data-demo="${NAME}"] .qb-slider { grid-template-columns: minmax(0, 1fr) 4.2em; gap: 0 .5em; }
  [data-demo="${NAME}"] .qb-slider-name { grid-column: 1 / -1; }
}
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

interface Scenario { key: string; label: string; frac: number; lfc: number; bg: number }
const SCENARIOS: Scenario[] = [
  { key: "none", label: "no global change", frac: 0, lfc: -2, bg: 1 },
  { key: "loss", label: "global loss (EZH2 inhibitor)", frac: 0.9, lfc: -2, bg: 1 },
  { key: "lossfrip", label: "loss + noisier libraries", frac: 0.9, lfc: -2, bg: 2 },
  { key: "frip", label: "noisier treated libraries only", frac: 0, lfc: -2, bg: 2 },
  { key: "gain", label: "global gain", frac: 0.5, lfc: 1, bg: 1 },
];

const sgn = (v: number, d = 2): string => (v < 0 ? "−" : v > 0 ? "+" : "") + Math.abs(v).toFixed(d);

function rich(svg: SVGSVGElement, x: number, y: number, parts: [string, string][], cls = "sm", anchor = "start"): void {
  const t = svgEl("text", { x, y, class: cls, "text-anchor": anchor }, svg);
  for (const [s, c] of parts) {
    const sp = svgEl("tspan", {}, t);
    sp.textContent = s;
    if (c) { sp.style.fill = c; sp.style.fontWeight = "600"; }
  }
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const p: Params = { ...DEFAULTS, seed: Number(ds.seed ?? 1) || 1 };
  let method: Method = (METHODS as readonly string[]).includes(ds.method ?? "") ? ds.method as Method : "mor";
  const sc0 = SCENARIOS.find(s => s.key === ds.scenario) ?? SCENARIOS[1];
  p.fracChanged = sc0.frac; p.globalLfc = sc0.lfc; p.bgTreated = sc0.bg;
  let sim: Sim = simulate(p);

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "gn-controls", "", host);
  const r1 = htmlEl("div", "gn-row", "", ctr);
  htmlEl("span", "gn-lab", "Scenario:", r1);
  const scBtns = presets(r1, SCENARIOS.map(sc => ({ label: sc.label, apply: () => {
    p.fracChanged = sc.frac; p.globalLfc = sc.lfc; p.bgTreated = sc.bg;
    sFrac.set(sc.frac * 100); sLfc.set(sc.lfc); sBg.set(Math.log2(sc.bg));
    resim();
  } })));
  scBtns[SCENARIOS.indexOf(sc0)].classList.add("active");
  const clearScenario = () => scBtns.forEach(b => b.classList.remove("active"));
  const r2 = htmlEl("div", "gn-row", "", ctr);
  htmlEl("span", "gn-lab", "Normalise by:", r2);
  const mBtns = presets(r2, METHODS.map(m => ({ label: METHOD_LABEL[m], apply: () => { method = m; schedule(false); } })));
  mBtns[METHODS.indexOf(method)].classList.add("active");

  const sl = htmlEl("div", "gn-sliders", "", ctr);
  const sFrac: Slider = slider(sl, "peaks that change globally", 0, 100, 5, p.fracChanged * 100,
    v => { p.fracChanged = v / 100; clearScenario(); schedule(true); }, v => `${v}%`);
  const sLfc: Slider = slider(sl, "their log2 fold change", -3, 2, 0.25, p.globalLfc,
    v => { p.globalLfc = v; clearScenario(); schedule(true); }, v => sgn(v));
  const sUp: Slider = slider(sl, `truly up peaks (×${2 ** UP_LFC})`, 0, 300, 10, p.nUp,
    v => { p.nUp = v; schedule(true); });
  const sBg: Slider = slider(sl, "extra background in treated (technical)", -1, 2, 0.25, Math.log2(p.bgTreated),
    v => { p.bgTreated = 2 ** v; clearScenario(); schedule(true); }, v => `×${+(2 ** v).toFixed(2)}`);
  const sCv: Slider = slider(sl, "replicate variation (CV)", 0.05, 0.5, 0.05, p.repCv,
    v => { p.repCv = v; schedule(true); });
  void sUp; void sCv;
  const r3 = htmlEl("div", "gn-row", "", ctr);
  button(r3, "draw again", () => { p.seed++; resim(); });
  htmlEl("span", "qb-note",
    "Simulated: 2000 peaks, 2000 background bins and a spike-in; 3 control and 3 treated libraries of about 2 million reads. Test: moderated t on log counts, BH 5%.", r3);

  // ---------------- plots ----------------
  const plots = htmlEl("div", "gn-plots", "", host);
  const W = 560, H = 470;
  const svgA = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "MA plot of simulated peaks" }, plots);
  const svgB = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Size factors and calls" }, plots);

  let pending = 0, needSim = false;
  function schedule(resimulate: boolean): void {
    needSim = needSim || resimulate;
    if (pending) return;
    pending = requestAnimationFrame(() => { pending = 0; if (needSim) { sim = simulate(p); needSim = false; } draw(); });
  }
  function resim(): void { schedule(true); }

  function draw(): void {
    const sf = sizeFactors(sim, method);
    const te = test(sim, sf);
    const ta = tally(sim, te);
    drawMA(te, ta.medLfc);
    drawRight(sf, ta);
  }

  // ---------------- (a) MA plot ----------------
  const aL = 58, aR = W - 14, aT = 84, aB = H - 52, YL = 4.5;
  function drawMA(te: ReturnType<typeof test>, med: number): void {
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = `(a) MA plot, normalised by ${METHOD_LABEL[method]}`;
    rich(svgA, 8, 50, [["truth: ", ""], ["● unchanged", C.unch], ["  ● global change", C.glob], [`  ● truly up ×${2 ** UP_LFC}`, C.up]]);
    rich(svgA, 8, 72, [["solid = called (BH 5%),  ", ""], ["○ false call", C.bad], ["  ·  median log2 FC ", ""], [sgn(med), C.fg]]);
    let a0 = Infinity, a1 = -Infinity;
    for (const a of te.A) { if (a < a0) a0 = a; if (a > a1) a1 = a; }
    a0 = Math.floor(a0); a1 = Math.ceil(a1);
    const sx = scale(a0, a1, aL, aR), sy = scale(-YL, YL, aB, aT);
    for (let v = -4; v <= 4; v += 2) {
      svgEl("line", { x1: aL, x2: aR, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: aL - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svgA).textContent = v < 0 ? `−${-v}` : String(v);
    }
    for (let v = Math.ceil(a0 / 2) * 2; v <= a1; v += 2) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: aB, y2: aB + 5, class: "axis" }, svgA);
      svgEl("text", { x: sx(v), y: aB + 21, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(v);
    }
    svgEl("line", { x1: aL, x2: aR, y1: aB, y2: aB, class: "axis" }, svgA);
    svgEl("text", { x: (aL + aR) / 2, y: H - 8, "text-anchor": "middle" }, svgA).textContent = "mean log2 normalised count";
    const ym = (aT + aB) / 2;
    svgEl("text", { x: 16, y: ym, "text-anchor": "middle", transform: `rotate(-90 16 ${ym})` }, svgA).textContent = "log2 FC, treated / control";
    // where the truth is
    const truthLines: [number, string, string][] = p.nUp > 0 ? [[UP_LFC, C.up, `true ×${2 ** UP_LFC}`]] : [];
    if (p.fracChanged > 0 && p.globalLfc !== 0) truthLines.push([p.globalLfc, C.glob, "true global"]);
    for (const [v, c, lab] of truthLines) {
      svgEl("line", { x1: aL, x2: aR, y1: sy(v), y2: sy(v), stroke: c, "stroke-width": 1.6, "stroke-dasharray": "6 4" }, svgA);
      const t = svgEl("text", { x: aR - 4, y: sy(v) - 5, "text-anchor": "end", class: "sm" }, svgA);
      t.textContent = lab; t.style.fill = c;
    }
    const clip = (v: number) => Math.max(-YL, Math.min(YL, v));
    const col = (i: number) => (sim.kind[i] === 2 ? C.up : sim.kind[i] === 1 ? C.glob : C.unch);
    const faint = svgEl("g", { "fill-opacity": 0.22 }, svgA);
    const solid = svgEl("g", { "fill-opacity": 0.85 }, svgA);
    const bad = svgEl("g", { fill: "none", stroke: C.bad, "stroke-width": 1.3 }, svgA);
    // draw unchanged first so the coloured ones are on top
    const order = Array.from({ length: sim.p.nPeaks }, (_, i) => i).sort((a, b) => sim.kind[a] - sim.kind[b]);
    for (const i of order) {
      const cx = sx(te.A[i]).toFixed(1), cy = sy(clip(te.lfc[i])).toFixed(1);
      const called = te.call[i] !== 0;
      svgEl("circle", { cx, cy, r: called ? 2.8 : 2.4, fill: col(i) }, called ? solid : faint);
      if (called && Math.sign(sim.trueLfc[i]) !== te.call[i]) svgEl("circle", { cx, cy, r: 4.3 }, bad);
    }
  }

  // ---------------- (b) size factors, FRiP, calls ----------------
  function drawRight(sf: number[], ta: ReturnType<typeof tally>): void {
    svgB.replaceChildren();
    const nS = sim.nS;
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) size factors and calls";
    rich(svgB, 8, 50, [["bars: ", ""], [METHOD_LABEL[method], C.fg], ["   ▬ true (per-cell amount)", C.bad]]);
    const bL = 70, bR = W - 14, bT = 72, bB = 240;
    const smax = Math.max(2, ...sf, ...sim.trueSf) * 1.05;
    const sy = scale(0, smax, bB, bT);
    for (let v = 0; v <= smax; v += smax > 3 ? 1 : 0.5) {
      svgEl("line", { x1: bL, x2: bR, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: bL - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svgB).textContent = String(v);
    }
    const ym = (bT + bB) / 2;
    svgEl("text", { x: 16, y: ym, "text-anchor": "middle", transform: `rotate(-90 16 ${ym})`, class: "sm" }, svgB).textContent = "size factor";
    const slot = (bR - bL) / nS;
    for (let j = 0; j < nS; j++) {
      const x0 = bL + j * slot, bw = slot * 0.56, xc = x0 + slot / 2;
      const c = sim.cond[j] ? C.trt : C.ctl;
      svgEl("rect", { x: xc - bw / 2, width: bw, y: sy(sf[j]), height: bB - sy(sf[j]), fill: c, "fill-opacity": 0.75 }, svgB);
      svgEl("line", { x1: xc - bw / 2 - 6, x2: xc + bw / 2 + 6, y1: sy(sim.trueSf[j]), y2: sy(sim.trueSf[j]), stroke: C.bad, "stroke-width": 3.5 }, svgB);
      const lab = svgEl("text", { x: xc, y: bB + 20, "text-anchor": "middle", "font-weight": 600 }, svgB);
      lab.textContent = (sim.cond[j] ? "T" : "C") + ((j % sim.p.reps) + 1); lab.style.fill = c;
      svgEl("text", { x: xc, y: bB + 40, "text-anchor": "middle", class: "sm" }, svgB).textContent = `${Math.round(100 * sim.frip[j])}%`;
    }
    svgEl("text", { x: bL - 6, y: bB + 40, "text-anchor": "end", class: "sm" }, svgB).textContent = "FRiP";

    // calls
    const tT = 318;
    svgEl("text", { x: 8, y: tT - 12, "font-weight": 600 }, svgB).textContent = `calls (${ta.changed} peaks truly changed)`;
    const cols = [150, 285, 420];
    const head = ["", "true", "false"];
    head.forEach((h, k) => { if (h) svgEl("text", { x: cols[k], y: tT + 14, "text-anchor": "middle", class: "sm" }, svgB).textContent = h; });
    const rows: [string, number, number][] = [["called up", ta.upTrue, ta.upFalse], ["called down", ta.downTrue, ta.downFalse]];
    rows.forEach(([lab, tr, fa], k) => {
      const y = tT + 46 + k * 34;
      svgEl("text", { x: 8, y }, svgB).textContent = lab;
      const a = svgEl("text", { x: cols[1], y, "text-anchor": "middle", class: "num" }, svgB);
      a.textContent = String(tr);
      const b = svgEl("text", { x: cols[2], y, "text-anchor": "middle", class: "num" }, svgB);
      b.textContent = String(fa); if (fa > 0) b.style.fill = C.bad;
    });
    svgEl("text", { x: 8, y: tT + 114 }, svgB).textContent = "missed";
    svgEl("text", { x: cols[1], y: tT + 114, "text-anchor": "middle", class: "num" }, svgB).textContent = String(ta.missed);
    svgEl("text", { x: 8, y: tT + 140, class: "sm" }, svgB).textContent = "missed = truly changed but not called";
  }

  draw();
}

autoMount(NAME, mount);
