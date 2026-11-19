// polycombmemory: read–write feedback, replication and memory of a Polycomb (H3K27me3) domain.
// Simulated: 150 cells, each with a domain of 60 nucleosomes, marked or unmarked.
// (a) kymographs (time × nucleosome) of four example cells, with the nucleation pulse
//     (e.g. cold at FLC) shaded and each replication marked;
// (b) fraction of the domain marked over time: 40 single cells and the mean of all 150;
// (c) histogram of the fraction marked per cell at the end (all-or-none or graded?).
// Controls: feedback strength, neighbourhood range, turnover, replication interval, pulse
// length, "start fully marked", presets, new cells.
//
// Auto-mounts on <div data-demo="polycombmemory">. data-* options (all optional):
//   data-preset="cold"  start preset: strong | weak | cold | long | nodiv
//   data-seed="1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount, Slider } from "../lib/ui";
import { simulate, meanFrac, fracSilenced, DEFAULTS, Params, Result } from "./model";

const NAME = "polycombmemory";
const C = { mark: "#2a6f97", unm: "#e9eef3", cold: "#9ecae1", fg: "#1b2a3a", cell: "#7f8c99", div: "#a8324a" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { display: block; }
[data-demo="${NAME}"] .pm-controls { display: flex; flex-direction: column; gap: .35em; margin-bottom: .4em; }
[data-demo="${NAME}"] .pm-row { display: flex; flex-wrap: wrap; gap: .35em .8em; align-items: center; }
[data-demo="${NAME}"] .pm-lab { color: var(--qb-muted); min-width: 4.5em; }
[data-demo="${NAME}"] .pm-sliders { display: grid; grid-template-columns: repeat(auto-fill, minmax(440px, 1fr)); gap: .25em 1.6em; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 16em minmax(5em, 1fr) 4.6em; }
[data-demo="${NAME}"] .qb-slider-name { white-space: normal; }
[data-demo="${NAME}"] .pm-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .pm-plots svg { flex: 1 1 460px; min-width: 0; max-width: 640px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 17px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #555; }
[data-demo="${NAME}"] svg image { image-rendering: pixelated; image-rendering: crisp-edges; }
.reveal [data-demo="${NAME}"] { font-size: 18px; }
@media (max-width: 760px) {
  [data-demo="${NAME}"] .pm-plots svg { max-width: 100%; flex-basis: 100%; }
  [data-demo="${NAME}"] .pm-sliders { grid-template-columns: 1fr; }
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

interface Preset { key: string; label: string; set: Partial<Params> }
const PRESETS: Preset[] = [
  { key: "strong", label: "strong feedback, start silenced", set: { feedback: 1.5, startMarked: true, pulse: 0, cycle: 22 } },
  { key: "weak", label: "weak feedback, start silenced", set: { feedback: 0.6, startMarked: true, pulse: 0, cycle: 22 } },
  { key: "nodiv", label: "weak feedback, no division", set: { feedback: 0.6, startMarked: true, pulse: 0, cycle: 0 } },
  { key: "cold", label: "2 days of cold", set: { feedback: 1.5, startMarked: false, pulse: 48, cycle: 22 } },
  { key: "long", label: "1 week of cold", set: { feedback: 1.5, startMarked: false, pulse: 168, cycle: 22 } },
];

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const p: Params = { ...DEFAULTS, seed: Number(ds.seed ?? 1) || 1 };
  const pr0 = PRESETS.find(x => x.key === ds.preset) ?? PRESETS[3];
  Object.assign(p, pr0.set);
  let res: Result = simulate(p);

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "pm-controls", "", host);
  const r1 = htmlEl("div", "pm-row", "", ctr);
  htmlEl("span", "pm-lab", "Presets:", r1);
  const pBtns = presets(r1, PRESETS.map(pr => ({ label: pr.label, apply: () => { Object.assign(p, pr.set); syncControls(); resim(); } })));
  pBtns[PRESETS.indexOf(pr0)].classList.add("active");
  const clear = () => pBtns.forEach(b => b.classList.remove("active"));
  const sl = htmlEl("div", "pm-sliders", "", ctr);
  // sliders update their value while dragging and re-simulate on release (~0.3 s per run)
  const mk = (label: string, min: number, max: number, step: number, key: keyof Params, fmt: (v: number) => string): Slider => {
    const s = slider(sl, label, min, max, step, p[key] as number,
      v => { (p as unknown as Record<string, number>)[key] = key === "range" && v >= 30 ? p.n : v; }, fmt);
    s.input.addEventListener("change", () => { clear(); resim(); });
    return s;
  };
  const sFb = mk("read–write feedback (per h)", 0, 3, 0.1, "feedback", v => v.toFixed(1));
  const sRange = mk("recruitment range (± nucleosomes)", 1, 30, 1, "range", v => (v >= 30 ? "all" : `±${v}`));
  const sOff = mk("turnover, mark lost (per h)", 0.02, 0.3, 0.01, "turnover", v => v.toFixed(2));
  const sCyc = mk("replication every", 0, 48, 2, "cycle", v => (v === 0 ? "never" : `${v} h`));
  const sPulse = mk("nucleation pulse (cold)", 0, 168, 6, "pulse", v => `${v} h`);
  const r3 = htmlEl("div", "pm-row", "", ctr);
  const cb = checkbox(r3, "start fully marked", p.startMarked, v => { p.startMarked = v; clear(); resim(); });
  button(r3, "new cells", () => { p.seed++; resim(); });
  const note = htmlEl("span", "qb-note", "", r3);
  function syncControls(): void {
    sFb.set(p.feedback); sRange.set(p.range); sOff.set(p.turnover); sCyc.set(p.cycle); sPulse.set(p.pulse); cb.checked = p.startMarked;
  }

  const plots = htmlEl("div", "pm-plots", "", host);
  const W = 560, H = 470;
  const svgA = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Kymographs of four simulated cells" }, plots);
  const svgB = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Fraction marked over time, and per cell at the end" }, plots);

  function resim(): void {
    note.textContent = "simulating…";
    setTimeout(() => { res = simulate(p); draw(); }, 10);
  }

  const tL = 72, tR = W - 22;
  const tx = () => scale(0, p.tEnd, tL, tR);

  function timeAxis(svg: SVGSVGElement, y: number, label: boolean): void {
    const sx = tx();
    svgEl("line", { x1: tL, x2: tR, y1: y, y2: y, class: "axis" }, svg);
    for (let t = 0; t <= p.tEnd; t += 96) {
      svgEl("line", { x1: sx(t), x2: sx(t), y1: y, y2: y + 5, class: "axis" }, svg);
      if (label) svgEl("text", { x: sx(t), y: y + 21, "text-anchor": "middle", class: "sm" }, svg).textContent = String(t);
    }
    if (label) svgEl("text", { x: tR, y: y + 40, "text-anchor": "end", class: "sm" }, svg).textContent = "time (h)";
  }

  function pulseBand(svg: SVGSVGElement, y0: number, y1: number, withLabel: boolean): void {
    if (p.pulse <= 0) return;
    const sx = tx();
    svgEl("rect", { x: sx(0), width: sx(p.pulse) - sx(0), y: y0, height: y1 - y0, fill: C.cold, "fill-opacity": 0.35 }, svg);
    if (withLabel) {
      const t = svgEl("text", { x: sx(0) + 4, y: y0 + 16, class: "sm" }, svg);
      t.textContent = "cold"; t.style.fill = "#2a6f97";
    }
  }

  function draw(): void {
    note.textContent = `Simulated: ${p.cells} cells × ${p.n} nucleosomes; marking ∝ feedback × (marked fraction nearby)²; each cell followed along one lineage.`;
    drawKymo();
    drawRight();
  }

  // ---------------- (a) kymographs ----------------
  function drawKymo(): void {
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) four cells: nucleosome × time";
    const leg = svgEl("text", { x: 8, y: 50, class: "sm" }, svgA);
    ([["■ H3K27me3", C.mark], ["   ▼ replication", C.div], ["   ┃ nucleation region", "#2a6f97"]] as [string, string][])
      .forEach(([s, c]) => { const sp = svgEl("tspan", {}, leg); sp.textContent = s; sp.style.fill = c; sp.style.fontWeight = "600"; });
    const N = p.n, nRec = res.times.length;
    const sx = tx();
    const top = 72, stripH = 72, gap = 14;
    const canvas = document.createElement("canvas");
    canvas.width = nRec; canvas.height = N;
    const ctx = canvas.getContext("2d")!;
    const hex = (h: string) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
    const cm = hex(C.mark), cu = hex(C.unm);
    res.kymo.forEach((km, e) => {
      const img = ctx.createImageData(nRec, N);
      for (let k = 0; k < nRec; k++) for (let i = 0; i < N; i++) {
        const o = 4 * (i * nRec + k), c = km[k * N + i] ? cm : cu;
        img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
      }
      ctx.putImageData(img, 0, 0);
      const y0 = top + e * (stripH + gap);
      svgEl("image", { x: tL, y: y0, width: tR - tL, height: stripH, preserveAspectRatio: "none", href: canvas.toDataURL() }, svgA);
      pulseBand(svgA, y0, y0 + stripH, false);
      svgEl("rect", { x: tL, y: y0, width: tR - tL, height: stripH, fill: "none", stroke: "#999", "stroke-width": 0.8 }, svgA);
      svgEl("text", { x: tL - 8, y: y0 + stripH / 2 - 4, "text-anchor": "end" }, svgA).textContent = `cell ${e + 1}`;
      for (const t of res.divisions[e]) {
        const x = sx(t);
        svgEl("path", { d: `M${x - 4},${y0 - 7} L${x + 4},${y0 - 7} L${x},${y0 - 1} Z`, fill: C.div }, svgA);
      }
      const n0 = Math.floor((N - p.nucWidth) / 2), yy = (i: number) => y0 + stripH * i / N;
      svgEl("line", { x1: tL - 3, x2: tL - 3, y1: yy(n0), y2: yy(n0 + p.nucWidth), stroke: "#2a6f97", "stroke-width": 3 }, svgA);
      const fe = res.final[res.examples[e]];
      svgEl("text", { x: tL - 8, y: y0 + stripH / 2 + 18, "text-anchor": "end", class: "sm" }, svgA).textContent = `${Math.round(100 * fe)}%`;
    });
    const yb = top + 4 * (stripH + gap) - gap + 18;
    timeAxis(svgA, yb, true);
    svgEl("text", { x: 8, y: yb + 40, class: "sm" }, svgA).textContent = `rows: ${N} nucleosomes; %: marked, last 24 h`;
  }

  // ---------------- (b) fraction marked over time, (c) end histogram ----------------
  function drawRight(): void {
    svgB.replaceChildren();
    const mf = meanFrac(res), sil = fracSilenced(res);
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) fraction of the domain marked";
    const leg = svgEl("text", { x: 8, y: 50, class: "sm" }, svgB);
    ([["— 40 single cells", C.cell], ["   — mean of all cells", C.mark]] as [string, string][])
      .forEach(([s, c]) => { const sp = svgEl("tspan", {}, leg); sp.textContent = s; sp.style.fill = c; sp.style.fontWeight = "600"; });
    const yT = 66, yB = 236;
    const sx = tx(), sy = scale(0, 1, yB, yT);
    pulseBand(svgB, yT, yB, true);
    for (const v of [0, 0.5, 1]) {
      svgEl("line", { x1: tL, x2: tR, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: tL - 6, y: sy(v) + 5, "text-anchor": "end", class: "sm" }, svgB).textContent = `${v * 100}%`;
    }
    timeAxis(svgB, yB, true);
    const path = (f: ArrayLike<number>) => {
      let d = "";
      for (let k = 0; k < f.length; k += 2) d += `${k ? "L" : "M"}${sx(res.times[k]).toFixed(1)},${sy(f[k]).toFixed(1)}`;
      return d;
    };
    const g = svgEl("g", { fill: "none", stroke: C.cell, "stroke-width": 1, "stroke-opacity": 0.45 }, svgB);
    for (let c = 0; c < Math.min(40, res.frac.length); c++) svgEl("path", { d: path(res.frac[c]) }, g);
    svgEl("path", { d: path(mf), fill: "none", stroke: C.mark, "stroke-width": 3.5 }, svgB);

    // (c) histogram
    const hT = 318;
    svgEl("text", { x: 8, y: hT - 10, class: "ttl" }, svgB).textContent = `(c) cells, last 24 h`;
    const big = svgEl("text", { x: W - 12, y: hT - 10, "text-anchor": "end", "font-weight": 700 }, svgB);
    big.textContent = `${Math.round(100 * sil)}% silenced (≥ 50% marked)`;
    big.style.fill = C.mark;
    const nb = 10, cnt = new Array(nb).fill(0);
    for (const v of res.final) cnt[Math.min(nb - 1, Math.floor(v * nb))]++;
    const hL = tL, hR = tR - 8, hTop = hT + 12, hB = H - 40;
    const hx = scale(0, 1, hL, hR), hy = scale(0, Math.max(...cnt), hB, hTop);
    cnt.forEach((n, b) => {
      if (n > 0) svgEl("rect", { x: hx(b / nb) + 1, width: hx(1 / nb) - hx(0) - 2, y: hy(n), height: hB - hy(n), fill: b >= nb / 2 ? C.mark : C.cell, "fill-opacity": 0.8 }, svgB);
      if (n > 0) svgEl("text", { x: hx((b + 0.5) / nb), y: hy(n) - 4, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(n);
    });
    svgEl("line", { x1: hL, x2: hR, y1: hB, y2: hB, class: "axis" }, svgB);
    for (const v of [0, 0.5, 1]) svgEl("text", { x: hx(v), y: hB + 20, "text-anchor": "middle", class: "sm" }, svgB).textContent = `${v * 100}%`;
    svgEl("text", { x: hR, y: hB + 38, "text-anchor": "end", class: "sm" }, svgB).textContent = "fraction of the domain marked, per cell";
  }

  draw();
}

autoMount(NAME, mount);
