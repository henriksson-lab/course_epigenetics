// motif2x2: motif enrichment as a 2×2 table, and why the threshold and the background matter.
// Simulated peaks (1000 target, 8000 other peaks, 300 bp, i.i.d. bases at each peak's GC).
// A PU.1-like motif is planted in some target peaks; a GC box (SP1-like) is planted nowhere.
// (a) the two 2×2 tables with fractions, odds ratio and the exact one-sided Fisher /
//     hypergeometric p-value, against all other peaks or a GC-matched background;
// (b) GC content of target, all other peaks and the GC-matched background;
// (c) odds ratio against the match threshold, for both motifs and both backgrounds.
//
// Auto-mounts on <div data-demo="motif2x2">. data-* options (all optional):
//   data-gc="0.08"     target GC minus background GC        data-bound="0.25"  fraction bound
//   data-p="5e-5"      per-window match threshold           data-bg="all" | "gc"
//   data-motifs="879"  motifs tested (for Bonferroni)        data-seed="1"
import { svgEl, htmlEl, slider, button, presets, scale, autoMount } from "../lib/ui";
import { simulate, table, scoreForP, MOTIFS, DEFAULTS, Data, Params, Table, meanOf } from "./model";

const NAME = "motif2x2";
const C = { t: "#e67e22", bg: "#5b6b7a", gcm: "#2a6f97", m: ["#2e8b57", "#a8324a"], fg: "#1b2a3a", ok: "#2e8b57", bad: "#c0392b" };
const BASE_COL = ["#2e8b57", "#2a6f97", "#d4a017", "#c0392b"];   // A C G T, as in common logos
const STYLE = `
[data-demo="${NAME}"].qb-demo { display: block; }
[data-demo="${NAME}"] .mo-controls { display: flex; flex-direction: column; gap: .35em; margin-bottom: .4em; }
[data-demo="${NAME}"] .mo-row { display: flex; flex-wrap: wrap; gap: .35em .6em; align-items: center; }
[data-demo="${NAME}"] .mo-lab { color: var(--qb-muted); min-width: 7.5em; }
[data-demo="${NAME}"] .mo-sliders { display: grid; grid-template-columns: repeat(auto-fill, minmax(440px, 1fr)); gap: .25em 1.6em; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 16em minmax(5em, 1fr) 4.6em; }
[data-demo="${NAME}"] .qb-slider-name { white-space: normal; }
[data-demo="${NAME}"] .mo-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .mo-plots svg { flex: 1 1 460px; min-width: 0; max-width: 640px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 17px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #555; }
[data-demo="${NAME}"] svg .num { font-size: 19px; font-family: ui-monospace, monospace; }
.reveal [data-demo="${NAME}"] { font-size: 18px; }
@media (max-width: 760px) {
  [data-demo="${NAME}"] .mo-plots svg { max-width: 100%; flex-basis: 100%; }
  [data-demo="${NAME}"] .mo-sliders { grid-template-columns: 1fr; }
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

/** 3.2e-5 style, or 10^-x for tiny values; input is log10. */
function fmtP(l10: number): string {
  if (l10 > -3) return (10 ** l10).toPrecision(2).replace(/\.?0+$/, "");
  const e = Math.floor(l10), m = 10 ** (l10 - e);
  const ms = m.toFixed(1);
  return ms === "10.0" ? `1 × 10${sup(e + 1)}` : `${ms} × 10${sup(e)}`;
}
function sup(n: number): string {
  const map: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
  return String(n).split("").map(ch => map[ch] ?? ch).join("");
}
function fmtThr(l10: number): string {
  const e = Math.floor(l10 + 1e-9), m = Math.round(10 ** (l10 - e) * 10) / 10;
  return m === 1 ? `1e${e}` : `${m}e${e}`;
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (v: string | undefined, d: number) => { const x = Number(v); return v !== undefined && Number.isFinite(x) ? x : d; };
  const p: Params = { ...DEFAULTS, seed: num(ds.seed, 1) || 1, gcShift: num(ds.gc, DEFAULTS.gcShift), fracBound: num(ds.bound, DEFAULTS.fracBound) };
  const nMotifs = Math.max(1, Math.round(num(ds.motifs, 879)));
  let l10thr = Math.max(-6, Math.min(-3, Math.log10(num(ds.p, 5e-5))));
  let gcBg = ds.bg === "gc";
  let data: Data = simulate(p);

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "mo-controls", "", host);
  const r1 = htmlEl("div", "mo-row", "", ctr);
  htmlEl("span", "mo-lab", "Background:", r1);
  const bgBtns = presets(r1, [
    { label: "all other peaks", apply: () => { gcBg = false; draw(); } },
    { label: "GC-matched peaks", apply: () => { gcBg = true; draw(); } },
  ]);
  bgBtns[gcBg ? 1 : 0].classList.add("active");
  const sl = htmlEl("div", "mo-sliders", "", ctr);
  slider(sl, "match threshold, p per window", 3, 6, 0.1, Math.round(-l10thr * 10) / 10,
    v => { l10thr = -v; draw(); }, v => fmtThr(-v));
  slider(sl, "target peaks with a PU.1 site", 0, 60, 5, Math.round(p.fracBound * 100),
    v => { p.fracBound = v / 100; draw(); }, v => `${v}%`);
  const sGc = slider(sl, "target GC − background GC", -0.04, 0.12, 0.01, p.gcShift,
    v => { p.gcShift = v; }, v => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(2));
  sGc.input.addEventListener("change", () => resim());   // re-simulating takes ~0.4 s: on release
  const r3 = htmlEl("div", "mo-row", "", ctr);
  button(r3, "draw again", () => { p.seed++; resim(); });
  const note = htmlEl("span", "qb-note", "", r3);

  const plots = htmlEl("div", "mo-plots", "", host);
  const W = 560, H = 470;
  const svgA = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "2×2 tables for two motifs" }, plots);
  const svgB = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "GC content and odds ratio against threshold" }, plots);

  function resim(): void {
    note.textContent = "simulating…";
    setTimeout(() => { data = simulate(p); draw(); }, 10);
  }

  function draw(): void {
    const L = p.len, w = MOTIFS[0].prob.length;
    const thrP = 10 ** l10thr;
    note.textContent = `Simulated: ${p.nTarget} target and ${p.nPool} other peaks of ${L} bp. ` +
      `Chance matches expected per peak at this threshold: ${(2 * (L - w + 1) * thrP).toPrecision(2)} (2 strands × ${L - w + 1} windows × p).`;
    const tabs = MOTIFS.map((_, m) => table(data, m, scoreForP(data.tails[m], thrP), p.fracBound, gcBg));
    drawTables(tabs);
    drawRight();
  }

  // ---------------- (a) the two 2×2 tables ----------------
  function drawTables(tabs: Table[]): void {
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = `(a) 2×2 tables, background: ${gcBg ? "GC-matched peaks" : "all other peaks"}`;
    const bonf = Math.log10(0.05 / nMotifs);
    tabs.forEach((t, m) => {
      const y0 = 44 + m * 212;
      const mo = MOTIFS[m];
      const ttl = svgEl("text", { x: 8, y: y0 + 22, "font-weight": 700 }, svgA);
      ttl.textContent = mo.name; ttl.style.fill = C.m[m];
      svgEl("text", { x: 8, y: y0 + 44, class: "sm" }, svgA).textContent = mo.note;
      logo(svgA, mo.prob, W - 196, y0 + 2, 188, 50);
      const cx = [210, 320, 440], ry = [y0 + 98, y0 + 126];
      ["match", "no match", "% match"].forEach((h, k) =>
        svgEl("text", { x: cx[k], y: y0 + 72, "text-anchor": "end", class: "sm" }, svgA).textContent = h);
      const rows: [string, number, number, string][] = [["target", t.a, t.b, C.t], ["background", t.c, t.d, gcBg ? C.gcm : C.bg]];
      rows.forEach(([lab, x, y, col], k) => {
        const lt = svgEl("text", { x: 8, y: ry[k] }, svgA); lt.textContent = lab; lt.style.fill = col; lt.style.fontWeight = "600";
        svgEl("text", { x: cx[0], y: ry[k], "text-anchor": "end", class: "num" }, svgA).textContent = String(x);
        svgEl("text", { x: cx[1], y: ry[k], "text-anchor": "end", class: "num" }, svgA).textContent = String(y);
        svgEl("text", { x: cx[2], y: ry[k], "text-anchor": "end", class: "num" }, svgA).textContent = (100 * x / (x + y)).toFixed(1) + "%";
      });
      svgEl("line", { x1: 120, x2: cx[2] + 4, y1: y0 + 79, y2: y0 + 79, class: "grid" }, svgA);
      const orS = Number.isFinite(t.or) ? t.or.toFixed(2) : "∞";
      const pass = t.log10p < bonf;
      const line = svgEl("text", { x: 8, y: y0 + 160 }, svgA);
      const add = (s: string, bold = false, col = "") => { const sp = svgEl("tspan", {}, line); sp.textContent = s; if (bold) sp.style.fontWeight = "700"; if (col) sp.style.fill = col; };
      add("OR "); add(orS, true); add("   p "); add(fmtP(t.log10p), true); add(`   expected a = ${t.expected.toFixed(0)}`);
      const l2 = svgEl("text", { x: 8, y: y0 + 184, class: "sm" }, svgA);
      l2.textContent = pass ? `✓ passes Bonferroni for ${nMotifs} motifs (p < ${fmtP(bonf)})` : `✗ fails Bonferroni for ${nMotifs} motifs (p < ${fmtP(bonf)})`;
      l2.style.fill = pass ? C.ok : C.bad;
    });
  }

  // ---------------- (b) GC histogram, (c) OR vs threshold ----------------
  function drawRight(): void {
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) GC content of the peaks";
    const hL = 50, hR = W - 14, hT = 62, hB = 176;
    const g0 = 0.25, g1 = 0.75, bw = 0.01, nb = Math.round((g1 - g0) / bw);
    const hist = (v: ArrayLike<number>, idx?: ArrayLike<number>) => {
      const h = new Float64Array(nb);
      const n = idx ? idx.length : v.length;
      for (let k = 0; k < n; k++) { const g = v[idx ? idx[k] : k]; const b = Math.floor((g - g0) / bw); if (b >= 0 && b < nb) h[b]++; }
      return h.map(x => x / n);
    };
    const hT_ = hist(data.gcT), hP = hist(data.gcP), hM = hist(data.gcP, data.matched);
    const ymax = Math.max(...hT_, ...hP, ...hM) * 1.05;
    const sx = scale(g0, g1, hL, hR), sy = scale(0, ymax, hB, hT);
    const mT = meanOf(data.gcT), mP = meanOf(data.gcP), mM = meanOf(data.gcP, data.matched);
    const leg = svgEl("text", { x: 8, y: 50, class: "sm" }, svgB);
    ([[`target ${mT.toFixed(2)}`, C.t], [`   all other peaks ${mP.toFixed(2)}`, C.bg], [`   GC-matched ${mM.toFixed(2)}`, C.gcm]] as [string, string][])
      .forEach(([s, c]) => { const sp = svgEl("tspan", {}, leg); sp.textContent = s; sp.style.fill = c; sp.style.fontWeight = "600"; });
    const area = (h: Float64Array, col: string, op: number) => {
      let d = `M${sx(g0)},${sy(0)}`;
      h.forEach((v, b) => { d += ` L${sx(g0 + b * bw).toFixed(1)},${sy(v).toFixed(1)} L${sx(g0 + (b + 1) * bw).toFixed(1)},${sy(v).toFixed(1)}`; });
      d += ` L${sx(g1)},${sy(0)} Z`;
      svgEl("path", { d, fill: col, "fill-opacity": op, stroke: col, "stroke-width": 1.5 }, svgB);
    };
    area(hP, C.bg, 0.25);
    area(hT_, C.t, 0.35);
    if (gcBg) area(hM, C.gcm, 0.0);
    svgEl("line", { x1: hL, x2: hR, y1: hB, y2: hB, class: "axis" }, svgB);
    for (let g = 0.3; g <= 0.71; g += 0.1) svgEl("text", { x: sx(g), y: hB + 19, "text-anchor": "middle", class: "sm" }, svgB).textContent = g.toFixed(1);
    svgEl("text", { x: hR, y: hB + 19, "text-anchor": "end", class: "sm" }, svgB).textContent = "GC";

    // (c) OR against threshold
    const cT = 236;
    svgEl("text", { x: 8, y: cT, class: "ttl" }, svgB).textContent = "(c) odds ratio against the threshold";
    const leg2 = svgEl("text", { x: 8, y: cT + 24, class: "sm" }, svgB);
    ([[`— ${MOTIFS[0].name}`, C.m[0]], [`   — ${MOTIFS[1].name}`, C.m[1]]] as [string, string][])
      .forEach(([s, c]) => { const sp = svgEl("tspan", {}, leg2); sp.textContent = s; sp.style.fill = c; sp.style.fontWeight = "600"; });
    svgEl("text", { x: 8, y: cT + 46, class: "sm" }, svgB).textContent =
      `solid: ${gcBg ? "GC-matched" : "all other peaks"}; dashed: ${gcBg ? "all other peaks" : "GC-matched"}; ● current threshold`;
    const cL = 50, cR = W - 30, cTop = cT + 62, cB = H - 40;
    const ORMIN = 0.5, ORMAX = 8;
    const sxx = scale(3, 6, cL, cR), syy = scale(Math.log2(ORMIN), Math.log2(ORMAX), cB, cTop);
    for (const v of [0.5, 1, 2, 4, 8]) {
      svgEl("line", { x1: cL, x2: cR, y1: syy(Math.log2(v)), y2: syy(Math.log2(v)), class: v === 1 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: cL - 6, y: syy(Math.log2(v)) + 5, "text-anchor": "end", class: "sm" }, svgB).textContent = String(v);
    }
    for (let e = 3; e <= 6; e++) svgEl("text", { x: sxx(e), y: cB + 19, "text-anchor": "middle", class: "sm" }, svgB).textContent = `1e−${e}`;
    svgEl("text", { x: cR, y: cB + 36, "text-anchor": "end", class: "sm" }, svgB).textContent = "threshold p per window (stricter →)";
    const cur = svgEl("line", { x1: sxx(-l10thr), x2: sxx(-l10thr), y1: cTop, y2: cB, stroke: "#555", "stroke-width": 1.4, "stroke-dasharray": "3 3" }, svgB);
    void cur;
    for (let m = 0; m < MOTIFS.length; m++) {
      for (const gm of [!gcBg, gcBg]) {
        let d = "", pen = false;
        for (let x = 3; x <= 6.001; x += 0.1) {
          const t = table(data, m, scoreForP(data.tails[m], 10 ** -x), p.fracBound, gm);
          if (t.a < 3 || t.c < 3 || !Number.isFinite(t.or) || t.or <= 0) { pen = false; continue; }
          const y = syy(Math.max(Math.log2(ORMIN), Math.min(Math.log2(ORMAX), Math.log2(t.or))));
          d += `${pen ? "L" : "M"}${sxx(x).toFixed(1)},${y.toFixed(1)} `;
          pen = true;
        }
        if (d) svgEl("path", { d, fill: "none", stroke: C.m[m], "stroke-width": gm === gcBg ? 3 : 1.8, "stroke-dasharray": gm === gcBg ? "" : "6 4" }, svgB);
      }
      const t = table(data, m, scoreForP(data.tails[m], 10 ** l10thr), p.fracBound, gcBg);
      if (t.a > 0 && t.c > 0) svgEl("circle", { cx: sxx(-l10thr), cy: syy(Math.max(Math.log2(ORMIN), Math.min(Math.log2(ORMAX), Math.log2(t.or)))), r: 6, fill: C.m[m] }, svgB);
    }
    svgEl("text", { x: 14, y: (cTop + cB) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${(cTop + cB) / 2})`, class: "sm" }, svgB).textContent = "odds ratio";
  }

  draw();
}

/** A small sequence logo (information content, letters stacked by probability). */
function logo(svg: SVGSVGElement, prob: number[][], x0: number, y0: number, w: number, h: number): void {
  const g = svgEl("g", {}, svg);
  const cw = w / prob.length;
  svgEl("line", { x1: x0, x2: x0 + w, y1: y0 + h, y2: y0 + h, class: "axis" }, g);
  prob.forEach((row, k) => {
    const ic = 2 + row.reduce((a, q) => a + (q > 0 ? q * Math.log2(q) : 0), 0);
    const order = [0, 1, 2, 3].sort((a, b) => row[a] - row[b]);
    let y = y0 + h;
    for (const b of order) {
      const hh = row[b] * ic / 2 * h;
      if (hh < 1) continue;
      // letters are ~0.72 em tall for capital glyphs: scale a 100-px glyph to height hh
      const t = svgEl("text", { x: 0, y: 0, "font-weight": 700, "text-anchor": "middle", "font-family": "Helvetica, Arial, sans-serif",
        transform: `translate(${x0 + cw * (k + 0.5)},${y}) scale(${cw / 70},${hh / 72})` }, g);
      t.textContent = "ACGT"[b];
      t.style.fill = BASE_COL[b];
      t.style.fontSize = "100px";
      y -= hh;
    }
  });
}

autoMount(NAME, mount);
