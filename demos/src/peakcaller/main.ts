// A Poisson peak caller in the style of MACS2/3, on simulated ChIP-seq (model.ts).
// Top track: ChIP pileup of extended fragments (brown), input pileup (grey, with an input),
// and lambda_local (blue); a strip below it says which estimate set lambda (BG, d, 1k, 10k).
// Bottom track: -log10 p or -log10 q (BH over the scored positions of this 20 kb), the
// cutoff, and the called peaks (green = at a true site, crimson = false).
// A draggable cursor gives the numbers of the by-hand calculation at one base.
//
// Auto-mounts on <div data-demo="peakcaller">. data-* options (all optional):
//   data-depth="20"  data-input="100"  data-signal="30"  data-bg="1"  data-bump="1"
//   data-control="off"  data-local="off"  data-score="p"  data-cut="2"  data-seed="1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount, draggable, svgPoint } from "../lib/ui";
import * as M from "./model";

const NAME = "peakcaller";
const COL = {
  chip: "#8c564b", input: "#9aa5b1", lam: "#2a6f97", score: "#5b6b7a", cut: "#c0392b",
  tp: "#2e8b57", fp: "#a8324a", ink: "#1b2a3a", bump: "#f3e6d8",
};
const SRC_COL: Record<M.Source, string> = { BG: "#c9d3dd", d: "#e67e22", "1k": "#6b4fa0", "10k": "#2a6f97" };
const SITE_NAMES = ["A", "B", "C", "D"];
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .pc-controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(22em, 100%), 1fr)); gap: .35em 1.4em; width: 100%; }
[data-demo="${NAME}"] .pc-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .pc-row .qb-presets { display: inline-flex; }
[data-demo="${NAME}"] svg.pc-main { width: 100%; height: auto; max-width: 1000px; }
[data-demo="${NAME}"] svg text { font-size: 16px; }
[data-demo="${NAME}"] svg .ttl { font-size: 17px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #5b6b7a; }
[data-demo="${NAME}"] .pc-out { display: flex; flex-wrap: wrap; gap: .6em; width: 100%; }
[data-demo="${NAME}"] .pc-out .qb-readout { flex: 1 1 22em; margin: 0; }
[data-demo="${NAME}"] table.pc-sites { border-collapse: collapse; font-size: .9em; flex: 1 1 22em; }
[data-demo="${NAME}"] table.pc-sites td, [data-demo="${NAME}"] table.pc-sites th { white-space: nowrap; padding: .1em .55em; text-align: right; border-bottom: 1px solid #e3e3e3; }
[data-demo="${NAME}"] table.pc-sites th { font-weight: 600; }
[data-demo="${NAME}"] .pc-yes { color: ${COL.tp}; font-weight: 600; }
[data-demo="${NAME}"] .pc-no { color: ${COL.fp}; font-weight: 600; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
.reveal [data-demo="${NAME}"] { font-size: 18px; }
.reveal [data-demo="${NAME}"] .pc-out { display: none; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const fmt = (v: number, d = 2): string => (Number.isFinite(v) ? v.toFixed(d) : "–");
const fmtScore = (v: number): string => (v === Infinity ? "∞" : v >= 100 ? v.toFixed(0) : v.toFixed(1));
const SUP: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const sup = (v: number): string => String(v).split("").map(c => SUP[c] ?? c).join("");
const on = (s: string | undefined, dflt: boolean): boolean => (s === undefined ? dflt : s === "on" || s === "1" || s === "true");

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, d: number, lo: number, hi: number) => {
    const v = Number(s); return Number.isFinite(v) && s !== undefined ? Math.min(hi, Math.max(lo, v)) : d;
  };
  const DEF: M.Params = {
    depth: num(ds.depth, 20, 5, 80), inputPct: num(ds.input, 100, 10, 100), signal: num(ds.signal, 30, 5, 100),
    bg: num(ds.bg, 1, 0.5, 3), bump: num(ds.bump, 1, 1, 10), control: on(ds.control, true), local: on(ds.local, true),
  };
  const st = { p: { ...DEF }, useQ: ds.score !== "p", cut: num(ds.cut, 2, 1, 20), seed: num(ds.seed, 1, 1, 1e6), cursor: M.SITES[0] };
  let pools = M.makePools(st.seed);

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "pc-controls", "", host);
  const sl = {
    depth: slider(ctr, "ChIP reads (M)", 5, 80, 5, st.p.depth, v => { st.p.depth = v; update(); }),
    input: slider(ctr, "input (% of ChIP)", 10, 100, 10, st.p.inputPct, v => { st.p.inputPct = v; update(); }, v => `${v}%`),
    signal: slider(ctr, "signal at A per 10 M", 5, 100, 5, st.p.signal, v => { st.p.signal = v; update(); }),
    bg: slider(ctr, "background (×)", 0.5, 3, 0.25, st.p.bg, v => { st.p.bg = v; update(); }, v => `×${v}`),
    bump: slider(ctr, "bump in 2 kb", 1, 10, 1, st.p.bump, v => { st.p.bump = v; update(); }, v => (v === 1 ? "off" : `×${v}`)),
    cut: slider(ctr, "cutoff −log10", 1, 20, 1, st.cut, v => { st.cut = v; update(); }),
  };
  const row = htmlEl("div", "pc-row", "", host);
  const cbCtl = checkbox(row, "use input control", st.p.control, v => { st.p.control = v; update(); });
  const cbLoc = checkbox(row, "local λ windows", st.p.local, v => { st.p.local = v; update(); });
  const scoreBtns = presets(row, [
    { label: "score: p", apply: () => { st.useQ = false; update(); } },
    { label: "score: q (BH)", apply: () => { st.useQ = true; update(); } },
  ]);
  scoreBtns[st.useQ ? 1 : 0].classList.add("active");
  button(row, "new sample", () => { st.seed++; pools = M.makePools(st.seed); update(); });

  const row2 = htmlEl("div", "pc-row", "", host);
  htmlEl("span", "", "scenarios:", row2);
  const setAll = (p: Partial<M.Params>, cut?: number) => {
    Object.assign(st.p, DEF, p);
    if (cut !== undefined) st.cut = cut;
    sl.depth.set(st.p.depth); sl.input.set(st.p.inputPct); sl.signal.set(st.p.signal);
    sl.bg.set(st.p.bg); sl.bump.set(st.p.bump); sl.cut.set(st.cut);
    cbCtl.checked = st.p.control; cbLoc.checked = st.p.local;
    update();
  };
  presets(row2, [
    { label: "start", apply: () => setAll({}, 2) },
    { label: "open-chromatin bump", apply: () => setAll({ bump: 8 }, 2) },
    { label: "… without input", apply: () => setAll({ bump: 8, control: false }, 2) },
    { label: "shallow input", apply: () => setAll({ inputPct: 10 }, 2) },
    { label: "deep", apply: () => setAll({ depth: 80 }, 2) },
  ]);

  htmlEl("p", "qb-note",
    "Simulated data: 20 kb of a 2.7-Gb genome, fragments d = 200 bp, four true binding sites A–D (▲). " +
    "q-values here are BH over the 2000 scored positions of this window; MACS does it over the whole genome. " +
    "Drag the cursor to read the numbers at one base.", host);

  // ---------------- SVG ----------------
  const W = 800, H = 506, L = 64, R = W - 16;
  const aT = 76, aB = 262, stripY = 270, stripH = 10, bT = 312, bB = 436, pkY = 444, pkH = 12;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, class: "pc-main", role: "img",
    "aria-label": "Simulated ChIP pileup, local lambda, -log10 p and called peaks" }, host);
  const gx = scale(0, M.SHOW, L, R);
  const out = htmlEl("div", "pc-out", "", host);
  const read = htmlEl("div", "qb-readout", "", out);
  const tbl = htmlEl("table", "pc-sites", "", out);

  const gDyn = svgEl("g", {}, svg);
  const gCur = svgEl("g", {}, svg);
  const curLine = svgEl("line", { y1: aT - 8, y2: pkY + pkH, stroke: COL.ink, "stroke-width": 1.5, "stroke-dasharray": "4 3" }, gCur);
  const curHandle = svgEl("rect", { y: aT - 22, width: 22, height: 18, rx: 4, fill: "#fff", stroke: COL.ink, "stroke-width": 2 }, gCur);
  const hit = svgEl("rect", { x: L, y: aT - 22, width: R - L, height: pkY + pkH - aT + 22, fill: "transparent" }, gCur);
  const toWorld = (ev: PointerEvent): [number, number] => { const [x] = svgPoint(svg, ev); return [gx.invert(x), 0]; };
  const setCursor = (x: number) => { st.cursor = Math.max(0, Math.min(M.SHOW - 1, x)); drawCursor(); };
  draggable(hit, toWorld, x => setCursor(x));
  hit.addEventListener("pointerdown", ev => setCursor(toWorld(ev)[0]));

  let res: M.Result;

  function sqrtAxis(max: number) {
    const top = Math.sqrt(max);
    const y = (v: number) => aB - (aB - aT) * Math.sqrt(Math.max(0, v)) / top;
    return y;
  }

  function update(): void {
    const { chip, input } = M.draw(st.p, pools);
    res = M.call(st.p, chip, input);
    const score = st.useQ ? res.qScore : res.pScore;
    const pk = M.peaks(res, score, st.cut);
    gDyn.replaceChildren();
    const n = res.pos.length;

    // title
    const t = svgEl("text", { x: 4, y: 20, class: "ttl" }, gDyn);
    t.textContent = `Simulated ChIP: ${st.p.depth} M reads` +
      (st.p.control ? `, input ${Math.round(st.p.depth * st.p.inputPct / 100 * 10) / 10} M` : ", no input") +
      (res.scale < 1 ? `  (ChIP pileup × ${fmt(res.scale, 2)})` : "");

    // ---- track A ----
    let ymax = 2;
    for (let i = 0; i < n; i++) ymax = Math.max(ymax, res.k[i], res.lam[i], st.p.control ? res.inputPile[i] : 0);
    const ticks = [0, 2, 5, 10, 20, 50, 100, 200, 500, 1000].filter(v => v <= ymax * 1.05);
    ymax = Math.max(ymax, ticks[ticks.length - 1]) * 1.04;
    const ya = sqrtAxis(ymax);
    if (st.p.bump > 1) {
      svgEl("rect", { x: gx(M.BUMP.from), width: gx(M.BUMP.to) - gx(M.BUMP.from), y: aT, height: aB - aT, fill: COL.bump }, gDyn);
      const bt = svgEl("text", { x: gx((M.BUMP.from + M.BUMP.to) / 2), y: aT + 16, "text-anchor": "middle", class: "sm" }, gDyn);
      bt.textContent = `bump ×${st.p.bump}`;
    }
    for (const v of ticks) {
      svgEl("line", { x1: L, x2: R, y1: ya(v), y2: ya(v), class: v === 0 ? "axis" : "grid" }, gDyn);
      svgEl("text", { x: L - 6, y: ya(v) + 5, "text-anchor": "end", class: "sm" }, gDyn).textContent = String(v);
    }
    const yl = svgEl("text", { x: 14, y: (aT + aB) / 2, "text-anchor": "middle", class: "sm", transform: `rotate(-90 14 ${(aT + aB) / 2})` }, gDyn);
    yl.textContent = "fragments (√ axis)";
    const area = (vals: Float64Array, y: (v: number) => number, y0: number) => {
      let d = `M${gx(0)},${y0}`;
      for (let i = 0; i < n; i++) d += `L${gx(res.pos[i] - M.STEP / 2).toFixed(1)},${y(vals[i]).toFixed(1)}L${gx(res.pos[i] + M.STEP / 2).toFixed(1)},${y(vals[i]).toFixed(1)}`;
      return d + `L${gx(M.SHOW)},${y0}Z`;
    };
    const line = (vals: Float64Array, y: (v: number) => number) => {
      let d = "";
      for (let i = 0; i < n; i++) d += `${i ? "L" : "M"}${gx(res.pos[i] - M.STEP / 2).toFixed(1)},${y(vals[i]).toFixed(1)}L${gx(res.pos[i] + M.STEP / 2).toFixed(1)},${y(vals[i]).toFixed(1)}`;
      return d;
    };
    if (st.p.control) svgEl("path", { d: area(res.inputPile, ya, aB), fill: COL.input, "fill-opacity": 0.55 }, gDyn);
    svgEl("path", { d: area(res.k, ya, aB), fill: COL.chip, "fill-opacity": 0.8 }, gDyn);
    svgEl("path", { d: line(res.lam, ya), fill: "none", stroke: COL.lam, "stroke-width": 2.6 }, gDyn);
    // legend
    const leg: [string, string, boolean][] = [["ChIP pileup k", COL.chip, true], ["input pileup", COL.input, st.p.control], ["λ local", COL.lam, true]];
    let lx = L + 4;
    for (const [lab, c, show] of leg) {
      if (!show) continue;
      svgEl("rect", { x: lx, y: 32, width: 16, height: 12, fill: c }, gDyn);
      const tt = svgEl("text", { x: lx + 21, y: 43 }, gDyn); tt.textContent = lab;
      lx += 30 + lab.length * 8.4;
    }
    // sites
    M.SITES.forEach((s, j) => {
      svgEl("path", { d: `M${gx(s)},${aB - 1} l-7,12 h14 z`, fill: COL.ink }, gDyn);
      const tt = svgEl("text", { x: gx(s) + 10, y: aB + 12, "font-weight": 600 }, gDyn);
      tt.textContent = SITE_NAMES[j];
    });
    // lambda source strip
    for (let i = 0; i < n;) {
      let j = i;
      while (j < n && res.src[j] === res.src[i]) j++;
      svgEl("rect", { x: gx(res.pos[i] - M.STEP / 2), width: gx(res.pos[j - 1] + M.STEP / 2) - gx(res.pos[i] - M.STEP / 2), y: stripY + 6, height: stripH, fill: SRC_COL[res.src[i]] }, gDyn);
      i = j;
    }
    let sx = R - 4 * 62 - 70;
    const st0 = svgEl("text", { x: sx, y: 43, class: "sm" }, gDyn); st0.textContent = "λ set by:";
    sx += 70;
    for (const s of ["BG", "d", "1k", "10k"] as M.Source[]) {
      svgEl("rect", { x: sx, y: 32, width: 14, height: 12, fill: SRC_COL[s] }, gDyn);
      const tt = svgEl("text", { x: sx + 18, y: 43, class: "sm" }, gDyn);
      tt.textContent = s === "BG" ? "λBG" : `λ${s}`;
      sx += 62;
    }
    const sl0 = svgEl("text", { x: L - 6, y: stripY + 15, "text-anchor": "end", class: "sm" }, gDyn); sl0.textContent = "λ by";

    // ---- track B ----
    let smax = Math.max(st.cut * 1.5, 5);
    for (let i = 0; i < n; i++) if (Number.isFinite(score[i])) smax = Math.max(smax, score[i]);
    const yb = (v: number) => bB - (bB - bT) * Math.sqrt(Math.min(v, smax)) / Math.sqrt(smax);
    const sticks = [0, 2, 10, 50, 200].filter(v => v <= smax);
    for (const v of sticks) {
      svgEl("line", { x1: L, x2: R, y1: yb(v), y2: yb(v), class: v === 0 ? "axis" : "grid" }, gDyn);
      svgEl("text", { x: L - 6, y: yb(v) + 5, "text-anchor": "end", class: "sm" }, gDyn).textContent = String(v);
    }
    const yl2 = svgEl("text", { x: 14, y: (bT + bB) / 2, "text-anchor": "middle", class: "sm", transform: `rotate(-90 14 ${(bT + bB) / 2})` }, gDyn);
    yl2.textContent = st.useQ ? "−log10 q (√)" : "−log10 p (√)";
    svgEl("path", { d: area(score, yb, bB), fill: COL.score, "fill-opacity": 0.75 }, gDyn);
    svgEl("line", { x1: L, x2: R, y1: yb(st.cut), y2: yb(st.cut), stroke: COL.cut, "stroke-width": 2, "stroke-dasharray": "7 4" }, gDyn);
    const ct = svgEl("text", { x: R - 4, y: yb(st.cut) - 6, "text-anchor": "end", "font-weight": 600 }, gDyn);
    ct.style.fill = COL.cut;
    ct.textContent = `cutoff ${st.useQ ? "q" : "p"} < 10${sup(-st.cut)}`;
    for (const p of pk) {
      svgEl("rect", { x: gx(p.from), width: Math.max(3, gx(p.to) - gx(p.from)), y: pkY, height: pkH, fill: p.site >= 0 ? COL.tp : COL.fp }, gDyn);
    }
    const pkLab = svgEl("text", { x: 4, y: pkY + 11, class: "sm" }, gDyn); pkLab.textContent = "peaks";
    // x axis
    for (let kb = 0; kb <= 20; kb += 2) {
      svgEl("line", { x1: gx(kb * 1000), x2: gx(kb * 1000), y1: pkY + pkH + 2, y2: pkY + pkH + 7, class: "axis" }, gDyn);
      svgEl("text", { x: gx(kb * 1000), y: pkY + pkH + 24, "text-anchor": "middle", class: "sm" }, gDyn).textContent = `${kb}`;
    }
    svgEl("text", { x: R, y: H - 4, "text-anchor": "end", class: "sm" }, gDyn).textContent = "position (kb)";
    const nTrue = pk.filter(p => p.site >= 0).length, nFalse = pk.length - nTrue;
    const sum = svgEl("text", { x: L + 4, y: bT - 6 }, gDyn);
    const a = svgEl("tspan", {}, sum); a.textContent = `${pk.length} peaks: `;
    const b = svgEl("tspan", { "font-weight": 600 }, sum); b.style.fill = COL.tp; b.textContent = `${nTrue} at true sites`;
    const c = svgEl("tspan", {}, sum); c.textContent = ", ";
    const e = svgEl("tspan", { "font-weight": 600 }, sum); e.style.fill = nFalse ? COL.fp : "#5b6b7a"; e.textContent = `${nFalse} false`;

    // ---- site table ----
    tbl.replaceChildren();
    const hr = htmlEl("tr", "", "", tbl);
    for (const h of ["site", "strength", "k", "λ", "−log10 p", "−log10 q", "called"]) htmlEl("th", "", h, hr);
    M.SITES.forEach((s, j) => {
      const i = M.indexAt(s);
      const tr = htmlEl("tr", "", "", tbl);
      const called = pk.some(p => p.site === j);
      for (const v of [SITE_NAMES[j], `${Math.round(M.SITE_REL[j] * 100)}%`, String(res.k[i]),
        `${fmt(res.lam[i])} ${res.src[i]}`, fmtScore(res.pScore[i]), fmtScore(res.qScore[i])]) htmlEl("td", "", v, tr);
      htmlEl("td", called ? "pc-yes" : "pc-no", called ? "yes" : "no", tr);
    });
    drawCursor();
  }

  function drawCursor(): void {
    const x = gx(st.cursor);
    curLine.setAttribute("x1", String(x)); curLine.setAttribute("x2", String(x));
    curHandle.setAttribute("x", String(x - 11));
    const i = M.indexAt(st.cursor);
    const P = res.parts;
    const lamParts = st.p.control
      ? `λBG ${fmt(P.BG)}  λd ${fmt(P.d[i])}  λ1k ${fmt(P["1k"][i])}  λ10k ${fmt(P["10k"][i])}`
      : `λBG ${fmt(P.BG)}  λ10k ${fmt(P["10k"][i])} (from the ChIP)`;
    const kTxt = res.scale < 1 ? `k = ⌈${fmt(res.rawChip[i], 0)} × ${fmt(res.scale, 2)}⌉ = ${res.k[i]}` : `k = ${res.k[i]}`;
    read.textContent =
      `cursor at ${(res.pos[i] / 1000).toFixed(2)} kb   ChIP pileup ${kTxt}\n` +
      `${lamParts}${st.p.local ? "" : "  (local windows off)"}\n` +
      `λ local = ${fmt(res.lam[i])} (${res.src[i]})   fold enrichment (k+1)/(λ+1) = ${fmt((res.k[i] + 1) / (res.lam[i] + 1), 1)}\n` +
      `p = P(X ≥ ${res.k[i]} | Poisson ${fmt(res.lam[i])}):  −log10 p = ${fmtScore(res.pScore[i])}   −log10 q = ${fmtScore(res.qScore[i])}`;
  }

  update();
}

autoMount(NAME, mount);
