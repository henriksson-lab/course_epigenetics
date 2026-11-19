// One convolutional filter reading DNA (model.ts).
// Top: the input sequence, one-hot encoded (rows A, C, G, T); click a cell to mutate that base.
// The filter window (a slider or "play" moves it along) is outlined on the grid and its score
// (the sum of the filter weights of the bases present) is written above it.
// Middle: the activation track, the filter score at every position (both strands optional),
// with the -bias line: ReLU(score + bias) keeps only what is above it; global max pooling
// keeps the largest. Bottom: in-silico mutagenesis, the change in P(bound) for every
// substitution (click a cell to make it). Below: the filter as a logo and a heatmap (click a
// cell, then +/-), presets, and the dense layer: z = w * maxpool + b, P(bound) = sigmoid(z).
//
// Auto-mounts on <div data-demo="cnnfilter">. data-* options (all optional):
//   data-example="ctcf" | "b2"   which example to start with      data-seq="ACGT..." own sequence
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount, rng } from "../lib/ui";
import * as M from "./model";

const NAME = "cnnfilter";
const COL = {
  ink: "#1b2a3a", muted: "#5b6b7a", grid: "#e3e7eb", light: "#f4f6f8",
  fwd: "#8c564b", rev: "#2a6f97", pos: "#8c564b", neg: "#2a6f97", up: "#2e8b57", down: "#a8324a",
  win: "#e67e22", sel: "#c0392b",
};
const BASECOL = ["#2e8b57", "#2a6f97", "#d4a017", "#c0392b"];     // A C G T
const CTCF_SEQ = "CCTGCACAGGCCGAAGAACCACTAGAGGGCAGCAGCACACGGGCTCCATA";  // hg38 chr5:43105871-43105920
const B2_SEQ = "TCCAGACTGGA";
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .cf-row { display: flex; flex-wrap: wrap; gap: .35em 1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .cf-row input[type=text] { font: 15px ui-monospace, monospace; padding: .2em .4em; flex: 1 1 18em; min-width: 10em; }
[data-demo="${NAME}"] .cf-scroll { width: 100%; overflow-x: auto; }
[data-demo="${NAME}"] .cf-scroll svg { width: 100%; min-width: 760px; height: auto; display: block; }
[data-demo="${NAME}"] .cf-bottom { display: flex; flex-wrap: wrap; gap: .8em 1.6em; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .cf-filter { flex: 1 1 600px; min-width: 0; max-width: 660px; display: flex; flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .cf-filter svg { width: 100%; height: auto; }
[data-demo="${NAME}"] .cf-net { flex: 1 1 320px; min-width: 0; display: flex; flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .cf-label { font-weight: 600; }
[data-demo="${NAME}"] .cf-p { font-size: 1.5em; font-weight: 700; }
.qb-demo[data-demo="${NAME}"] .qb-readout { font-family: inherit; line-height: 1.5; }
[data-demo="${NAME}"] .qb-btn:disabled { opacity: .45; cursor: default; }
[data-demo="${NAME}"] svg text { font-size: 16px; }
[data-demo="${NAME}"] svg .ttl { font-size: 18px; font-weight: 700; }
[data-demo="${NAME}"] svg .sm { font-size: 15px; fill: ${COL.muted}; }
[data-demo="${NAME}"] svg .cell { cursor: pointer; }
.reveal [data-demo="${NAME}"] { font-size: 17px; gap: .3em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const f1 = (v: number): string => (Math.abs(v) < 0.05 ? "0" : v.toFixed(1));
const f2 = (v: number): string => v.toFixed(2);
const f3 = (v: number): string => v.toFixed(3);
const sgn = (v: number, d = 3): string => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(d);

/** Blend white -> colour by t in [0,1]. */
function tint(hex: string, t: number): string {
  const c = [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
  const m = c.map(v => Math.round(255 + (v - 255) * Math.max(0, Math.min(1, t))));
  return `rgb(${m[0]},${m[1]},${m[2]})`;
}

interface Example { filter: () => M.Filter; seq: string; bias: number; w: number; b: number; both: boolean }

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const CTCF_W = M.logOdds(M.CTCF_COUNTS);
  const EX: Record<string, Example> = {
    ctcf: { filter: () => CTCF_W.map(r => r.slice()), seq: CTCF_SEQ, bias: -10, w: 0.6, b: -3, both: true },
    b2: { filter: () => M.B2_FILTER.map(r => r.slice()), seq: B2_SEQ, bias: -5, w: 1.2, b: -2, both: false },
  };
  const ex0 = EX[ds.example ?? "ctcf"] ?? EX.ctcf;
  const userSeq = ds.seq ? M.cleanSeq(ds.seq) : "";
  const st = {
    filter: ex0.filter(), seq: userSeq.length >= 4 ? userSeq : ex0.seq,
    bias: ex0.bias, w: ex0.w, b: ex0.b, both: ex0.both,
    pos: 0, sel: null as null | [number, number], randSeed: 1, shufSeed: 1,
    filterName: ds.example === "b2" ? "B2 filter (CCAG)" : "CTCF PWM (JASPAR MA0139.2, log2 odds)",
  };
  const undo: string[] = [];
  let lastMsg = "";
  let timer = 0;

  // ---------------- sequence controls ----------------
  const r1 = htmlEl("div", "cf-row", "", host);
  htmlEl("span", "cf-label", "Sequence:", r1);
  const seqPre = presets(r1, [
    { label: "CTCF peak (chr5)", apply: () => setSeq(CTCF_SEQ, "") },
    { label: "reverse complement", apply: () => setSeq(M.revComp(st.seq), "reverse complement") },
    { label: "shuffle bases", apply: () => { const r = rng(1000 + st.shufSeed++); setSeq(M.shuffleSeq(st.seq, r.unif), "shuffled (same bases, new order)"); } },
    { label: "B2: TCCAGACTGGA", apply: () => setSeq(B2_SEQ, "") },
  ]);
  const r1b = htmlEl("div", "cf-row", "", host);
  const inp = htmlEl("input", "", "", r1b);
  inp.type = "text"; inp.spellcheck = false; inp.setAttribute("aria-label", "DNA sequence");
  inp.addEventListener("change", () => {
    const s = M.cleanSeq(inp.value);
    if (s.length < 4 || s.length > 80) { lastMsg = "Type 4–80 bases (A, C, G, T)."; inp.value = st.seq; draw(); }
    else setSeq(s, "typed sequence");
    inp.blur();
  });
  const bUndo = button(r1b, "undo", () => { const s = undo.pop(); if (s) { st.seq = s; lastMsg = "undone"; afterSeq(); } });
  const r1c = htmlEl("div", "cf-row", "", host);
  const sPos = slider(r1c, "filter position", 0, 1, 1, 0, v => { stop(); st.pos = v; draw(); }, v => String(v + 1));
  sPos.input.parentElement!.style.flex = "1 1 22em";
  const bPlay = button(r1c, "▶ slide along", () => play());

  // ---------------- the sequence pipeline (SVG) ----------------
  const scroll = htmlEl("div", "cf-scroll", "", host);
  // drawn at the width available (760-1100), so text keeps its nominal size
  const W = Math.max(760, Math.min(1100, host.clientWidth || 1100)), ml = 74, mr = 16;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} 520`, role: "img", "aria-label": "One-hot sequence, filter activation track and in-silico mutagenesis map" }, scroll);
  const gSeq = svgEl("g", {}, svg);

  // ---------------- filter + network ----------------
  const bottom = htmlEl("div", "cf-bottom", "", host);
  const fbox = htmlEl("div", "cf-filter", "", bottom);
  const fRow = htmlEl("div", "cf-row", "", fbox);
  htmlEl("span", "cf-label", "Filter:", fRow);
  const filtPre = presets(fRow, [
    { label: "CTCF PWM", apply: () => setFilter(CTCF_W.map(r => r.slice()), "CTCF PWM (JASPAR MA0139.2, log2 odds)", -10) },
    { label: "GC filter", apply: () => setFilter(M.gcFilter(8), "GC filter (+1 for C/G, −1 for A/T)", -4) },
    { label: "random", apply: () => { const r = rng(77 + st.randSeed++); setFilter(M.randomFilter(15, r.norm), "random filter (untrained)", -5); } },
    { label: "B2 (CCAG)", apply: () => setFilter(M.B2_FILTER.map(r => r.slice()), "B2 filter (CCAG)", -5) },
  ]);
  const fsvg = svgEl("svg", { viewBox: "0 0 640 300", role: "img", "aria-label": "Filter weights as a logo and a heatmap" }, fbox);
  const gF = svgEl("g", {}, fsvg);
  const eRow = htmlEl("div", "cf-row", "", fbox);
  const selLab = htmlEl("span", "", "", eRow);
  selLab.style.flex = "1 1 100%";
  const eBtns = [-1, -0.5, 0.5, 1].map(dv => button(eRow, (dv > 0 ? "+" : "−") + Math.abs(dv), () => editSel(v => v + dv)));
  eBtns.push(button(eRow, "set 0", () => editSel(() => 0)));

  const net = htmlEl("div", "cf-net", "", bottom);
  htmlEl("div", "cf-label", "Bias, ReLU, max pooling, dense layer", net);
  const sBias = slider(net, "filter bias", -30, 5, 0.5, st.bias, v => { st.bias = v; recompute(); });
  checkbox(net, "scan both strands (add the reverse-complement filter)", st.both, v => { st.both = v; recompute(); });
  const sW = slider(net, "dense weight w", -2, 2, 0.05, st.w, v => { st.w = v; recompute(); });
  const sB = slider(net, "dense bias b", -10, 5, 0.1, st.b, v => { st.b = v; recompute(); });
  const out = htmlEl("div", "qb-readout", "", net);
  const exRow = htmlEl("div", "cf-row", "", net);
  htmlEl("span", "", "Examples:", exRow);
  presets(exRow, [
    { label: "CTCF", apply: () => loadExample(EX.ctcf, "CTCF PWM (JASPAR MA0139.2, log2 odds)") },
    { label: "by-hand B2–B4", apply: () => loadExample(EX.b2, "B2 filter (CCAG)") },
  ]);
  const note = htmlEl("div", "qb-note", "", host);
  note.textContent = "The CTCF peak sequence is 50 bp of hg38 (chr5:43,105,871–43,105,920) around the summit of a GM12878 CTCF ChIP-seq peak. " +
    "The filter is not trained: the CTCF preset is the JASPAR PWM as log2 odds (pseudocount 1), and the bias and dense layer are set by hand.";

  // ---------------- state changes ----------------
  function setSeq(s: string, msg: string): void {
    if (s === st.seq) return;
    undo.push(st.seq); st.seq = s; lastMsg = msg; afterSeq();
  }
  function afterSeq(): void { stop(); inp.value = st.seq; clampPos(); recompute(); }
  function setFilter(f: M.Filter, name: string, bias: number): void {
    stop(); st.filter = f; st.filterName = name; st.sel = null;
    st.bias = bias; sBias.set(bias); clampPos(); lastMsg = ""; recompute();
  }
  function loadExample(e: Example, name: string): void {
    stop();
    if (st.seq !== e.seq) undo.push(st.seq);
    st.filter = e.filter(); st.filterName = name; st.seq = e.seq; st.sel = null;
    st.bias = e.bias; st.w = e.w; st.b = e.b; st.both = e.both;
    sBias.set(st.bias); sW.set(st.w); sB.set(st.b);
    (net.querySelector("input[type=checkbox]") as HTMLInputElement).checked = st.both;
    seqPre.forEach(b => b.classList.remove("active")); filtPre.forEach(b => b.classList.remove("active"));
    inp.value = st.seq; lastMsg = ""; clampPos(); recompute();
  }
  function editSel(fn: (v: number) => number): void {
    if (!st.sel) return;
    const [b, j] = st.sel;
    st.filter[b][j] = Math.round(fn(st.filter[b][j]) * 100) / 100;
    st.filterName = st.filterName.endsWith("(edited)") ? st.filterName : st.filterName + " (edited)";
    filtPre.forEach(x => x.classList.remove("active"));
    recompute();
  }
  function mutate(i: number, b: number): void {
    const ref = st.seq[i];
    if (M.BASES[b] === ref) return;
    const p0 = M.forward(st.seq, netNow()).p;
    undo.push(st.seq);
    st.seq = st.seq.slice(0, i) + M.BASES[b] + st.seq.slice(i + 1);
    const p1 = M.forward(st.seq, netNow()).p;
    lastMsg = `${ref}→${M.BASES[b]} at position ${i + 1}: P(bound) ${f3(p0)} → ${f3(p1)} (Δ ${sgn(p1 - p0)})`;
    seqPre.forEach(x => x.classList.remove("active"));
    inp.value = st.seq; recompute();
  }
  const k = (): number => st.filter[0].length;
  function clampPos(): void {
    const n = Math.max(0, st.seq.length - k());
    sPos.input.max = String(n);
    if (st.pos > n) st.pos = n;
    sPos.set(st.pos);
  }
  const netNow = (): M.Net => ({ filter: st.filter, bias: st.bias, w: st.w, b: st.b, both: st.both });
  function stop(): void { if (timer) { clearInterval(timer); timer = 0; bPlay.textContent = "▶ slide along"; } }
  function play(): void {
    if (timer) { stop(); return; }
    const n = st.seq.length - k();
    if (n < 0) return;
    if (st.pos >= n) st.pos = 0;
    bPlay.textContent = "■ stop";
    timer = window.setInterval(() => {
      st.pos++; sPos.set(st.pos); draw();
      if (st.pos >= n) stop();
    }, 160);
    draw();
  }

  let fw!: M.Forward, ismMap!: number[][];
  function recompute(): void { fw = M.forward(st.seq, netNow()); ismMap = M.ism(st.seq, netNow()); draw(); }

  // ---------------- drawing ----------------
  function draw(): void { drawSeq(); drawFilter(); drawReadout(); bUndo.disabled = undo.length === 0; }

  function drawSeq(): void {
    gSeq.replaceChildren();
    const L = st.seq.length, K = k();
    const cw = Math.min(26, (W - ml - mr) / L);
    const x = (i: number): number => ml + i * cw;
    const ch = 22, gy = 66;                 // one-hot grid top
    const ty0 = 266, tyH = 150;             // activation track
    const iy = ty0 + tyH + 110;            // ISM grid top
    const H = iy + 4 * ch + 34;
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const t = (xx: number, yy: number, s: string, attrs: Record<string, string | number> = {}) => {
      const e = svgEl("text", { x: xx, y: yy, ...attrs }, gSeq); e.textContent = s;
      if (attrs["font-size"]) e.style.fontSize = attrs["font-size"] + "px";
      return e;
    };
    // --- one-hot grid
    t(4, 20, "Input: one-hot encoded sequence", { class: "ttl" });
    t(W - mr, 20, `${L} bp · click a cell to change that base`, { class: "sm", "text-anchor": "end" });
    const fsz = Math.min(16, Math.max(13, cw * 0.8));
    for (let i = 0; i < L; i++) {
      const b = M.baseIndex(st.seq[i]);
      const lt = t(x(i) + cw / 2, gy - 8, st.seq[i], { "text-anchor": "middle", "font-weight": 700, "font-size": fsz, "font-family": "ui-monospace, monospace" });
      lt.style.fill = BASECOL[b];
      for (let r = 0; r < 4; r++) {
        const on = r === b;
        const c = svgEl("rect", {
          x: x(i) + 1, y: gy + r * ch + 1, width: cw - 2, height: ch - 2, rx: 2, class: "cell",
          fill: on ? BASECOL[r] : COL.light, stroke: on ? "none" : "#dde3e8",
        }, gSeq);
        c.addEventListener("click", () => mutate(i, r));
        const ttl = svgEl("title", {}, c); ttl.textContent = on ? `${M.BASES[r]} at ${i + 1}` : `change ${st.seq[i]}→${M.BASES[r]} at ${i + 1}`;
        if (cw >= 16) {
          const v = t(x(i) + cw / 2, gy + r * ch + ch / 2 + 5, on ? "1" : "0", { "text-anchor": "middle", "font-size": 13, "pointer-events": "none" });
          v.style.fill = on ? "#fff" : "#b3bcc5";
        }
      }
    }
    M.BASES.forEach((b, r) => {
      const e = t(ml - 10, gy + r * ch + ch / 2 + 6, b, { "text-anchor": "end", "font-weight": 700 }); e.style.fill = BASECOL[r];
    });
    // positions ruler
    for (let i = 0; i < L; i += L > 30 ? 10 : 5) t(x(i) + cw / 2, gy + 4 * ch + 18, String(i + 1), { class: "sm", "text-anchor": "middle" });
    // --- windows: best match (dashed) and the slider window (solid)
    const n = L - K + 1;
    if (n >= 1) {
      if (fw.argmax >= 0 && fw.pooled > 0) {
        svgEl("rect", {
          x: x(fw.argmax) - 1, y: gy - 2, width: K * cw + 2, height: 4 * ch + 4, fill: "none", rx: 3,
          stroke: fw.strand === "+" ? COL.fwd : COL.rev, "stroke-width": 2.5, "stroke-dasharray": "6 4", "pointer-events": "none",
        }, gSeq);
      }
      const p = Math.min(st.pos, n - 1);
      svgEl("rect", { x: x(p) - 2, y: gy - 30, width: K * cw + 4, height: 4 * ch + 32, fill: "none", rx: 4, stroke: COL.win, "stroke-width": 3, "pointer-events": "none" }, gSeq);
      const sc = fw.fwd[p];
      const lab = `window ${p + 1}–${p + K}: score ${f2(sc)}` + (fw.rev ? ` (other strand ${f2(fw.rev[p])})` : "");
      const lx = Math.min(W - mr, Math.max(4, x(p)));
      const anchor = x(p) + 360 > W ? "end" : "start";
      const tl = t(anchor === "end" ? Math.min(W - mr, x(p) + K * cw) : lx, gy + 4 * ch + 38, lab, { "font-weight": 700, "text-anchor": anchor });
      tl.style.fill = COL.win;
    } else {
      t(ml, gy + 4 * ch + 40, `The sequence (${L} bp) is shorter than the filter (${K} bp).`, { "font-weight": 700 }).style.fill = COL.sel;
    }

    // --- activation track
    t(4, ty0 - 30, "Filter score at each position (bar under the middle of its window)", { class: "ttl" });
    const [lo, hi] = M.scoreRange(st.filter);
    const R = Math.max(1, hi, -st.bias + 1, ...fw.fwd.map(Math.abs).filter(v => v <= hi)) * 1.08;
    const ylo = -0.6 * R;
    const y = scale(ylo, R, ty0 + tyH, ty0);
    svgEl("line", { x1: ml, x2: W - mr, y1: y(0), y2: y(0), stroke: COL.muted, "stroke-width": 1 }, gSeq);
    t(ml - 8, y(0) + 5, "0", { class: "sm", "text-anchor": "end" });
    t(ml - 8, y(R / 1.08) + 5, f1(R / 1.08), { class: "sm", "text-anchor": "end" });
    const thr = -st.bias;
    if (thr <= R && thr >= ylo) {
      svgEl("line", { x1: ml, x2: W - mr, y1: y(thr), y2: y(thr), stroke: COL.sel, "stroke-width": 2, "stroke-dasharray": "7 5" }, gSeq);
      // short label, on the side away from the max-pool bar
      const right = fw.argmax < 0 || x(fw.argmax + K / 2) < W - mr - 160;
      const e = t(right ? W - mr : ml + 4, y(thr) - 6, `−bias = ${f1(thr)}`,
        { "text-anchor": right ? "end" : "start", "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" });
      e.style.fill = COL.sel;
    }
    const tracks: [number[], string, number][] = fw.rev ? [[fw.fwd, COL.fwd, -0.22], [fw.rev, COL.rev, 0.22]] : [[fw.fwd, COL.fwd, 0]];
    const bw = cw * (fw.rev ? 0.42 : 0.7);
    for (const [arr, col, off] of tracks) {
      arr.forEach((s, i) => {
        const cx = x(i + K / 2) + off * cw;                // bar under the middle of its window
        const sv = Math.max(ylo, Math.min(R, s));
        const top = Math.max(sv, 0), bot = Math.min(sv, 0);
        svgEl("rect", { x: cx - bw / 2, y: y(top), width: bw, height: Math.max(0.5, y(bot) - y(top)), fill: col, "fill-opacity": 0.28 }, gSeq);
        if (s > thr) {                                      // the part that passes ReLU
          const b0 = Math.max(thr, 0);
          svgEl("rect", { x: cx - bw / 2, y: y(Math.min(R, s)), width: bw, height: Math.max(1, y(b0) - y(Math.min(R, s))), fill: col }, gSeq);
        }
        if (s < ylo) svgEl("path", { d: `M${cx - 4},${y(ylo) - 6} l4,7 l4,-7 z`, fill: col }, gSeq);
      });
    }
    // slider window bar marker
    if (n >= 1) {
      const p = Math.min(st.pos, n - 1), cx = x(p + K / 2);
      svgEl("rect", { x: cx - cw * 0.5, y: ty0 - 4, width: cw, height: tyH + 8, fill: "none", stroke: COL.win, "stroke-width": 2.5, rx: 3 }, gSeq);
    }
    // max pool mark
    if (fw.argmax >= 0) {
      const arr = fw.strand === "+" ? fw.fwd : fw.rev!;
      const s = arr[fw.argmax];
      const cx = x(fw.argmax + K / 2) + (fw.rev ? (fw.strand === "+" ? -0.22 : 0.22) * cw : 0);
      const cy = y(Math.max(ylo, Math.min(R, s)));
      svgEl("circle", { cx, cy, r: 6, fill: "#fff", stroke: COL.ink, "stroke-width": 2.5 }, gSeq);
      const txt = fw.pooled > 0 ? `max pool = ${f2(fw.pooled)} (${fw.strand} strand)` : "max pool = 0: no position passes the bias";
      const right = cx < W - 330;
      t(right ? cx + 12 : cx - 12, Math.max(ty0 - 8, cy - 10), txt, { "font-weight": 700, "text-anchor": right ? "start" : "end", stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" });
    }
    if (fw.rev) {
      const lg = t(4, ty0 + tyH + 26, "", { class: "sm" });
      lg.innerHTML = `<tspan fill="${COL.fwd}" font-weight="700">■ filter (+ strand)</tspan>   <tspan fill="${COL.rev}" font-weight="700">■ reverse-complement filter (− strand)</tspan>`;
    }
    t(4, ty0 + tyH + (fw.rev ? 48 : 26), `red line: −bias; light: score; dark: the part that passes ReLU; ▼ = below ${f1(ylo)} (lowest possible ${f1(lo)})`, { class: "sm" });

    // --- ISM
    t(4, iy - 14, "In-silico mutagenesis: change in P(bound) for each substitution", { class: "ttl" });
    let mx = 0.02;
    for (const r of ismMap) for (const v of r) mx = Math.max(mx, Math.abs(v));
    t(4, iy + 4 * ch + 22, `red = lower P(bound), green = higher (darkest ±${f2(mx)}); dot = the base in the sequence`, { class: "sm" });
    for (let i = 0; i < L; i++) {
      for (let r = 0; r < 4; r++) {
        const ref = M.BASES[r] === st.seq[i];
        const v = ismMap[r][i];
        const c = svgEl("rect", {
          x: x(i) + 1, y: iy + r * ch + 1, width: cw - 2, height: ch - 2, rx: 2, class: "cell",
          fill: ref ? "#fff" : v < 0 ? tint(COL.down, -v / mx) : tint(COL.up, v / mx),
          stroke: ref ? "#c5ccd3" : "none",
        }, gSeq);
        if (ref) svgEl("circle", { cx: x(i) + cw / 2, cy: iy + r * ch + ch / 2, r: 2.5, fill: COL.muted, "pointer-events": "none" }, gSeq);
        c.addEventListener("click", () => mutate(i, r));
        const tt = svgEl("title", {}, c);
        tt.textContent = ref ? `${st.seq[i]} at ${i + 1} (reference)` : `${st.seq[i]}→${M.BASES[r]} at ${i + 1}: Δ ${sgn(v)}`;
      }
    }
    M.BASES.forEach((b, r) => {
      const e = t(ml - 10, iy + r * ch + ch / 2 + 6, "→" + b, { "text-anchor": "end", "font-weight": 700 }); e.style.fill = BASECOL[r];
    });
  }

  function drawFilter(): void {
    gF.replaceChildren();
    const K = k(), fm = 40;
    const cw = Math.min(46, (640 - fm - 10) / K);
    const x = (j: number): number => fm + j * cw;
    const [, hi] = M.scoreRange(st.filter);
    let mxAbs = 0.5;
    for (const r of st.filter) for (const v of r) mxAbs = Math.max(mxAbs, Math.abs(v));
    const t = (xx: number, yy: number, s: string, attrs: Record<string, string | number> = {}) => {
      const e = svgEl("text", { x: xx, y: yy, ...attrs }, gF); e.textContent = s;
      if (attrs["font-size"]) e.style.fontSize = attrs["font-size"] + "px";
      return e;
    };
    t(2, 18, st.filterName, { class: "ttl" });
    t(2, 40, `width ${K} · consensus ${M.consensus(st.filter)} · best score ${f2(hi)}`, { class: "sm" });
    // logo of the positive weights (letter height ∝ weight)
    const ly = 140, lh = 86;
    let colMax = 0.5;
    for (let j = 0; j < K; j++) colMax = Math.max(colMax, st.filter.reduce((s, r) => s + Math.max(0, r[j]), 0));
    for (let j = 0; j < K; j++) {
      const order = [0, 1, 2, 3].filter(r => st.filter[r][j] > 0).sort((a, b) => st.filter[a][j] - st.filter[b][j]);
      let yb = ly;
      for (const r of order) {
        const h = lh * st.filter[r][j] / colMax;
        if (h < 1) continue;
        // capital letters are about 0.72 em tall: font size 100 -> 72 px glyph
        const e = svgEl("text", {
          x: 0, y: 0, "font-size": 100, "font-weight": 700, "font-family": "Helvetica, Arial, sans-serif", "text-anchor": "middle",
          transform: `translate(${x(j) + cw / 2},${yb}) scale(${(cw * 0.95) / 70},${h / 72})`,
        }, gF);
        e.textContent = M.BASES[r]; e.style.fill = BASECOL[r]; e.style.fontSize = "100px";
        yb -= h;
      }
    }
    svgEl("line", { x1: fm, x2: x(K), y1: ly, y2: ly, stroke: COL.muted }, gF);
    t(fm - 6, ly - lh + 8, "logo", { class: "sm", "text-anchor": "end", transform: `rotate(-90 ${fm - 22} ${ly - lh / 2})`, x: fm - 22, y: ly - lh / 2 + 5 });
    // heatmap
    const hy = 152, hh = 30;
    const used = new Set<string>();
    const n = st.seq.length - K + 1;
    if (n >= 1) {
      const p = Math.min(st.pos, n - 1);
      for (let j = 0; j < K; j++) used.add(`${M.baseIndex(st.seq[p + j])},${j}`);
    }
    for (let r = 0; r < 4; r++) {
      const e = t(fm - 8, hy + r * hh + hh / 2 + 6, M.BASES[r], { "text-anchor": "end", "font-weight": 700 }); e.style.fill = BASECOL[r];
      for (let j = 0; j < K; j++) {
        const v = st.filter[r][j];
        const c = svgEl("rect", {
          x: x(j) + 1, y: hy + r * hh + 1, width: cw - 2, height: hh - 2, rx: 2, class: "cell",
          fill: v >= 0 ? tint(COL.pos, v / mxAbs) : tint(COL.neg, -v / mxAbs),
        }, gF);
        c.addEventListener("click", () => { st.sel = [r, j]; draw(); });
        if (cw >= 30) {
          const tv = t(x(j) + cw / 2, hy + r * hh + hh / 2 + 5, f1(v), { "text-anchor": "middle", "font-size": 15, "pointer-events": "none" });
          tv.style.fill = Math.abs(v) / mxAbs > 0.55 ? "#fff" : COL.ink;
        }
        if (used.has(`${r},${j}`))
          svgEl("rect", { x: x(j) + 1.5, y: hy + r * hh + 1.5, width: cw - 3, height: hh - 3, rx: 2, fill: "none", stroke: COL.win, "stroke-width": 3, "pointer-events": "none" }, gF);
        if (st.sel && st.sel[0] === r && st.sel[1] === j)
          svgEl("rect", { x: x(j) - 1, y: hy + r * hh - 1, width: cw + 2, height: hh + 2, rx: 3, fill: "none", stroke: COL.ink, "stroke-width": 3, "stroke-dasharray": "4 3", "pointer-events": "none" }, gF);
      }
    }
    for (let j = 0; j < K; j += K > 10 ? 2 : 1) t(x(j) + cw / 2, hy + 4 * hh + 18, String(j + 1), { class: "sm", "text-anchor": "middle" });
    t(fm, hy + 4 * hh + 42, "orange outline: the weights the current window adds up · click a cell to edit it", { class: "sm" });
    // edit row
    const has = !!st.sel;
    eBtns.forEach(b => (b.disabled = !has));
    selLab.textContent = has ? `weight of ${M.BASES[st.sel![0]]} at column ${st.sel![1] + 1}: ${f2(st.filter[st.sel![0]][st.sel![1]])}` : "Click a weight in the heatmap to edit it.";
  }

  function drawReadout(): void {
    const lines: string[] = [];
    const K = k();
    if (fw.argmax < 0) lines.push(`Sequence shorter than the filter: nothing to scan.`);
    else {
      const s = (fw.strand === "+" ? fw.fwd : fw.rev!)[fw.argmax];
      lines.push(`best window: ${fw.argmax + 1}–${fw.argmax + K} (${fw.strand} strand), score ${f2(s)}`);
      lines.push(`max pool = ReLU(${f2(s)} ${st.bias < 0 ? "−" : "+"} ${f1(Math.abs(st.bias))}) = <b>${f2(fw.pooled)}</b>`);
    }
    lines.push(`z = w × max pool + b = ${f2(st.w)} × ${f2(fw.pooled)} ${st.b < 0 ? "−" : "+"} ${f2(Math.abs(st.b))} = <b>${f2(fw.z)}</b>`);
    lines.push(`<span class="cf-p">P(bound) = σ(z) = ${f3(fw.p)}</span>`);
    if (lastMsg) lines.push(`<span style="color:${COL.sel}">${lastMsg}</span>`);
    out.innerHTML = lines.join("<br>");
  }

  inp.value = st.seq;
  clampPos();
  fw = M.forward(st.seq, netNow());
  if (fw.argmax >= 0) { st.pos = fw.argmax; sPos.set(st.pos); }
  recompute();
}

autoMount(NAME, mount);
