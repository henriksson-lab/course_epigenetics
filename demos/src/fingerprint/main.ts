// Fingerprint plot and library complexity on a simulated ChIP-seq library (model.ts).
// Left: the fingerprint (bins sorted by read count; cumulative fraction of reads against the
// fraction of bins), with an input of the same depth and the diagonal for comparison.
// Right: how many positions have 1, 2, 3, 4, 5+ reads, and NRF, PBC1, PBC2, FRiP.
//
// Auto-mounts on <div data-demo="fingerprint">. data-* options (all optional):
//   data-f="30" (% reads from enriched bins)  data-bins="0.5" (% of bins enriched)
//   data-complexity="10"  data-dedup="on"  data-input="off"  data-seed="1"
import { svgEl, htmlEl, slider, checkbox, presets, scale, autoMount } from "../lib/ui";
import * as M from "./model";

const NAME = "fingerprint";
const BINS = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40];
const CPLX = [0.3, 0.5, 1, 2, 3, 5, 10, 30, 100];
const COL = { chip: "#e67e22", input: "#5b6b7a", diag: "#b8c2cc", good: "#2e8b57", mild: "#d4a017", bad: "#a8324a", bar: "#2a6f97", ink: "#1b2a3a" };
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .fp-controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(22em, 100%), 1fr)); gap: .35em 1.4em; width: 100%; }
[data-demo="${NAME}"] .fp-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .fp-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .fp-plots svg { flex: 1 1 480px; min-width: 0; max-width: 560px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 17px; }
[data-demo="${NAME}"] svg .ttl { font-size: 18px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #5b6b7a; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
.reveal [data-demo="${NAME}"] { font-size: 19px; }
@media (max-width: 760px) { [data-demo="${NAME}"] .fp-plots svg { max-width: 100%; flex-basis: 100%; } }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const nearest = (arr: number[], v: number) => arr.reduce((b, x, i) => (Math.abs(Math.log(x / v)) < Math.abs(Math.log(arr[b] / v)) ? i : b), 0);
const pct = (v: number, d = 0) => `${(100 * v).toFixed(d)}%`;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const numv = (s: string | undefined, d: number) => { const v = Number(s); return s !== undefined && Number.isFinite(v) ? v : d; };
  const seed = numv(ds.seed, 1);
  const st = {
    f: Math.min(80, Math.max(0, numv(ds.f, 30))), bi: nearest(BINS, numv(ds.bins, 0.5)), ci: nearest(CPLX, numv(ds.complexity, 10)),
    dedup: ds.dedup === "on", input: ds.input !== "off",
  };
  const genome = M.makeGenome(seed);
  const inputLib = M.simulate({ f: 0, enrichedPct: 0, complexity: 30 }, genome, seed + 500);

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "fp-controls", "", host);
  const sF = slider(ctr, "enriched reads", 0, 80, 1, st.f, v => { st.f = v; schedule(); }, v => `${v}%`);
  const sB = slider(ctr, "enriched bins", 0, BINS.length - 1, 1, st.bi, v => { st.bi = v; schedule(); }, v => `${BINS[v]}%`);
  const sC = slider(ctr, "complexity", 0, CPLX.length - 1, 1, st.ci, v => { st.ci = v; schedule(); }, v => `×${CPLX[v]}`);
  const row = htmlEl("div", "fp-row", "", host);
  checkbox(row, "remove duplicates first", st.dedup, v => { st.dedup = v; schedule(); });
  checkbox(row, "show input", st.input, v => { st.input = v; schedule(); });
  const row2 = htmlEl("div", "fp-row", "", host);
  htmlEl("span", "", "samples:", row2);
  const setAll = (f: number, bins: number, cx: number) => {
    st.f = f; st.bi = nearest(BINS, bins); st.ci = nearest(CPLX, cx);
    sF.set(st.f); sB.set(st.bi); sC.set(st.ci); schedule();
  };
  presets(row2, [
    { label: "input-like", apply: () => setAll(0, 0.5, 30) },
    { label: "failed ChIP", apply: () => setAll(3, 0.5, 10) },
    { label: "TF ChIP", apply: () => setAll(30, 0.5, 10) },
    { label: "broad mark", apply: () => setAll(50, 20, 10) },
    { label: "bottlenecked TF", apply: () => setAll(30, 0.5, 0.5) },
  ]);
  htmlEl("p", "qb-note",
    `Simulated library: ${M.READS / 1000}k reads in ${M.NBINS / 1000}k bins of 1 kb (8% of bins unmappable). The input has the same depth and a complex library.`, host);

  const plots = htmlEl("div", "fp-plots", "", host);
  const W = 560, H = 470;
  const svgL = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Fingerprint plot" }, plots);
  const svgR = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Reads per position and library complexity" }, plots);

  let pending = 0;
  function schedule(): void {
    if (pending) return;
    pending = window.setTimeout(() => { pending = 0; draw(); }, 20);
  }

  const L = 66, R = W - 16, T = 84, B = H - 62;
  const x = scale(0, 1, L, R), y = scale(0, 1, B, T);

  function draw(): void {
    const lib = M.simulate({ f: st.f / 100, enrichedPct: BINS[st.bi], complexity: CPLX[st.ci] }, genome, seed);
    const counts = st.dedup ? lib.countsDedup : lib.countsAll;
    const inCounts = st.dedup ? inputLib.countsDedup : inputLib.countsAll;
    const fp = M.fingerprint(counts), fpIn = M.fingerprint(inCounts);
    const top10 = M.topShare(counts, 0.1), top1 = M.topShare(counts, 0.01);
    const fripS = M.frip(counts, lib.enriched), fripIn = M.frip(inCounts, lib.enriched);

    // ---- fingerprint ----
    svgL.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgL).textContent = "Fingerprint (1-kb bins)";
    const sub = svgEl("text", { x: 8, y: 50 }, svgL);
    const a1 = svgEl("tspan", {}, sub); a1.textContent = "top 10% of bins: ";
    const a2 = svgEl("tspan", { "font-weight": 600 }, sub); a2.style.fill = COL.chip; a2.textContent = pct(top10);
    const a3 = svgEl("tspan", {}, sub); a3.textContent = " of reads; top 1%: ";
    const a4 = svgEl("tspan", { "font-weight": 600 }, sub); a4.style.fill = COL.chip; a4.textContent = pct(top1);
    for (const v of [0, 0.25, 0.5, 0.75, 1]) {
      svgEl("line", { x1: L, x2: R, y1: y(v), y2: y(v), class: v === 0 ? "axis" : "grid" }, svgL);
      svgEl("line", { x1: x(v), x2: x(v), y1: T, y2: B, class: v === 0 ? "axis" : "grid" }, svgL);
      svgEl("text", { x: L - 6, y: y(v) + 5, "text-anchor": "end", class: "sm" }, svgL).textContent = String(v);
      svgEl("text", { x: x(v), y: B + 22, "text-anchor": "middle", class: "sm" }, svgL).textContent = String(v);
    }
    svgEl("text", { x: (L + R) / 2, y: H - 16, "text-anchor": "middle", class: "sm" }, svgL).textContent = "fraction of bins (ranked by read count)";
    const yl = svgEl("text", { x: 16, y: (T + B) / 2, "text-anchor": "middle", class: "sm", transform: `rotate(-90 16 ${(T + B) / 2})` }, svgL);
    yl.textContent = "cumulative fraction of reads";
    svgEl("line", { x1: x(0), y1: y(0), x2: x(1), y2: y(1), stroke: COL.diag, "stroke-width": 2, "stroke-dasharray": "6 5" }, svgL);
    const path = (c: { x: number[]; y: number[] }) => c.x.map((v, i) => `${i ? "L" : "M"}${x(v).toFixed(1)},${y(c.y[i]).toFixed(1)}`).join("");
    if (st.input) svgEl("path", { d: path(fpIn), fill: "none", stroke: COL.input, "stroke-width": 2.5, "stroke-dasharray": "8 4" }, svgL);
    svgEl("path", { d: path(fp), fill: "none", stroke: COL.chip, "stroke-width": 3.2 }, svgL);
    // top-10% marker
    const y90 = fp.y[Math.round(0.9 * (fp.x.length - 1))];
    svgEl("line", { x1: x(0.9), x2: x(0.9), y1: y(y90), y2: y(1), stroke: COL.chip, "stroke-width": 2, "stroke-dasharray": "3 3" }, svgL);
    // legend
    const lg = (yy: number, c: string, dash: string, lab: string) => {
      svgEl("line", { x1: L + 12, x2: L + 44, y1: yy, y2: yy, stroke: c, "stroke-width": 3, "stroke-dasharray": dash }, svgL);
      svgEl("text", { x: L + 52, y: yy + 5 }, svgL).textContent = lab;
    };
    lg(T + 16, COL.chip, "", "sample" + (st.dedup ? " (duplicates removed)" : ""));
    if (st.input) lg(T + 40, COL.input, "8 4", "input");
    lg(T + (st.input ? 64 : 40), COL.diag, "6 5", "perfectly uniform");

    // ---- complexity ----
    svgR.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgR).textContent = "Reads per position (before duplicate removal)";
    const sub2 = svgEl("text", { x: 8, y: 50, class: "sm" }, svgR);
    sub2.textContent = `${lib.reads.toLocaleString("en-US").replace(/,/g, " ")} reads at ${lib.distinct.toLocaleString("en-US").replace(/,/g, " ")} distinct positions (5' end + strand)`;
    const bx = 70, bR = W - 20, bT = 76, bB = 250;
    const tot = lib.distinct;
    const fr = lib.perPos.map(v => v / tot);
    const yb = scale(0, 1, bB, bT);
    for (const v of [0, 0.5, 1]) {
      svgEl("line", { x1: bx, x2: bR, y1: yb(v), y2: yb(v), class: v === 0 ? "axis" : "grid" }, svgR);
      svgEl("text", { x: bx - 6, y: yb(v) + 5, "text-anchor": "end", class: "sm" }, svgR).textContent = pct(v);
    }
    const bw = (bR - bx) / 5;
    ["1", "2", "3", "4", "5+"].forEach((lab, i) => {
      const h = fr[i];
      svgEl("rect", { x: bx + i * bw + 8, width: bw - 16, y: yb(h), height: bB - yb(h), fill: i === 0 ? COL.bar : COL.bad, "fill-opacity": 0.8 }, svgR);
      svgEl("text", { x: bx + (i + 0.5) * bw, y: bB + 20, "text-anchor": "middle" }, svgR).textContent = lab;
      svgEl("text", { x: bx + (i + 0.5) * bw, y: yb(h) - 6, "text-anchor": "middle", class: "sm" }, svgR).textContent = pct(h, h < 0.1 ? 1 : 0);
    });
    svgEl("text", { x: bx + 2.5 * bw, y: bB + 40, "text-anchor": "middle", class: "sm" }, svgR).textContent = "reads at the same position";
    const yl2 = svgEl("text", { x: 16, y: (bT + bB) / 2, "text-anchor": "middle", class: "sm", transform: `rotate(-90 16 ${(bT + bB) / 2})` }, svgR);
    yl2.textContent = "positions";

    // metrics table
    const cls = (c: string) => (c === "ideal" || c === "none" ? COL.good : c === "compliant" || c === "mild" ? COL.mild : COL.bad);
    const rowsY = 320;
    const rowT = (i: number, name: string, val: string, note: string, col: string) => {
      const yy = rowsY + i * 27;
      svgEl("text", { x: 20, y: yy, "font-weight": 600 }, svgR).textContent = name;
      const v = svgEl("text", { x: 190, y: yy, "text-anchor": "end", "font-weight": 600 }, svgR);
      v.style.fill = col; v.textContent = val;
      svgEl("text", { x: 206, y: yy, class: "sm" }, svgR).textContent = note;
    };
    rowT(0, "NRF", lib.nrf.toFixed(2), `${M.nrfClass(lib.nrf)}  (distinct / reads)`, cls(M.nrfClass(lib.nrf)));
    rowT(1, "PBC1", lib.pbc1.toFixed(2), `${M.pbc1Class(lib.pbc1)} bottlenecking  (M1 / distinct)`, cls(M.pbc1Class(lib.pbc1)));
    rowT(2, "PBC2", Number.isFinite(lib.pbc2) ? lib.pbc2.toFixed(1) : "∞", `${M.pbc2Class(lib.pbc2)} bottlenecking  (M1 / M2)`, cls(M.pbc2Class(lib.pbc2)));
    rowT(3, "duplicates", pct(1 - lib.nrf), "of reads (markdup would remove)", COL.ink);
    rowT(4, "FRiP", pct(fripS, 1), `reads in enriched bins${st.dedup ? " (deduplicated)" : ""}; input ${pct(fripIn, 1)}`, COL.ink);
  }

  draw();
}

autoMount(NAME, mount);
