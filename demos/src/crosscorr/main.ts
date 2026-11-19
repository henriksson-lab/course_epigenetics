// Strand cross-correlation on simulated single-end ChIP-seq (model.ts).
// Left: + strand and − strand read 5' ends summed around all binding sites (or around the
// left edges of all repeat stretches); the − strand is drawn shifted left by k.
// Right: CC(k) for k = 0..600, with the current shift, the fragment-length peak, the
// read-length ("phantom") peak and the minimum; NSC and RSC as in ENCODE.
//
// Auto-mounts on <div data-demo="crosscorr">. data-* options (all optional):
//   data-signal="15" (% of fragments from sites)  data-frag="200"  data-readlen="36"
//   data-shift="0"  data-align="sites|repeats"  data-seed="1"
import { svgEl, htmlEl, slider, presets, scale, autoMount } from "../lib/ui";
import * as M from "./model";

const NAME = "crosscorr";
const COL = { plus: "#2a6f97", minus: "#c0392b", cc: "#1b2a3a", frag: "#2e8b57", read: "#e67e22", min: "#5b6b7a", k: "#6b4fa0" };
const READS = 200000;
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .xc-controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(22em, 100%), 1fr)); gap: .35em 1.4em; width: 100%; }
[data-demo="${NAME}"] .xc-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .xc-plots { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; }
[data-demo="${NAME}"] .xc-plots svg { flex: 1 1 480px; min-width: 0; max-width: 560px; height: auto; }
[data-demo="${NAME}"] .xc-shift .qb-slider-name { font-weight: 600; color: ${COL.k}; }
[data-demo="${NAME}"] svg text { font-size: 17px; }
[data-demo="${NAME}"] svg .ttl { font-size: 18px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #5b6b7a; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
[data-demo="${NAME}"] .xc-busy { opacity: .55; transition: opacity .1s; }
.reveal [data-demo="${NAME}"] { font-size: 19px; }
@media (max-width: 760px) { [data-demo="${NAME}"] .xc-plots svg { max-width: 100%; flex-basis: 100%; } }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

function rich(svg: SVGElement, x: number, y: number, parts: [string, string][], cls = ""): void {
  const t = svgEl("text", { x, y, class: cls }, svg);
  for (const [s, c] of parts) {
    const sp = svgEl("tspan", {}, t);
    sp.textContent = s;
    if (c) { sp.style.fill = c; sp.style.fontWeight = "600"; }
  }
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, d: number, lo: number, hi: number) => {
    const v = Number(s); return s !== undefined && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
  };
  const st = {
    signal: num(ds.signal, 15, 0, 40), frag: num(ds.frag, 200, 80, 300), readLen: num(ds.readlen, 36, 25, 100),
    k: num(ds.shift, 0, 0, M.KMAX), align: ds.align === "repeats" ? "repeats" : "sites", seed: num(ds.seed, 1, 1, 1e6),
  };
  const genome = M.makeGenome(st.seed);
  const repeatAnchors = genome.repeatStarts.subarray(0, 20000);

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "xc-controls", "", host);
  const sk = slider(ctr, "shift − strand, k", 0, M.KMAX, 1, st.k, v => { st.k = v; drawShift(); });
  sk.input.closest("label")!.classList.add("xc-shift");
  const sSig = slider(ctr, "reads from sites", 0, 40, 1, st.signal, v => { st.signal = v; schedule(); }, v => `${v}%`);
  const sFrag = slider(ctr, "fragment length", 80, 300, 5, st.frag, v => { st.frag = v; schedule(); });
  const sRead = slider(ctr, "read length L", 25, 100, 1, st.readLen, v => { st.readLen = v; schedule(); });
  const row = htmlEl("div", "xc-row", "", host);
  htmlEl("span", "", "left panel, aligned on:", row);
  const alignBtns = presets(row, [
    { label: "binding sites", apply: () => { st.align = "sites"; drawLeft(); drawShift(); } },
    { label: "repeat edges", apply: () => { st.align = "repeats"; drawLeft(); drawShift(); } },
  ]);
  alignBtns[st.align === "sites" ? 0 : 1].classList.add("active");
  const row2 = htmlEl("div", "xc-row", "", host);
  htmlEl("span", "", "samples:", row2);
  const setAll = (sig: number, frag: number, L: number) => {
    st.signal = sig; st.frag = frag; st.readLen = L;
    sSig.set(sig); sFrag.set(frag); sRead.set(L);
    schedule();
  };
  presets(row2, [
    { label: "good TF ChIP", apply: () => setAll(15, 200, 36) },
    { label: "weak ChIP", apply: () => setAll(4, 200, 36) },
    { label: "input", apply: () => setAll(0, 200, 36) },
    { label: "100-bp reads", apply: () => setAll(15, 200, 100) },
    { label: "short fragments", apply: () => setAll(15, 95, 36) },
  ]);
  htmlEl("p", "qb-note",
    `Simulated single-end reads: ${READS / 1000}k reads on a 1-Mb genome with ${M.NSITES} binding sites; short repeats make about 19% of start positions unmappable for 36-bp reads (9% for 100-bp reads).`, host);

  const plots = htmlEl("div", "xc-plots", "", host);
  const W = 560, H = 430;
  const svgL = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Read 5' ends on the two strands around binding sites" }, plots);
  const svgR = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Strand cross-correlation curve" }, plots);
  const gL = svgEl("g", {}, svgL), gLs = svgEl("g", {}, svgL);
  const gR = svgEl("g", {}, svgR), gRs = svgEl("g", {}, svgR);

  let rd: M.Reads, cc: Float64Array, q: M.QC;
  let aggP: Float64Array, aggM: Float64Array, view = { half: 0, bin: 0, extra: 0 };

  function compute(): void {
    rd = M.simulate({ reads: READS, signal: st.signal / 100, frag: st.frag, readLen: st.readLen }, genome, st.seed);
    cc = M.crossCorr(rd);
    q = M.qc(cc, st.readLen);
  }

  let pending = 0;
  function schedule(): void {
    if (pending) return;
    plots.classList.add("xc-busy");
    pending = window.setTimeout(() => { pending = 0; compute(); drawLeft(); drawRight(); drawShift(); plots.classList.remove("xc-busy"); }, 30);
  }

  // ---------------- left panel ----------------
  const lL = 56, lR = W - 14, lT = 92, lB = H - 60;
  function drawLeft(): void {
    gL.replaceChildren();
    const sites = st.align === "sites";
    view = sites ? { half: 400, bin: 10, extra: M.KMAX } : { half: 60, bin: 1, extra: M.KMAX };
    // + strand: [-half, half); − strand: [-half, half + extra) so it can be shifted left by up to KMAX
    const anchors = sites ? genome.sites : repeatAnchors;
    aggP = M.aggregate(rd.fwd, anchors, view.half, view.bin);
    // − strand over [-half, half + extra) around each anchor: a window of half-width
    // half + extra/2 centred extra/2 to the right
    aggM = M.aggregate(rd.rev, Array.from(anchors, a => a + view.extra / 2), view.half + view.extra / 2, view.bin);
    const t = svgEl("text", { x: 8, y: 24, class: "ttl" }, gL);
    t.textContent = sites ? `Read 5' ends around ${M.NSITES} binding sites` : `Read 5' ends at ${anchors.length / 1000}k repeat left edges`;
    const sub = svgEl("text", { x: 8, y: 48, class: "sm" }, gL);
    sub.textContent = sites ? "summed over all sites, 10-bp bins" : "reads that start in a repeat are often unmappable";
  }

  function drawShift(): void { drawShiftL(); drawShiftR(); }

  function drawShiftL(): void {
    gLs.replaceChildren();
    if (!aggP) return;
    const { half, bin } = view;
    const x = scale(-half, half, lL, lR);
    // − strand bins cover [-half, half + extra); shifted: position − k
    const nP = aggP.length;
    let ymax = 1;
    for (let i = 0; i < nP; i++) ymax = Math.max(ymax, aggP[i]);
    for (let i = 0; i < aggM.length; i++) ymax = Math.max(ymax, aggM[i]);
    const y = scale(0, ymax * 1.08, lB, lT);
    const k = st.k;
    // axis
    svgEl("line", { x1: lL, x2: lR, y1: lB, y2: lB, class: "axis" }, gLs);
    svgEl("line", { x1: x(0), x2: x(0), y1: lT - 6, y2: lB, class: "grid", "stroke-dasharray": "3 3" }, gLs);
    const ticks = st.align === "sites" ? [-400, -200, 0, 200, 400] : [-60, -30, 0, 30, 60];
    for (const v of ticks) {
      svgEl("line", { x1: x(v), x2: x(v), y1: lB, y2: lB + 5, class: "axis" }, gLs);
      svgEl("text", { x: x(v), y: lB + 22, "text-anchor": "middle", class: "sm" }, gLs).textContent = v < 0 ? `−${-v}` : String(v);
    }
    svgEl("text", { x: (lL + lR) / 2, y: H - 14, "text-anchor": "middle", class: "sm" }, gLs).textContent =
      st.align === "sites" ? "position relative to the site (bp)" : "position relative to the repeat's left edge (bp)";
    const yl = svgEl("text", { x: 16, y: (lT + lB) / 2, "text-anchor": "middle", class: "sm", transform: `rotate(-90 16 ${(lT + lB) / 2})` }, gLs);
    yl.textContent = "read 5' ends";
    // + strand
    let dP = "";
    for (let i = 0; i < nP; i++) {
      const a = -half + i * bin;
      dP += `${i ? "L" : "M"}${x(a).toFixed(1)},${y(aggP[i]).toFixed(1)}L${x(a + bin).toFixed(1)},${y(aggP[i]).toFixed(1)}`;
    }
    svgEl("path", { d: dP, fill: "none", stroke: COL.plus, "stroke-width": 2.4 }, gLs);
    // − strand, shifted left by k
    let dM = "", started = false;
    for (let i = 0; i < aggM.length; i++) {
      const a = -half + i * bin - k;
      if (a < -half || a + bin > half) continue;
      dM += `${started ? "L" : "M"}${x(a).toFixed(1)},${y(aggM[i]).toFixed(1)}L${x(a + bin).toFixed(1)},${y(aggM[i]).toFixed(1)}`;
      started = true;
    }
    svgEl("path", { d: dM, fill: "none", stroke: COL.minus, "stroke-width": 2.4, "stroke-opacity": 0.9 }, gLs);
    rich(gLs, lL, 74, [["+ strand", COL.plus], ["   ", ""], [`− strand, shifted left by k = ${k}`, COL.minus]]);
  }

  // ---------------- right panel ----------------
  const rL = 70, rR = W - 14, rT = 92, rB = H - 60;
  const xk = scale(0, M.KMAX, rL, rR);
  let yc = scale(0, 1, rB, rT);
  function drawRight(): void {
    gR.replaceChildren();
    let lo = Infinity, hi = -Infinity;
    for (const v of cc) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const pad = (hi - lo) * 0.08 || 0.01;
    yc = scale(Math.max(0, lo - pad), hi + pad, rB, rT);
    const t = svgEl("text", { x: 8, y: 24, class: "ttl" }, gR);
    t.textContent = "Strand cross-correlation CC(k)";
    const pass = (ok: boolean) => (ok ? COL.frag : COL.minus);
    rich(gR, 8, 50, [["NSC = ", ""], [q.nsc.toFixed(2), pass(q.nsc >= 1.05)], ["   RSC = ", ""], [q.rsc.toFixed(2), pass(q.rsc >= 0.8)],
      [`   fragment ≈ ${q.frag} bp`, ""]]);
    rich(gR, 8, 74, [["ENCODE flags NSC < 1.05, RSC < 0.8", ""]], "sm");
    // grid
    const step = niceStep((hi - lo + 2 * pad) / 4);
    for (let v = Math.ceil(Math.max(0, lo - pad) / step) * step; v <= hi + pad; v += step) {
      svgEl("line", { x1: rL, x2: rR, y1: yc(v), y2: yc(v), class: "grid" }, gR);
      svgEl("text", { x: rL - 6, y: yc(v) + 5, "text-anchor": "end", class: "sm" }, gR).textContent = v.toFixed(step < 0.01 ? 3 : 2);
    }
    svgEl("line", { x1: rL, x2: rR, y1: rB, y2: rB, class: "axis" }, gR);
    for (const v of [0, 100, 200, 300, 400, 500, 600]) {
      svgEl("line", { x1: xk(v), x2: xk(v), y1: rB, y2: rB + 5, class: "axis" }, gR);
      svgEl("text", { x: xk(v), y: rB + 22, "text-anchor": "middle", class: "sm" }, gR).textContent = String(v);
    }
    svgEl("text", { x: (rL + rR) / 2, y: H - 14, "text-anchor": "middle", class: "sm" }, gR).textContent = "shift k of the − strand (bp)";
    // min line
    svgEl("line", { x1: rL, x2: rR, y1: yc(q.ccMin), y2: yc(q.ccMin), stroke: COL.min, "stroke-width": 1.5, "stroke-dasharray": "6 4" }, gR);
    svgEl("text", { x: rR - 2, y: yc(q.ccMin) + 18, "text-anchor": "end", class: "sm" }, gR).textContent = `min ${q.ccMin.toFixed(3)}`;
    // curve
    let d = "";
    for (let k = 0; k < cc.length; k++) d += `${k ? "L" : "M"}${xk(k).toFixed(1)},${yc(cc[k]).toFixed(1)}`;
    svgEl("path", { d, fill: "none", stroke: COL.cc, "stroke-width": 2 }, gR);
    // read-length and fragment markers
    const mark = (k: number, v: number, col: string, lab: string, dy: number) => {
      svgEl("circle", { cx: xk(k), cy: yc(v), r: 6, fill: "none", stroke: col, "stroke-width": 2.5 }, gR);
      const right = xk(k) > rL + 0.62 * (rR - rL);
      const tx = svgEl("text", { x: xk(k) + (right ? -9 : 9), y: yc(v) + dy, "font-weight": 600, "text-anchor": right ? "end" : "start" }, gR);
      tx.style.fill = col; tx.textContent = lab;
    };
    mark(st.readLen, q.ccRead, COL.read, `L = ${st.readLen}: ${q.ccRead.toFixed(3)}`, q.frag > st.readLen && Math.abs(yc(q.ccRead) - yc(q.ccFrag)) < 24 ? 24 : -9);
    mark(q.frag, q.ccFrag, COL.frag, `fragment ${q.frag}: ${q.ccFrag.toFixed(3)}`, -9);
  }

  function drawShiftR(): void {
    gRs.replaceChildren();
    if (!cc) return;
    const k = st.k;
    svgEl("line", { x1: xk(k), x2: xk(k), y1: rT - 4, y2: rB, stroke: COL.k, "stroke-width": 2 }, gRs);
    svgEl("circle", { cx: xk(k), cy: yc(cc[k]), r: 5, fill: COL.k }, gRs);
    const tx = svgEl("text", { x: rR, y: 74, "text-anchor": "end", "font-weight": 600 }, gRs);
    tx.style.fill = COL.k; tx.textContent = `CC(${k}) = ${cc[k].toFixed(3)}`;
  }

  compute(); drawLeft(); drawRight(); drawShift();
}

function niceStep(x: number): number {
  const p = 10 ** Math.floor(Math.log10(x));
  const m = x / p;
  return (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p;
}

autoMount(NAME, mount);
