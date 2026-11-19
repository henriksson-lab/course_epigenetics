// Is an overlap more than chance? The answer depends on the null. Two simulated interval sets on
// a toy chromosome with a gap (unassembled arm), a centromere, blacklisted repeats and two
// GC-rich, gene-dense regions where both sets prefer to be. Shuffle B ×1/×10/×100 to build the
// null histogram of overlap counts; the red line is the observed count. Toggles restrict the
// shuffle: keep out of gaps/blacklist (bedtools shuffle -excl), and keep each interval in its own
// GC class (a matched background). Enrichment = observed / null mean, p = (k + 1) / (N + 1).
// Model: model.ts.
//
// Auto-mounts on <div data-demo="shufflenull">. data-* options (all optional):
//   data-scenario="assoc|none"   data-excl="on"   data-gc="on"   data-seed="3"
import { svgEl, htmlEl, button, presets, checkbox, scale, autoMount } from "../lib/ui";
import { makeData, overlapCount, shuffler, stats, Data, Iv, Scenario, NullOpts,
  L, GAPS, CENTRO, BLACK, GC, EXCL, NA, NB, inAny } from "./model";

const NAME = "shufflenull";
const COL = { A: "#2a6f97", B: "#e67e22", obs: "#c0392b", null: "#9aa7b3", gc: "#2e8b57", text: "#1b2a3a", muted: "#5b6b7a" };
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .sn-row { display: flex; flex-wrap: wrap; gap: .35em 1em; align-items: center; }
[data-demo="${NAME}"] .sn-row .qb-presets { display: inline-flex; }
[data-demo="${NAME}"] .qb-readout.sn-readout { font-family: inherit; width: 100%; box-sizing: border-box; font-size: .9em; line-height: 1.5; }
[data-demo="${NAME}"] .sn-wrap { width: 100%; max-width: 900px; overflow-x: auto; }
[data-demo="${NAME}"] svg { width: 100%; min-width: 600px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 16.5px; }
[data-demo="${NAME}"] svg .ttl { font-size: 17.5px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #5b6b7a; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
@media (max-width: 760px) { [data-demo="${NAME}"] svg text { font-size: 20px; } [data-demo="${NAME}"] svg .ttl { font-size: 21px; } [data-demo="${NAME}"] svg .sm { font-size: 19.5px; } }
.reveal [data-demo="${NAME}"] { font-size: 18px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

function niceStep(x: number): number {
  const p = 10 ** Math.floor(Math.log10(x));
  const m = x / p;
  return (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p;
}

interface NullRun { counts: number[]; next: () => { b: Iv[]; count: number }; last: Iv[] | null }

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const seed = Number(ds.seed ?? 3) || 3;
  const st = {
    sc: (ds.scenario === "none" ? "none" : "assoc") as Scenario,
    excl: ds.excl === "on", gc: ds.gc === "on",
  };
  const data: Record<Scenario, Data> = { assoc: makeData("assoc", seed), none: makeData("none", seed) };
  const runs = new Map<string, NullRun>();
  const key = () => `${st.sc}|${st.excl}|${st.gc}`;
  const run = (): NullRun => {
    const k = key();
    let r = runs.get(k);
    if (!r) {
      const opt: NullOpts = { excl: st.excl, gc: st.gc };
      const s = 100 + (st.sc === "assoc" ? 0 : 1) + 2 * Number(st.excl) + 4 * Number(st.gc);
      r = { counts: [], next: shuffler(data[st.sc], opt, s), last: null };
      runs.set(k, r);
    }
    return r;
  };

  // ---------------- controls ----------------
  const r1 = htmlEl("div", "sn-row", "", host);
  htmlEl("span", "", "the data:", r1);
  const scItems: [Scenario, string][] = [["assoc", "B is placed near A (real association)"], ["none", "B only shares the GC preference"]];
  const scBtns = presets(r1, scItems.map(([s, l]) => ({ label: l, apply: () => { st.sc = s; draw(); } })));
  const r2 = htmlEl("div", "sn-row", "", host);
  htmlEl("span", "", "null:", r2);
  checkbox(r2, "keep out of gaps, centromere and blacklist", st.excl, v => { st.excl = v; draw(); });
  checkbox(r2, "match GC region", st.gc, v => { st.gc = v; draw(); });
  const r3 = htmlEl("div", "sn-row", "", host);
  const shuf = (n: number) => { const r = run(); for (let i = 0; i < n; i++) { const o = r.next(); r.counts.push(o.count); r.last = o.b; } draw(); };
  button(r3, "shuffle B ×1", () => shuf(1));
  button(r3, "×10", () => shuf(10));
  button(r3, "×100", () => shuf(100));
  button(r3, "clear shuffles", () => { runs.delete(key()); draw(); });
  const syncSc = () => scBtns.forEach((b, i) => b.classList.toggle("active", scItems[i][0] === st.sc));
  scBtns.forEach(b => b.addEventListener("click", syncSc));
  syncSc();
  const readout = htmlEl("div", "qb-readout sn-readout", "", host);

  // ---------------- drawing ----------------
  const W = 800, X0 = 160, X1 = W - 16;
  const yAnn = 72, yA = 112, yB = 140, yS = 168, rowH = 18;
  const hT = 248, hB = 438, H = hB + 50;
  htmlEl("p", "qb-note", "Hatched: gaps (runs of N) · grey ticks: blacklisted repeats · dark bar: centromere · red bars in the histogram: shuffles at least as large as the observed count.", host);
  const wrap = htmlEl("div", "sn-wrap", "", host);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
    "aria-label": "Two interval sets on a toy chromosome and the null distribution of their overlap from shuffles" }, wrap);
  htmlEl("p", "qb-note",
    `Simulated data on a toy 10 Mb chromosome: ${NA} A intervals (5–15 kb), ${NB} B intervals (10–30 kb). ` +
    "Neither set has intervals in the gaps or the blacklist; both prefer the GC-rich regions. " +
    "Each setting keeps its own list of shuffles (fixed seeds).", host);
  const sx = scale(0, L, X0, X1);
  const layer = svgEl("g", {}, svg);
  const defs = svgEl("defs", {}, svg);
  const hid = `sn-hatch-${Math.random().toString(36).slice(2, 8)}`;
  const pat = svgEl("pattern", { id: hid, width: 7, height: 7, patternUnits: "userSpaceOnUse", patternTransform: "rotate(45)" }, defs);
  svgEl("rect", { width: 7, height: 7, fill: "#fff" }, pat);
  svgEl("line", { x1: 0, y1: 0, x2: 0, y2: 7, stroke: "#5b6b7a", "stroke-width": 3 }, pat);

  const ivRow = (ivs: Iv[], y: number, col: string, op: (i: number) => number, outline?: (i: number) => boolean) => {
    const g = svgEl("g", {}, layer);
    ivs.forEach((iv, i) => {
      const x = sx(iv[0]), w = Math.max(2.2, sx(iv[1]) - sx(iv[0]));
      const attrs: Record<string, string | number> = { x, y, width: w, height: rowH, fill: col, "fill-opacity": op(i) };
      if (outline && outline(i)) Object.assign(attrs, { stroke: COL.obs, "stroke-width": 1.6 });
      svgEl("rect", attrs, g);
    });
  };

  function draw(): void {
    const d = data[st.sc];
    const r = run();
    const obs = overlapCount(d.A, d.B);
    const s = stats(obs, r.counts);
    layer.replaceChildren();

    // ---- chromosome ----
    svgEl("text", { x: 8, y: 20, class: "ttl" }, layer).textContent = "toy chromosome (10 Mb)";
    for (const g of GC) {
      svgEl("rect", { x: sx(g[0]), y: yAnn - 24, width: sx(g[1]) - sx(g[0]), height: yS + rowH + 4 - (yAnn - 24), fill: COL.gc, "fill-opacity": 0.1 }, layer);
      const t = svgEl("text", { x: (sx(g[0]) + sx(g[1])) / 2, y: yAnn - 30, "text-anchor": "middle", class: "sm" }, layer);
      t.textContent = "GC-rich, gene-dense"; t.style.fill = COL.gc; t.style.fontWeight = "600";
    }
    svgEl("rect", { x: X0, y: yAnn, width: X1 - X0, height: 16, rx: 8, fill: "#eaf1f7", stroke: COL.text, "stroke-width": 1 }, layer);
    for (const g of GAPS) svgEl("rect", { x: sx(g[0]), y: yAnn, width: Math.max(3, sx(g[1]) - sx(g[0])), height: 16, fill: `url(#${hid})` }, layer);
    for (const g of BLACK) svgEl("rect", { x: sx(g[0]), y: yAnn - 3, width: Math.max(3, sx(g[1]) - sx(g[0])), height: 22, fill: COL.muted }, layer);
    svgEl("rect", { x: sx(CENTRO[0]), y: yAnn + 3, width: sx(CENTRO[1]) - sx(CENTRO[0]), height: 10, fill: COL.text }, layer);
    const lab = (x: number, txt: string, anchor = "middle") => {
      const t = svgEl("text", { x, y: yAnn - 8, "text-anchor": anchor, class: "sm" }, layer); t.textContent = txt;
    };
    lab(sx(GAPS[0][0]), "gap (N)", "start");
    svgEl("text", { x: X0 - 10, y: yAnn + 13, "text-anchor": "end", class: "sm" }, layer).textContent = "regions";

    // ---- interval rows ----
    const hitA = d.A.map(a => d.B.some(b => a[0] < b[1] && b[0] < a[1]));
    ivRow(d.A, yA, COL.A, i => (hitA[i] ? 1 : 0.35));
    ivRow(d.B, yB, COL.B, () => 0.9);
    const rowLab = (y: number, txt: string, col: string) => {
      const t = svgEl("text", { x: X0 - 10, y: y + 14, "text-anchor": "end" }, layer);
      t.textContent = txt; t.style.fill = col; t.style.fontWeight = "600";
    };
    rowLab(yA, "A", COL.A); rowLab(yB, "B", COL.B); rowLab(yS, "shuffled B", COL.B);
    let lastNote = "";
    if (r.last) {
      const bad = r.last.map(iv => inAny(iv, EXCL));
      ivRow(r.last, yS, COL.B, () => 0.5, i => bad[i]);
      const nBad = bad.filter(x => x).length;
      lastNote = `last shuffle: ${overlapCount(d.A, r.last)} A intervals overlap the shuffled B` +
        (nBad ? `; ${nBad} shuffled B intervals (red outline) fell in a gap or the blacklist, where no peak can be` : "");
    } else {
      svgEl("text", { x: X0, y: yS + 14, class: "sm" }, layer).textContent = "press “shuffle B” to place B at random";
    }

    // ---- histogram ----
    const lo = 10, hi = 110;
    const hx = scale(lo - 0.5, hi + 0.5, X0, X1);
    const cnt = new Map<number, number>();
    for (const c of r.counts) cnt.set(c, (cnt.get(c) ?? 0) + 1);
    const ymax = Math.max(1, ...cnt.values());
    const ystep = niceStep(Math.max(1, ymax / 4));
    const ytop = Math.max(ystep, Math.ceil(ymax / ystep) * ystep);
    const hy = scale(0, ytop, hB, hT);
    svgEl("text", { x: 8, y: hT - 22, class: "ttl" }, layer).textContent =
      `null distribution: ${r.counts.length} shuffle${r.counts.length === 1 ? "" : "s"} of B`;
    for (let v = 0; v <= ytop + 1e-9; v += ystep) {
      svgEl("line", { x1: X0, x2: X1, y1: hy(v), y2: hy(v), class: v === 0 ? "axis" : "grid" }, layer);
      svgEl("text", { x: X0 - 8, y: hy(v) + 5, "text-anchor": "end", class: "sm" }, layer).textContent = String(+v.toFixed(3));
    }
    const yl = svgEl("text", { x: 40, y: (hT + hB) / 2, "text-anchor": "middle", transform: `rotate(-90 40 ${(hT + hB) / 2})`, class: "sm" }, layer);
    yl.textContent = "shuffles";
    const bw = hx(1) - hx(0);
    for (const [c, n] of cnt) {
      const cc = Math.max(lo, Math.min(hi, c));
      svgEl("rect", { x: hx(cc) - bw / 2 + 0.5, y: hy(n), width: Math.max(1, bw - 1), height: hB - hy(n), fill: c >= obs ? COL.obs : COL.null, "fill-opacity": c >= obs ? 0.6 : 0.85 }, layer);
    }
    if (r.counts.length) {
      const mx = hx(Math.max(lo, Math.min(hi, s.mean)));
      svgEl("line", { x1: mx, x2: mx, y1: hT - 4, y2: hB, stroke: COL.text, "stroke-width": 1.5, "stroke-dasharray": "4 4" }, layer);
      svgEl("text", { x: mx, y: hT - 6, "text-anchor": "middle", class: "sm" }, layer).textContent = `null mean ${s.mean.toFixed(1)}`;
    }
    const ox = hx(obs);
    svgEl("line", { x1: ox, x2: ox, y1: hT - 4, y2: hB, stroke: COL.obs, "stroke-width": 3 }, layer);
    const right = !r.counts.length || obs >= s.mean;
    const ot = svgEl("text", { x: right ? ox + 6 : ox - 6, y: hT + 14, "text-anchor": right ? "start" : "end" }, layer);
    ot.textContent = `observed ${obs}`; ot.style.fill = COL.obs; ot.style.fontWeight = "600";
    for (let v = 10; v <= hi; v += 10) {
      svgEl("line", { x1: hx(v), x2: hx(v), y1: hB, y2: hB + 5, class: "axis" }, layer);
      svgEl("text", { x: hx(v), y: hB + 22, "text-anchor": "middle", class: "sm" }, layer).textContent = String(v);
    }
    svgEl("text", { x: (X0 + X1) / 2, y: hB + 44, "text-anchor": "middle" }, layer).textContent =
      `number of A intervals (of ${NA}) that overlap a B interval`;

    // ---- readout ----
    const nullName = !st.excl && !st.gc ? "anywhere on the chromosome (naive)"
      : [st.excl ? "outside gaps, centromere and blacklist" : "", st.gc ? "within its own GC class" : ""].filter(x => x).join(", ");
    readout.innerHTML = `observed: <b style="color:${COL.obs}">${obs}</b> of ${NA} A intervals overlap a B interval · null: each B interval shuffled ${nullName}<br>` +
      (s.n === 0 ? "no shuffles yet: press “shuffle B”."
        : `N = ${s.n} shuffles · null mean ${s.mean.toFixed(1)} · <b>enrichment ${obs} / ${s.mean.toFixed(1)} = ${s.enrich.toFixed(2)}</b> · ` +
          `k = ${s.k} shuffle${s.k === 1 ? "" : "s"} ≥ ${obs} · <b>p = (${s.k} + 1) / (${s.n} + 1) = ${s.p < 0.001 ? s.p.toExponential(1) : s.p.toFixed(3)}</b>`) +
      `<br><span style="color:${COL.muted}">Dark A intervals overlap an observed B interval.${lastNote ? " " + lastNote[0].toUpperCase() + lastNote.slice(1) + "." : ""}</span>`;
  }

  draw();
}

autoMount(NAME, mount);
