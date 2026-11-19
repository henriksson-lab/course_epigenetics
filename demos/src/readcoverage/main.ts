// From reads to a coverage track. Simulated single-end ChIP-seq around one binding site:
// fragments that carry the protein contain the site; each is read from one end, so + strand
// reads sit left of the site and − strand reads right of it. Extending each read to the
// fragment length rebuilds the fragments and merges the two humps into one peak on the site.
// PCR copies make spikes of identical reads; removing duplicates keeps one read per start and
// strand. The coverage track is shown as raw reads per bin, CPM, or log2 fold change over a
// (simulated) input with a pseudocount. Model and its assumptions: model.ts.
//
// Auto-mounts on <div data-demo="readcoverage">. data-* options (all optional):
//   data-frag="200" data-read="36" data-bin="10" data-extend="on" data-norm="raw|cpm|fc" data-seed="1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount } from "../lib/ui";
import { makeDraws, simulate, track, summarize, binCounts, DEFAULTS, Params, Norm, Read, WIN, SITE, FOOT } from "./model";

const NAME = "readcoverage";
const COL = { plus: "#2a6f97", minus: "#e67e22", pcr: "#c0392b", cov: "#5b6b7a", input: "#8c564b", site: "#1b2a3a", frag: "#9aa7b3" };
const BINS = [1, 10, 25, 50, 100];
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .rc-ctrl { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: .35em 1.4em; width: 100%; }
[data-demo="${NAME}"] .rc-row { display: flex; flex-wrap: wrap; gap: .35em .9em; align-items: center; }
[data-demo="${NAME}"] .rc-row .qb-presets { display: inline-flex; }
[data-demo="${NAME}"] .rc-hidden { display: none !important; }
[data-demo="${NAME}"] .rc-wrap { width: 100%; max-width: 900px; overflow-x: auto; }
[data-demo="${NAME}"] svg { width: 100%; min-width: 600px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 16.5px; }
[data-demo="${NAME}"] svg .ttl { font-size: 17.5px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #5b6b7a; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
[data-demo="${NAME}"] .qb-readout.rc-readout { font-family: inherit; width: 100%; box-sizing: border-box; font-size: .88em; line-height: 1.5; }
[data-demo="${NAME}"] .rc-legend { display: flex; flex-wrap: wrap; gap: .2em 1.1em; font-size: .88em; color: #5b6b7a; }
[data-demo="${NAME}"] .rc-legend i { display: inline-block; width: 1.1em; height: .65em; margin-right: .35em; }
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
const fmt = (v: number): string => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2));
const signed = (v: number): string => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v);
const mega = (v: number): string => (v / 1e6).toFixed(2) + " M";

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const p0: Params = {
    ...DEFAULTS,
    fragLen: Number(ds.frag ?? DEFAULTS.fragLen) || DEFAULTS.fragLen,
    readLen: Number(ds.read ?? DEFAULTS.readLen) || DEFAULTS.readLen,
    bin: BINS.includes(Number(ds.bin)) ? Number(ds.bin) : DEFAULTS.bin,
    extend: ds.extend === "on",
    norm: (["raw", "cpm", "fc"].includes(ds.norm ?? "") ? ds.norm : "raw") as Norm,
  };
  const seed = Number(ds.seed ?? 1) || 1;
  const draws = makeDraws(seed);
  const p: Params = { ...p0 };
  let showFrag = false;

  // ---------------- controls ----------------
  const ctrl = htmlEl("div", "rc-ctrl", "", host);
  const sF = slider(ctrl, "fragment length", 80, 400, 10, p.fragLen, v => { p.fragLen = v; draw(); }, v => `${v} bp`);
  const sR = slider(ctrl, "read length", 25, 150, 1, p.readLen, v => { p.readLen = v; draw(); }, v => `${v} bp`);
  const row1 = htmlEl("div", "rc-row", "", host);
  const cExt = checkbox(row1, "extend reads to the fragment length", p.extend, v => { p.extend = v; draw(); });
  const cFrag = checkbox(row1, "show the true fragments (simulation only)", showFrag, v => { showFrag = v; draw(); });
  const row2 = htmlEl("div", "rc-row", "", host);
  const cDup = checkbox(row2, "add PCR duplicates", p.dups, v => { p.dups = v; draw(); });
  const cDedup = checkbox(row2, "remove duplicates", p.dedup, v => { p.dedup = v; draw(); });
  const row3 = htmlEl("div", "rc-row", "", host);
  htmlEl("span", "", "track:", row3);
  const normItems: [Norm, string][] = [["raw", "raw counts"], ["cpm", "CPM"], ["fc", "log2 fold change over input"]];
  const nBtns = presets(row3, normItems.map(([n, l]) => ({ label: l, apply: () => { p.norm = n; draw(); } })));
  const pcWrap = htmlEl("span", "rc-row", "", row3);
  htmlEl("span", "", "pseudocount:", pcWrap);
  const PCS = [0.1, 1, 5];
  const pcBtns = presets(pcWrap, PCS.map(v => ({ label: String(v), apply: () => { p.pseudo = v; draw(); } })));
  const row4 = htmlEl("div", "rc-row", "", host);
  htmlEl("span", "", "bin size:", row4);
  const bBtns = presets(row4, BINS.map(b => ({ label: `${b} bp`, apply: () => { p.bin = b; draw(); } })));
  button(row4, "reset", () => {
    Object.assign(p, p0); showFrag = false;
    sF.set(p.fragLen); sR.set(p.readLen); cExt.checked = p.extend; cFrag.checked = false;
    cDup.checked = p.dups; cDedup.checked = p.dedup; syncButtons(); draw();
  });
  function syncButtons(): void {
    nBtns.forEach((b, i) => b.classList.toggle("active", normItems[i][0] === p.norm));
    pcBtns.forEach((b, i) => b.classList.toggle("active", PCS[i] === p.pseudo));
    bBtns.forEach((b, i) => b.classList.toggle("active", BINS[i] === p.bin));
  }
  syncButtons();
  [...nBtns, ...pcBtns, ...bBtns].forEach(b => b.addEventListener("click", syncButtons));

  // ---------------- drawing ----------------
  const W = 800, X0 = 70, X1 = W - 20;
  const rT = 30, rB = 250;            // reads panel
  const cT = 318, cB = 488;           // coverage panel
  const H = cB + 48;
  const readout = htmlEl("div", "qb-readout rc-readout", "", host);
  const legend = htmlEl("div", "rc-legend", "", host);
  const wrap = htmlEl("div", "rc-wrap", "", host);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
    "aria-label": "Simulated single-end reads around a binding site and the coverage track made from them" }, wrap);
  htmlEl("p", "qb-note",
    "Simulated data: 190 ChIP fragments (150 of them carry the bound protein) and an input sequenced at half the depth. " +
    "The library sizes stand for the whole genome; the window shows 1200 bp of it.", host);
  const sx = scale(0, WIN, X0, X1);
  const clipId = `rc-clip-${Math.random().toString(36).slice(2, 8)}`;
  const defs = svgEl("defs", {}, svg);
  const cp = svgEl("clipPath", { id: clipId }, defs);
  svgEl("rect", { x: X0, y: 0, width: X1 - X0, height: H }, cp);
  const layer = svgEl("g", {}, svg);

  function pile(reads: Read[]): number[] {
    // greedy rows over the drawn span; reads sorted by left end
    const idx = reads.map((_, i) => i);
    const lo = (r: Read) => (showFrag ? Math.min(r.start, r.fragStart) : r.start);
    const hi = (r: Read) => (showFrag ? Math.max(r.end, r.fragEnd) : r.end);
    idx.sort((a, b) => lo(reads[a]) - lo(reads[b]) || Number(reads[a].pcr) - Number(reads[b].pcr));
    const rowEnd: number[] = [], row = new Array(reads.length).fill(0);
    for (const i of idx) {
      let k = rowEnd.findIndex(e => e + 3 <= lo(reads[i]));
      if (k < 0) { k = rowEnd.length; rowEnd.push(0); }
      rowEnd[k] = hi(reads[i]); row[i] = k;
    }
    return row;
  }

  function draw(): void {
    pcWrap.classList.toggle("rc-hidden", p.norm !== "fc");
    const s = simulate(draws, p);
    const t = track(s, p);
    const m = summarize(s, p, t);
    layer.replaceChildren();

    // ---- the numbers (HTML, so they reflow on a narrow screen) ----
    const removed = s.removedPcr + s.removedIndep;
    const unitShort = p.norm === "raw" ? "reads" : p.norm === "cpm" ? "CPM" : "log2 FC";
    const span = (txt: string, c: string) => `<span style="color:${c};font-weight:600">${txt}</span>`;
    readout.innerHTML =
      `ChIP reads in the window: <b>${m.kept}</b> (${span("+ " + m.nPlus, COL.plus)} / ${span("− " + m.nMinus, COL.minus)})` +
      (p.dups ? ` · ${span("PCR copies added: " + s.pcrCopies, COL.pcr)}` : "") +
      ` · library: ${mega(s.chipLib)} reads<br>` +
      (p.dedup ? `removed as duplicates: <b>${removed}</b> (${s.removedPcr} PCR copies, ${s.removedIndep} independent fragments that happen to share a start)<br>` : "") +
      `read hump centres: + ${signed(m.cPlus)} bp, − ${signed(m.cMinus)} bp<br>` +
      `<b>highest bin: ${fmt(m.maxVal)} ${unitShort} at ${signed(m.maxPos)} bp</b> · bin with the site: ${fmt(m.siteVal)} · peak width at half height: ${m.halfWidth} bp`;
    const lg: [string, string, number][] = [["+ strand read", COL.plus, 1], ["− strand read", COL.minus, 1], ["coverage of all reads", COL.cov, 0.55]];
    if (p.extend) lg.push(["extension", COL.cov, 0.3]);
    if (p.dups) lg.push(["PCR copy", COL.pcr, 1]);
    if (p.dedup) lg.push(["removed (faint)", COL.cov, 0.15]);
    if (p.norm !== "fc") lg.push(["input, scaled to the ChIP depth (dashed)", COL.input, 1]);
    legend.innerHTML = lg.map(([txt, c, o]) =>
      `<span><i style="background:${c};opacity:${o}"></i>${txt}</span>`).join("");

    // ---- reads panel ----
    svgEl("text", { x: 8, y: rT - 6, class: "ttl" }, layer).textContent = "reads";
    const gReads = svgEl("g", { "clip-path": `url(#${clipId})` }, layer);
    const reads = s.chip;
    const rows = pile(reads);
    const nRows = Math.max(1, ...rows.map(r => r + 1));
    const rh = Math.min(7, (rB - rT) / nRows);
    const hBar = Math.max(1.2, rh * 0.75);
    reads.forEach((r, i) => {
      const y = rT + rows[i] * rh;
      const g = svgEl("g", r.removed ? { opacity: 0.15 } : {}, gReads);
      const col = r.pcr ? COL.pcr : r.plus ? COL.plus : COL.minus;
      if (showFrag) svgEl("line", { x1: sx(r.fragStart), x2: sx(r.fragEnd), y1: y + hBar / 2, y2: y + hBar / 2, stroke: COL.frag, "stroke-width": 1 }, g);
      const rl = Math.min(p.readLen, r.fragEnd - r.fragStart);
      const r0 = r.plus ? r.start : r.end - rl, r1 = r.plus ? r.start + rl : r.end;
      if (r.end - r.start > rl) // the extension, lighter
        svgEl("rect", { x: sx(r.start), y, width: Math.max(0.5, sx(r.end) - sx(r.start)), height: hBar, fill: col, "fill-opacity": 0.3 }, g);
      svgEl("rect", { x: sx(r0), y, width: Math.max(0.8, sx(r1) - sx(r0)), height: hBar, fill: col }, g);
    });
    // site marker across both panels
    svgEl("rect", { x: sx(SITE - FOOT), y: rT - 2, width: sx(SITE + FOOT) - sx(SITE - FOOT), height: cB - rT + 2, fill: "#1b2a3a", "fill-opacity": 0.06 }, layer);
    svgEl("line", { x1: sx(SITE), x2: sx(SITE), y1: rT - 2, y2: cB, stroke: COL.site, "stroke-width": 1.3, "stroke-dasharray": "5 4" }, layer);
    const sl = svgEl("text", { x: sx(SITE), y: rB + 17, "text-anchor": "middle", class: "sm" }, layer);
    sl.textContent = "binding site"; sl.style.fontWeight = "600";

    // ---- coverage panel ----
    const v = t.values;
    const fc = p.norm === "fc";
    const inputVals = fc ? [] : t.inputScaled.map(x => (p.norm === "cpm" ? x * 1e6 / s.chipLib : x));
    const plusV = binCounts(s.chip, p.bin, "plus"), minusV = binCounts(s.chip, p.bin, "minus");
    const k = p.norm === "cpm" ? 1e6 / s.chipLib : 1;
    let lo = 0, hi = Math.max(...v, ...inputVals, 1e-9);
    if (fc) { lo = Math.min(-1, ...v); hi = Math.max(1, ...v); }
    const step = niceStep((hi - lo) / 4);
    lo = Math.floor(lo / step) * step; hi = Math.max(lo + step, Math.ceil(hi / step) * step);
    const sy = scale(lo, hi, cB, cT);
    svgEl("text", { x: 8, y: cT - 16, class: "ttl" }, layer).textContent =
      `coverage: ${t.unit}, ${p.bin} bp bins${fc ? `, pseudocount ${p.pseudo}` : ""}`;
    for (let y = lo; y <= hi + 1e-9; y += step) {
      svgEl("line", { x1: X0, x2: X1, y1: sy(y), y2: sy(y), class: Math.abs(y) < 1e-9 ? "axis" : "grid" }, layer);
      svgEl("text", { x: X0 - 6, y: sy(y) + 5, "text-anchor": "end", class: "sm" }, layer).textContent = String(+y.toFixed(6));
    }
    const gCov = svgEl("g", { "clip-path": `url(#${clipId})` }, layer);
    const y0 = sy(fc ? 0 : lo);
    let d = `M${sx(0)},${y0}`;
    v.forEach((val, i) => {
      const a = sx(i * p.bin), b = sx(Math.min(WIN, (i + 1) * p.bin));
      d += `L${a.toFixed(1)},${sy(val).toFixed(1)}L${b.toFixed(1)},${sy(val).toFixed(1)}`;
    });
    d += `L${sx(WIN)},${y0}Z`;
    svgEl("path", { d, fill: COL.cov, "fill-opacity": 0.55, stroke: COL.cov, "stroke-width": 1 }, gCov);
    const stepLine = (vals: number[], col: string, w: number, dash = "") => {
      let dd = "";
      vals.forEach((val, i) => {
        const a = sx(i * p.bin), b = sx(Math.min(WIN, (i + 1) * p.bin));
        dd += `${i ? "L" : "M"}${a.toFixed(1)},${sy(val).toFixed(1)}L${b.toFixed(1)},${sy(val).toFixed(1)}`;
      });
      svgEl("path", { d: dd, fill: "none", stroke: col, "stroke-width": w, ...(dash ? { "stroke-dasharray": dash } : {}) }, gCov);
    };
    if (!fc) {
      stepLine(plusV.map(x => x * k), COL.plus, 1.8);
      stepLine(minusV.map(x => x * k), COL.minus, 1.8);
      stepLine(inputVals, COL.input, 2, "6 4");
    }
    // x axis
    svgEl("line", { x1: X0, x2: X1, y1: cB, y2: cB, class: "axis" }, layer);
    for (let x = -600; x <= 600; x += 200) {
      svgEl("line", { x1: sx(SITE + x), x2: sx(SITE + x), y1: cB, y2: cB + 5, class: "axis" }, layer);
      svgEl("text", { x: sx(SITE + x), y: cB + 22, "text-anchor": "middle", class: "sm" }, layer).textContent = signed(x);
    }
    svgEl("text", { x: (X0 + X1) / 2, y: cB + 42, "text-anchor": "middle" }, layer).textContent = "position relative to the binding site (bp)";
  }

  draw();
}

autoMount(NAME, mount);
