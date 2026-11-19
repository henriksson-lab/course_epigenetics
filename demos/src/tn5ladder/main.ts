// Where the ATAC-seq nucleosome ladder comes from. A simulated chromatin fibre with nucleosomes
// (linker length, nucleosome occupancy) and a nucleosome-depleted region (NDR) in the middle.
// Tn5 inserts into accessible DNA (linkers, the NDR, empty nucleosome slots; rarely inside a
// nucleosome), the DNA is cut at every insertion, each end gets adapter A or B, and only A–B
// fragments are sequenced. The histogram collects the sequenced fragment lengths from many fibres:
// nucleosome-free, mono- and di-nucleosome peaks. More or less Tn5 shows over- and
// under-tagmentation. Model: model.ts.
//
// Auto-mounts on <div data-demo="tn5ladder">. data-* options (all optional):
//   data-linker="45" data-ndr="300" data-occ="0.9" data-tn5="1" data-helix="on" data-log="on"
//   data-fibres="2000" data-seed="1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount, rng } from "../lib/ui";
import { histogram, fibre, peaks, DEFAULTS, Params, Hist, Fibre, MIN_LEN, MAX_LEN, FIB } from "./model";

const NAME = "tn5ladder";
const COL = {
  A: "#2a6f97", B: "#e67e22", seq: "#2e8b57", lost: "#9aa7b3", nuc: "#2a9d8f", dna: "#1b2a3a",
  text: "#1b2a3a", muted: "#5b6b7a", ins: "#c0392b", bar: "#5b6b7a",
};
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .tl-ctrl { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: .35em 1.4em; width: 100%; }
[data-demo="${NAME}"] .tl-row { display: flex; flex-wrap: wrap; gap: .35em 1em; align-items: center; }
[data-demo="${NAME}"] .tl-row .qb-presets { display: inline-flex; }
[data-demo="${NAME}"] .qb-readout.tl-readout { font-family: inherit; width: 100%; box-sizing: border-box; font-size: .9em; line-height: 1.5; }
[data-demo="${NAME}"] .tl-wrap { width: 100%; max-width: 900px; overflow-x: auto; }
[data-demo="${NAME}"] svg { width: 100%; min-width: 600px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 16.5px; }
[data-demo="${NAME}"] svg .ttl { font-size: 17.5px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15.5px; fill: #5b6b7a; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
[data-demo="${NAME}"] .tl-legend { display: flex; flex-wrap: wrap; gap: .2em 1.1em; font-size: .88em; color: #5b6b7a; align-items: center; }
[data-demo="${NAME}"] .tl-legend i { display: inline-block; width: 1.4em; height: .8em; vertical-align: middle; box-sizing: border-box; }
[data-demo="${NAME}"] .tl-legend i.tl-sq { width: .45em; height: .9em; }
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
const pct = (a: number, b: number) => (b ? (100 * a / b).toFixed(1) : "–") + "%";
const tn5Fmt = (l2: number) => "×" + +(2 ** l2).toFixed(2);
const thin = (n: number) => n.toLocaleString("en-US").replace(/,/g, " ");

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (v: string | undefined, d: number) => (v !== undefined && Number.isFinite(Number(v)) ? Number(v) : d);
  const p0: Params = {
    linker: num(ds.linker, DEFAULTS.linker), ndr: num(ds.ndr, DEFAULTS.ndr), occ: num(ds.occ, DEFAULTS.occ),
    tn5: num(ds.tn5, DEFAULTS.tn5), helix: ds.helix === "on",
  };
  const NF = Math.max(200, Math.round(num(ds.fibres, 2000)));
  const seed = num(ds.seed, 1);
  const p: Params = { ...p0 };
  let logY = ds.log === "on";
  let example = 0;

  // ---------------- controls ----------------
  const ctrl = htmlEl("div", "tl-ctrl", "", host);
  const sL = slider(ctrl, "linker length", 10, 80, 1, p.linker, v => { p.linker = v; schedule(); }, v => `${v} bp`);
  const sN = slider(ctrl, "NDR width", 0, 800, 10, p.ndr, v => { p.ndr = v; schedule(); }, v => `${v} bp`);
  const sO = slider(ctrl, "nucleosome occupancy", 0.5, 1, 0.05, p.occ, v => { p.occ = v; schedule(); }, v => `${Math.round(v * 100)}%`);
  const sT = slider(ctrl, "Tn5 amount", -2, 3, 0.25, Math.log2(p.tn5), v => { p.tn5 = 2 ** v; schedule(); }, tn5Fmt);
  const row = htmlEl("div", "tl-row", "", host);
  const cH = checkbox(row, "10.5 bp helical preference", p.helix, v => { p.helix = v; schedule(); });
  const cLog = checkbox(row, "log scale", logY, v => { logY = v; draw(); });
  button(row, "another fibre", () => { example++; draw(); });
  const row2 = htmlEl("div", "tl-row", "", host);
  htmlEl("span", "", "Tn5 per nucleus:", row2);
  const tPresets: [string, number][] = [["too little (×¼)", -2], ["standard (×1)", 0], ["too much (×8)", 3]];
  presets(row2, tPresets.map(([l, v]) => ({ label: l, apply: () => { p.tn5 = 2 ** v; sT.set(v); schedule(); } })));
  button(row2, "reset", () => {
    Object.assign(p, p0); logY = ds.log === "on"; example = 0;
    sL.set(p.linker); sN.set(p.ndr); sO.set(p.occ); sT.set(Math.log2(p.tn5)); cH.checked = p.helix; cLog.checked = logY;
    schedule();
  });
  const readout = htmlEl("div", "qb-readout tl-readout", "", host);

  // ---------------- drawing ----------------
  const W = 800, X0 = 80, X1 = W - 20;
  const yDNA = 74, yF1 = 112, yF2 = 132;
  const hT = 212, hB = 450, H = hB + 52;
  const legend = htmlEl("div", "tl-legend", "", host);
  legend.innerHTML = [
    [`<b style="color:${COL.ins}">▼</b>`, "Tn5 insertion"],
    [`<i style="background:${COL.nuc};border-radius:4px"></i>`, "nucleosome"],
    [`<i style="border:1px dashed ${COL.muted};border-radius:4px"></i>`, "empty nucleosome slot"],
    [`<i class="tl-sq" style="background:${COL.A}"></i>`, "adapter A"],
    [`<i class="tl-sq" style="background:${COL.B}"></i>`, "adapter B"],
    [`<i style="background:${COL.seq};height:.3em"></i>`, "A–B fragment: sequenced"],
    [`<i style="border-top:2px dashed ${COL.lost};height:0"></i>`, "A–A or B–B (or &lt; 35 bp): lost"],
  ].map(([sym, txt]) => `<span>${sym} ${txt}</span>`).join("");
  const wrap = htmlEl("div", "tl-wrap", "", host);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
    "aria-label": "A chromatin fibre cut by Tn5 and the fragment-length histogram of many such fibres" }, wrap);
  htmlEl("p", "qb-note",
    `Simulated: ${thin(NF)} fibres of ${FIB / 1000} kb (more when there is less Tn5, so the histogram is not just noise), each with an NDR in the middle; nucleosome core 147 bp; ` +
    `Tn5 cuts inside a nucleosome at 1/10 of the rate in free DNA. Fragments < ${MIN_LEN} bp are lost in the clean-up; ` +
    `fragments > ${MAX_LEN} bp are not shown. The 9-bp target duplication and Tn5's sequence bias are ignored.`, host);
  const layer = svgEl("g", {}, svg);
  const V0 = 1300, V1 = 3700;              // the part of the example fibre that is drawn
  const fx = scale(V0, V1, X0, X1);
  // the same sequencing depth at every Tn5 amount: with less Tn5, more nuclei (fibres) are needed
  const nFib = () => Math.round(NF / Math.min(1, p.tn5));
  let hist: Hist = histogram(p, nFib(), seed);

  let pending = 0;
  function schedule(): void {
    if (pending) return;
    pending = requestAnimationFrame(() => { pending = 0; hist = histogram(p, nFib(), seed); draw(); });
  }

  function drawFibre(fb: Fibre): void {
    svgEl("text", { x: 8, y: 20, class: "ttl" }, layer).textContent = `one simulated fibre (the middle ${(V1 - V0) / 1000} kb)`;
    const clip = (a: number, b: number): [number, number] | null => (b <= V0 || a >= V1 ? null : [Math.max(a, V0), Math.min(b, V1)]);
    // NDR
    const nd = clip(fb.ndr[0], fb.ndr[1]);
    if (nd && nd[1] > nd[0]) {
      svgEl("rect", { x: fx(nd[0]), y: yDNA - 30, width: fx(nd[1]) - fx(nd[0]), height: 60, fill: "#2e8b57", "fill-opacity": 0.1 }, layer);
      const t = svgEl("text", { x: (fx(nd[0]) + fx(nd[1])) / 2, y: yDNA - 34, "text-anchor": "middle", class: "sm" }, layer);
      t.textContent = "NDR"; t.style.fill = "#2e8b57"; t.style.fontWeight = "600";
    }
    svgEl("line", { x1: X0, x2: X1, y1: yDNA, y2: yDNA, stroke: COL.dna, "stroke-width": 2 }, layer);
    for (const n of fb.nucs) {
      const c = clip(n.start, n.end);
      if (!c) continue;
      if (n.present)
        svgEl("rect", { x: fx(c[0]), y: yDNA - 15, width: fx(c[1]) - fx(c[0]), height: 30, rx: 12, fill: COL.nuc, "fill-opacity": 0.75, stroke: COL.text, "stroke-width": 1 }, layer);
      else
        svgEl("rect", { x: fx(c[0]), y: yDNA - 15, width: fx(c[1]) - fx(c[0]), height: 30, rx: 12, fill: "none", stroke: COL.muted, "stroke-width": 1, "stroke-dasharray": "3 3" }, layer);
    }
    for (const x of fb.ins) {
      if (x < V0 || x > V1) continue;
      svgEl("path", { d: `M${fx(x)},${yDNA - 4}l-4,-11h8z`, fill: COL.ins }, layer);
    }
    // fragments, alternating rows
    fb.frags.forEach((f, i) => {
      const c = clip(f.start, f.end);
      if (!c) return;
      const y = i % 2 ? yF2 : yF1;
      const len = f.end - f.start;
      const seq = f.ok && len >= MIN_LEN && len <= MAX_LEN;
      svgEl("line", { x1: fx(c[0]) + 1, x2: fx(c[1]) - 1, y1: y, y2: y, stroke: seq ? COL.seq : COL.lost, "stroke-width": seq ? 5 : 3,
        ...(f.ok ? {} : { "stroke-dasharray": "4 3" }) }, layer);
      if (f.start >= V0) svgEl("rect", { x: fx(f.start) + 1, y: y - 6, width: 5, height: 12, fill: f.a === "A" ? COL.A : COL.B }, layer);
      if (f.end <= V1) svgEl("rect", { x: fx(f.end) - 6, y: y - 6, width: 5, height: 12, fill: f.b === "A" ? COL.A : COL.B }, layer);
    });
  }

  function drawHist(h: Hist): void {
    const hx = scale(0, MAX_LEN, X0, X1);
    const c = h.counts;
    const cmax = Math.max(1, ...c);
    let hy: (v: number) => number;
    const ticks: number[] = [];
    if (logY) {
      const top = 10 ** Math.ceil(Math.log10(cmax));
      const ly = scale(0, Math.log10(top), hB, hT);
      hy = v => ly(Math.log10(Math.max(1, v)));
      for (let v = 1; v <= top; v *= 10) ticks.push(v);
    } else {
      const step = niceStep(cmax / 4);
      const top = Math.ceil(cmax / step) * step;
      hy = scale(0, top, hB, hT);
      for (let v = 0; v <= top + 1e-9; v += step) ticks.push(v);
    }
    svgEl("text", { x: 8, y: hT - 42, class: "ttl" }, layer).textContent =
      `sequenced fragment lengths (${thin(h.seq)} fragments, ${thin(nFib())} fibres)`;
    // class bands
    const bands: [number, number, string][] = [[0, 100, "nucleosome-free"], [180, 250, "mono"], [315, 473, "di"]];
    for (const [a, b, lab] of bands) {
      svgEl("rect", { x: hx(a), y: hT - 6, width: hx(b) - hx(a), height: hB - hT + 6, fill: "#2a6f97", "fill-opacity": 0.07 }, layer);
      svgEl("text", { x: (hx(a) + hx(b)) / 2, y: hT - 13, "text-anchor": "middle", class: "sm" }, layer).textContent = lab;
    }
    for (const v of ticks) {
      svgEl("line", { x1: X0, x2: X1, y1: hy(v), y2: hy(v), class: v === ticks[0] ? "axis" : "grid" }, layer);
      svgEl("text", { x: X0 - 6, y: hy(v) + 5, "text-anchor": "end", class: "sm" }, layer).textContent = thin(v);
    }
    let d = `M${hx(0)},${hB}`;
    c.forEach((v, i) => {
      const a = hx(i * h.bin), b = hx((i + 1) * h.bin);
      const y = v > 0 ? hy(v) : hB;
      d += `L${a.toFixed(1)},${y.toFixed(1)}L${b.toFixed(1)},${y.toFixed(1)}`;
    });
    d += `L${hx(MAX_LEN)},${hB}Z`;
    svgEl("path", { d, fill: COL.bar, "fill-opacity": 0.6, stroke: COL.bar, "stroke-width": 1 }, layer);
    svgEl("rect", { x: hx(0), y: hT - 6, width: hx(MIN_LEN) - hx(0), height: hB - hT + 6, fill: "#fff", "fill-opacity": 0.6 }, layer);
    for (let x = 0; x <= MAX_LEN; x += 100) {
      svgEl("line", { x1: hx(x), x2: hx(x), y1: hB, y2: hB + 5, class: "axis" }, layer);
      svgEl("text", { x: hx(x), y: hB + 22, "text-anchor": "middle", class: "sm" }, layer).textContent = String(x);
    }
    svgEl("text", { x: (X0 + X1) / 2, y: hB + 44, "text-anchor": "middle" }, layer).textContent = `fragment length (bp), ${h.bin} bp bins`;
    const yl = svgEl("text", { x: 22, y: (hT + hB) / 2, "text-anchor": "middle", transform: `rotate(-90 22 ${(hT + hB) / 2})`, class: "sm" }, layer);
    yl.textContent = logY ? "fragments (log scale)" : "fragments";
  }

  function draw(): void {
    layer.replaceChildren();
    const fb = fibre(p, rng(seed * 7919 + 17 + example).unif);
    drawFibre(fb);
    drawHist(hist);
    const h = hist;
    const ab = h.total - h.sameAdapter;
    const pk = peaks(h).slice(0, 4);
    readout.innerHTML =
      `Tn5 insertions: ${h.insPerKb.toFixed(1)} per kb · fragments between two insertions: ${thin(h.total)} · ` +
      `<b>A–A or B–B (not amplified): ${pct(h.sameAdapter, h.total)}</b> · A–B but &lt; ${MIN_LEN} bp: ${pct(h.tooShort, ab)} of A–B · ` +
      `A–B &gt; ${MAX_LEN} bp: ${pct(h.tooLong, ab)}<br>` +
      `of the sequenced fragments: <b>nucleosome-free (&lt; 100 bp) ${pct(h.nfr, h.seq)}</b> · <b>mono (180–250) ${pct(h.mono, h.seq)}</b> · ` +
      `di (315–473) ${pct(h.di, h.seq)} · peaks of the ladder near ${pk.join(", ")} bp`;
  }

  draw();
}

autoMount(NAME, mount);
