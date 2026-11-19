// Leakage through overlapping windows and homologous copies (model.ts).
// Top: the simulated genome, 10 chromosomes x 30 regions; each region is cut into 3 windows
// (staggered bars; their overlap set by a slider), coloured by the split shown (orange = test,
// grey = training); a small square under each region gives its label (green open, crimson
// closed); purple regions are homologous copies of a region on another chromosome (line).
// Bottom: 1-nearest-neighbour accuracy, 15-NN ROC AUC and the GC-score AUC (no training) for
// a random window split, a region split and a chromosome split, with the chance line, and how
// often the nearest training window is a shifted copy of the test window or a homologous copy.
//
// Auto-mounts on <div data-demo="leakagesim">. data-* options (all optional):
//   data-overlap="0.7"  data-dup="0"  data-signal="0.08"  data-seed="1"  data-show="random"
import { svgEl, htmlEl, slider, button, presets, scale, autoMount } from "../lib/ui";
import * as M from "./model";

const NAME = "leakagesim";
const COL = {
  ink: "#1b2a3a", muted: "#5b6b7a", grid: "#e3e7eb", test: "#e67e22", train: "#b8c2cc",
  open: "#2e8b57", closed: "#a8324a", dup: "#6b4fa0", acc: "#2a6f97", knn: "#6b4fa0", gc: "#9aa5b1", chance: "#c0392b",
};
const SPLIT_LABEL: Record<M.SplitName, string> = { random: "random windows", region: "whole regions", chromosome: "whole chromosomes" };
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .ls-controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(19em, 1fr)); gap: .3em 1.4em; width: 100%; }
[data-demo="${NAME}"] .ls-row { display: flex; flex-wrap: wrap; gap: .35em 1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .ls-scroll { width: 100%; overflow-x: auto; }
[data-demo="${NAME}"] .ls-scroll svg { width: 100%; min-width: 720px; height: auto; display: block; }
[data-demo="${NAME}"] svg.ls-chart { width: 100%; height: auto; max-width: 1100px; }
[data-demo="${NAME}"] svg text { font-size: 16px; }
[data-demo="${NAME}"] svg .ttl { font-size: 18px; font-weight: 700; }
[data-demo="${NAME}"] svg .sm { font-size: 15px; fill: ${COL.muted}; }
[data-demo="${NAME}"] svg .num { font-size: 17px; font-weight: 700; }
[data-demo="${NAME}"] .ls-busy { color: ${COL.muted}; font-style: italic; }
.reveal [data-demo="${NAME}"] { font-size: 17px; gap: .25em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const f2 = (v: number): string => v.toFixed(2);
const pct = (v: number): string => Math.round(100 * v) + "%";
const num = (s: string | undefined, d: number, lo: number, hi: number): number => {
  const v = Number(s);
  return s !== undefined && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
};

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const DEF: M.Params = {
    overlap: num(ds.overlap, 0.7, 0, 0.9), dup: num(ds.dup, 0, 0, 0.5), signal: num(ds.signal, 0.08, 0, 0.12),
    seed: num(ds.seed, 1, 1, 1e6),
  };
  const p: M.Params = { ...DEF };
  let show: M.SplitName = (M.SPLITS as string[]).includes(ds.show ?? "") ? ds.show as M.SplitName : "random";
  let cur: ReturnType<typeof M.runAll> | null = null;
  let pending = 0;

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "ls-controls", "", host);
  const sl = {
    overlap: slider(ctr, "window overlap", 0, 0.9, 0.05, p.overlap, v => { p.overlap = v; schedule(); }, pct),
    dup: slider(ctr, "regions that are copies", 0, 0.5, 0.05, p.dup, v => { p.dup = v; schedule(); }, pct),
    signal: slider(ctr, "true signal (GC open − closed)", 0, 0.12, 0.01, p.signal, v => { p.signal = v; schedule(); }, f2),
  };
  const row = htmlEl("div", "ls-row", "", host);
  const pre = presets(row, [
    { label: "unit R1: 70% overlap", apply: () => setAll({ overlap: 0.7, dup: 0, signal: 0.08 }) },
    { label: "no signal at all", apply: () => setAll({ signal: 0 }) },
    { label: "no overlap", apply: () => setAll({ overlap: 0 }) },
    { label: "30% copies, no signal", apply: () => setAll({ overlap: 0, dup: 0.3, signal: 0 }) },
  ]);
  button(row, "new random genome", () => { p.seed++; schedule(); });
  button(row, "reset", () => { pre.forEach(b => b.classList.remove("active")); setAll({ ...DEF }); });
  const row2 = htmlEl("div", "ls-row", "", host);
  htmlEl("span", "", "Show test set of:", row2);
  const showBtns = presets(row2, M.SPLITS.map(s => ({ label: SPLIT_LABEL[s], apply: () => { show = s; drawMap(); } })));
  const busy = htmlEl("span", "ls-busy", "", row2);

  function setAll(q: Partial<M.Params>): void {
    Object.assign(p, q);
    sl.overlap.set(p.overlap); sl.dup.set(p.dup); sl.signal.set(p.signal);
    schedule();
  }

  // ---------------- genome map ----------------
  const scroll = htmlEl("div", "ls-scroll", "", host);
  // drawn at the width available (720-1100), so text keeps its nominal size
  const W = Math.max(720, Math.min(1100, host.clientWidth || 1100)), ml = 70, mr = 10, rowH = 32;
  const top = W < 1000 ? 80 : 58;
  const H = top + M.NCHROM * rowH + 34;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Simulated genome: regions, windows and the test set of the chosen split" }, scroll);
  const gMap = svgEl("g", {}, svg);

  // ---------------- results chart ----------------
  const hw = host.clientWidth;
  const narrow = hw > 0 && hw < 640, medium = !narrow && hw > 0 && hw < 1000;
  const CW = narrow ? 420 : medium ? 780 : 1100, CH = narrow ? 480 : 380;
  const chart = svgEl("svg", { class: "ls-chart", viewBox: `0 0 ${CW} ${CH}`, role: "img", "aria-label": "Test accuracy and AUC for three ways of splitting" }, host);
  const gChart = svgEl("g", {}, chart);
  const note = htmlEl("div", "qb-note", "", host);
  note.textContent = `Simulated data: ${M.NCHROM * M.PER_CHROM} regions with random sequence; ${M.NWIN} windows of ${M.WINLEN} bp per region; ` +
    `features = counts of the 256 4-mers; k-nearest neighbours as in the unit's R1. Copies differ from their source at ${100 * M.MUT}% of bases.`;

  function schedule(): void {
    busy.textContent = "computing…";
    const id = ++pending;
    window.setTimeout(() => {
      if (id !== pending) return;
      cur = M.runAll(p);
      busy.textContent = "";
      drawMap(); drawChart();
    }, 60);
  }

  // ---------------- drawing ----------------
  function drawMap(): void {
    if (!cur) return;
    showBtns.forEach((b, i) => b.classList.toggle("active", M.SPLITS[i] === show));
    gMap.replaceChildren();
    const { data } = cur, test = cur.tests[show];
    const t = (x: number, y: number, s: string, attrs: Record<string, string | number> = {}) => {
      const e = svgEl("text", { x, y, ...attrs }, gMap); e.textContent = s; return e;
    };
    const nTest = test.filter(Boolean).length;
    t(4, 20, `Test set of the "${SPLIT_LABEL[show]}" split: ${nTest} of ${test.length} windows`, { class: "ttl" });
    // legend
    let lx = 4, ly = 44;
    const item = (draw: (x: number) => void, txt: string, w: number) => {
      if (lx + w > W) { lx = 4; ly += 22; }
      draw(lx); t(lx + 22, ly, txt, { class: "sm" }); lx += w;
    };
    item(x => svgEl("rect", { x, y: ly - 10, width: 16, height: 6, fill: COL.test }, gMap), "test window", 120);
    item(x => svgEl("rect", { x, y: ly - 10, width: 16, height: 6, fill: COL.train }, gMap), "training window", 150);
    item(x => svgEl("rect", { x: x + 4, y: ly - 11, width: 9, height: 9, fill: COL.open }, gMap), "open", 70);
    item(x => svgEl("rect", { x: x + 4, y: ly - 11, width: 9, height: 9, fill: COL.closed }, gMap), "closed", 82);
    item(x => svgEl("rect", { x, y: ly - 13, width: 16, height: 14, fill: "none", stroke: COL.dup, "stroke-width": 2, rx: 2 }, gMap), "homologous copy (dark line: one in test, one in training)", 400);
    const bw = (W - ml - mr) / M.PER_CHROM;
    const len = M.WINLEN + (M.NWIN - 1) * data.step;
    const bx = (r: number): number => ml + (r % M.PER_CHROM) * bw;
    const by = (r: number): number => top + Math.floor(r / M.PER_CHROM) * rowH;
    for (let c = 0; c < M.NCHROM; c++) {
      const y = top + c * rowH;
      const e = t(ml - 10, y + 17, `chr${c + 1}`, { "text-anchor": "end", "font-weight": M.TEST_CHROMS.includes(c) && show === "chromosome" ? 700 : 400 });
      if (M.TEST_CHROMS.includes(c) && show === "chromosome") e.style.fill = COL.test;
      svgEl("line", { x1: ml, x2: W - mr, y1: y + 12, y2: y + 12, stroke: COL.grid, "stroke-width": 1 }, gMap);
    }
    const links = svgEl("g", {}, gMap);
    data.regions.forEach((rg, r) => {
      const x0 = bx(r) + 3, w0 = bw - 6, y0 = by(r);
      const sx = scale(0, len, x0, x0 + w0);
      for (let w = 0; w < M.NWIN; w++) {
        const wi = r * M.NWIN + w;
        svgEl("rect", {
          x: sx(w * data.step), y: y0 + 3 + w * 6, width: Math.max(1, sx(w * data.step + M.WINLEN) - sx(w * data.step)), height: 5,
          fill: test[wi] ? COL.test : COL.train,
        }, gMap);
      }
      svgEl("rect", { x: x0 + w0 / 2 - 4.5, y: y0 + 22, width: 9, height: 9, fill: rg.label ? COL.open : COL.closed }, gMap);
      if (rg.dupOf >= 0) {
        svgEl("rect", { x: x0 - 2, y: y0, width: w0 + 4, height: 32, fill: "none", stroke: COL.dup, "stroke-width": 2, rx: 3 }, gMap);
        const s = rg.dupOf;
        const xa = bx(s) + bw / 2, ya = by(s) + 12, xb = x0 + w0 / 2, yb = y0 + 12;
        const mx = (xa + xb) / 2 + (yb - ya) * 0.25;
        // a pair "leaks" when one copy is in the test set and the other in training
        const inTest = (q: number): boolean => test[q * M.NWIN] || test[q * M.NWIN + 1] || test[q * M.NWIN + 2];
        const inTrain = (q: number): boolean => !test[q * M.NWIN] || !test[q * M.NWIN + 1] || !test[q * M.NWIN + 2];
        const leak = (inTest(r) && inTrain(s)) || (inTest(s) && inTrain(r));
        svgEl("path", {
          d: `M${xa},${ya} Q${mx},${(ya + yb) / 2} ${xb},${yb}`, fill: "none", stroke: COL.dup,
          "stroke-width": leak ? 2 : 1.2, "stroke-opacity": leak ? 0.85 : 0.18,
        }, links);
      }
    });
    t(4, H - 8, `Each region: ${M.NWIN} windows of ${M.WINLEN} bp, ${data.step} bp apart (overlap ${pct(p.overlap)}); square = label of the region.`, { class: "sm" });
  }

  function drawChart(): void {
    if (!cur) return;
    gChart.replaceChildren();
    const t = (x: number, y: number, s: string, attrs: Record<string, string | number> = {}) => {
      const e = svgEl("text", { x, y, ...attrs }, gChart); e.textContent = s; return e;
    };
    t(4, 22, narrow ? "Test performance by split" : "How good does the model look? Test performance by split", { class: "ttl" });
    if (narrow) t(4, CH - 4, "same / copy: where the nearest training window is from", { class: "sm", style: "font-size:13px" });
    const ml2 = 56, mr2 = 10, yTop = narrow ? 100 : 84, yBot = CH - (narrow ? 112 : 98);
    const y = scale(0, 1, yBot, yTop);
    // legend
    const items: [string, string][] = [["1-NN accuracy", COL.acc], [`${M.KNN_K}-NN ROC AUC`, COL.knn], ["GC score AUC (no training)", COL.gc], ["chance", COL.chance]];
    let lx = 4, lyy = 50;
    for (const [s, c] of items) {
      if (s === "chance") svgEl("line", { x1: lx, x2: lx + 16, y1: lyy - 5, y2: lyy - 5, stroke: c, "stroke-width": 2, "stroke-dasharray": "5 3" }, gChart);
      else svgEl("rect", { x: lx, y: lyy - 12, width: 14, height: 14, fill: c }, gChart);
      t(lx + 20, lyy, s, { class: "sm" });
      const w = s.length * 8 + 40;
      lx += w;
      if (narrow && lx > CW - 140) { lx = 4; lyy += 22; }
    }
    for (const v of [0, 0.25, 0.5, 0.75, 1]) {
      svgEl("line", { x1: ml2, x2: CW - mr2, y1: y(v), y2: y(v), stroke: COL.grid }, gChart);
      t(ml2 - 8, y(v) + 5, String(v), { class: "sm", "text-anchor": "end" });
    }
    svgEl("line", { x1: ml2, x2: CW - mr2, y1: y(0.5), y2: y(0.5), stroke: COL.chance, "stroke-width": 2, "stroke-dasharray": "7 5" }, gChart);
    const gw = (CW - ml2 - mr2) / 3;
    M.SPLITS.forEach((s, g) => {
      const r = cur!.res[s];
      const gx = ml2 + g * gw;
      const vals = [r.acc1, r.aucKnn, r.aucGc];
      const bwid = Math.min(70, gw / 4.2);
      const x0 = gx + gw / 2 - 1.5 * bwid - 6;
      vals.forEach((v, i) => {
        const x = x0 + i * (bwid + 6);
        svgEl("rect", { x, y: y(v), width: bwid, height: yBot - y(v), fill: [COL.acc, COL.knn, COL.gc][i], "fill-opacity": s === show ? 1 : 0.85 }, gChart);
        t(x + bwid / 2, y(v) - 6, f2(v), { class: "num", "text-anchor": "middle", ...(narrow ? { style: "font-size:15px" } : {}) });
      });
      const cx = gx + gw / 2;
      const e = t(cx, yBot + 24, narrow ? s : SPLIT_LABEL[s], { "text-anchor": "middle", "font-weight": 700 });
      if (s === show) e.style.fill = COL.test;
      t(cx, yBot + 45, narrow ? `n = ${r.nTest}` : `${r.nTest} test windows`, { class: "sm", "text-anchor": "middle" });
      if (narrow) {
        t(cx, yBot + 66, `same ${pct(r.nnSame)}`, { class: "sm", "text-anchor": "middle" });
        t(cx, yBot + 86, `copy ${pct(r.nnHomolog)}`, { class: "sm", "text-anchor": "middle" });
      } else {
        t(cx, yBot + 66, `nearest training window:`, { class: "sm", "text-anchor": "middle" });
        t(cx, yBot + 86, `same region ${pct(r.nnSame)} · copy ${pct(r.nnHomolog)}`, { class: "sm", "text-anchor": "middle" });
      }
    });
  }

  cur = M.runAll(p);
  drawMap(); drawChart();
}

autoMount(NAME, mount);
