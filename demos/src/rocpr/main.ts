// Thresholds, ROC and precision-recall under class imbalance (model.ts).
// Top left: simulated classifier scores of open (positive) and closed (negative) windows,
// two normal curves drawn as expected numbers of windows (so a rare class is a small curve),
// with a draggable threshold and draggable handles for the separation (peak of the open
// curve) and the spread (its shoulder). Top right: the confusion matrix at the threshold,
// as expected counts in 10 000 test windows, and precision / recall / FPR / accuracy.
// Bottom: ROC and PR curves with the current point (draggable along the curve), ROC AUC,
// average precision and the PR baseline (fraction open). Optional: a simulated finite
// test set (empirical, tie-aware curves) and a pinned copy of the curves for comparison.
//
// Auto-mounts on <div data-demo="rocpr">. data-* options (all optional):
//   data-d="1.53"  separation of the means (in sd of the closed windows)
//   data-sd="1"    spread of both classes      data-pi="0.46"  fraction of open windows
//   data-t="0.77"  threshold                   data-seed="1"   simulated test set
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount, draggable, svgPoint, rng } from "../lib/ui";
import * as M from "./model";

const NAME = "rocpr";
const COL = {
  pos: "#2e8b57", neg: "#a8324a", ink: "#1b2a3a", muted: "#5b6b7a", grid: "#e3e7eb",
  thr: "#1b2a3a", pin: "#9aa5b1", emp: "#2a6f97", base: "#e67e22", light: "#f4f6f8",
};
const NWIN = 10000;      // confusion matrix: expected counts in this many test windows
const NSIM = 2000;       // simulated finite test set
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .rp-controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(17em, 1fr)); gap: .3em 1.4em; width: 100%; }
[data-demo="${NAME}"] .rp-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] svg.rp-main { width: 100%; height: auto; max-width: 1100px; }
[data-demo="${NAME}"] svg text { font-size: 17px; }
[data-demo="${NAME}"] svg .ttl { font-size: 18px; font-weight: 700; }
[data-demo="${NAME}"] svg .sm { font-size: 16px; fill: ${COL.muted}; }
[data-demo="${NAME}"] svg .big { font-size: 20px; font-weight: 700; }
[data-demo="${NAME}"] svg .num { font-size: 22px; font-weight: 700; }
.reveal [data-demo="${NAME}"] { font-size: 18px; gap: .25em; }
.reveal [data-demo="${NAME}"] .rp-controls { gap: 0 1em; }
.reveal [data-demo="${NAME}"] svg.rp-main { max-height: 430px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const f2 = (v: number): string => v.toFixed(2);
const f3 = (v: number): string => v.toFixed(3);
const pct = (v: number): string => (v < 0.1 ? (100 * v).toFixed(1) : (100 * v).toFixed(0)) + "%";
const num = (s: string | undefined, d: number, lo: number, hi: number): number => {
  const v = Number(s);
  return s !== undefined && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
};

interface Box { x: number; y: number; w: number; h: number }

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const DEF = {
    d: num(ds.d, 1.53, 0, 4), sd: num(ds.sd, 1, 0.5, 2), pi: num(ds.pi, 0.46, 0.01, 0.5), t: num(ds.t, 0.77, -4, 7),
  };
  const st = { ...DEF, seed: num(ds.seed, 1, 1, 1e6), sim: false, sameHeight: false };
  let pinned: { roc: [number, number][]; pr: [number, number][]; label: string } | null = null;
  let sim: { c: M.EmpPoint[]; auc: number; ap: number; n1: number } | null = null;

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "rp-controls", "", host);
  const sl = {
    d: slider(ctr, "separation of means", 0, 4, 0.01, st.d, v => { st.d = v; changed(); }),
    sd: slider(ctr, "spread (sd, both)", 0.5, 2, 0.01, st.sd, v => { st.sd = v; changed(); }),
    pi: slider(ctr, "fraction open", -2, Math.log10(0.5), 0.001, Math.log10(st.pi),
      v => { st.pi = 10 ** v; changed(); }, v => pct(10 ** v)),
    t: slider(ctr, "threshold", -4, 7, 0.01, st.t, v => { st.t = v; draw(); }),
  };
  const row = htmlEl("div", "rp-row", "", host);
  const pre = presets(row, [
    { label: "unit test set: 46% open", apply: () => setPi(0.46) },
    { label: "subsampled: 7.7% open", apply: () => setPi(0.077) },
    { label: "genome-wide: 1% open", apply: () => setPi(0.01) },
    { label: "useless model", apply: () => { st.d = 0; sl.d.set(0); changed(); } },
  ]);
  const row2 = htmlEl("div", "rp-row", "", host);
  const bPin = button(row2, "pin these curves", () => {
    if (pinned) { pinned = null; } else {
      const c = M.curve(M3());
      pinned = { roc: c.map(r => [r.fpr, r.tpr]), pr: c.map(r => [r.tpr, r.precision]), label: `pinned: ${pct(st.pi)} open` };
    }
    draw();
  });
  checkbox(row2, `simulated test set (${NSIM} windows)`, st.sim, v => { st.sim = v; resim(); draw(); });
  const bNew = button(row2, "new sample", () => { st.seed++; resim(); draw(); });
  checkbox(row2, "draw both classes at the same height", st.sameHeight, v => { st.sameHeight = v; draw(); });
  button(row2, "reset", () => {
    Object.assign(st, DEF); pinned = null;
    sl.d.set(st.d); sl.sd.set(st.sd); sl.pi.set(Math.log10(st.pi)); sl.t.set(st.t);
    pre.forEach(b => b.classList.remove("active"));
    changed();
  });

  function setPi(p: number): void { st.pi = p; sl.pi.set(Math.log10(p)); changed(); }
  const M3 = (): M.Binormal => ({ d: st.d, sd: st.sd, pi: st.pi });

  // ---------------- layout ----------------
  // layout by the width available, so that text keeps about its nominal size
  const hw = host.clientWidth;
  const narrow = hw > 0 && hw < 640, medium = !narrow && hw > 0 && hw < 1000;
  const W = narrow ? 460 : medium ? 760 : 1100;
  let L: { dist: Box; cm: Box; roc: Box; pr: Box; H: number };
  if (narrow) {
    L = {
      dist: { x: 0, y: 0, w: 460, h: 310 }, cm: { x: 0, y: 320, w: 460, h: 300 },
      roc: { x: 0, y: 630, w: 460, h: 460 }, pr: { x: 0, y: 1100, w: 460, h: 460 }, H: 1565,
    };
  } else if (medium) {
    L = {
      dist: { x: 0, y: 0, w: 760, h: 300 }, cm: { x: 0, y: 306, w: 760, h: 196 },
      roc: { x: 0, y: 512, w: 380, h: 430 }, pr: { x: 380, y: 512, w: 380, h: 430 }, H: 945,
    };
  } else {
    L = {
      dist: { x: 0, y: 0, w: 690, h: 310 }, cm: { x: 710, y: 0, w: 390, h: 310 },
      roc: { x: 0, y: 330, w: 540, h: 500 }, pr: { x: 560, y: 330, w: 540, h: 500 }, H: 835,
    };
  }
  const svg = svgEl("svg", {
    class: "rp-main", viewBox: `0 0 ${W} ${L.H}`, role: "img",
    "aria-label": "Score distributions with a threshold, the confusion matrix, and ROC and precision-recall curves",
  }, host);
  svgEl("rect", { x: 0, y: 0, width: W, height: L.H, fill: "#fff" }, svg);
  const note = htmlEl("div", "qb-note", "", host);
  note.textContent = "Simulated scores: each class is a normal distribution (binormal model). " +
    "Curves and the confusion matrix are exact expectations for an infinitely large test set; " +
    `tick "simulated test set" to see what ${NSIM} windows give.`;

  // ----- distribution panel -----
  const D = L.dist;
  const dm = { l: 14, r: 14, t: 66, b: 52 };
  const X0 = -5, X1 = 8;
  const dx = scale(X0, X1, D.x + dm.l, D.x + D.w - dm.r);
  const dBase = D.y + D.h - dm.b, dTop = D.y + dm.t;
  const gDist = svgEl("g", {}, svg);
  const gDistTop = svgEl("g", {}, svg);
  const thrLine = svgEl("line", { stroke: COL.thr, "stroke-width": 3 }, gDistTop);
  const thrHit = svgEl("rect", { width: 28, fill: "transparent" }, gDistTop);
  const thrHandle = svgEl("circle", { r: 11, class: "handle", stroke: COL.thr }, gDistTop);
  const peakHandle = svgEl("circle", { r: 10, class: "handle", stroke: COL.pos }, gDistTop);
  const sdHandle = svgEl("rect", { width: 16, height: 16, rx: 3, class: "handle", stroke: COL.pos }, gDistTop);
  const toW = (ev: PointerEvent): [number, number] => svgPoint(svg, ev);
  const clampT = (v: number): number => Math.min(7, Math.max(-4, Math.round(v * 100) / 100));
  const setT = (v: number): void => { st.t = clampT(v); sl.t.set(st.t); draw(); };
  for (const h of [thrHandle, thrHit, thrLine]) draggable(h, toW, x => setT(dx.invert(x)));
  draggable(peakHandle, toW, x => {
    st.d = Math.min(4, Math.max(0, Math.round(dx.invert(x) * 100) / 100)); sl.d.set(st.d); changed();
  });
  draggable(sdHandle, toW, x => {
    st.sd = Math.min(2, Math.max(0.5, Math.round((dx.invert(x) - st.d) * 100) / 100)); sl.sd.set(st.sd); changed();
  });

  // ----- confusion matrix panel -----
  const gCm = svgEl("g", {}, svg);

  // ----- ROC / PR panels -----
  const pm = { l: 62, r: 16, t: 58, b: 56 };
  function frame(B: Box) {
    const side = Math.min(B.w - pm.l - pm.r, B.h - pm.t - pm.b);
    const x0 = B.x + pm.l, y0 = B.y + pm.t;
    return { sx: scale(0, 1, x0, x0 + side), sy: scale(0, 1, y0 + side, y0), side, x0, y0 };
  }
  const fR = frame(L.roc), fP = frame(L.pr);
  const gRoc = svgEl("g", {}, svg), gPr = svgEl("g", {}, svg);
  const rocDot = svgEl("circle", { r: 10, class: "handle", stroke: COL.thr }, svg);
  const prDot = svgEl("circle", { r: 10, class: "handle", stroke: COL.thr }, svg);
  draggable(rocDot, toW, (x, y) => setT(M.nearestThreshold(M3(), fR.sx.invert(x), fR.sy.invert(y), "roc")));
  draggable(prDot, toW, (x, y) => setT(M.nearestThreshold(M3(), fP.sx.invert(x), fP.sy.invert(y), "pr")));

  function axes(g: SVGGElement, f: ReturnType<typeof frame>, xl: string, yl: string, title: string, B: Box): void {
    svgEl("text", { x: B.x + 4, y: B.y + 22, class: "ttl" }, g).textContent = title;
    for (const v of [0, 0.25, 0.5, 0.75, 1]) {
      svgEl("line", { x1: f.sx(0), x2: f.sx(1), y1: f.sy(v), y2: f.sy(v), stroke: COL.grid }, g);
      svgEl("line", { x1: f.sx(v), x2: f.sx(v), y1: f.sy(0), y2: f.sy(1), stroke: COL.grid }, g);
      svgEl("text", { x: f.sx(v), y: f.sy(0) + 22, "text-anchor": "middle", class: "sm" }, g).textContent = String(v);
      svgEl("text", { x: f.sx(0) - 8, y: f.sy(v) + 5, "text-anchor": "end", class: "sm" }, g).textContent = String(v);
    }
    svgEl("rect", { x: f.sx(0), y: f.sy(1), width: f.side, height: f.side, fill: "none", stroke: COL.muted }, g);
    svgEl("text", { x: f.sx(0.5), y: f.sy(0) + 46, "text-anchor": "middle" }, g).textContent = xl;
    const yc = f.sy(0.5), xc = B.x + 16;
    svgEl("text", { x: xc, y: yc, "text-anchor": "middle", transform: `rotate(-90 ${xc} ${yc})` }, g).textContent = yl;
  }
  const poly = (pts: [number, number][], f: ReturnType<typeof frame>, attrs: Record<string, string | number>, g: Element) =>
    svgEl("polyline", { points: pts.map(([a, b]) => `${f.sx(a).toFixed(1)},${f.sy(b).toFixed(1)}`).join(" "), fill: "none", ...attrs }, g);
  // empirical curves as steps (ties: diagonal), starting at (0,0) / (0, precision at first threshold)
  const empRoc = (c: M.EmpPoint[]): [number, number][] => [[0, 0], ...c.map(p => [p.fpr, p.tpr] as [number, number])];
  const empPr = (c: M.EmpPoint[]): [number, number][] => {
    const out: [number, number][] = [[0, c[0].precision]];
    let r0 = 0;
    for (const p of c) { out.push([r0, p.precision], [p.tpr, p.precision]); r0 = p.tpr; }
    return out;
  };

  function resim(): void {
    if (!st.sim) { sim = null; bNew.disabled = true; return; }
    bNew.disabled = false;
    // a fixed stream per seed (independent of the model), so moving a slider shifts the same windows
    const { norm } = rng(st.seed * 7919 + 17);
    const s = M.simulate(M3(), NSIM, norm);
    const c = M.empiricalCurve(s.score, s.label);
    sim = { c, auc: M.aucTrapezoid(c), ap: M.averagePrecision(c), n1: s.label.reduce((x, y) => x + y, 0) };
  }

  // ---------------- drawing ----------------
  function drawDist(): void {
    gDist.replaceChildren();
    const m = M3();
    svgEl("text", { x: D.x + 4, y: D.y + 22, class: "ttl" }, gDist).textContent = "Scores of the test windows (simulated)";
    // densities as expected numbers of windows per unit score
    const wPos = st.sameHeight ? 0.5 : m.pi, wNeg = st.sameHeight ? 0.5 : 1 - m.pi;
    const peak = Math.max(wPos, wNeg) * M.dnorm(0, 0, m.sd);
    const dy = scale(0, peak * 1.08, dBase, dTop);
    const N = 260;
    const xs = Array.from({ length: N + 1 }, (_, i) => X0 + (X1 - X0) * i / N);
    const area = (mu: number, w: number, from: number, to: number, fill: string, op: number) => {
      const a = Math.max(X0, from), b = Math.min(X1, to);
      if (b <= a) return;
      const pts: string[] = [`${dx(a)},${dBase}`];
      const k = Math.max(2, Math.ceil((b - a) / (X1 - X0) * N));
      for (let i = 0; i <= k; i++) { const x = a + (b - a) * i / k; pts.push(`${dx(x).toFixed(1)},${dy(w * M.dnorm(x, mu, m.sd)).toFixed(1)}`); }
      pts.push(`${dx(b)},${dBase}`);
      svgEl("polygon", { points: pts.join(" "), fill, "fill-opacity": op }, gDist);
    };
    // closed: TN left, FP right; open: FN left, TP right
    area(0, wNeg, X0, st.t, COL.neg, 0.10); area(0, wNeg, st.t, X1, COL.neg, 0.38);
    area(m.d, wPos, X0, st.t, COL.pos, 0.10); area(m.d, wPos, st.t, X1, COL.pos, 0.42);
    const line = (mu: number, w: number, col: string) => svgEl("polyline", {
      points: xs.map(x => `${dx(x).toFixed(1)},${dy(w * M.dnorm(x, mu, m.sd)).toFixed(1)}`).join(" "),
      fill: "none", stroke: col, "stroke-width": 2.5,
    }, gDist);
    line(0, wNeg, COL.neg); line(m.d, wPos, COL.pos);
    // axis
    svgEl("line", { x1: dx(X0), x2: dx(X1), y1: dBase, y2: dBase, stroke: COL.ink, "stroke-width": 1.5 }, gDist);
    for (let v = -4; v <= 8; v += 2) {
      svgEl("line", { x1: dx(v), x2: dx(v), y1: dBase, y2: dBase + 5, stroke: COL.ink }, gDist);
      svgEl("text", { x: dx(v), y: dBase + 22, "text-anchor": "middle", class: "sm" }, gDist).textContent = String(v);
    }
    svgEl("text", { x: dx((X0 + X1) / 2), y: dBase + 45, "text-anchor": "middle" }, gDist).textContent = "score (higher = more likely open)";
    // class labels next to their curves
    const lab = (x: number, y: number, txt: string, col: string, anchor = "middle") => {
      const t = svgEl("text", { x, y, "text-anchor": anchor, "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, gDist);
      t.style.fill = col; t.textContent = txt;
    };
    const hiY = (mu: number, w: number) => Math.min(dBase - 8, dy(w * M.dnorm(mu + 0.7 * m.sd, mu, m.sd)));
    lab(dx(-1.2 * m.sd) - 10, hiY(0, wNeg), `closed (${pct(1 - m.pi)})`, COL.neg, "end");
    const px = dx(m.d + 1.2 * m.sd) + 10;
    lab(Math.min(px, dx(X1) - 4), hiY(m.d, wPos), `open (${pct(m.pi)})`, COL.pos, px > dx(X1) - 130 ? "end" : "start");
    // handles: threshold, peak of open curve (separation), shoulder (spread)
    const xt = dx(Math.max(X0, Math.min(X1, st.t)));
    setAttrsN(thrLine, { x1: xt, x2: xt, y1: dTop - 6, y2: dBase });
    setAttrsN(thrHit, { x: xt - 14, y: dTop - 6, height: dBase - dTop + 6 });
    setAttrsN(thrHandle, { cx: xt, cy: dTop - 6 });
    const yp = dy(wPos * M.dnorm(m.d, m.d, m.sd));
    setAttrsN(peakHandle, { cx: dx(m.d), cy: Math.min(dBase - 4, yp) });
    const ys = dy(wPos * M.dnorm(m.d + m.sd, m.d, m.sd));
    setAttrsN(sdHandle, { x: dx(m.d + m.sd) - 8, y: Math.min(dBase - 12, ys) - 8 });
    // arrow "called open"
    const ya = dTop - 1;
    if (xt < dx(X1) - 220) {
      svgEl("text", { x: xt + 16, y: ya, class: "sm" }, gDist).textContent = `threshold ${f2(st.t)}: called open →`;
    } else {
      svgEl("text", { x: xt - 16, y: ya, "text-anchor": "end", class: "sm" }, gDist).textContent = `threshold ${f2(st.t)}`;
    }
  }

  function setAttrsN(e: Element, a: Record<string, number>): void {
    for (const [k, v] of Object.entries(a)) e.setAttribute(k, v.toFixed(1));
  }

  function drawCm(r: M.Rates): void {
    gCm.replaceChildren();
    const B = L.cm, m = M3();
    const P = m.pi * NWIN, Nn = (1 - m.pi) * NWIN;
    const TP = r.tpr * P, FN = P - TP, FP = r.fpr * Nn, TN = Nn - FP;
    svgEl("text", { x: B.x + 4, y: B.y + 22, class: "ttl" }, gCm).textContent = `Confusion matrix (${NWIN.toLocaleString("en")} windows)`;
    const cx = B.x + 96, cy = B.y + 66, cw = 128, ch = 50;
    const head = (x: number, y: number, t: string, anchor = "middle") =>
      (svgEl("text", { x, y, "text-anchor": anchor, class: "sm" }, gCm).textContent = t);
    head(cx + cw / 2, cy - 10, "called open"); head(cx + 1.5 * cw, cy - 10, "called closed");
    head(cx - 8, cy + ch / 2 + 6, "open", "end"); head(cx - 8, cy + 1.5 * ch + 6, "closed", "end");
    const cell = (i: number, j: number, v: number, name: string, col: string, op: number) => {
      svgEl("rect", { x: cx + j * cw, y: cy + i * ch, width: cw - 3, height: ch - 3, rx: 4, fill: col, "fill-opacity": op }, gCm);
      const t = svgEl("text", { x: cx + j * cw + cw / 2 - 1, y: cy + i * ch + ch / 2 + 7, "text-anchor": "middle", class: "big" }, gCm);
      t.textContent = `${name} ${Math.round(v).toLocaleString("en")}`;
    };
    cell(0, 0, TP, "TP", COL.pos, 0.42); cell(0, 1, FN, "FN", COL.pos, 0.10);
    cell(1, 0, FP, "FP", COL.neg, 0.38); cell(1, 1, TN, "TN", COL.neg, 0.10);
    const lines: [string, string, string][] = [
      ["precision", f3(r.precision), "TP/(TP+FP)"],
      ["recall = TPR", f3(r.tpr), "TP/(TP+FN)"],
      ["FPR", f3(r.fpr), "FP/(FP+TN)"],
      ["accuracy", f3(r.accuracy), `all closed: ${f3(1 - m.pi)}`],
    ];
    // metrics below the matrix, or beside it when the panel is wide
    const beside = B.w >= 700;
    const x0 = beside ? B.x + 366 : B.x + 8;
    let y = beside ? cy + 6 : cy + 2 * ch + 30;
    for (const [a, b, c] of lines) {
      svgEl("text", { x: x0, y, "font-weight": 600 }, gCm).textContent = a;
      svgEl("text", { x: x0 + 184, y, "text-anchor": "end", "font-weight": 700 }, gCm).textContent = b;
      svgEl("text", { x: x0 + 196, y, class: "sm" }, gCm).textContent = c;
      y += 25;
    }
    if (r.accuracy < 1 - m.pi + 1e-9) {
      const t = svgEl("text", { x: x0, y: y + 2, "font-weight": 700 }, gCm);
      t.style.fill = "#c0392b"; t.textContent = "worse than calling everything closed!";
    }
  }

  function drawCurves(r: M.Rates): void {
    const m = M3(), c = M.curve(m);
    // ROC
    gRoc.replaceChildren();
    axes(gRoc, fR, "false positive rate (FPR)", "true positive rate (recall)", "ROC curve", L.roc);
    poly([[0, 0], [1, 1]], fR, { stroke: COL.muted, "stroke-width": 1.5, "stroke-dasharray": "6 5" }, gRoc);
    if (sim) poly(empRoc(sim.c), fR, { stroke: COL.emp, "stroke-width": 1.6 }, gRoc);
    poly(c.map(p => [p.fpr, p.tpr]), fR, { stroke: COL.ink, "stroke-width": 3 }, gRoc);
    if (pinned) {
      poly(pinned.roc, fR, { stroke: COL.pin, "stroke-width": 4, "stroke-dasharray": "2 7", "stroke-linecap": "round" }, gRoc);
      const t = svgEl("text", { x: fR.sx(0.97), y: fR.sy(0.05), "text-anchor": "end", "font-weight": 600 }, gRoc);
      t.style.fill = COL.muted; t.textContent = `··· ${pinned.label}`;
    }
    const tx = fR.sx(1), ty = L.roc.y + 22;
    const auc = svgEl("text", { x: tx, y: ty, "text-anchor": "end", class: "num" }, gRoc);
    auc.textContent = `AUC = ${f3(M.rocAuc(m))}`;
    if (sim) {
      const t = svgEl("text", { x: tx, y: ty + 24, "text-anchor": "end", "font-weight": 600 }, gRoc);
      t.style.fill = COL.emp; t.textContent = `sample of ${NSIM}: ${f3(sim.auc)}`;
    }
    setAttrsN(rocDot, { cx: fR.sx(r.fpr), cy: fR.sy(r.tpr) });
    // PR
    gPr.replaceChildren();
    axes(gPr, fP, "recall (TPR)", "precision", "Precision–recall", L.pr);
    const yb = fP.sy(m.pi);
    svgEl("line", { x1: fP.sx(0), x2: fP.sx(1), y1: yb, y2: yb, stroke: COL.base, "stroke-width": 2.5, "stroke-dasharray": "8 5" }, gPr);
    const ly = L.pr.y + 46, lx = L.pr.x + 4;
    svgEl("line", { x1: lx, x2: lx + 30, y1: ly - 6, y2: ly - 6, stroke: COL.base, "stroke-width": 2.5, "stroke-dasharray": "8 5" }, gPr);
    const bl = svgEl("text", { x: lx + 38, y: ly, "font-weight": 600 }, gPr);
    bl.style.fill = COL.base; bl.textContent = `random: ${f3(m.pi)}`;
    if (sim) poly(empPr(sim.c), fP, { stroke: COL.emp, "stroke-width": 1.6 }, gPr);
    poly(c.map(p => [p.tpr, p.precision]), fP, { stroke: COL.ink, "stroke-width": 3 }, gPr);
    const ap = svgEl("text", { x: fP.sx(1), y: L.pr.y + 22, "text-anchor": "end", class: "num" }, gPr);
    ap.textContent = `AP = ${f3(M.avgPrecision(m))}`;
    if (sim) {
      const t = svgEl("text", { x: fP.sx(1), y: L.pr.y + 46, "text-anchor": "end", "font-weight": 600 }, gPr);
      t.style.fill = COL.emp; t.textContent = `sample (${sim.n1} open): ${f3(sim.ap)}`;
    }
    if (pinned) {
      const t = svgEl("text", { x: fP.sx(0) + 8, y: fP.sy(1) + 24, "font-weight": 600, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, gPr);
      t.style.fill = COL.muted; t.textContent = `··· ${pinned.label}`;
    }
    if (pinned) poly(pinned.pr, fP, { stroke: COL.pin, "stroke-width": 4, "stroke-dasharray": "2 7", "stroke-linecap": "round" }, gPr);
    setAttrsN(prDot, { cx: fP.sx(r.tpr), cy: fP.sy(r.precision) });
    bPin.textContent = pinned ? "unpin curves" : "pin these curves";
  }

  function changed(): void { resim(); draw(); }

  function draw(): void {
    const r = M.rates(M3(), st.t);
    drawDist(); drawCm(r); drawCurves(r);
  }

  resim();
  draw();
}

autoMount(NAME, mount);
