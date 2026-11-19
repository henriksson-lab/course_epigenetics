// A/B compartments: paint the true A/B pattern of a toy chromosome, simulate its Hi-C map
// (distance decay × compartment preference × noise), then follow the recipe: observed/expected,
// correlation of rows, first eigenvector by power iteration, sign chosen with GC content.
//
// <div data-demo="compartments"
//      data-n="48"           number of bins (24..64)
//      data-pattern="blocks" first pattern: blocks | fine | island | random
//      data-seed="4"></div>
//
// The model and the recipe are in model.ts (no DOM).
import { svgEl, htmlEl, button, presets, checkbox, slider, rng, autoMount, svgPoint, Slider } from "../lib/ui";
import { analyse, normals, pattern, gcTrack, Mat, Result } from "./model";

const NAME = "compartments";
const C_A = "#2e8b57", C_B = "#a8324a", C_TEXT = "#1b2a3a", C_MUTED = "#5b6b7a";

const STYLE = `
[data-demo="${NAME}"].qb-demo { display: block; }
[data-demo="${NAME}"] .cp-svg { width: 100%; height: auto; display: block; touch-action: none; }
[data-demo="${NAME}"] .cp-bottom { display: flex; flex-wrap: wrap; gap: .5em 1.4em; align-items: flex-start; margin-top: .3em; }
[data-demo="${NAME}"] .cp-col { display: flex; flex-direction: column; gap: .4em; flex: 1 1 340px; min-width: 0; }
[data-demo="${NAME}"] .cp-row { display: flex; flex-wrap: wrap; gap: .35em .9em; align-items: center; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 9.6em 1fr 2.8em; }
[data-demo="${NAME}"] .qb-readout { font-size: 15px; min-height: 7.2em; }
[data-demo="${NAME}"] image { image-rendering: pixelated; image-rendering: crisp-edges; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; margin: 0; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

type RGB = [number, number, number];
const REDS: RGB[] = [[255, 255, 255], [254, 224, 210], [252, 146, 114], [222, 45, 38], [120, 8, 16]];
function ramp(stops: RGB[], t: number): RGB {
  t = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const k = Math.min(stops.length - 2, Math.floor(t)), f = t - k;
  const a = stops[k], b = stops[k + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}
function diverging(t: number): RGB {           // -1 blue, 0 white, +1 red
  t = Math.max(-1, Math.min(1, t));
  const c: RGB = t >= 0 ? [178, 24, 43] : [33, 102, 172], f = Math.abs(t);
  return [255 + f * (c[0] - 255), 255 + f * (c[1] - 255), 255 + f * (c[2] - 255)];
}

/** The shared CSS sets every svg text to the foreground colour; keep the colours given as fill. */
function keepTextFill(root: Element): void {
  root.querySelectorAll<SVGTextElement>("text[fill]").forEach(t => { t.style.fill = t.getAttribute("fill") ?? ""; });
}

function mount(host: HTMLElement): void {
  injectStyle();
  const n = Math.max(24, Math.min(64, Number(host.dataset.n ?? 48) || 48));
  const seed0 = Number(host.dataset.seed ?? 4) || 4;
  const pat0 = host.dataset.pattern ?? "blocks";

  // ---- state ----
  let isA = pattern(pat0, n, seed0);
  let ratio = 2, noise = 0.25, alpha = 1;
  let skipOE = false, orient = true;
  let noiseSeed = seed0, startSeed = seed0;
  let z: Mat = normals(n, rng(noiseSeed).unif);
  let gcNoise: number[] = [];
  let v0: number[] = [];
  let res: Result;
  const drawGcNoise = () => { const r = rng(noiseSeed + 7); gcNoise = Array.from({ length: n }, () => r.norm()); };
  const drawStart = () => { const r = rng(startSeed + 31); v0 = Array.from({ length: n }, () => r.norm()); };
  const gc = () => gcTrack(isA, gcNoise);

  // ---- SVG ----
  const W = 1000, P = 300, gap = 46, top = 34;
  const px = [0, P + gap, 2 * (P + gap)].map(x => x + 4);   // 4 .. 996
  const LW = 140, TW = W - LW - 6, bw = TW / n;
  const yTrue = top + P + 30, hTrue = 30;
  const yGC = yTrue + hTrue + 12, hGC = 46;
  const yE1 = yGC + hGC + 16, hE1 = 92;
  const yCall = yE1 + hE1 + 10, hCall = 26;
  const H = yCall + hCall + 8;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, class: "cp-svg", role: "img",
    "aria-label": "Simulated contact map, observed/expected, correlation matrix, and the first eigenvector with A/B calls" }, host);
  const titles = ["observed (log scale)", "observed / expected", "correlation (of O/E)"];
  const images: SVGImageElement[] = [];
  const titleEls: SVGTextElement[] = [];
  const off = document.createElement("canvas");
  off.width = n; off.height = n;
  const offCtx = off.getContext("2d")!;
  for (let k = 0; k < 3; k++) {
    titleEls.push(svgEl("text", { x: px[k] + P / 2, y: top - 10, "text-anchor": "middle", "font-size": 21, "font-weight": 600, fill: C_TEXT }, svg));
    titleEls[k].textContent = titles[k];
    images.push(svgEl("image", { x: px[k], y: top, width: P, height: P, preserveAspectRatio: "none" }, svg));
    svgEl("rect", { x: px[k], y: top, width: P, height: P, fill: "none", stroke: "#999" }, svg);
    if (k < 2) svgEl("text", { x: px[k] + P + gap / 2, y: top + P / 2 + 8, "text-anchor": "middle", "font-size": 26, fill: C_MUTED }, svg).textContent = "→";
  }
  const gTracks = svgEl("g", {}, svg);
  const gTrue = svgEl("g", { style: "cursor: pointer" }, svg);

  // ---- controls ----
  const bottom = htmlEl("div", "cp-bottom", "", host);
  const colL = htmlEl("div", "cp-col", "", bottom);
  const colR = htmlEl("div", "cp-col", "", bottom);
  const patRow = htmlEl("div", "cp-row", "", colL);
  htmlEl("span", "qb-note", "pattern:", patRow);
  const patBtns = presets(patRow, ["blocks", "fine", "island", "random"].map(k => ({
    label: k === "island" ? "one A island" : k, apply: () => { isA = pattern(k, n, seed0); update(); } })));
  const markPat = (k: string | null) => patBtns.forEach((b, i) => b.classList.toggle("active", ["blocks", "fine", "island", "random"][i] === k));
  const sR: Slider = slider(colL, "A/B preference r", 1, 4, 0.1, ratio, v => { ratio = v; update(); }, v => v.toFixed(1));
  const sN: Slider = slider(colL, "noise", 0, 1.5, 0.05, noise, v => { noise = v; update(); }, v => v.toFixed(2));
  const sA: Slider = slider(colL, "decay exponent", 0, 2, 0.1, alpha, v => { alpha = v; update(); }, v => v.toFixed(1));
  const r2 = htmlEl("div", "cp-row", "", colL);
  const skipBox = checkbox(r2, "skip O/E", skipOE, v => { skipOE = v; update(); });
  const orientBox = checkbox(r2, "orient E1 with GC", orient, v => { orient = v; update(); });
  const r3 = htmlEl("div", "cp-row", "", colL);
  button(r3, "New start vector", () => { startSeed++; drawStart(); update(); });
  button(r3, "New noise", () => { noiseSeed++; z = normals(n, rng(noiseSeed).unif); drawGcNoise(); update(); });
  button(r3, "Reset", () => reset());
  const readout = htmlEl("div", "qb-readout", "", colR);
  htmlEl("div", "qb-note", "Click or drag on the top track to paint bins A or B. Simulated data: contacts fall as "
    + "distance^−exponent, same-type pairs contact r times more than A–B pairs, times random noise.", colR);

  function reset(): void {
    isA = pattern(pat0, n, seed0); markPat(pat0);
    ratio = 2; noise = 0.25; alpha = 1; sR.set(ratio); sN.set(noise); sA.set(alpha);
    skipOE = false; orient = true; skipBox.checked = false; orientBox.checked = true;
    noiseSeed = seed0; startSeed = seed0;
    z = normals(n, rng(noiseSeed).unif); drawGcNoise(); drawStart();
    update();
  }

  // ---- painting the true pattern ----
  let painting: boolean | null = null;
  const binAt = (ev: PointerEvent) => {
    const [x] = svgPoint(svg, ev);
    return Math.floor((x - LW) / bw);
  };
  gTrue.addEventListener("pointerdown", ev => {
    const b = binAt(ev);
    if (b < 0 || b >= n) return;
    ev.preventDefault();
    (ev.target as Element).setPointerCapture?.(ev.pointerId);
    painting = !isA[b];
    isA[b] = painting; markPat(null); update();
  });
  gTrue.addEventListener("pointermove", ev => {
    if (painting === null) return;
    const b = binAt(ev);
    if (b >= 0 && b < n && isA[b] !== painting) { isA[b] = painting; update(); }
  });
  const stop = () => { painting = null; };
  gTrue.addEventListener("pointerup", stop);
  gTrue.addEventListener("pointercancel", stop);

  // ---- drawing ----
  function matrixImage(M: Mat, color: (v: number) => RGB): string {
    const img = offCtx.createImageData(n, n);
    for (let i = 0; i < n * n; i++) {
      const c = color(M[i]);
      img.data[4 * i] = c[0]; img.data[4 * i + 1] = c[1]; img.data[4 * i + 2] = c[2]; img.data[4 * i + 3] = 255;
    }
    offCtx.putImageData(img, 0, 0);
    return off.toDataURL();
  }

  function update(): void {
    const g = gc();
    res = analyse({ isA, ratio, alpha, noise, z }, g, v0, { skipOE, orient });
    // matrices
    let lmax = -Infinity, lmin = Infinity;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const l = Math.log10(res.O[i * n + j]);
      if (l > lmax) lmax = l;
      if (l < lmin) lmin = l;
    }
    const lo = Math.max(lmin, lmax - 2.5);
    images[0].setAttribute("href", matrixImage(res.O, v => ramp(REDS, (Math.log10(v) - lo) / (lmax - lo))));
    images[1].setAttribute("href", matrixImage(res.OE, v => diverging(Math.log2(v) / 1.5)));
    images[2].setAttribute("href", matrixImage(res.C, v => diverging(v)));
    titleEls[2].textContent = skipOE ? "correlation (of observed)" : "correlation (of O/E)";
    images[1].setAttribute("opacity", skipOE ? "0.25" : "1");
    drawTracks(g);
    drawReadout();
    keepTextFill(svg);
  }

  function drawTracks(g: number[]): void {
    gTracks.replaceChildren(); gTrue.replaceChildren();
    const x = (i: number) => LW + i * bw;
    const lab = (y: number, t: string, sub = "") => {
      svgEl("text", { x: LW - 10, y, "text-anchor": "end", "font-size": 19, fill: C_TEXT }, gTracks).textContent = t;
      if (sub) svgEl("text", { x: LW - 10, y: y + 20, "text-anchor": "end", "font-size": 18, fill: C_MUTED }, gTracks).textContent = sub;
    };
    // true pattern (paintable)
    lab(yTrue + 12, "true type", "click to paint");
    svgEl("rect", { x: LW, y: yTrue - 4, width: TW, height: hTrue + 8, fill: "transparent" }, gTrue);
    for (let i = 0; i < n; i++) {
      svgEl("rect", { x: x(i) + 0.5, y: yTrue, width: bw - 1, height: hTrue, fill: isA[i] ? C_A : C_B }, gTrue);
      if (bw >= 14) svgEl("text", { x: x(i) + bw / 2, y: yTrue + 21, "text-anchor": "middle", "font-size": 18, "font-weight": 700, fill: "#fff", "pointer-events": "none" }, gTrue)
        .textContent = isA[i] ? "A" : "B";
    }
    // GC
    lab(yGC + 30, "GC");
    const gmin = 0.32, gmax = 0.48;
    for (let i = 0; i < n; i++) {
      const h = hGC * Math.max(0, Math.min(1, (g[i] - gmin) / (gmax - gmin)));
      svgEl("rect", { x: x(i) + 1, y: yGC + hGC - h, width: bw - 2, height: h, fill: "#8a96a3" }, gTracks);
    }
    svgEl("line", { x1: LW, x2: LW + TW, y1: yGC + hGC, y2: yGC + hGC, stroke: "#bbb" }, gTracks);
    // E1
    lab(yE1 + hE1 / 2 - 4, "E1", orient ? "oriented" : "as solved");
    const mid = yE1 + hE1 / 2;
    let vmax = 0;
    for (let i = 0; i < n; i++) vmax = Math.max(vmax, Math.abs(res.E1[i]));
    for (let i = 0; i < n; i++) {
      const h = (hE1 / 2 - 2) * res.E1[i] / (vmax || 1);
      svgEl("rect", { x: x(i) + 1, y: Math.min(mid, mid - h), width: bw - 2, height: Math.abs(h), fill: res.E1[i] > 0 ? C_A : C_B }, gTracks);
    }
    svgEl("line", { x1: LW, x2: LW + TW, y1: mid, y2: mid, stroke: C_TEXT, "stroke-width": 1.2 }, gTracks);
    // calls
    lab(yCall + 19, "call");
    for (let i = 0; i < n; i++) {
      const a = res.callA[i];
      svgEl("rect", { x: x(i) + 0.5, y: yCall, width: bw - 1, height: hCall, fill: a ? C_A : C_B, opacity: 0.85 }, gTracks);
      if (a !== isA[i]) svgEl("circle", { cx: x(i) + bw / 2, cy: yCall + hCall / 2, r: 5, fill: "#111" }, gTracks);
    }
  }

  function drawReadout(): void {
    const pct = (v: number) => `${Math.round(100 * v)}%`;
    const nOK = Math.round(res.agree * n);
    const lines = [
      `power iteration: ${res.e1.iters} steps${res.e1.converged ? "" : " (not converged)"}`,
      `λ1 = ${res.e1.lambda.toFixed(1)} (${pct(res.e1.lambda / n)} of all variance), λ2 = ${res.lambda2.toFixed(1)}`,
      `corr(E1, GC) = ${res.rGC >= 0 ? "+" : "−"}${Math.abs(res.rGC).toFixed(2)}` + (orient ? (res.flipped ? "   (solver's sign was flipped)" : "   (solver's sign kept)") : "   (sign not checked)"),
      `calls = true type in ${nOK} of ${n} bins (${pct(res.agree)})` + (nOK === 0 ? ": all swapped!" : "") + " (● wrong)",
      `strength (AA + BB)/(AB + BA) from O/E: ${isFinite(res.strength) ? res.strength.toFixed(2) : "–"}`,
    ];
    readout.textContent = lines.join("\n");
  }

  markPat(pat0);
  reset();
}

autoMount(NAME, mount);
