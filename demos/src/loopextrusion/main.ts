// Loop extrusion: cohesin extrudes loops on a 200-bin chain until it meets a CTCF site that
// points into the loop. Many simulated cells are averaged into a Hi-C-like contact map (a
// triangle above the chain), and its distance decay P(s) is plotted.
//
// <div data-demo="loopextrusion"
//      data-maxdist="100"   largest distance (bins) shown in the triangle
//      data-preset="domains" first preset: domains | flip | delete | none | nocohesin
//      data-seed="3"></div>
//
// The model is in sim.ts. Every change of a control restarts the simulation from the same
// seed, so the same settings always give the same map.
import { svgEl, htmlEl, button, presets, checkbox, slider, rng, autoMount, svgPoint, Slider } from "../lib/ui";
import { Population, Params, Site, Dir, contactDecay, siteArray, sitePairs } from "./sim";

const NAME = "loopextrusion";
const N = 200, KB = 5;                 // bins, kb per bin
const CELLS = 30, SNAP_EVERY = 5, SNAP_ROUNDS = 40;
const C_CTCF = "#8c564b", C_COH = "#2a9d8f", C_TEXT = "#1b2a3a", C_MUTED = "#5b6b7a";

const STYLE = `
[data-demo="${NAME}"].qb-demo { display: block; }
[data-demo="${NAME}"] .le-fig { position: relative; width: 100%; }
[data-demo="${NAME}"] canvas { display: block; width: 100%; height: auto; cursor: crosshair; }
[data-demo="${NAME}"] .le-track { display: block; width: 100%; height: auto; cursor: pointer; }
[data-demo="${NAME}"] .le-bottom { display: flex; flex-wrap: wrap; gap: .6em 1.4em; align-items: flex-start; margin-top: .3em; }
[data-demo="${NAME}"] .le-col { display: flex; flex-direction: column; gap: .4em; flex: 1 1 330px; min-width: 0; }
[data-demo="${NAME}"] .le-ps { flex: 1 1 330px; min-width: 0; max-width: 460px; }
[data-demo="${NAME}"] .le-ps svg { width: 100%; height: auto; }
[data-demo="${NAME}"] .le-row { display: flex; flex-wrap: wrap; gap: .35em .8em; align-items: center; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 10.5em 1fr 2.6em; }
[data-demo="${NAME}"] .qb-readout { min-height: 9.5em; font-size: 15px; }
[data-demo="${NAME}"] .le-status { font-size: .85em; color: var(--qb-muted); min-height: 1.3em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; margin: 0; }
.reveal [data-demo="${NAME}"] .qb-readout { min-height: 0; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

interface Setup { label: string; key: string; sites: Site[]; load: number; life: number; pStop: number; oneSided: boolean }
const S = (pos: number, dir: Dir): Site => ({ pos, dir });
const BASE_SITES = [S(30, 1), S(80, -1), S(82, 1), S(130, -1), S(160, 1)];
const PRESETS: Setup[] = [
  { label: "Two domains + stripe", key: "domains", sites: BASE_SITES, load: 20, life: 60, pStop: 0.9, oneSided: false },
  { label: "Flip one site", key: "flip", sites: [S(30, 1), S(80, 1), S(82, 1), S(130, -1), S(160, 1)], load: 20, life: 60, pStop: 0.9, oneSided: false },
  { label: "Delete one site", key: "delete", sites: [S(30, 1), S(82, 1), S(130, -1), S(160, 1)], load: 20, life: 60, pStop: 0.9, oneSided: false },
  { label: "No CTCF", key: "none", sites: [], load: 20, life: 60, pStop: 0.9, oneSided: false },
  { label: "No cohesin", key: "nocohesin", sites: BASE_SITES, load: 0, life: 60, pStop: 0.9, oneSided: false },
];

// colour maps
const REDS: [number, number, number][] = [[255, 255, 255], [254, 224, 210], [252, 146, 114], [222, 45, 38], [140, 10, 20], [60, 0, 5]];
function ramp(stops: [number, number, number][], t: number): [number, number, number] {
  t = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const k = Math.min(stops.length - 2, Math.floor(t)), f = t - k;
  const a = stops[k], b = stops[k + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}
function diverging(t: number): [number, number, number] {   // t in [-1, 1]: blue - white - red
  t = Math.max(-1, Math.min(1, t));
  const w: [number, number, number] = [255, 255, 255], r: [number, number, number] = [178, 24, 43], b: [number, number, number] = [33, 102, 172];
  const c = t >= 0 ? r : b, f = Math.abs(t);
  return [w[0] + f * (c[0] - w[0]), w[1] + f * (c[1] - w[1]), w[2] + f * (c[2] - w[2])];
}

const fmtKb = (bins: number) => {
  const kb = bins * KB;
  return kb >= 1000 ? (kb / 1000).toFixed(kb % 1000 === 0 ? 0 : 2) + " Mb" : Math.round(kb) + " kb";
};

/** The shared CSS sets every svg text to the foreground colour; keep the colours given as fill. */
function keepTextFill(root: Element): void {
  root.querySelectorAll<SVGTextElement>("text[fill]").forEach(t => { t.style.fill = t.getAttribute("fill") ?? ""; });
}

function mount(host: HTMLElement): void {
  injectStyle();
  const maxD = Math.max(20, Math.min(N - 1, Number(host.dataset.maxdist ?? 110) || 110));
  const seed = Number(host.dataset.seed ?? 3) || 3;
  const startKey = host.dataset.preset ?? "domains";
  const start = PRESETS.find(p => p.key === startKey) ?? PRESETS[0];

  // ---- state ----
  let sites: Site[] = start.sites.map(s => ({ ...s }));
  let load = start.load, life = start.life, pStop = start.pStop, oneSided = start.oneSided;
  let showOE = false, playing = true;
  let pop: Population | null = null;
  let rounds = 0;                       // snapshot rounds done
  let mean: Float64Array = new Float64Array(N * N);
  let ps: Float64Array = new Float64Array(N);
  let prevPs: Float64Array | null = null, prevLabel = "";
  let pick: [number, number] | null = null;
  let frame = 0;

  // ---- DOM: the map (canvas) and the chain (svg) share the x scale: W units for N bins ----
  const W = 1000, bw = W / N;
  const mapH = Math.ceil(maxD * bw / 2) + 34;
  const fig = htmlEl("div", "le-fig", "", host);
  const canvas = htmlEl("canvas", "", "", fig);
  canvas.width = W; canvas.height = mapH;
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", "Simulated contact map, shown as a triangle above the chain");
  const ctx = canvas.getContext("2d")!;
  const off = document.createElement("canvas");
  off.width = N; off.height = N;
  const offCtx = off.getContext("2d")!;
  const img = offCtx.createImageData(N, N);

  const TH = 118;
  const track = svgEl("svg", { viewBox: `0 0 ${W} ${TH}`, class: "le-track", role: "img",
    "aria-label": "The chain with CTCF sites and the cohesin loops of one cell" }, fig);
  const gArcs = svgEl("g", {}, track);
  const gChain = svgEl("g", {}, track);
  const gSites = svgEl("g", {}, track);
  const gLabels = svgEl("g", {}, track);

  const bottom = htmlEl("div", "le-bottom", "", host);
  const colL = htmlEl("div", "le-col", "", bottom);
  const psBox = htmlEl("div", "le-ps", "", bottom);
  const colR = htmlEl("div", "le-col", "", bottom);

  const presetBtns = presets(colL, PRESETS.map(p => ({ label: p.label, apply: () => applyPreset(p) })));
  const markPreset = (p: Setup | null) => presetBtns.forEach((b, i) => b.classList.toggle("active", PRESETS[i] === p));
  const onSlide = () => {
    load = sLoad.get(); life = sLife.get(); pStop = sStop.get();
    markPreset(null); restart();
  };
  const sLoad: Slider = slider(colL, "loads per 100 steps", 0, 80, 5, load, onSlide, v => String(v));
  const sLife: Slider = slider(colL, "lifetime (steps)", 5, 150, 5, life, onSlide, v => String(v));
  const sStop: Slider = slider(colL, "CTCF stop prob.", 0, 1, 0.05, pStop, onSlide, v => v.toFixed(2));
  const r2 = htmlEl("div", "le-row", "", colL);
  const oneBox = checkbox(r2, "one-sided extrusion", oneSided, v => { oneSided = v; markPreset(null); restart(); });
  const oeBox = checkbox(r2, "show obs/exp", showOE, v => { showOE = v; drawMap(); });
  const r3 = htmlEl("div", "le-row", "", colL);
  const playBtn = button(r3, "Pause", () => { playing = !playing; playBtn.textContent = playing ? "Pause" : "Play"; });
  button(r3, "Reset", () => { applyPreset(start); });

  const psSvg = svgEl("svg", { viewBox: "0 0 440 330", role: "img", "aria-label": "Contact probability against distance, log-log" }, psBox);

  const readout = htmlEl("div", "qb-readout", "", colR);
  const status = htmlEl("div", "le-status", "", colR);
  htmlEl("div", "qb-note", "Click the chain: an empty spot adds a → site; a → site turns into ←; a ← site is removed. "
    + "Click the map to read one pixel. Simulated data: one step moves a free leg by one 5-kb bin.", colR);

  // ---- simulation ----
  function params(): Params {
    return { n: N, sites, loadRate: load / 100, lifetime: life, pStop, oneSided, alpha: 1 };
  }
  function restart(): void {
    if (pop && rounds >= SNAP_ROUNDS) { prevPs = ps.slice(); }
    pop = new Population(params(), CELLS, rng(seed).unif);
    pop.run(Math.min(400, 2 * life + 20), 0);    // reach the steady state
    rounds = 0;
    siteArr = siteArray(N, sites);
    drawTrack();
  }
  function applyPreset(p: Setup): void {
    sites = p.sites.map(s => ({ ...s }));
    load = p.load; life = p.life; pStop = p.pStop; oneSided = p.oneSided;
    sLoad.set(load); sLife.set(life); sStop.set(pStop); oneBox.checked = oneSided;
    markPreset(p);
    prevLabel = "";
    restart();
  }

  function tick(): void {
    frame++;
    if (pop && rounds < SNAP_ROUNDS) {
      const t0 = performance.now();
      while (rounds < SNAP_ROUNDS && performance.now() - t0 < 14) { pop.run(SNAP_EVERY, SNAP_EVERY); rounds++; }
      mean = pop.mean();
      ps = contactDecay(N, mean);
      drawMap(); drawPs(); drawReadout(); drawTrack();
    } else if (pop && playing && frame % 3 === 0) {
      // after the map is done, keep showing one cell (cell 1) in motion; the map is not changed
      pop.cells[0].step(pop.p, siteArr, pop.rnd);
      drawTrack();
    }
    requestAnimationFrame(tick);
  }
  let siteArr: Int8Array = siteArray(N, []);

  // ---- drawing: map ----
  const base = mapH - 2;
  function drawMap(): void {
    const d = img.data;
    const lo = -2.1, hi = -0.5;
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      const k = 4 * (i * N + j);
      let c: [number, number, number] = [255, 255, 255];
      if (j > i) {
        const v = mean[i * N + j];
        if (showOE) c = diverging(ps[j - i] > 0 ? Math.log2(v / ps[j - i]) / 1.6 : 0);
        else c = ramp(REDS, v > 0 ? (Math.log10(v) - lo) / (hi - lo) : 0);
      } else if (j === i) c = showOE ? [200, 200, 200] : [60, 0, 5];
      d[k] = c[0]; d[k + 1] = c[1]; d[k + 2] = c[2]; d[k + 3] = 255;
    }
    offCtx.putImageData(img, 0, 0);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, W, mapH);
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, W, mapH);
    ctx.save();
    ctx.beginPath();
    ctx.moveTo(0, base); ctx.lineTo(W, base); ctx.lineTo(W - maxD * bw / 2, base - maxD * bw / 2);
    ctx.lineTo(maxD * bw / 2, base - maxD * bw / 2); ctx.closePath();
    ctx.clip();
    ctx.imageSmoothingEnabled = false;
    // matrix column u = j, row v = i  ->  x = bw/2 (u + v), y = base - bw/2 (u - v)
    ctx.setTransform(bw / 2, -bw / 2, bw / 2, bw / 2, 0, base);
    ctx.drawImage(off, 0, 0);
    ctx.restore();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    // picked pixel
    if (pick) {
      const [i, j] = pick;
      const cx = (i + j + 1) * bw / 2, cy = base - (j - i) * bw / 2;
      ctx.strokeStyle = C_TEXT; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(cx, cy, 9, 0, 2 * Math.PI); ctx.stroke();
    }
    // labels and colour key
    ctx.fillStyle = C_TEXT; ctx.font = "20px system-ui, sans-serif"; ctx.textBaseline = "top";
    ctx.fillText(showOE ? "observed / expected (simulated)" : "contacts, log scale (simulated)", 8, 6);
    const kx = W - 250, ky = 10;
    for (let k = 0; k < 120; k++) {
      const t = k / 119;
      const c = showOE ? diverging(2 * t - 1) : ramp(REDS, t);
      ctx.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
      ctx.fillRect(kx + k, ky, 1.2, 14);
    }
    ctx.strokeStyle = "#999"; ctx.lineWidth = 1; ctx.strokeRect(kx, ky, 120, 14);
    ctx.fillStyle = C_MUTED; ctx.font = "17px system-ui, sans-serif"; ctx.textAlign = "right";
    ctx.fillText(showOE ? "½×" : "low", kx - 6, ky - 1);
    ctx.textAlign = "left";
    ctx.fillText(showOE ? "3×" : "high", kx + 126, ky - 1);
    ctx.textAlign = "left";
  }

  // ---- drawing: chain track ----
  const yChain = 70;
  function drawTrack(): void {
    gArcs.replaceChildren(); gChain.replaceChildren(); gSites.replaceChildren(); gLabels.replaceChildren();
    const x = (i: number) => (i + 0.5) * bw;
    svgEl("line", { x1: x(0), x2: x(N - 1), y1: yChain, y2: yChain, stroke: C_TEXT, "stroke-width": 3 }, gChain);
    // ruler
    for (let b = 0; b <= N; b += 40) {
      svgEl("line", { x1: b * bw, x2: b * bw, y1: yChain + 6, y2: yChain + 12, stroke: C_MUTED, "stroke-width": 1.5 }, gLabels);
      svgEl("text", { x: Math.min(W - 30, Math.max(18, b * bw)), y: yChain + 34, "text-anchor": "middle", "font-size": 19, fill: C_MUTED }, gLabels)
        .textContent = b === 0 ? "0" : fmtKb(b);
    }
    // loops of cell 1: an arc below the chain from leg to leg
    if (pop) {
      for (const c of pop.cells[0].coh) {
        const x1 = x(c.l), x2 = x(c.r), h = Math.min(46, 6 + (c.r - c.l) * 0.7);
        svgEl("path", { d: `M${x1},${yChain - 3} C${x1},${yChain - 3 - h * 1.3} ${x2},${yChain - 3 - h * 1.3} ${x2},${yChain - 3}`,
          fill: "none", stroke: C_COH, "stroke-width": 3 }, gArcs);
        for (const [xx, st] of [[x1, c.lStall], [x2, c.rStall]] as const)
          svgEl("circle", { cx: xx, cy: yChain - 3, r: 5, fill: st ? C_COH : "#fff", stroke: C_COH, "stroke-width": 2.5 }, gArcs);
      }
    }
    // CTCF sites: a triangle pointing in the motif direction
    for (const s of sites) {
      const cx = x(s.pos), y = yChain + 14;
      const d = s.dir === 1 ? `M${cx - 7},${y - 8} L${cx + 9},${y} L${cx - 7},${y + 8} Z` : `M${cx + 7},${y - 8} L${cx - 9},${y} L${cx + 7},${y + 8} Z`;
      svgEl("path", { d, fill: C_CTCF }, gSites);
    }
    svgEl("text", { x: 4, y: 20, "font-size": 19, fill: C_MUTED }, gLabels).textContent =
      pop && pop.cells[0].coh.length ? "one cell: cohesin loops (● = stopped at CTCF)" : "one cell: no cohesin loaded";
    svgEl("text", { x: W - 4, y: 20, "font-size": 19, fill: C_CTCF, "text-anchor": "end" }, gLabels).textContent = "▶ ◀ CTCF motif direction";
    keepTextFill(gLabels);
  }
  track.addEventListener("click", ev => {
    const [px] = svgPoint(track, ev);
    const b = Math.max(0, Math.min(N - 1, Math.floor(px / bw)));
    let best = -1, bd = 4;
    sites.forEach((s, k) => { const dd = Math.abs(s.pos - b); if (dd < bd) { bd = dd; best = k; } });
    if (best < 0) sites = [...sites, { pos: b, dir: 1 as Dir }].sort((a, c) => a.pos - c.pos);
    else if (sites[best].dir === 1) sites = sites.map((s, k) => (k === best ? { pos: s.pos, dir: -1 as Dir } : s));
    else sites = sites.filter((_, k) => k !== best);
    markPreset(null);
    restart();
  });
  canvas.addEventListener("click", ev => {
    const r = canvas.getBoundingClientRect();
    const px = (ev.clientX - r.left) * W / r.width, py = (ev.clientY - r.top) * mapH / r.height;
    const sum = 2 * px / bw - 1, diff = 2 * (base - py) / bw;   // i + j, j - i
    const j = Math.round((sum + diff) / 2), i = Math.round((sum - diff) / 2);
    if (i >= 0 && j < N && j - i >= 1 && j - i <= maxD) { pick = [i, j]; drawMap(); drawReadout(); }
  });

  // ---- drawing: P(s) ----
  function drawPs(): void {
    psSvg.replaceChildren();
    const L = 72, R = 425, T = 34, B = 270;
    const sMax = Math.min(N - 1, 150);
    const lx = (s: number) => L + (R - L) * Math.log10(s) / Math.log10(sMax);
    const yLo = -3, yHi = 0.3;
    const ly = (v: number) => B - (B - T) * (Math.log10(v) - yLo) / (yHi - yLo);
    svgEl("rect", { x: L, y: T, width: R - L, height: B - T, fill: "none", stroke: "#ccc" }, psSvg);
    for (const s of [1, 2, 5, 10, 20, 50, 100]) {
      svgEl("line", { x1: lx(s), x2: lx(s), y1: T, y2: B, stroke: "#eee" }, psSvg);
      if (s === 1 || s === 10 || s === 100)
        svgEl("text", { x: lx(s), y: B + 22, "text-anchor": "middle", "font-size": 18, fill: C_MUTED }, psSvg).textContent = fmtKb(s);
    }
    for (const e of [-3, -2, -1, 0]) {
      svgEl("line", { x1: L, x2: R, y1: ly(10 ** e), y2: ly(10 ** e), stroke: "#eee" }, psSvg);
      svgEl("text", { x: L - 6, y: ly(10 ** e) + 6, "text-anchor": "end", "font-size": 18, fill: C_MUTED }, psSvg).textContent = e === 0 ? "1" : `10${"⁻"}${"⁰¹²³"[-e]}`;
    }
    svgEl("text", { x: (L + R) / 2, y: B + 50, "text-anchor": "middle", "font-size": 19, fill: C_TEXT }, psSvg).textContent = "distance s";
    svgEl("text", { x: 6, y: 22, "font-size": 19, fill: C_TEXT }, psSvg).textContent = "P(s): mean contact at distance s";
    const line = (f: (s: number) => number, attrs: Record<string, string | number>) => {
      let d = "";
      for (let s = 1; s <= sMax; s++) d += (s === 1 ? "M" : "L") + lx(s).toFixed(1) + "," + ly(Math.max(1e-4, f(s))).toFixed(1);
      svgEl("path", { d, fill: "none", ...attrs }, psSvg);
    };
    line(s => 1 / s, { stroke: "#999", "stroke-width": 2.5, "stroke-dasharray": "7 6" });
    if (prevPs) { const pp = prevPs; line(s => pp[s], { stroke: "#e9a1a0", "stroke-width": 3 }); }
    line(s => ps[s], { stroke: "#c0392b", "stroke-width": 3.5 });
    const lg = (y: number, col: string, dash: string, txt: string) => {
      svgEl("line", { x1: R - 150, x2: R - 120, y1: y, y2: y, stroke: col, "stroke-width": 3, "stroke-dasharray": dash }, psSvg);
      svgEl("text", { x: R - 112, y: y + 6, "font-size": 17, fill: C_TEXT }, psSvg).textContent = txt;
    };
    lg(T + 18, "#c0392b", "", "now");
    lg(T + 40, "#999", "7 6", "s⁻¹, no cohesin");
    if (prevPs) lg(T + 62, "#e9a1a0", "", "previous");
    keepTextFill(psSvg);
  }

  function drawReadout(): void {
    if (!pop) return;
    const nS = Math.max(1, pop.nSnap);
    const nc = pop.sumCoh / nS;
    const lines: string[] = [];
    lines.push(`cohesins per cell: ${nc.toFixed(1)}` + (nc > 0.05 ? `  (one per ${fmtKb(N / nc)})` : ""));
    if (pop.nLoop > 0) {
      lines.push(`mean loop size: ${fmtKb(pop.sumLoop / pop.nLoop)}`);
      lines.push(`legs stopped at CTCF: ${Math.round(100 * pop.sumStalledLegs / pop.sumLegs)}%`);
    } else lines.push("no loops");
    lines.push(`P(s) at 250 kb: ${(ps[50] * 50).toFixed(1)}× the no-cohesin line`);
    const pairs = sitePairs(N, sites, mean, ps, 5, maxD).slice(0, 3);
    if (pairs.length) {
      lines.push("strongest pixels between two CTCF sites (obs/exp):");
      const arrow = (s: Site) => (s.dir === 1 ? "▶" : "◀");
      for (const p of pairs)
        lines.push(`  ${fmtKb(p.a.pos)} ${arrow(p.a)} – ${fmtKb(p.b.pos)} ${arrow(p.b)}  ${p.kind.padEnd(10)} ${p.oe.toFixed(1)}`);
    }
    if (pick) {
      const [i, j] = pick, v = mean[i * N + j], e = ps[j - i];
      lines.push(`clicked pixel ${fmtKb(i)} – ${fmtKb(j)}: obs/exp ${e > 0 ? (v / e).toFixed(2) : "–"}`);
    }
    readout.textContent = lines.join("\n");
    status.textContent = rounds < SNAP_ROUNDS ? `averaging ${rounds * CELLS} of ${SNAP_ROUNDS * CELLS} snapshots…`
      : `map = average of ${SNAP_ROUNDS * CELLS} snapshots (${CELLS} cells × ${SNAP_ROUNDS} times)`;
  }

  applyPreset(start);
  requestAnimationFrame(tick);
}

autoMount(NAME, mount);
