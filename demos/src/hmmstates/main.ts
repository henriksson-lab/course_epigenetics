// Chromatin states with a hidden Markov model, ChromHMM style: three hidden states along
// 200-bp bins emit four binary marks. Generate a sequence, then decode it with Viterbi, with
// the forward–backward posterior, and with the emissions alone (each bin on its own).
//
// <div data-demo="hmmstates"
//      data-n="150"     number of bins (40..300)
//      data-seed="2"></div>
//
// The data are generated with the "true" self-transition; the decoder has its own
// self-transition slider (1/3 = no preference for staying = emissions only). Both use the
// emission table. Any change re-simulates with the same random numbers. Model code: hmm.ts.
import { svgEl, htmlEl, button, presets, slider, rng, autoMount, Slider } from "../lib/ui";
import { Model, stayMatrix, simulate, forwardBackward, viterbi, emissionOnly, posteriorArgmax, segments, agreement } from "./hmm";

const NAME = "hmmstates";
const MARKS = ["H3K4me3", "H3K27ac", "H3K36me3", "H3K27me3"];
const MARK_COL = ["#c0392b", "#e67e22", "#2e8b57", "#2a6f97"];
const STATES = ["TSS", "Tx", "Pc"];
const STATE_COL = ["#c0392b", "#2e8b57", "#8a96a3"];
const C_TEXT = "#1b2a3a", C_MUTED = "#5b6b7a";

const EMIT_PRESETS: { label: string; E: number[][] }[] = [
  { label: "clear marks", E: [[0.85, 0.75, 0.05, 0.05], [0.05, 0.10, 0.70, 0.03], [0.10, 0.02, 0.03, 0.65]] },
  { label: "noisy marks", E: [[0.55, 0.45, 0.15, 0.15], [0.15, 0.20, 0.45, 0.10], [0.20, 0.10, 0.10, 0.45]] },
  { label: "sparse marks", E: [[0.40, 0.30, 0.02, 0.02], [0.02, 0.05, 0.25, 0.02], [0.05, 0.01, 0.02, 0.25]] },
];

const STYLE = `
[data-demo="${NAME}"].qb-demo { display: block; }
[data-demo="${NAME}"] .hs-svg { width: 100%; height: auto; display: block; }
[data-demo="${NAME}"] .hs-bottom { display: flex; flex-wrap: wrap; gap: .5em 1.4em; align-items: flex-start; margin-top: .3em; }
[data-demo="${NAME}"] .hs-col { display: flex; flex-direction: column; gap: .4em; flex: 1 1 340px; min-width: 0; }
[data-demo="${NAME}"] .hs-row { display: flex; flex-wrap: wrap; gap: .35em .6em; align-items: center; }
[data-demo="${NAME}"] .hs-head { font-weight: 600; margin-top: .2em; }
[data-demo="${NAME}"] .qb-slider { grid-template-columns: 11.2em 1fr 3em; }
[data-demo="${NAME}"] .qb-readout { font-size: 15px; min-height: 8.2em; }
[data-demo="${NAME}"] table.hs-em { border-collapse: collapse; font-size: 15px; width: 100%; }
[data-demo="${NAME}"] table.hs-em th { font-weight: 600; padding: 0 .2em .1em; text-align: center; white-space: nowrap; }
[data-demo="${NAME}"] table.hs-em td { padding: .05em .2em; text-align: center; }
[data-demo="${NAME}"] table.hs-em input[type=range] { width: 100%; min-width: 3.2em; margin: 0; }
[data-demo="${NAME}"] table.hs-em output { display: block; font-family: ui-monospace, monospace; font-size: 15px; }
[data-demo="${NAME}"] .hs-sw { display: inline-block; width: 15px; height: 13px; border-radius: 2px; margin-right: .25em; vertical-align: -1px; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; margin: 0; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

/** The shared CSS sets every svg text to the foreground colour; keep the colours given as fill. */
function keepTextFill(root: Element): void {
  root.querySelectorAll<SVGTextElement>("text[fill]").forEach(t => { t.style.fill = t.getAttribute("fill") ?? ""; });
}

function mount(host: HTMLElement): void {
  injectStyle();
  const T = Math.max(40, Math.min(300, Number(host.dataset.n ?? 150) || 150));
  const seed0 = Number(host.dataset.seed ?? 2) || 2;
  const K = 3, M = 4;

  // ---- state ----
  let E = EMIT_PRESETS[0].E.map(r => r.slice());
  let stayTrue = 0.9, stayDec = 0.9;
  let seed = seed0;
  let us: number[] = [];
  const drawUs = () => { const r = rng(seed); us = Array.from({ length: T * (1 + M) }, () => r.unif()); };

  // ---- SVG tracks ----
  // a slide is ~1150 px wide, a page ~800 px: a narrower viewBox on a page keeps the text large
  const W = host.closest(".reveal") ? 1100 : 820, LW = 170, TW = W - LW - 4, bw = TW / T;
  const hMark = 19, gMark = 3;
  const yMarks = 30;
  const yTrue = yMarks + M * (hMark + gMark) + 14, hS = 24;
  const yVit = yTrue + hS + 18;
  const yPost = yVit + hS + 18, hPost = 62;
  const yEmi = yPost + hPost + 14;
  const H = yEmi + hS + 14;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, class: "hs-svg", role: "img",
    "aria-label": "Binary marks per bin, the true states, the Viterbi path, posterior probabilities and emissions-only decoding" }, host);

  // ---- controls ----
  const bottom = htmlEl("div", "hs-bottom", "", host);
  const colL = htmlEl("div", "hs-col", "", bottom);
  const colR = htmlEl("div", "hs-col", "", bottom);

  htmlEl("div", "hs-head", "Generate (simulated data)", colL);
  const sTrue: Slider = slider(colL, "true self-transition", 0.5, 0.99, 0.01, stayTrue, v => { stayTrue = v; update(); }, v => v.toFixed(2));
  const gRow = htmlEl("div", "hs-row", "", colL);
  button(gRow, "New sequence", () => { seed++; drawUs(); update(); });
  button(gRow, "Reset", () => reset());

  htmlEl("div", "hs-head", "Decode", colL);
  const sDec: Slider = slider(colL, "decoder self-transition", 0.34, 0.995, 0.005, stayDec, v => { stayDec = v; decPresets.forEach(b => b.classList.remove("active")); update(); }, v => v.toFixed(3));
  const dRow = htmlEl("div", "hs-row", "", colL);
  const decPresets = presets(dRow, [
    { label: "1/3 (emissions only)", apply: () => setDec(1 / 3) },
    { label: "= true", apply: () => setDec(stayTrue) },
    { label: "0.99", apply: () => setDec(0.99) },
  ]);
  function setDec(v: number): void { stayDec = Math.max(1 / 3, v); sDec.set(Math.round(stayDec * 1000) / 1000); update(); }
  const readout = htmlEl("div", "qb-readout", "", colL);

  htmlEl("div", "hs-head", "Emissions: P(mark present | state)", colR);
  const eRow = htmlEl("div", "hs-row", "", colR);
  const emBtns = presets(eRow, EMIT_PRESETS.map(p => ({ label: p.label, apply: () => setE(p.E) })));
  const table = htmlEl("table", "hs-em", "", colR);
  const hr = htmlEl("tr", "", "", table);
  htmlEl("th", "", "", hr);
  MARKS.forEach((mk, k) => { const th = htmlEl("th", "", mk.replace("H3", ""), hr); th.style.color = MARK_COL[k]; th.title = mk; });
  const cells: { input: HTMLInputElement; out: HTMLOutputElement }[][] = [];
  for (let s = 0; s < K; s++) {
    const tr = htmlEl("tr", "", "", table);
    const th = htmlEl("th", "", "", tr);
    const sw = htmlEl("span", "hs-sw", "", th); sw.style.background = STATE_COL[s];
    th.append(STATES[s]);
    th.style.textAlign = "left";
    cells.push([]);
    for (let m = 0; m < M; m++) {
      const td = htmlEl("td", "", "", tr);
      const input = htmlEl("input", "", "", td);
      Object.assign(input, { type: "range", min: "0.01", max: "0.99", step: "0.01", value: String(E[s][m]) });
      input.style.accentColor = MARK_COL[m];
      input.setAttribute("aria-label", `P(${MARKS[m]} | ${STATES[s]})`);
      const out = htmlEl("output", "", E[s][m].toFixed(2), td);
      input.addEventListener("input", () => { E[s][m] = Number(input.value); out.value = E[s][m].toFixed(2); emBtns.forEach(b => b.classList.remove("active")); update(); });
      input.addEventListener("change", () => input.blur());
      input.addEventListener("pointerup", () => input.blur());
      cells[s].push({ input, out });
    }
  }
  htmlEl("div", "qb-note", "States: TSS = active promoter, Tx = transcribed gene body, Pc = Polycomb-repressed. "
    + "Bins of 200 bp; a black tick marks a wrongly decoded bin. Leaving a state costs a factor (1 − stay)/2 per switch.", colR);

  function setE(Ev: number[][]): void {
    E = Ev.map(r => r.slice());
    for (let s = 0; s < K; s++) for (let m = 0; m < M; m++) { cells[s][m].input.value = String(E[s][m]); cells[s][m].out.value = E[s][m].toFixed(2); }
    update();
  }
  function reset(): void {
    seed = seed0; drawUs();
    stayTrue = 0.9; sTrue.set(0.9);
    emBtns.forEach((b, i) => b.classList.toggle("active", i === 0));
    decPresets.forEach((b, i) => b.classList.toggle("active", i === 1));
    stayDec = 0.9; sDec.set(0.9);
    setE(EMIT_PRESETS[0].E);
  }

  // ---- compute + draw ----
  function update(): void {
    const pi = new Array(K).fill(1 / K);
    const truth: Model = { pi, A: stayMatrix(K, stayTrue), E };
    const dec: Model = { pi, A: stayMatrix(K, stayDec), E };
    const { s, o } = simulate(truth, us, T);
    const fb = forwardBackward(dec, o);
    const vit = viterbi(dec, o);
    const pmax = posteriorArgmax(fb.post);
    const emi = emissionOnly(dec, o);
    draw(s, o, vit, fb.post, emi);
    const pc = (a: number[]) => `${Math.round(100 * agreement(a, s))}%`;
    const row = (name: string, a: number[] | null) =>
      `${name.padEnd(11)}${(a ? pc(a) : "–").padStart(6)}${String(segments(a ?? s)).padStart(10)}`;
    readout.textContent = [
      `${"".padEnd(11)}${"correct".padStart(8)}${"segments".padStart(10)}`,
      row("true", null), row("Viterbi", vit), row("posterior", pmax), row("emissions", emi),
      `mean run: true ${(1 / (1 - stayTrue)).toFixed(1)}, decoder ${(1 / (1 - stayDec)).toFixed(1)} bins`,
      `log P(marks | decoder) = ${fb.logLik.toFixed(1).replace("-", "−")}`,
    ].join("\n");
  }

  function draw(s: number[], o: number[][], vit: number[], post: number[][], emi: number[]): void {
    svg.replaceChildren();
    const x = (t: number) => LW + t * bw;
    const lab = (y: number, txt: string, col = C_TEXT, size = 18) =>
      svgEl("text", { x: LW - 10, y, "text-anchor": "end", "font-size": size, fill: col }, svg).textContent = txt;
    // ruler
    for (let t = 0; t <= T; t += 25) {
      svgEl("line", { x1: x(t), x2: x(t), y1: 20, y2: 26, stroke: C_MUTED }, svg);
      svgEl("text", { x: Math.min(W - 26, x(t)), y: 17, "text-anchor": "middle", "font-size": 16, fill: C_MUTED }, svg).textContent = `${(t * 0.2).toFixed(0)} kb`;
    }
    // marks
    for (let m = 0; m < M; m++) {
      const y = yMarks + m * (hMark + gMark);
      lab(y + hMark - 3, MARKS[m], MARK_COL[m], 17);
      svgEl("rect", { x: LW, y, width: TW, height: hMark, fill: "#f4f6f8" }, svg);
      for (let t = 0; t < T; t++) if (o[t][m])
        svgEl("rect", { x: x(t), y, width: Math.max(1, bw - 0.6), height: hMark, fill: MARK_COL[m] }, svg);
    }
    // state strips
    const strip = (path: number[], y: number, wrong: boolean) => {
      let a = 0;
      for (let t = 1; t <= T; t++) {
        if (t === T || path[t] !== path[a]) {
          svgEl("rect", { x: x(a), y, width: x(t) - x(a), height: hS, fill: STATE_COL[path[a]] }, svg);
          a = t;
        }
      }
      if (wrong) for (let t = 0; t < T; t++) if (path[t] !== s[t])
        svgEl("rect", { x: x(t), y: y + hS + 2, width: Math.max(1.5, bw - 0.6), height: 6, fill: "#111" }, svg);
    };
    lab(yTrue + 18, "true state");
    strip(s, yTrue, false);
    lab(yVit + 18, "Viterbi");
    strip(vit, yVit, true);
    lab(yPost + 26, "posterior");
    lab(yPost + 46, "P(state | marks)", C_MUTED, 16);
    // stacked posterior
    for (let k = 0; k < K; k++) {
      let d = "";
      for (let t = 0; t < T; t++) {
        let below = 0;
        for (let q = 0; q < k; q++) below += post[t][q];
        const y0 = yPost + hPost * (1 - below), y1 = yPost + hPost * (1 - below - post[t][k]);
        d += `M${x(t).toFixed(1)},${y0.toFixed(1)}H${(x(t) + bw).toFixed(1)}V${y1.toFixed(1)}H${x(t).toFixed(1)}Z`;
      }
      svgEl("path", { d, fill: STATE_COL[k] }, svg);
    }
    svgEl("rect", { x: LW, y: yPost, width: TW, height: hPost, fill: "none", stroke: "#bbb" }, svg);
    lab(yEmi + 18, "emissions only");
    strip(emi, yEmi, true);
    keepTextFill(svg);
  }

  reset();
}

autoMount(NAME, mount);
