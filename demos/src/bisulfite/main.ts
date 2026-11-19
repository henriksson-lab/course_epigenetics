// Bisulfite sequencing of a short double-stranded sequence.
// (a) the genome, both strands, with lollipops on CpGs (click to cycle unmethylated ->
//     methylated -> 50% of molecules); non-CpG Cs are unmethylated; optional C/T SNP at a CpG.
// (b) the reads each strand gives after conversion and PCR (OT, OB; + CTOT, CTOB for a
//     non-directional library), 5'->3' as sequenced; converted bases in red.
// (c) three-letter alignment: reads converted in silico, against the C->T and G->A genomes.
// (d) methylation calls: each read compared with the original genome (Bismark XM letters).
// (e) simulated reads (coverage, conversion efficiency): pileup and per-CpG beta.
// Auto-mounts on <div data-demo="bisulfite">. Optional data-*:
//   data-seq="TACG..."  top strand        data-meth="1 0 1"  level per CpG (0, 0.5, 1)
//   data-snp="3"        CpG number with the SNP option        data-panels="abcde"
//   data-coverage="20"  data-eff="99.5"   data-seed="1"
//   data-genotype="ref|het|hom"  (SNP at start)   data-directional="off" (non-directional at start)
import { svgEl, htmlEl, slider, button, presets, checkbox, autoMount, numList } from "../lib/ui";
import * as M from "./model";

const NAME = "bisulfite";
const W = 800;
const C = { fg: "#1b2a3a", muted: "#5b6b7a", light: "#eaf1f7", cell: "#f4f6f8", cellLine: "#dfe6ec",
  red: "#c0392b", blue: "#2a6f97", orange: "#e67e22", green: "#2e8b57" };
const MONO = "ui-monospace, 'Courier New', monospace";
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .bs-ctrl { display: grid; grid-template-columns: repeat(auto-fit, minmax(270px, 1fr)); gap: .35em 1.4em; width: 100%; }
[data-demo="${NAME}"] .bs-row { display: flex; flex-wrap: wrap; gap: .4em; align-items: center; }
[data-demo="${NAME}"] .bs-seq { font: 16px ${MONO}; width: 15em; max-width: 100%; padding: .15em .3em; border: 1px solid #bbb; border-radius: 4px; text-transform: uppercase; }
[data-demo="${NAME}"] .bs-scroll { width: 100%; overflow-x: auto; }
[data-demo="${NAME}"] svg { width: 100%; min-width: 760px; max-width: 1000px; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 16px; }
[data-demo="${NAME}"] svg .ttl { font-size: 18px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 15px; fill: ${C.muted}; }
[data-demo="${NAME}"] svg .base { font-family: ${MONO}; font-weight: 700; font-size: 17px; }
[data-demo="${NAME}"] svg .xm { font-family: ${MONO}; font-size: 15px; font-weight: 700; }
[data-demo="${NAME}"] svg .lolli { cursor: pointer; }
[data-demo="${NAME}"] .qb-note { margin: 0; }
.reveal [data-demo="${NAME}"] { font-size: 19px; gap: .3em; }
.reveal [data-demo="${NAME}"] .qb-note { font-size: 16px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const PRESETS: { label: string; seq: string; meth: number[] }[] = [
  { label: "default", seq: "TACGTTCGCTAAGCGTTCAG", meth: [1, 0, 1] },
  { label: "B1 (by hand)", seq: "ACGTCCGATCAG", meth: [1, 0] },
  { label: "CpG-rich", seq: "GCGCTCGACGTTCGCGAC", meth: [0, 0, 0, 0, 0, 0] },
];

interface State {
  seq: string;
  levels: Map<number, number>;
  snpCpg: number;            // 1-based CpG number for the SNP
  genotype: M.Genotype;
  directional: boolean;
  coverage: number;
  eff: number;               // fraction
  seed: number;
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const panels = (ds.panels ?? "abcde").toLowerCase();
  const seq0 = M.cleanSeq(ds.seq ?? PRESETS[0].seq) || PRESETS[0].seq;
  const st: State = {
    seq: seq0, levels: new Map(), snpCpg: Math.max(1, Number(ds.snp ?? 3) || 3),
    genotype: (["ref", "het", "hom"].includes(ds.genotype ?? "") ? ds.genotype : "ref") as M.Genotype,
    directional: ds.directional !== "off", coverage: Math.round(Number(ds.coverage ?? 20)) || 20,
    eff: (Number(ds.eff ?? 99.5) || 99.5) / 100, seed: Number(ds.seed ?? 1) || 1,
  };
  const setLevels = (vals: number[]) => {
    st.levels = new Map(M.cpgs(st.seq).map((p, k) => [p, vals[k] ?? 1]));
  };
  setLevels(numList(ds.meth, ds.seq ? [] : PRESETS[0].meth));

  // ---------------- controls ----------------
  const ctr = htmlEl("div", "bs-ctrl", "", host);
  const c1 = htmlEl("div", "", "", ctr);
  const r1 = htmlEl("div", "bs-row", "", c1);
  htmlEl("span", "", "top strand 5'→3':", r1);
  const inp = htmlEl("input", "bs-seq", "", r1);
  inp.value = st.seq; inp.maxLength = 30; inp.spellcheck = false;
  const apply = () => {
    const s = M.cleanSeq(inp.value).slice(0, 30);
    if (s.length < 6) { inp.value = st.seq; return; }
    const old = st.levels;
    st.seq = s; inp.value = s;
    st.levels = new Map(M.cpgs(s).map(p => [p, old.get(p) ?? 1]));
    redraw();
  };
  inp.addEventListener("change", () => { apply(); inp.blur(); });
  button(r1, "apply", apply);
  const r2 = htmlEl("div", "bs-row", "", c1);
  htmlEl("span", "", "examples:", r2);
  presets(r2, PRESETS.map(p => ({ label: p.label, apply: () => { st.seq = p.seq; inp.value = p.seq; setLevels(p.meth); redraw(); } })));

  const c2 = htmlEl("div", "", "", ctr);
  checkbox(c2, "non-directional library (all four strands)", !st.directional, v => { st.directional = !v; redraw(); });
  const r3 = htmlEl("div", "bs-row", "", c2);
  const snpLab = htmlEl("span", "", "", r3);
  const snpBtns = presets(r3, [
    { label: "no SNP", apply: () => { st.genotype = "ref"; redraw(); } },
    { label: "C/T", apply: () => { st.genotype = "het"; redraw(); } },
    { label: "T/T", apply: () => { st.genotype = "hom"; redraw(); } },
  ]);
  snpBtns[["ref", "het", "hom"].indexOf(st.genotype)].classList.add("active");

  const c3 = htmlEl("div", "", "", ctr);
  slider(c3, "reads", 1, 100, 1, st.coverage, v => { st.coverage = v; drawSim(); });
  slider(c3, "conversion", 90, 100, 0.1, st.eff * 100, v => { st.eff = v / 100; drawSim(); }, v => v.toFixed(1) + "%");
  button(c3, "draw reads again", () => { st.seed++; drawSim(); });

  htmlEl("p", "qb-note", "Simulated molecules and reads. Every read covers the whole sequence (real reads are longer than this toy genome and start anywhere).", host);

  const wrapMain = htmlEl("div", "bs-scroll", "", host), wrapSim = htmlEl("div", "bs-scroll", "", host);
  const svgMain = svgEl("svg", { viewBox: `0 0 ${W} 600`, role: "img", "aria-label": "Bisulfite conversion, read types and alignment" }, wrapMain);
  const svgSim = svgEl("svg", { viewBox: `0 0 ${W} 420`, role: "img", "aria-label": "Simulated reads and per-CpG methylation" }, wrapSim);
  if (!/[abcd]/.test(panels)) wrapMain.style.display = "none";
  if (!panels.includes("e")) wrapSim.style.display = "none";

  // ---------------- geometry ----------------
  const LBL = 128, X0 = 136;
  let cw = 30;
  const cx = (i: number) => X0 + i * cw + cw / 2;
  const snpPos = (): number => {
    const cg = M.cpgs(st.seq);
    if (st.genotype === "ref" || !cg.length) return -1;
    return cg[Math.min(st.snpCpg, cg.length) - 1];
  };

  function baseRow(svg: SVGSVGElement, y: number, s: string, colour: (i: number) => string, left: string, right: string,
    fill = C.cell): void {
    if (left) svgEl("text", { x: LBL, y: y + 20, "text-anchor": "end" }, svg).textContent = left;
    for (let i = 0; i < s.length; i++) {
      svgEl("rect", { x: X0 + i * cw + 1, y, width: cw - 2, height: 28, fill, stroke: C.cellLine }, svg);
      const t = svgEl("text", { x: cx(i), y: y + 20, "text-anchor": "middle", class: "base" }, svg);
      t.textContent = s[i]; t.style.fill = colour(i);
    }
    if (right) svgEl("text", { x: X0 + s.length * cw + 14, y: y + 20, class: "sm" }, svg).textContent = right;
  }

  function lollipop(svg: SVGSVGElement, x: number, yBase: number, dir: -1 | 1, level: number, onClick?: () => void): void {
    const g = svgEl("g", onClick ? { class: "lolli" } : {}, svg);
    const yc = yBase + dir * 18;
    svgEl("line", { x1: x, x2: x, y1: yBase, y2: yc, stroke: C.fg, "stroke-width": 2 }, g);
    svgEl("circle", { cx: x, cy: yc, r: 7.5, fill: level >= 1 ? C.fg : "#fff", stroke: C.fg, "stroke-width": 2 }, g);
    if (level > 0 && level < 1) svgEl("path", { d: `M${x},${yc - 7.5} A7.5,7.5 0 0 0 ${x},${yc + 7.5} Z`, fill: C.fg }, g);
    if (onClick) {
      svgEl("rect", { x: x - Math.max(10, cw / 2), y: Math.min(yBase, yc - 12), width: Math.max(20, cw), height: 32, fill: "transparent" }, g);
      g.addEventListener("click", onClick);
    }
  }

  // ---------------- (a)–(d) ----------------
  function drawMain(): void {
    const seq = st.seq, L = seq.length;
    cw = Math.min(28, Math.floor(390 / L));
    const svg = svgMain;
    svg.replaceChildren();
    const sp = snpPos();
    const cg = M.cpgs(seq);
    const types: M.ReadType[] = st.directional ? ["OT", "OB"] : M.READ_TYPES;
    // one molecule: a 50% CpG is methylated on this molecule
    const meth = new Set(cg.filter(p => (st.levels.get(p) ?? 0) > 0));
    const reads = M.readsFromMolecule({ top: seq, meth });
    let y = 0;

    if (panels.includes("a")) {
      svgEl("text", { x: 8, y: y + 22, class: "ttl" }, svg).textContent =
        "(a) Genome: click a lollipop to change its methylation";
      y += 64;
      const yTop = y, yBot = y + 30;
      // context labels for non-CpG Cs
      for (let i = 0; i < L; i++) {
        const tc = M.topContext(seq, i), bc = M.bottomContext(seq, i);
        if (tc && tc !== "CG") svgEl("text", { x: cx(i), y: yTop - 8, "text-anchor": "middle", class: "sm" }, svg).textContent = tc;
        if (bc && bc !== "CG") svgEl("text", { x: cx(i), y: yBot + 28 + 20, "text-anchor": "middle", class: "sm" }, svg).textContent = bc;
      }
      baseRow(svg, yTop, seq, () => C.fg, "top  5'", "3'");
      baseRow(svg, yBot, M.comp(seq), () => C.fg, "bottom  3'", "5'");
      for (const p of cg) {
        const lev = st.levels.get(p) ?? 0;
        const cyc = () => { st.levels.set(p, lev === 0 ? 1 : lev === 1 ? 0.5 : 0); redraw(); };
        lollipop(svg, cx(p), yTop, -1, lev, cyc);
        lollipop(svg, cx(p + 1), yBot + 28, 1, lev, cyc);
      }
      if (sp >= 0) {
        for (const yy of [yTop, yBot])
          svgEl("rect", { x: X0 + sp * cw, y: yy - 1, width: cw, height: 30, fill: "none", stroke: C.red, "stroke-width": 2.5 }, svg);
      }
      const tx = X0 + L * cw + 40;
      const nTopC = [...seq].filter((_, i) => M.topContext(seq, i) && M.topContext(seq, i) !== "CG").length;
      const nBotC = [...seq].filter((_, i) => M.bottomContext(seq, i) && M.bottomContext(seq, i) !== "CG").length;
      svgEl("text", { x: tx, y: yTop + 6, class: "sm" }, svg).textContent = `${cg.length} CpGs (symmetric)`;
      svgEl("text", { x: tx, y: yTop + 28, class: "sm" }, svg).textContent = `non-CpG Cs: ${nTopC} top, ${nBotC} bottom,`;
      svgEl("text", { x: tx, y: yTop + 48, class: "sm" }, svg).textContent = "all unmethylated";
      if (sp >= 0) {
        const t = svgEl("text", { x: tx, y: yTop + 70, class: "sm" }, svg);
        t.textContent = `SNP at CpG ${cg.indexOf(sp) + 1}: sample is ${st.genotype === "het" ? "C/T" : "T/T"}`;
        t.style.fill = C.red;
      }
      y = yBot + 28 + 34;
    }

    if (panels.includes("b")) {
      y += 10;
      svgEl("text", { x: 8, y: y + 20, class: "ttl" }, svg).textContent =
        "(b) Reads after bisulfite and PCR (5'→3', as sequenced)";
      if (sp >= 0) svgEl("text", { x: W - 8, y: y + 20, class: "sm", "text-anchor": "end" }, svg).textContent = "(reference allele)";
      y += 34;
      const orig: Record<M.ReadType, string> = { OT: seq, OB: M.revcomp(seq), CTOT: M.revcomp(seq), CTOB: seq };
      for (const t of types) {
        const r = reads[t];
        const nC = [...r].filter(b => b === "C").length, nG = [...r].filter(b => b === "G").length;
        baseRow(svg, y, r, i => (r[i] !== orig[t][i] ? C.red : C.fg), `${t}  5'`,
          t === "OT" || t === "OB" ? `3'   C: ${nC} of ${L}` : `3'   G: ${nG} of ${L}`);
        y += 34;
      }
      svgEl("text", { x: 8, y: y + 12, class: "sm" }, svg).textContent = types.length === 2
        ? "Red: converted bases. Directional library: read 1 is OT or OB. OB runs 5'→3' on the bottom strand."
        : "Red: converted bases. CTOT and CTOB are PCR copies (reverse complements) of OT and OB.";
      y += 26;
    }

    if (panels.includes("c")) {
      y += 12;
      svgEl("text", { x: 8, y: y + 20, class: "ttl" }, svg).textContent =
        "(c) Three-letter alignment of reads converted in silico";
      y += 30;
      svgEl("text", { x: 8, y: y + 16, class: "sm" }, svg).textContent =
        "Reads shown on the + strand, as a genome browser does; red = mismatch with the genome above.";
      y += 30;
      const gCT = M.ct(seq), gGA = M.ga(seq);
      const groups: [string, string, M.ReadType[]][] = [["genome C→T", gCT, ["OT", "CTOT"]], ["genome G→A", gGA, ["OB", "CTOB"]]];
      for (const [lab, g, ts] of groups) {
        baseRow(svg, y, g, () => C.fg, lab, "", C.light);
        y += 34;
        for (const t of ts) {
          if (!types.includes(t)) continue;
          const al = M.align(reads[t], seq, st.directional);
          const own = al.find(a => a.type === t)!;
          const others = al.filter(a => a !== own).map(a => a.mismatches);
          const best = Math.min(...others);
          const conv = own.readConv;
          baseRow(svg, y, own.shown, i => (own.shown[i] !== g[i] ? C.red : C.fg), `${t} read, ${conv}`,
            "", "#fff");
          const tx = X0 + L * cw + 14;
          const ok = own.mismatches < best;
          const t1 = svgEl("text", { x: tx, y: y + 12, class: "sm" }, svg);
          t1.textContent = `${own.mismatches} mismatches here (${own.strand} strand)`;
          t1.style.fill = ok ? C.green : C.red;
          svgEl("text", { x: tx, y: y + 30, class: "sm" }, svg).textContent =
            `${ok ? "unique;" : "tie: discarded;"} next best ${best}, plain ${M.plainMismatches(reads[t], seq)}`;
          y += 38;
        }
      }
      svgEl("text", { x: 8, y: y + 10, class: "sm" }, svg).textContent =
        "next best = other converted genome or strand; plain = ordinary alignment, no conversion.";
      y += 24;
    }

    if (panels.includes("d")) {
      y += 14;
      svgEl("text", { x: 8, y: y + 20, class: "ttl" }, svg).textContent =
        "(d) Calls: original read against original genome (XM letters)";
      y += 34;
      baseRow(svg, y, seq, () => C.fg, "genome (+)", "", C.light);
      y += 34;
      for (const t of types) {
        const r = M.refOrientation(t, reads[t]);
        const calls = M.callRead(t, reads[t], seq);
        baseRow(svg, y, r, i => (r[i] !== seq[i] ? C.red : C.fg), `${t}`, "", "#fff");
        for (const c of calls) {
          const tt = svgEl("text", { x: cx(c.pos), y: y + 46, "text-anchor": "middle", class: "xm" }, svg);
          tt.textContent = c.letter;
          tt.style.fill = c.meth ? C.fg : C.muted;
        }
        const xm = [...seq].map((_, i) => calls.find(c => c.pos === i)?.letter ?? ".").join("");
        svgEl("text", { x: X0 + L * cw + 14, y: y + 20, class: "sm" }, svg).textContent = `XM:Z:${xm}`;
        y += 56;
      }
      svgEl("text", { x: 8, y: y + 6, class: "sm" }, svg).textContent =
        "Z/z = methylated/unmethylated CpG, X/x CHG, H/h CHH; XM runs along the + strand, as in the BAM.";
      y += 22;
      svgEl("text", { x: 8, y: y + 6, class: "sm" }, svg).textContent =
        "Top-strand reads (OT, CTOT) are called at the genome's Cs, bottom-strand reads (OB, CTOB) at its Gs.";
      y += 20;
    }
    svg.setAttribute("viewBox", `0 0 ${W} ${y + 10}`);
  }

  // ---------------- (e) simulated reads ----------------
  function drawSim(): void {
    const seq = st.seq, L = seq.length;
    cw = Math.min(28, Math.floor(390 / L));
    const svg = svgSim;
    svg.replaceChildren();
    const sp = snpPos();
    const reads = M.simulateReads({ seq, levels: st.levels, snpPos: sp, genotype: st.genotype, directional: st.directional,
      coverage: st.coverage, eff: st.eff, seed: st.seed });
    const topR = reads.filter(r => r.type === "OT" || r.type === "CTOT");
    const botR = reads.filter(r => r.type === "OB" || r.type === "CTOB");
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent =
      `(e) ${reads.length} simulated reads (● C kept = methylated, ○ T = unmethylated)`;
    const HALF = 140, yMid = 40 + HALF + 22;
    const maxSide = Math.max(topR.length, botR.length, 1);
    const rh = Math.min(14, HALF / maxSide);
    const sz = Math.max(2, Math.min(rh - 1, cw * 0.55));
    const small = rh < 6;
    const drawRead = (rd: M.SimRead, yy: number) => {
      if (!small) svgEl("line", { x1: X0 + 2, x2: X0 + L * cw - 2, y1: yy, y2: yy, stroke: "#c3ccd5", "stroke-width": Math.min(2, rh / 4) }, svg);
      for (const c of rd.calls) {
        if (c.meth === null) continue;
        const nonCpG = c.context !== "CG";
        const col = c.meth ? (nonCpG ? C.orange : C.fg) : (small ? "#b4c1ce" : "#fff");
        const w = small ? cw * 0.62 : sz, h = small ? Math.max(1.5, rh - 0.8) : sz;
        svgEl("rect", { x: cx(c.pos) - w / 2, y: yy - h / 2, width: w, height: h, rx: small ? 0 : sz / 2,
          fill: col, stroke: c.meth || small ? "none" : C.muted, "stroke-width": 1.3 }, svg);
      }
      if (sp >= 0 && (rd.type === "OB" || rd.type === "CTOB") && rd.ref[sp] === "T") {
        if (rh >= 9) {
          const t = svgEl("text", { x: cx(sp), y: yy + 5, "text-anchor": "middle", class: "xm" }, svg);
          t.textContent = "T"; t.style.fill = C.red;
        } else svgEl("rect", { x: cx(sp) - cw * 0.31, y: yy - Math.max(1.5, rh - 0.8) / 2, width: cw * 0.62, height: Math.max(1.5, rh - 0.8), fill: C.red }, svg);
      }
    };
    topR.forEach((rd, k) => drawRead(rd, yMid - 22 - rh / 2 - k * rh));
    botR.forEach((rd, k) => drawRead(rd, yMid + 22 + rh / 2 + k * rh));
    // genome in the middle
    for (let i = 0; i < L; i++) {
      const t = svgEl("text", { x: cx(i), y: yMid + 6, "text-anchor": "middle", class: "base" }, svg);
      t.textContent = seq[i];
      t.style.fill = i === sp ? C.red : C.muted;
    }
    svgEl("text", { x: LBL, y: yMid - 22 - Math.min(HALF, topR.length * rh) / 2 + 5, "text-anchor": "end" }, svg).textContent = `top: ${topR.length} reads`;
    svgEl("text", { x: LBL, y: yMid + 22 + Math.min(HALF, botR.length * rh) / 2 + 5, "text-anchor": "end" }, svg).textContent = `bottom: ${botR.length}`;
    svgEl("text", { x: LBL, y: yMid + 6, "text-anchor": "end", class: "sm" }, svg).textContent = "genome";

    // per-CpG summary below
    const cgs = M.summarizeCpgs(seq, st.levels, reads);
    const yb0 = yMid + 22 + HALF + 34, bh = 70;
    svgEl("text", { x: LBL, y: yb0 + 4, "text-anchor": "end" }, svg).textContent = "β per CpG";
    svgEl("text", { x: LBL, y: yb0 + 24, "text-anchor": "end", class: "sm" }, svg).textContent = "both strands";
    svgEl("line", { x1: X0, x2: X0 + L * cw, y1: yb0 + bh, y2: yb0 + bh, stroke: C.fg }, svg);
    svgEl("line", { x1: X0, x2: X0 + L * cw, y1: yb0, y2: yb0, stroke: "#e3e3e3", "stroke-dasharray": "3 3" }, svg);
    svgEl("text", { x: X0 + L * cw + 6, y: yb0 + 5, class: "sm" }, svg).textContent = "1";
    svgEl("text", { x: X0 + L * cw + 6, y: yb0 + bh + 5, class: "sm" }, svg).textContent = "0";
    for (const c of cgs) {
      const x = (cx(c.pos) + cx(c.pos + 1)) / 2, w = Math.min(2 * cw - 8, 34);
      const n = c.nTop + c.nBot, m = c.mTop + c.mBot;
      if (n > 0) {
        svgEl("rect", { x: x - w / 2, y: yb0 + bh * (1 - c.beta), width: w, height: bh * c.beta, fill: C.blue, "fill-opacity": 0.8 }, svg);
        const [lo, hi] = wilson(m, n);
        svgEl("line", { x1: x, x2: x, y1: yb0 + bh * (1 - hi), y2: yb0 + bh * (1 - lo), stroke: C.fg, "stroke-width": 2 }, svg);
        for (const v of [lo, hi]) svgEl("line", { x1: x - 5, x2: x + 5, y1: yb0 + bh * (1 - v), y2: yb0 + bh * (1 - v), stroke: C.fg, "stroke-width": 2 }, svg);
      }
      svgEl("line", { x1: x - w / 2 - 4, x2: x + w / 2 + 4, y1: yb0 + bh * (1 - c.level), y2: yb0 + bh * (1 - c.level), stroke: C.green, "stroke-width": 3, "stroke-dasharray": "6 3" }, svg);
      svgEl("text", { x, y: yb0 + bh + 20, "text-anchor": "middle", class: "sm" }, svg).textContent = n ? c.beta.toFixed(2) : "–";
      svgEl("text", { x, y: yb0 + bh + 38, "text-anchor": "middle", class: "sm" }, svg).textContent = `${m}/${n}`;
    }
    const nc = M.nonCpg(reads);
    const tx = X0 + L * cw + 30;
    let ty = yb0 + 4;
    const line = (s: string, col = C.muted) => { const t = svgEl("text", { x: tx, y: ty, class: "sm" }, svg); t.textContent = s; t.style.fill = col; ty += 21; };
    line("bar: measured β, 95% CI");
    line("green dashes: true level", C.green);
    line(`orange: non-CpG C kept, ${nc.m} of ${nc.n}`, nc.m ? C.orange : C.muted);
    line(`non-conversion ${nc.n ? (100 * nc.m / nc.n).toFixed(1) : "–"}% (true ${(100 * (1 - st.eff)).toFixed(1)}%)`);
    if (sp >= 0) line("red: T allele in bottom reads", C.red);
    svg.setAttribute("viewBox", `0 0 ${W} ${yb0 + bh + 50}`);
  }

  function redraw(): void {
    const cg = M.cpgs(st.seq);
    snpLab.textContent = cg.length ? `SNP at CpG ${Math.min(st.snpCpg, cg.length)}:` : "SNP (no CpG):";
    drawMain();
    drawSim();
  }
  redraw();
}

/** Wilson 95% interval for m of n. */
function wilson(m: number, n: number): [number, number] {
  const z = 1.96, p = m / n, d = 1 + z * z / n;
  const c = (p + z * z / (2 * n)) / d, h = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d;
  return [Math.max(0, c - h), Math.min(1, c + h)];
}

autoMount(NAME, mount);
