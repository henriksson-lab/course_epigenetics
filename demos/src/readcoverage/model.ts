// Model for the readcoverage demo (no DOM; testable with node).
//
// Simulated single-end ChIP-seq around one binding site in a WIN bp window.
// - ChIP fragments that carry the bound protein contain the site: the protein protects
//   FOOT bp on each side, so a fragment of length L starts uniformly in
//   [SITE + FOOT - L, SITE - FOOT]. Background ChIP fragments and all input fragments
//   are placed uniformly (also a margin outside the window, so the edges are not empty).
// - Fragment lengths vary around the "fragment length" setting (SD 15%).
// - Each fragment is sequenced from one end, chosen at random: a + strand read covers the
//   first readLen bases from the left end, a − strand read the last readLen bases.
// - Extension (as deepTools --extendReads <len>): every read is extended from its 5' end, in
//   its own direction, to the *assumed* fragment length (the setting), not to the fragment's
//   true length, which single-end data cannot know.
// - PCR duplicates: a few fragments get extra identical copies (same 5' end, same strand).
//   Removing duplicates keeps one read per (strand, 5' position), as samtools markdup does
//   for single-end reads; this also removes independent fragments that happen to share a start.
// - Library sizes: CHIP_LIB and INPUT_LIB reads genome-wide. With duplicates added, the
//   genome-wide library grows by the same fraction as the reads in the window; removing
//   duplicates brings it back (genome-wide, chance collisions of independent starts are rare).
// - Coverage per bin = number of reads (or extended reads) overlapping the bin, as
//   bamCoverage counts. CPM = count × 1e6 / library. Fold change (bamCompare --operation log2
//   --scaleFactorsMethod readCount): input scaled to the ChIP's depth, then
//   log2((ChIP + pc) / (input_scaled + pc)).
import { rng } from "../lib/ui";

export const WIN = 1200, SITE = 600, FOOT = 18, MARGIN = 500;
export const CHIP_LIB = 2_000_000, INPUT_LIB = 1_000_000;
export const N_SITE = 150, N_BG = 40, N_INPUT = 20; // N_BG, N_INPUT over WIN + 2*MARGIN
export const LEN_CV = 0.15;

export type Norm = "raw" | "cpm" | "fc";

export interface Params {
  fragLen: number; readLen: number; extend: boolean;
  dups: boolean; dedup: boolean; norm: Norm; pseudo: number; bin: number;
}

export const DEFAULTS: Params = { fragLen: 200, readLen: 36, extend: false, dups: false, dedup: false, norm: "raw", pseudo: 1, bin: 10 };

/** The random draws of one simulated data set, independent of the settings
 *  (common random numbers: moving a slider changes the reads smoothly). */
export interface Draws {
  chip: { u: number; z: number; plus: boolean; site: boolean; copies: number }[];
  input: { u: number; z: number; plus: boolean }[];
}

export function makeDraws(seed: number): Draws {
  const r = rng(seed);
  const chip: Draws["chip"] = [];
  for (let i = 0; i < N_SITE + N_BG; i++) {
    const u = r.unif(), z = Math.max(-2.5, Math.min(2.5, r.norm())), plus = r.unif() < 0.5;
    // about 5% of fragments are over-amplified: 2-20 extra identical copies
    const a = r.unif(), b = r.unif();
    const copies = a < 0.05 ? Math.min(20, 2 + Math.floor(-Math.log(1 - b) * 5)) : 0;
    chip.push({ u, z, plus, site: i < N_SITE, copies });
  }
  const input: Draws["input"] = [];
  for (let i = 0; i < N_INPUT; i++) {
    const u = r.unif(), z = Math.max(-2.5, Math.min(2.5, r.norm())), plus = r.unif() < 0.5;
    input.push({ u, z, plus });
  }
  return { chip, input };
}

export interface Read {
  start: number; end: number;        // shown/counted interval, [start, end)
  plus: boolean;
  five: number;                      // 5' end position (+: start, −: end − 1)
  fragStart: number; fragEnd: number;// the true fragment (known only in a simulation)
  frag: number;                      // index of the fragment
  pcr: boolean;                      // an extra PCR copy
  removed: boolean;                  // removed by duplicate removal
}

function fragLength(fragLen: number, z: number): number {
  return Math.max(2 * FOOT + 10, Math.round(fragLen * (1 + LEN_CV * z)));
}

function makeRead(fs: number, L: number, plus: boolean, p: Params, frag: number, pcr: boolean): Read {
  const fe = fs + L;
  const r = Math.min(p.readLen, L);           // a read cannot be longer than its fragment
  const len = p.extend ? Math.max(r, p.fragLen) : r;
  return plus
    ? { start: fs, end: fs + len, plus, five: fs, fragStart: fs, fragEnd: fe, frag, pcr, removed: false }
    : { start: fe - len, end: fe, plus, five: fe - 1, fragStart: fs, fragEnd: fe, frag, pcr, removed: false };
}

export interface Sim {
  chip: Read[]; input: Read[];
  pcrCopies: number;           // extra PCR copies in the window data
  removedPcr: number;          // removed copies that were PCR duplicates
  removedIndep: number;        // removed reads that were independent fragments with the same start
  chipLib: number; inputLib: number;
}

/** Shift of read 5' ends outside the window is fine; reads are kept if they touch the window. */
export function simulate(d: Draws, p: Params): Sim {
  const chip: Read[] = [];
  let pcrCopies = 0, nOrig = 0;
  d.chip.forEach((f, i) => {
    const L = fragLength(p.fragLen, f.z);
    const fs = f.site
      ? Math.round(SITE + FOOT - L + f.u * (L - 2 * FOOT))
      : Math.round(-MARGIN + f.u * (WIN + 2 * MARGIN) - L / 2);
    const copies = p.dups ? f.copies : 0;
    const rd = makeRead(fs, L, f.plus, p, i, false);
    const keep = (x: Read) => x.end > 0 && x.start < WIN;
    if (keep(rd)) { chip.push(rd); nOrig++; }
    for (let c = 0; c < copies; c++) {
      const cp = makeRead(fs, L, f.plus, p, i, true);
      if (keep(cp)) { chip.push(cp); pcrCopies++; }
    }
  });
  let removedPcr = 0, removedIndep = 0;
  if (p.dedup) {
    // keep one read per (strand, 5' end); prefer an original read over a PCR copy
    const seen = new Set<string>();
    const order = chip.map((r, i) => i).sort((a, b) => Number(chip[a].pcr) - Number(chip[b].pcr) || a - b);
    for (const i of order) {
      const r = chip[i], key = (r.plus ? "+" : "-") + r.five;
      if (seen.has(key)) { r.removed = true; if (r.pcr) removedPcr++; else removedIndep++; }
      else seen.add(key);
    }
  }
  const input: Read[] = [];
  d.input.forEach((f, i) => {
    const L = fragLength(p.fragLen, f.z);
    const fs = Math.round(-MARGIN + f.u * (WIN + 2 * MARGIN) - L / 2);
    const rd = makeRead(fs, L, f.plus, p, i, false);
    if (rd.end > 0 && rd.start < WIN) input.push(rd);
  });
  const dupFrac = p.dups && !p.dedup && nOrig > 0 ? pcrCopies / nOrig : 0;
  return { chip, input, pcrCopies, removedPcr, removedIndep,
    chipLib: Math.round(CHIP_LIB * (1 + dupFrac)), inputLib: INPUT_LIB };
}

/** Number of reads overlapping each bin of width `bin` over [0, WIN). */
export function binCounts(reads: Read[], bin: number, which: "all" | "plus" | "minus" = "all"): number[] {
  const nb = Math.ceil(WIN / bin);
  const out = new Array(nb).fill(0);
  for (const r of reads) {
    if (r.removed) continue;
    if (which === "plus" && !r.plus) continue;
    if (which === "minus" && r.plus) continue;
    const b0 = Math.max(0, Math.floor(r.start / bin)), b1 = Math.min(nb - 1, Math.floor((r.end - 1) / bin));
    for (let b = b0; b <= b1; b++) out[b]++;
  }
  return out;
}

export interface Track { values: number[]; chipCounts: number[]; inputScaled: number[]; unit: string }

/** The track in the chosen unit. */
export function track(s: Sim, p: Params): Track {
  const chipCounts = binCounts(s.chip, p.bin);
  const inp = binCounts(s.input, p.bin);
  const scaleIn = s.chipLib / s.inputLib;
  const inputScaled = inp.map(v => v * scaleIn);
  let values: number[], unit: string;
  if (p.norm === "raw") { values = chipCounts.slice(); unit = "reads per bin"; }
  else if (p.norm === "cpm") { values = chipCounts.map(v => v * 1e6 / s.chipLib); unit = "CPM"; }
  else { values = chipCounts.map((v, i) => Math.log2((v + p.pseudo) / (inputScaled[i] + p.pseudo))); unit = "log2 fold change over input"; }
  return { values, chipCounts, inputScaled, unit };
}

export interface Summary {
  nPlus: number; nMinus: number; kept: number;
  maxVal: number; maxPos: number;    // bin centre of the highest bin, relative to the site (bp)
  siteVal: number;                   // value in the bin that contains the site
  halfWidth: number;                 // width (bp) of the contiguous region around the max above half of it
  cPlus: number; cMinus: number;     // centres of the + and − strand read humps, relative to the site (bp)
}

export function summarize(s: Sim, p: Params, t: Track): Summary {
  let nPlus = 0, nMinus = 0;
  for (const r of s.chip) if (!r.removed) { if (r.plus) nPlus++; else nMinus++; }
  const v = t.values;
  let iMax = 0;
  for (let i = 1; i < v.length; i++) if (v[i] > v[iMax]) iMax = i;
  const iSite = Math.floor(SITE / p.bin);
  // half-maximum width (above the track's median, so log2 tracks work too)
  const sorted = v.slice().sort((a, b) => a - b), base = sorted[Math.floor(sorted.length / 2)];
  const half = base + (v[iMax] - base) / 2;
  let lo = iMax, hi = iMax;
  while (lo > 0 && v[lo - 1] >= half) lo--;
  while (hi < v.length - 1 && v[hi + 1] >= half) hi++;
  // centres of the + and − strand read humps: mean read midpoint within ±300 bp of the site
  const centre = (plus: boolean) => {
    let sum = 0, n = 0;
    for (const r of s.chip) {
      if (r.removed || r.plus !== plus) continue;
      const m = (r.start + r.end) / 2 - SITE;
      if (Math.abs(m) <= 300) { sum += m; n++; }
    }
    return n ? sum / n : 0;
  };
  const cPlus = centre(true), cMinus = centre(false);
  return {
    nPlus, nMinus, kept: nPlus + nMinus,
    maxVal: v[iMax], maxPos: Math.round((iMax + 0.5) * p.bin - SITE),
    siteVal: v[iSite], halfWidth: (hi - lo + 1) * p.bin, cPlus: Math.round(cPlus), cMinus: Math.round(cMinus),
  };
}
