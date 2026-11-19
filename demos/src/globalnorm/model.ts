// globalnorm model (pure, no DOM). Simulated ChIP/ATAC libraries of two conditions.
//
// Per cell, each sample contains some amount of material in every peak, in every background
// bin, and (optionally) a fixed amount of spike-in chromatin. A library is sequenced to its
// own depth: reads are drawn in proportion to the material (Poisson), so the sequencer only
// ever sees *proportions*. Knowing the per-cell amounts, the true size factor of sample j is
//   s_j ∝ depth_j / (total material per cell)_j,
// i.e. dividing counts by s_j gives numbers proportional to the amount per cell.
//
// Truth: n_up peaks go up by UP_LFC in the treated condition; a fraction of the other peaks
// change by the global log2 fold change (e.g. -2: H3K27me3 down to 25% after EZH2 inhibition).
// Technical: the treated libraries can carry extra background (unenriched DNA, dead cells,
// a poorer IP), which lowers their FRiP without any change in the cells' chromatin.
import { rng } from "../lib/ui";
import { rpois, median, bh, tTwoSided } from "./stats";

export const UP_LFC = 2;
export const METHODS = ["total", "mor", "tmm", "bg", "spike"] as const;
export type Method = typeof METHODS[number];
export const METHOD_LABEL: Record<Method, string> = {
  total: "total reads", mor: "median of ratios (peaks)", tmm: "TMM (peaks)",
  bg: "background bins", spike: "spike-in",
};

export interface Params {
  nPeaks: number; nBins: number; reps: number;
  fracChanged: number;   // fraction of the non-up peaks that change globally
  globalLfc: number;     // log2 fold change of those peaks (treated / control)
  nUp: number;           // truly up peaks (log2 FC = UP_LFC)
  bgTreated: number;     // technical extra background in treated libraries (×)
  repCv: number;         // biological replicate CV per peak (log-normal)
  depth: number;         // mean reads per library
  frip0: number;         // FRiP of a control library
  spikeFrac: number;     // spike-in share of the control library's reads
  seed: number;
}

export const DEFAULTS: Params = {
  nPeaks: 2000, nBins: 2000, reps: 3, fracChanged: 0, globalLfc: -2, nUp: 100,
  bgTreated: 1, repCv: 0.15, depth: 2e6, frip0: 0.3, spikeFrac: 0.03, seed: 1,
};

export interface Sim {
  p: Params;
  nS: number;                // samples: first reps control, then reps treated
  cond: number[];            // 0 control, 1 treated
  peaks: Float64Array[];     // [sample][peak] counts
  bins: Float64Array[];      // [sample][bin] counts
  spike: number[];           // spike-in reads per sample
  trueLfc: Float64Array;     // per peak
  kind: Uint8Array;          // 0 unchanged, 1 global change, 2 truly up
  trueSf: number[];          // true size factors (geometric mean 1)
  frip: number[];            // reads in peaks / human reads
}

export function simulate(p: Params): Sim {
  const r = rng(p.seed * 9973 + 17);
  const nS = 2 * p.reps;
  const cond = Array.from({ length: nS }, (_, j) => (j < p.reps ? 0 : 1));
  // peak strengths: log-normal, scaled so that peaks carry frip0 of the control material
  const base = new Float64Array(p.nPeaks);
  for (let i = 0; i < p.nPeaks; i++) base[i] = Math.exp(1.0 * r.norm());
  const sb = base.reduce((a, b) => a + b, 0);
  for (let i = 0; i < p.nPeaks; i++) base[i] *= p.frip0 / sb;
  // background bins: nearly constant (mappability differences only), 1 - frip0 in total
  const bbase = new Float64Array(p.nBins);
  for (let k = 0; k < p.nBins; k++) bbase[k] = Math.exp(0.2 * r.norm());
  const sbb = bbase.reduce((a, b) => a + b, 0);
  for (let k = 0; k < p.nBins; k++) bbase[k] *= (1 - p.frip0) / sbb;
  const spikeMat = p.spikeFrac / (1 - p.spikeFrac);   // so that it is spikeFrac of a control library

  // truth: random order, first nUp are up, then a fraction of the rest change globally
  const order = Array.from({ length: p.nPeaks }, (_, i) => i);
  for (let i = order.length - 1; i > 0; i--) { const j = Math.floor(r.unif() * (i + 1)); [order[i], order[j]] = [order[j], order[i]]; }
  const trueLfc = new Float64Array(p.nPeaks), kind = new Uint8Array(p.nPeaks);
  const nUp = Math.min(p.nUp, p.nPeaks);
  const nGlob = Math.round(p.fracChanged * (p.nPeaks - nUp));
  for (let q = 0; q < p.nPeaks; q++) {
    const i = order[q];
    if (q < nUp) { trueLfc[i] = UP_LFC; kind[i] = 2; }
    else if (q < nUp + nGlob && p.globalLfc !== 0) { trueLfc[i] = p.globalLfc; kind[i] = 1; }
  }

  const peaks: Float64Array[] = [], bins: Float64Array[] = [], spike: number[] = [];
  const sfRaw: number[] = [], frip: number[] = [];
  const cv = p.repCv, s2 = Math.log(1 + cv * cv), sd = Math.sqrt(s2);
  for (let j = 0; j < nS; j++) {
    const t = cond[j];
    const depth = p.depth * Math.exp(0.25 * r.norm());
    const pm = new Float64Array(p.nPeaks), bm = new Float64Array(p.nBins);
    let tot = 0;
    for (let i = 0; i < p.nPeaks; i++) {
      pm[i] = base[i] * Math.pow(2, t * trueLfc[i]) * Math.exp(sd * r.norm() - s2 / 2);
      tot += pm[i];
    }
    const bgF = (t ? p.bgTreated : 1) * Math.exp(0.08 * r.norm());     // technical, per library
    for (let k = 0; k < p.nBins; k++) { bm[k] = bbase[k] * bgF * Math.exp(0.05 * r.norm()); tot += bm[k]; }
    const sp = spikeMat * Math.exp(0.05 * r.norm());                   // pipetting error ~5%
    tot += sp;
    const scaleJ = depth / tot;
    const pc = new Float64Array(p.nPeaks), bc = new Float64Array(p.nBins);
    let rip = 0, hum = 0;
    for (let i = 0; i < p.nPeaks; i++) { pc[i] = rpois(pm[i] * scaleJ, r.unif, r.norm); rip += pc[i]; }
    for (let k = 0; k < p.nBins; k++) { bc[k] = rpois(bm[k] * scaleJ, r.unif, r.norm); hum += bc[k]; }
    hum += rip;
    peaks.push(pc); bins.push(bc);
    spike.push(rpois(sp * scaleJ, r.unif, r.norm));
    sfRaw.push(scaleJ);
    frip.push(rip / hum);
  }
  return { p, nS, cond, peaks, bins, spike, trueLfc, kind, trueSf: geoNorm(sfRaw), frip };
}

/** Scale a vector to geometric mean 1. */
export function geoNorm(v: number[]): number[] {
  const g = Math.exp(v.reduce((a, x) => a + Math.log(x), 0) / v.length);
  return v.map(x => x / g);
}

/** DESeq2 median of ratios on a [sample][feature] matrix (features with a zero are skipped). */
export function medianOfRatios(m: Float64Array[]): number[] {
  const nS = m.length, n = m[0].length;
  const ratios: number[][] = m.map(() => []);
  for (let i = 0; i < n; i++) {
    let lg = 0, ok = true;
    for (let j = 0; j < nS; j++) { if (m[j][i] <= 0) { ok = false; break; } lg += Math.log(m[j][i]); }
    if (!ok) continue;
    lg /= nS;
    for (let j = 0; j < nS; j++) ratios[j].push(Math.log(m[j][i]) - lg);
  }
  return geoNorm(ratios.map(rr => Math.exp(median(rr))));
}

/** edgeR-style TMM (trim 30% of M, 5% of A, precision weights); returns effective-library size factors. */
export function tmm(m: Float64Array[]): number[] {
  const nS = m.length, n = m[0].length;
  const lib = m.map(c => c.reduce((a, b) => a + b, 0));
  // reference: the sample whose upper quartile / library size is closest to the mean
  const uq = m.map((c, j) => quantile(Array.from(c), 0.75) / lib[j]);
  const muq = uq.reduce((a, b) => a + b, 0) / nS;
  let ref = 0;
  uq.forEach((v, j) => { if (Math.abs(v - muq) < Math.abs(uq[ref] - muq)) ref = j; });
  const nf = m.map((c, j) => {
    if (j === ref) return 1;
    const rows: { M: number; A: number; w: number }[] = [];
    const Nj = lib[j], Nr = lib[ref];
    for (let i = 0; i < n; i++) {
      const a = c[i], b = m[ref][i];
      if (a <= 0 || b <= 0) continue;
      const M = Math.log2((a / Nj) / (b / Nr)), A = 0.5 * Math.log2((a / Nj) * (b / Nr));
      const v = (Nj - a) / (Nj * a) + (Nr - b) / (Nr * b);
      rows.push({ M, A, w: 1 / v });
    }
    const k = rows.length;
    const lowM = Math.floor(0.3 * k), lowA = Math.floor(0.05 * k);
    const byM = rows.map((_, i) => i).sort((x, y) => rows[x].M - rows[y].M);
    const byA = rows.map((_, i) => i).sort((x, y) => rows[x].A - rows[y].A);
    const keep = new Uint8Array(k);
    for (let q = lowM; q < k - lowM; q++) keep[byM[q]]++;
    for (let q = lowA; q < k - lowA; q++) keep[byA[q]]++;
    let sw = 0, swm = 0;
    for (let i = 0; i < k; i++) if (keep[i] === 2) { sw += rows[i].w; swm += rows[i].w * rows[i].M; }
    return Math.pow(2, swm / sw);
  });
  const nfn = geoNorm(nf);
  return geoNorm(lib.map((l, j) => l * nfn[j]));
}

export function quantile(v: number[], q: number): number {
  const s = v.slice().sort((a, b) => a - b);
  const h = (s.length - 1) * q, lo = Math.floor(h);
  return s[lo] + (h - lo) * ((s[Math.min(lo + 1, s.length - 1)]) - s[lo]);
}

export function sizeFactors(s: Sim, method: Method): number[] {
  switch (method) {
    case "total": return geoNorm(s.peaks.map((c, j) => c.reduce((a, b) => a + b, 0) + s.bins[j].reduce((a, b) => a + b, 0)));
    case "mor": return medianOfRatios(s.peaks);
    case "tmm": return tmm(s.peaks);
    case "bg": return medianOfRatios(s.bins);
    case "spike": return geoNorm(s.spike.map(x => Math.max(x, 1)));
  }
}

export interface Test {
  A: Float64Array;       // mean log2 normalised count
  lfc: Float64Array;     // treated - control, log2
  padj: number[];
  call: Int8Array;       // +1 up, -1 down, 0 not called
}

/**
 * limma-trend-like test on log2(normalised count + 0.5): per-peak variance pooled within
 * groups, moderated towards a trend in A (10 bins, prior df 4), t-test, BH.
 */
export function test(s: Sim, sf: number[], alpha = 0.05): Test {
  const n = s.p.nPeaks, nS = s.nS, R = s.p.reps;
  const A = new Float64Array(n), lfc = new Float64Array(n), s2 = new Float64Array(n);
  const y = new Float64Array(nS);
  for (let i = 0; i < n; i++) {
    let m0 = 0, m1 = 0;
    for (let j = 0; j < nS; j++) { y[j] = Math.log2(s.peaks[j][i] / sf[j] + 0.5); if (s.cond[j]) m1 += y[j]; else m0 += y[j]; }
    m0 /= R; m1 /= R;
    let ss = 0;
    for (let j = 0; j < nS; j++) ss += (y[j] - (s.cond[j] ? m1 : m0)) ** 2;
    A[i] = (m0 + m1) / 2; lfc[i] = m1 - m0; s2[i] = ss / (nS - 2);
  }
  const d = nS - 2, d0 = 4;
  // prior variance: median of s2 in 10 bins of A, scaled from the chi-square median to the mean
  const chiMed = d * Math.pow(1 - 2 / (9 * d), 3) / d;
  const ord = Array.from({ length: n }, (_, i) => i).sort((a, b) => A[a] - A[b]);
  const prior = new Float64Array(n);
  const NB = 10;
  for (let b = 0; b < NB; b++) {
    const lo = Math.floor(b * n / NB), hi = Math.floor((b + 1) * n / NB);
    const idx = ord.slice(lo, hi);
    const pv = median(idx.map(i => s2[i])) / chiMed;
    for (const i of idx) prior[i] = pv;
  }
  const pv: number[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const sm = (d0 * prior[i] + d * s2[i]) / (d0 + d);
    const t = lfc[i] / Math.sqrt(sm * (2 / R));
    pv[i] = tTwoSided(t, d0 + d);
  }
  const padj = bh(pv);
  const call = new Int8Array(n);
  for (let i = 0; i < n; i++) if (padj[i] < alpha) call[i] = lfc[i] > 0 ? 1 : -1;
  return { A, lfc, padj, call };
}

export interface Tally { upTrue: number; upFalse: number; downTrue: number; downFalse: number; missed: number; changed: number; medLfc: number }

export function tally(s: Sim, t: Test): Tally {
  const o: Tally = { upTrue: 0, upFalse: 0, downTrue: 0, downFalse: 0, missed: 0, changed: 0, medLfc: median(t.lfc) };
  for (let i = 0; i < s.p.nPeaks; i++) {
    const tr = s.trueLfc[i];
    if (tr !== 0) o.changed++;
    if (t.call[i] > 0) { if (tr > 0) o.upTrue++; else o.upFalse++; }
    else if (t.call[i] < 0) { if (tr < 0) o.downTrue++; else o.downFalse++; }
    else if (tr !== 0) o.missed++;
  }
  return o;
}
