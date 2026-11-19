// motif2x2 model (pure, no DOM).
// Simulated peak sequences: each peak is L bases drawn independently with its own GC content.
// Target peaks are GC-richer than the background pool by `gcShift`. A PU.1-like motif is
// planted in a fraction of the target peaks (and in 5% of the background: bound sites occur in
// other peaks too). A GC-box (SP1-like) motif is planted nowhere: any enrichment it shows is
// base composition. Peaks are scanned on both strands with integer log-odds PWMs (log2 vs equal
// base frequencies, ×100); a peak "has a match" when its best window reaches the score whose
// per-window p-value (under equal base frequencies, exact by dynamic programming) is ≤ pThr.
import { rng } from "../lib/ui";

export const BASES = "ACGT";
/** Probability matrices, rows = positions, columns = A C G T. Approximations of JASPAR-like logos. */
export const MOTIFS: { key: string; name: string; note: string; prob: number[][] }[] = [
  {
    key: "pu1", name: "PU.1-like (ETS)", note: "bound in target peaks", prob: [
      [0.70, 0.08, 0.12, 0.10], [0.72, 0.08, 0.08, 0.12], [0.68, 0.07, 0.17, 0.08], [0.15, 0.04, 0.77, 0.04],
      [0.85, 0.04, 0.07, 0.04], [0.03, 0.02, 0.93, 0.02], [0.02, 0.02, 0.94, 0.02], [0.94, 0.02, 0.02, 0.02],
      [0.88, 0.03, 0.05, 0.04], [0.08, 0.12, 0.62, 0.18]],
  },
  {
    key: "sp1", name: "GC box (SP1-like)", note: "not bound anywhere", prob: [
      [0.10, 0.10, 0.70, 0.10], [0.05, 0.05, 0.85, 0.05], [0.03, 0.03, 0.90, 0.04], [0.05, 0.10, 0.80, 0.05],
      [0.02, 0.90, 0.05, 0.03], [0.02, 0.03, 0.92, 0.03], [0.03, 0.03, 0.90, 0.04], [0.05, 0.10, 0.80, 0.05],
      [0.10, 0.10, 0.70, 0.10], [0.15, 0.40, 0.25, 0.20]],
  },
];

export interface Params {
  nTarget: number; nPool: number; len: number;
  gcBg: number;        // mean GC of the background pool
  gcShift: number;     // target mean GC − background mean GC
  fracBound: number;   // fraction of target peaks with a planted PU.1 site
  fracBoundBg: number; // fraction of background peaks with a planted site
  matchPerTarget: number; // GC-matched background peaks per target peak
  seed: number;
}
export const DEFAULTS: Params = {
  nTarget: 1000, nPool: 8000, len: 300, gcBg: 0.44, gcShift: 0.08, fracBound: 0.25, fracBoundBg: 0.05,
  matchPerTarget: 2, seed: 1,
};

export function intMatrix(prob: number[][]): number[][] {
  return prob.map(row => row.map(p => Math.round(100 * Math.log2(p / 0.25))));
}

/** Exact upper-tail distribution of the window score under equal base frequencies. */
export function scoreTail(M: number[][]): { min: number; tail: Float64Array } {
  let min = 0, max = 0;
  for (const row of M) { min += Math.min(...row); max += Math.max(...row); }
  let dist = new Float64Array(max - min + 1);
  // start: score 0 at offset -min
  let cur = new Map<number, number>([[0, 1]]);
  for (const row of M) {
    const nxt = new Map<number, number>();
    for (const [s, pr] of cur) for (const v of row) nxt.set(s + v, (nxt.get(s + v) ?? 0) + pr / 4);
    cur = nxt;
  }
  for (const [s, pr] of cur) dist[s - min] += pr;
  const tail = new Float64Array(dist.length);
  let acc = 0;
  for (let k = dist.length - 1; k >= 0; k--) { acc += dist[k]; tail[k] = acc; }
  dist = new Float64Array(0);
  return { min, tail };
}

/** Smallest integer score whose per-window P(S ≥ score) ≤ p. */
export function scoreForP(st: { min: number; tail: Float64Array }, p: number): number {
  for (let k = 0; k < st.tail.length; k++) if (st.tail[k] <= p) return st.min + k;
  return st.min + st.tail.length;
}

/** Best window score on either strand; M and Mrc are flat [position * 4 + base] integer matrices. */
function maxScore(seq: Uint8Array, L: number, M: Int32Array, Mrc: Int32Array): number {
  const w = M.length >> 2;
  let best = -1e9;
  for (let i = 0; i + w <= L; i++) {
    let f = 0, r = 0;
    for (let k = 0; k < w; k++) { const q = (k << 2) | seq[i + k]; f += M[q]; r += Mrc[q]; }
    if (f > best) best = f;
    if (r > best) best = r;
  }
  return best;
}

export interface Data {
  p: Params;
  gcT: Float64Array; gcP: Float64Array;            // realised GC per peak
  // per motif: best score without / with the planted site, per peak
  tScore: Int32Array[]; tScorePl: Int32Array[];     // target
  pScore: Int32Array[];                             // pool (background planted already applied)
  plantU: Float64Array;                             // target planted iff plantU < fracBound
  matched: Int32Array;                              // indices into the pool: GC-matched background
  tails: { min: number; tail: Float64Array }[];
}

export function simulate(p: Params): Data {
  const r = rng(p.seed * 7907 + 3);
  const L = p.len, nT = p.nTarget, nP = p.nPool;
  const Ms = MOTIFS.map(m => intMatrix(m.prob));
  const flat = Ms.map(M => Int32Array.from(M.flat()));
  const flatRc = Ms.map(M => Int32Array.from(M.slice().reverse().map(row => [row[3], row[2], row[1], row[0]]).flat()));
  const clamp = (g: number) => Math.min(0.85, Math.max(0.2, g));
  const gcT = new Float64Array(nT), gcP = new Float64Array(nP);
  const seq = new Uint8Array(L);
  const draw = (gc: number): number => {
    let n = 0;
    for (let i = 0; i < L; i++) {
      const u = r.unif();
      const b = u < gc / 2 ? 1 : u < gc ? 2 : u < gc + (1 - gc) / 2 ? 0 : 3;
      seq[i] = b; if (b === 1 || b === 2) n++;
    }
    return n / L;
  };
  const plant = (m: number[][]): void => {
    const w = m.length, pos = Math.floor(r.unif() * (L - w + 1)), rc = r.unif() < 0.5;
    for (let k = 0; k < w; k++) {
      const row = m[k];
      let u = r.unif(), b = 0;
      while (b < 3 && u > row[b]) { u -= row[b]; b++; }
      if (rc) seq[pos + w - 1 - k] = 3 - b; else seq[pos + k] = b;
    }
  };
  const nM = MOTIFS.length;
  const tScore = MOTIFS.map(() => new Int32Array(nT)), tScorePl = MOTIFS.map(() => new Int32Array(nT));
  const pScore = MOTIFS.map(() => new Int32Array(nP));
  const plantU = new Float64Array(nT);
  for (let i = 0; i < nT; i++) {
    gcT[i] = draw(clamp(p.gcBg + p.gcShift + 0.04 * r.norm()));
    for (let m = 0; m < nM; m++) tScore[m][i] = maxScore(seq, L, flat[m], flatRc[m]);
    plantU[i] = r.unif();
    plant(MOTIFS[0].prob);
    for (let m = 0; m < nM; m++) tScorePl[m][i] = maxScore(seq, L, flat[m], flatRc[m]);
  }
  for (let i = 0; i < nP; i++) {
    gcP[i] = draw(clamp(p.gcBg + 0.07 * r.norm()));
    if (r.unif() < p.fracBoundBg) plant(MOTIFS[0].prob);
    for (let m = 0; m < nM; m++) pScore[m][i] = maxScore(seq, L, flat[m], flatRc[m]);
  }
  // GC-matched background: for each target, matchPerTarget pool peaks from the same 1% GC bin
  const bins = new Map<number, number[]>();
  for (let i = 0; i < nP; i++) { const b = Math.round(gcP[i] * 100); if (!bins.has(b)) bins.set(b, []); bins.get(b)!.push(i); }
  const matched: number[] = [];
  for (let i = 0; i < nT; i++) {
    const b0 = Math.round(gcT[i] * 100);
    let pool: number[] | undefined;
    for (let d = 0; d < 30 && !pool; d++) pool = bins.get(b0 + d) ?? bins.get(b0 - d);
    if (!pool) continue;
    for (let k = 0; k < p.matchPerTarget; k++) matched.push(pool[Math.floor(r.unif() * pool.length)]);
  }
  return { p, gcT, gcP, tScore, tScorePl, pScore, plantU, matched: Int32Array.from(matched), tails: Ms.map(scoreTail) };
}

/** Best score of target peak i for motif m, given the bound fraction. */
export function targetScore(d: Data, m: number, i: number, fracBound: number): number {
  return d.plantU[i] < fracBound ? d.tScorePl[m][i] : d.tScore[m][i];
}

export interface Table { a: number; b: number; c: number; d: number; or: number; log10p: number; expected: number }

/** 2×2 table for motif m at score threshold thr, against the whole pool or the GC-matched set. */
export function table(d: Data, m: number, thr: number, fracBound: number, gcMatched: boolean): Table {
  let a = 0;
  const nT = d.p.nTarget;
  for (let i = 0; i < nT; i++) if (targetScore(d, m, i, fracBound) >= thr) a++;
  let c = 0, nB = 0;
  if (gcMatched) { for (const j of d.matched) { nB++; if (d.pScore[m][j] >= thr) c++; } }
  else { nB = d.p.nPool; for (let j = 0; j < nB; j++) if (d.pScore[m][j] >= thr) c++; }
  return fill(a, nT - a, c, nB - c);
}

export function fill(a: number, b: number, c: number, dd: number): Table {
  const N = a + b + c + dd, K = a + c, n = a + b;
  return { a, b, c, d: dd, or: (a * dd) / (b * c), log10p: hyperUpperLog10(a, K, N, n), expected: n * K / N };
}

let lf: Float64Array = new Float64Array([0, 0]);
function logFact(n: number): number {
  if (n >= lf.length) {
    const nl = new Float64Array(Math.max(n + 1, 2 * lf.length));
    nl.set(lf);
    for (let k = lf.length; k < nl.length; k++) nl[k] = nl[k - 1] + Math.log(k);
    lf = nl;
  }
  return lf[n];
}
const lchoose = (n: number, k: number) => logFact(n) - logFact(k) - logFact(n - k);

/**
 * log10 P(X ≥ a) for X ~ Hypergeometric(N total, K successes, n draws): one-sided Fisher's test
 * for enrichment, summed exactly in log space (log-sum-exp), so tiny p-values do not underflow.
 */
export function hyperUpperLog10(a: number, K: number, N: number, n: number): number {
  const hi = Math.min(K, n), lo = Math.max(a, 0, n - (N - K));
  if (lo > hi) return 0;
  const den = lchoose(N, n);
  let mx = -Infinity;
  const terms: number[] = [];
  for (let x = lo; x <= hi; x++) { const t = lchoose(K, x) + lchoose(N - K, n - x) - den; terms.push(t); if (t > mx) mx = t; }
  let s = 0;
  for (const t of terms) s += Math.exp(t - mx);
  return Math.min(0, (mx + Math.log(s)) / Math.LN10);
}

export function meanOf(v: ArrayLike<number>, idx?: ArrayLike<number>): number {
  let s = 0;
  const n = idx ? idx.length : v.length;
  for (let k = 0; k < n; k++) s += v[idx ? idx[k] : k];
  return s / n;
}
