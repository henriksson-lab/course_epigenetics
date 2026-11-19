// A/B compartments on a toy chromosome, and the eigenvector recipe that finds them.
// No DOM, so node tests can import it.
//
// Simulated contact matrix (bins i != j, s = |i - j|):
//   O_ij = s^-alpha × (sqrt(r) if bins i and j are the same type, 1/sqrt(r) if not) × noise_ij
// so same-type pairs contact r times more than A–B pairs at the same distance; noise_ij is
// log-normal, exp(sd·z_ij − sd²/2) with z_ij standard normal (mean 1). The recipe:
//   expected E(s) = mean of diagonal s; O/E = O / E(s);
//   C = Pearson correlation between the rows of O/E (or of O, to show why O/E is needed);
//   E1 = eigenvector of C with the largest eigenvalue (power iteration);
//   sign: flip E1 so that it correlates positively with GC content; E1 > 0 is called A.

import { rng } from "../lib/ui";

export type Mat = Float64Array;          // n*n, row-major

/**
 * Starting patterns: "blocks" and "fine" are run lengths alternating A, B, A, ... (summing
 * to 48, scaled to n bins; "fine" is irregular on purpose: a strictly periodic pattern lines up
 * with the diagonals, and dividing by the expected value then removes the signal itself).
 */
export function pattern(kind: string, n: number, seed: number): boolean[] {
  if (kind === "island") return Array.from({ length: n }, (_, i) => i >= Math.round(n * 0.42) && i < Math.round(n * 0.58));
  if (kind === "random") { const r = rng(seed + 100); return Array.from({ length: n }, () => r.unif() < 0.5); }
  const runs = kind === "fine" ? [3, 2, 4, 2, 3, 3, 2, 4, 3, 2, 3, 2, 4, 3, 2, 3, 3] : [6, 8, 4, 6, 10, 5, 3, 6];
  const out: boolean[] = [];
  runs.forEach((len, k) => { for (let q = 0; q < Math.round(len * n / 48); q++) out.push(k % 2 === 0); });
  while (out.length < n) out.push(out[out.length - 1]);
  return out.slice(0, n);
}

/** GC content per bin: 0.43 in A, 0.37 in B, plus noise (gcNoise: standard normal per bin). */
export function gcTrack(isA: boolean[], gcNoise: number[]): number[] {
  return isA.map((a, i) => 0.40 + (a ? 0.03 : -0.03) + 0.018 * gcNoise[i]);
}

export interface SimInput {
  isA: boolean[];       // true type of each bin
  ratio: number;        // r: same-type / different-type contact ratio (>= 1)
  alpha: number;        // distance-decay exponent
  noise: number;        // sd of log noise
  z: Mat;               // standard normal numbers, symmetric, n*n (fixed per seed)
}

/** The simulated (balanced) contact matrix. The diagonal holds E(1), the s = 1 level. */
export function simulate(inp: SimInput): Mat {
  const n = inp.isA.length, O = new Float64Array(n * n);
  const hi = Math.sqrt(inp.ratio), lo = 1 / hi, sd = inp.noise;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const s = j - i;
    const v = Math.pow(s, -inp.alpha) * (inp.isA[i] === inp.isA[j] ? hi : lo) * Math.exp(sd * inp.z[i * n + j] - sd * sd / 2);
    O[i * n + j] = v; O[j * n + i] = v;
  }
  const e = expected(O, n);
  for (let i = 0; i < n; i++) O[i * n + i] = e[1];
  return O;
}

/** Symmetric matrix of standard normal numbers from a uniform generator (Box–Muller). */
export function normals(n: number, unif: () => number): Mat {
  const z = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    const u = 1 - unif(), v = unif();
    const x = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    z[i * n + j] = x; z[j * n + i] = x;
  }
  return z;
}

/** E(s): mean of the off-diagonal s (s = 1..n-1); E[0] = 0. */
export function expected(M: Mat, n: number): Float64Array {
  const e = new Float64Array(n);
  for (let s = 1; s < n; s++) {
    let sum = 0;
    for (let i = 0; i + s < n; i++) sum += M[i * n + i + s];
    e[s] = sum / (n - s);
  }
  return e;
}

/** Observed / expected; the diagonal is set to 1. */
export function obsExp(M: Mat, n: number): Mat {
  const e = expected(M, n), R = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++)
    R[i * n + j] = i === j ? 1 : M[i * n + j] / e[Math.abs(i - j)];
  return R;
}

/** Pearson correlation between the rows of M. */
export function corrRows(M: Mat, n: number): Mat {
  const Z = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    let m = 0;
    for (let j = 0; j < n; j++) m += M[i * n + j];
    m /= n;
    let ss = 0;
    for (let j = 0; j < n; j++) { const d = M[i * n + j] - m; Z[i * n + j] = d; ss += d * d; }
    const sd = Math.sqrt(ss) || 1;
    for (let j = 0; j < n; j++) Z[i * n + j] /= sd;
  }
  const C = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    let s = 0;
    for (let k = 0; k < n; k++) s += Z[i * n + k] * Z[j * n + k];
    C[i * n + j] = s; C[j * n + i] = s;
  }
  return C;
}

export interface Eig { v: Float64Array; lambda: number; iters: number; converged: boolean }

/**
 * Power iteration: v <- C v / |C v| from the start vector v0, until v changes by less than
 * tol (or maxIt steps). For a symmetric positive semi-definite C (a correlation matrix) this
 * converges to the eigenvector with the largest eigenvalue, with an arbitrary sign that
 * depends on v0. lambda = v' C v.
 */
export function powerIteration(C: Mat, n: number, v0: ArrayLike<number>, tol = 1e-10, maxIt = 2000): Eig {
  let v = Float64Array.from(v0);
  normalise(v);
  const w = new Float64Array(n);
  let iters = 0, converged = false;
  while (iters < maxIt) {
    iters++;
    matVec(C, n, v, w);
    normalise(w);
    let d = 0;
    for (let i = 0; i < n; i++) d = Math.max(d, Math.abs(w[i] - v[i]));
    v = Float64Array.from(w);
    if (d < tol) { converged = true; break; }
  }
  matVec(C, n, v, w);
  let lambda = 0;
  for (let i = 0; i < n; i++) lambda += v[i] * w[i];
  return { v, lambda, iters, converged };
}

/** Second eigenvalue by deflation: power iteration on C - lambda1 v1 v1'. */
export function secondEigenvalue(C: Mat, n: number, e1: Eig, v0: ArrayLike<number>): number {
  const D = Float64Array.from(C);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) D[i * n + j] -= e1.lambda * e1.v[i] * e1.v[j];
  return powerIteration(D, n, v0, 1e-9, 3000).lambda;
}

function matVec(C: Mat, n: number, v: ArrayLike<number>, out: Float64Array): void {
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += C[i * n + j] * v[j];
    out[i] = s;
  }
}
function normalise(v: Float64Array): void {
  let s = 0;
  for (let i = 0; i < v.length; i++) s += v[i] * v[i];
  s = Math.sqrt(s) || 1;
  for (let i = 0; i < v.length; i++) v[i] /= s;
}

/** Pearson correlation of two vectors. */
export function pearson(a: ArrayLike<number>, b: ArrayLike<number>): number {
  const n = a.length;
  let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) {
    const x = a[i] - ma, y = b[i] - mb;
    sab += x * y; saa += x * x; sbb += y * y;
  }
  return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : 0;
}

/**
 * Saddle-style compartment strength from O/E and a calling: (mean O/E of A–A pairs + mean
 * of B–B pairs) / (2 × mean of A–B pairs), over pairs at least minS bins apart.
 */
export function strength(OE: Mat, n: number, callA: boolean[], minS = 2): number {
  let aa = 0, na = 0, bb = 0, nb = 0, ab = 0, nab = 0;
  for (let i = 0; i < n; i++) for (let j = i + minS; j < n; j++) {
    const v = OE[i * n + j];
    if (callA[i] && callA[j]) { aa += v; na++; }
    else if (!callA[i] && !callA[j]) { bb += v; nb++; }
    else { ab += v; nab++; }
  }
  if (!na || !nb || !nab) return NaN;
  return (aa / na + bb / nb) / (2 * ab / nab);
}

/** The full recipe. */
export interface Result {
  O: Mat; OE: Mat; C: Mat; e1: Eig; lambda2: number;
  E1: Float64Array;        // after orientation (if requested)
  flipped: boolean;        // true if the solver's sign was flipped
  rGC: number;             // corr(E1 as shown, GC)
  callA: boolean[];        // E1 > 0
  agree: number;           // fraction of bins whose call matches the true type
  strength: number;        // (AA + BB) / (AB + BA) from O/E with the calls
}

export function analyse(inp: SimInput, gc: ArrayLike<number>, v0: ArrayLike<number>,
  opts: { skipOE: boolean; orient: boolean }): Result {
  const n = inp.isA.length;
  const O = simulate(inp);
  const OE = obsExp(O, n);
  const C = corrRows(opts.skipOE ? O : OE, n);
  const e1 = powerIteration(C, n, v0);
  const lambda2 = secondEigenvalue(C, n, e1, Array.from(v0, (x, i) => x + 0.37 * Math.cos(i)));
  const r0 = pearson(e1.v, gc);
  const flipped = opts.orient && r0 < 0;
  const E1 = Float64Array.from(e1.v, x => (flipped ? -x : x));
  const callA = Array.from(E1, x => x > 0);
  let k = 0;
  for (let i = 0; i < n; i++) if (callA[i] === inp.isA[i]) k++;
  return { O, OE, C, e1, lambda2, E1, flipped, rGC: flipped ? -r0 : r0, callA, agree: k / n,
    strength: strength(OE, n, callA) };
}
