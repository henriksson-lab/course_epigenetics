// Simulated two-group bisulfite data and three per-CpG tests (no DOM; testable with node).
//
// For CpG i, group g has a true mean level mu_ig. Each replicate j of group g has its own
// true level p_ijg ~ Beta with mean mu_ig and overdispersion phi (Beta(a, b) with
// a + b = 1/phi - 1, so Var(X) = N p (1-p) [1 + (N-1) phi], as in the unit's Background).
// Reads: N_ijg ~ Poisson(coverage) (at least 1), methylated reads X ~ Binomial(N, p_ijg).
// A fraction of CpGs is truly differential: mu_i2 = mu_i1 +/- delta.
//
// Tests (two-sided p-values):
//  - fisher: Fisher's exact test on reads pooled per group (ignores replicates);
//  - bb:     beta-binomial Wald test, z = (mean beta_2 - mean beta_1) / SE, where the SE uses
//            the beta-binomial variance with one phi estimated from all CpGs (method of
//            moments on the replicates within groups), the limit of DSS's shrinkage;
//  - t:      Student's t-test on the per-sample beta values (2n - 2 df).
import { rng } from "../lib/ui";

type U = () => number;

export interface SimParams {
  cpgs: number;
  coverage: number;
  reps: number;     // replicates per group
  phi: number;      // overdispersion between replicates
  delta: number;    // true |delta beta| of differential CpGs
  frac: number;     // fraction of CpGs that are differential
  seed: number;
}

export interface SimData {
  mu1: Float64Array; mu2: Float64Array;
  isDiff: Uint8Array;
  // per CpG, per sample (group 1 samples 0..reps-1, then group 2): reads and methylated reads
  N: Int32Array[]; X: Int32Array[];
  reps: number;
}

function normal(u: U): number {
  const a = 1 - u(), b = u();
  return Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * b);
}

/** Gamma(shape a, scale 1), Marsaglia–Tsang. */
export function rgamma(a: number, u: U): number {
  if (a < 1) return rgamma(a + 1, u) * Math.pow(u(), 1 / a);
  const d = a - 1 / 3, c = 1 / Math.sqrt(9 * d);
  for (;;) {
    let x: number, v: number;
    do { x = normal(u); v = 1 + c * x; } while (v <= 0);
    v = v * v * v;
    const w = u();
    if (w < 1 - 0.0331 * x ** 4) return d * v;
    if (Math.log(w) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
  }
}

export function rbeta(a: number, b: number, u: U): number {
  const x = rgamma(a, u), y = rgamma(b, u);
  return x + y > 0 ? x / (x + y) : (u() < a / (a + b) ? 1 : 0);
}

export function rpois(lambda: number, u: U): number {
  let k = 0;
  while (lambda > 0) {
    const l = Math.min(lambda, 30);
    lambda -= l;
    const L = Math.exp(-l);
    let p = u();
    while (p > L) { k++; p *= u(); }
  }
  return k;
}

export function rbinom(n: number, p: number, u: U): number {
  if (p <= 0) return 0;
  if (p >= 1) return n;
  let k = 0;
  for (let i = 0; i < n; i++) if (u() < p) k++;
  return k;
}

/** True level of one replicate: Beta with mean mu and overdispersion phi (phi = 0: mu). */
export function sampleLevel(mu: number, phi: number, u: U): number {
  if (phi <= 0) return mu;
  const s = 1 / phi - 1;
  return rbeta(mu * s, (1 - mu) * s, u);
}

/** Group-1 mean levels: mostly high (as in most of the genome), some low, some intermediate. */
function drawMu(u: U): number {
  const r = u();
  if (r < 0.6) return 0.7 + 0.25 * u();
  if (r < 0.8) return 0.02 + 0.1 * u();
  return 0.15 + 0.6 * u();
}

export function simulate(P: SimParams): SimData {
  // CpG-level truth from its own stream, so it does not change with coverage or replicates.
  const rc = rng(P.seed);
  const mu1 = new Float64Array(P.cpgs), mu2 = new Float64Array(P.cpgs), isDiff = new Uint8Array(P.cpgs);
  for (let i = 0; i < P.cpgs; i++) {
    const m = drawMu(rc.unif), ud = rc.unif(), us = rc.unif();
    mu1[i] = m; mu2[i] = m;
    if (ud < P.frac && P.delta > 0) {
      isDiff[i] = 1;
      let s = us < 0.5 ? 1 : -1;
      if (m + s * P.delta > 0.99 || m + s * P.delta < 0.01) s = -s;
      mu2[i] = Math.min(0.99, Math.max(0.01, m + s * P.delta));
    }
  }
  const r = rng(P.seed + 7919);
  const S = 2 * P.reps;
  const N: Int32Array[] = [], X: Int32Array[] = [];
  for (let i = 0; i < P.cpgs; i++) {
    const n = new Int32Array(S), x = new Int32Array(S);
    for (let j = 0; j < S; j++) {
      const mu = j < P.reps ? mu1[i] : mu2[i];
      const p = sampleLevel(mu, P.phi, r.unif);
      n[j] = Math.max(1, rpois(P.coverage, r.unif));
      x[j] = rbinom(n[j], p, r.unif);
    }
    N.push(n); X.push(x);
  }
  return { mu1, mu2, isDiff, N, X, reps: P.reps };
}

// ---------------------------------------------------------------- distributions

let lf = new Float64Array([0]);
/** log(k!) from a growing table. */
export function logFact(k: number): number {
  if (k >= lf.length) {
    const m = Math.max(k + 1, 2 * lf.length);
    const t = new Float64Array(m);
    t.set(lf);
    for (let i = lf.length; i < m; i++) t[i] = t[i - 1] + Math.log(i);
    lf = t;
  }
  return lf[k];
}

/** Two-sided Fisher exact test of the 2x2 table [[a, b], [c, d]] (as R's fisher.test). */
export function fisher(a: number, b: number, c: number, d: number): number {
  const r1 = a + b, r2 = c + d, c1 = a + c, n = r1 + r2;
  const lo = Math.max(0, c1 - r2), hi = Math.min(r1, c1);
  const base = logFact(r1) + logFact(r2) + logFact(c1) + logFact(n - c1) - logFact(n);
  const lp = (k: number) => base - logFact(k) - logFact(r1 - k) - logFact(c1 - k) - logFact(r2 - c1 + k);
  const l0 = lp(a);
  let p = 0;
  for (let k = lo; k <= hi; k++) {
    const l = lp(k);
    if (l <= l0 + 1e-7) p += Math.exp(l);
  }
  return Math.min(1, p);
}

/** Standard normal upper tail. */
export function pnormUpper(z: number): number {
  // Complementary error function, W. J. Cody-style rational approximation (via erfc).
  return 0.5 * erfc(z / Math.SQRT2);
}

function erfc(x: number): number {
  // Numerical Recipes erfcc (fractional error < 1.2e-7), adequate for p-values here.
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
    t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}

function lgamma(x: number): number {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x;
  const tmp = x + 5.5 - (x + 0.5) * Math.log(x + 5.5);
  let ser = 1.000000000190015;
  for (const ci of c) ser += ci / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

/** Regularised incomplete beta I_x(a, b) (continued fraction, Numerical Recipes betacf). */
export function pbeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  const cf = (aa: number, bb: number, xx: number): number => {
    const qab = aa + bb, qap = aa + 1, qam = aa - 1;
    let c = 1, d = 1 - qab * xx / qap;
    if (Math.abs(d) < 1e-300) d = 1e-300;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 300; m++) {
      const m2 = 2 * m;
      let an = m * (bb - m) * xx / ((qam + m2) * (aa + m2));
      d = 1 + an * d; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = 1 + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d; h *= d * c;
      an = -(aa + m) * (qab + m) * xx / ((aa + m2) * (qap + m2));
      d = 1 + an * d; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = 1 + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 1e-12) break;
    }
    return h;
  };
  return x < (a + 1) / (a + b + 2) ? bt * cf(a, b, x) / a : 1 - bt * cf(b, a, 1 - x) / b;
}

/** Two-sided p-value of a t statistic with df degrees of freedom. */
export function ptTwoSided(t: number, df: number): number {
  return pbeta(df / (df + t * t), df / 2, 0.5);
}

/** Benjamini–Hochberg adjusted p-values (NaN stays NaN). */
export function bh(p: Float64Array): Float64Array {
  const idx: number[] = [];
  p.forEach((v, i) => { if (!Number.isNaN(v)) idx.push(i); });
  idx.sort((a, b) => p[a] - p[b]);
  const m = idx.length, q = new Float64Array(p.length).fill(NaN);
  let run = 1;
  for (let r = m - 1; r >= 0; r--) {
    run = Math.min(run, p[idx[r]] * m / (r + 1));
    q[idx[r]] = run;
  }
  return q;
}

// ---------------------------------------------------------------- tests

export function testFisher(d: SimData): Float64Array {
  const n = d.reps, out = new Float64Array(d.N.length);
  d.N.forEach((N, i) => {
    let m1 = 0, n1 = 0, m2 = 0, n2 = 0;
    for (let j = 0; j < n; j++) { m1 += d.X[i][j]; n1 += N[j]; }
    for (let j = n; j < 2 * n; j++) { m2 += d.X[i][j]; n2 += N[j]; }
    out[i] = fisher(m1, n1 - m1, m2, n2 - m2);
  });
  return out;
}

/** Method-of-moments estimate of one common phi from replicates within groups. */
export function estimatePhi(d: SimData): number {
  const n = d.reps;
  if (n < 2) return NaN;
  let num = 0, a = 0, b = 0;
  for (let i = 0; i < d.N.length; i++) {
    for (let g = 0; g < 2; g++) {
      let mb = 0;
      for (let j = g * n; j < (g + 1) * n; j++) mb += d.X[i][j] / d.N[i][j];
      mb /= n;
      if (mb < 0.1 || mb > 0.9) continue;   // stable part of the range
      let ss = 0, inv = 0;
      for (let j = g * n; j < (g + 1) * n; j++) { ss += (d.X[i][j] / d.N[i][j] - mb) ** 2; inv += 1 / d.N[i][j]; }
      // E[ss] = m(1-m)(1 - 1/n) [ sum 1/N_j + phi sum (1 - 1/N_j) ]
      num += ss / (mb * (1 - mb) * (1 - 1 / n)) - inv;
      a += n - inv;
      b++;
    }
  }
  return b > 0 && a > 0 ? Math.max(0, num / a) : 0;
}

export function testBB(d: SimData, phi: number): Float64Array {
  const n = d.reps, out = new Float64Array(d.N.length);
  if (n < 2) return out.fill(NaN);
  const clamp = (m: number) => Math.min(0.995, Math.max(0.005, m));
  d.N.forEach((N, i) => {
    const mean = [0, 0], v = [0, 0];
    for (let g = 0; g < 2; g++) {
      for (let j = g * n; j < (g + 1) * n; j++) mean[g] += d.X[i][j] / N[j];
      mean[g] /= n;
      const m = clamp(mean[g]);
      for (let j = g * n; j < (g + 1) * n; j++) v[g] += m * (1 - m) * (1 + (N[j] - 1) * phi) / N[j];
      v[g] /= n * n;
    }
    const z = (mean[1] - mean[0]) / Math.sqrt(v[0] + v[1]);
    out[i] = 2 * pnormUpper(Math.abs(z));
  });
  return out;
}

export function testT(d: SimData): Float64Array {
  const n = d.reps, out = new Float64Array(d.N.length);
  if (n < 2) return out.fill(NaN);
  d.N.forEach((N, i) => {
    const m = [0, 0];
    let ss = 0;
    for (let g = 0; g < 2; g++) {
      for (let j = g * n; j < (g + 1) * n; j++) m[g] += d.X[i][j] / N[j];
      m[g] /= n;
      for (let j = g * n; j < (g + 1) * n; j++) ss += (d.X[i][j] / N[j] - m[g]) ** 2;
    }
    const df = 2 * n - 2, s2 = ss / df;
    const diff = m[1] - m[0];
    if (s2 <= 0) { out[i] = diff === 0 ? 1 : 0; return; }
    out[i] = ptTwoSided(diff / Math.sqrt(s2 * 2 / n), df);
  });
  return out;
}

export interface Summary { calls: number; trueCalls: number; falseCalls: number; nDiff: number; nNull: number; nullBelow05: number }

/** Calls at BH-FDR < fdr, split by the truth. */
export function summarize(p: Float64Array, isDiff: Uint8Array, fdr = 0.05): Summary {
  const q = bh(p);
  let calls = 0, tc = 0, nd = 0, nb = 0;
  for (let i = 0; i < p.length; i++) {
    if (isDiff[i]) nd++;
    else if (p[i] < 0.05) nb++;
    if (q[i] < fdr) { calls++; if (isDiff[i]) tc++; }
  }
  return { calls, trueCalls: tc, falseCalls: calls - tc, nDiff: nd, nNull: p.length - nd, nullBelow05: nb };
}
