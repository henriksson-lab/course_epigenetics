// Simulated "epigenetic clock" data and an elastic net fitted by coordinate descent
// (no DOM; testable with node). The elastic net minimises, as glmnet documents,
//   (1 / 2n) sum_i (y_i - a - sum_j b_j x_ij)^2 + lambda [ alpha sum_j |b_j| + (1 - alpha)/2 sum_j b_j^2 ]
// on standardised features (mean 0, variance 1 with the 1/n convention), and the
// coefficients are reported on the original scale. (glmnet also scales y internally, which
// makes its ridge part differ by a factor sd(y) when alpha < 1; for y with sd 1 the two agree.)
import { rng } from "../lib/ui";

// ---------------------------------------------------------------- simulation

export interface ClockSim {
  ages: Float64Array;   // all people: training first, then test
  nTrainMax: number;
  nTest: number;
  pMax: number;
  kMax: number;         // the first kMax CpGs may be age-associated
  base: Float64Array;   // pMax baseline beta (at age 50)
  slope: Float64Array;  // pMax change in beta per decade (0 for noise CpGs)
  eps: Float64Array;    // people x pMax standard normal noise, row-major
}

/** Draw ages, CpG baselines, slopes and the noise once; the sliders then only rescale. */
export function simulate(seed: number, nTrainMax: number, nTest: number, pMax: number, kMax: number): ClockSim {
  const r = rng(seed);
  const n = nTrainMax + nTest;
  const ages = new Float64Array(n);
  for (let i = 0; i < n; i++) ages[i] = 20 + 60 * r.unif();
  const base = new Float64Array(pMax), slope = new Float64Array(pMax);
  for (let j = 0; j < pMax; j++) base[j] = 0.25 + 0.5 * r.unif();
  for (let j = 0; j < kMax; j++) slope[j] = (r.unif() < 0.6 ? 1 : -1) * (0.008 + 0.022 * r.unif());
  const eps = new Float64Array(n * pMax);
  for (let i = 0; i < eps.length; i++) eps[i] = r.norm();
  return { ages, nTrainMax, nTest, pMax, kMax, base, slope, eps };
}

export interface Data {
  X: Float64Array[];    // columns: X[j][i], i over people
  y: Float64Array;
  train: number[];
  test: number[];
  isAge: boolean[];     // true for CpGs whose mean changes with age
}

/**
 * Build the data for given slider values: nTrain training people (the first nTrain of the
 * training pool), the fixed test people, p CpGs of which the first k change with age,
 * noise SD sigma (in beta units) for every CpG.
 */
export function buildData(s: ClockSim, nTrain: number, p: number, k: number, sigma: number): Data {
  const people: number[] = [];
  for (let i = 0; i < nTrain; i++) people.push(i);
  for (let i = 0; i < s.nTest; i++) people.push(s.nTrainMax + i);
  const n = people.length;
  const X: Float64Array[] = [];
  const isAge: boolean[] = [];
  for (let j = 0; j < p; j++) {
    const col = new Float64Array(n);
    const sl = j < k ? s.slope[j] : 0;
    for (let ii = 0; ii < n; ii++) {
      const i = people[ii];
      col[ii] = s.base[j] + sl * (s.ages[i] - 50) / 10 + sigma * s.eps[i * s.pMax + j];
    }
    X.push(col);
    isAge.push(j < k);
  }
  const y = Float64Array.from(people, i => s.ages[i]);
  const train = Array.from({ length: nTrain }, (_, i) => i);
  const test = Array.from({ length: s.nTest }, (_, i) => nTrain + i);
  return { X, y, train, test, isAge };
}

/** Pearson correlation of column x with y over the rows in idx. */
export function corr(x: Float64Array, y: Float64Array, idx: number[]): number {
  let mx = 0, my = 0;
  for (const i of idx) { mx += x[i]; my += y[i]; }
  mx /= idx.length; my /= idx.length;
  let sxy = 0, sxx = 0, syy = 0;
  for (const i of idx) { const a = x[i] - mx, b = y[i] - my; sxy += a * b; sxx += a * a; syy += b * b; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

/** Indices of the m columns with the largest |correlation| with y over the rows in idx. */
export function topCorrelated(X: Float64Array[], y: Float64Array, idx: number[], m: number): number[] {
  const r = X.map(col => Math.abs(corr(col, y, idx)));
  return r.map((v, j) => [v, j]).sort((a, b) => b[0] - a[0]).slice(0, m).map(v => v[1]);
}

// ---------------------------------------------------------------- elastic net

export interface Fit {
  lambda: number;
  a0: number;            // intercept, original scale
  beta: Float64Array;    // coefficients, original scale (one per column passed in)
  nonzero: number;
  sweeps: number;
}

/** Standardised training design, kept between fits along a path. */
export interface Prepared {
  n: number;
  Z: Float64Array[];     // standardised training columns
  mean: Float64Array;
  sd: Float64Array;
  ymean: number;
  yc: Float64Array;      // centred training y
}

export function prepare(cols: Float64Array[], y: Float64Array, rows: number[]): Prepared {
  const n = rows.length, p = cols.length;
  const Z: Float64Array[] = [], mean = new Float64Array(p), sd = new Float64Array(p);
  for (let j = 0; j < p; j++) {
    const c = cols[j];
    let m = 0;
    for (const i of rows) m += c[i];
    m /= n;
    let v = 0;
    for (const i of rows) v += (c[i] - m) ** 2;
    const s = Math.sqrt(v / n);
    const z = new Float64Array(n);
    if (s > 0) rows.forEach((i, ii) => { z[ii] = (c[i] - m) / s; });
    Z.push(z); mean[j] = m; sd[j] = s;
  }
  let ym = 0;
  for (const i of rows) ym += y[i];
  ym /= n;
  const yc = Float64Array.from(rows, i => y[i] - ym);
  return { n, Z, mean, sd, ymean: ym, yc };
}

/** Smallest lambda at which all coefficients are 0 (for alpha > 0). */
export function lambdaMax(pr: Prepared, alpha: number): number {
  let m = 0;
  for (const z of pr.Z) {
    let s = 0;
    for (let i = 0; i < pr.n; i++) s += z[i] * pr.yc[i];
    m = Math.max(m, Math.abs(s) / pr.n);
  }
  return m / Math.max(alpha, 1e-3);
}

const soft = (z: number, g: number): number => (z > g ? z - g : z < -g ? z + g : 0);

/**
 * Coordinate descent at one lambda, starting from the standardised coefficients b0
 * (warm start; modified in place and returned in the Fit). Cycles over all coordinates,
 * then over the active set until it converges, then checks all coordinates again.
 */
export function enetFit(pr: Prepared, lambda: number, alpha: number, b0?: Float64Array,
  tol = 1e-6, maxSweeps = 1000): Fit & { b: Float64Array } {
  const { n, Z, yc } = pr, p = Z.length;
  const b = b0 ?? new Float64Array(p);
  const r = Float64Array.from(yc);
  for (let j = 0; j < p; j++) if (b[j] !== 0) { const z = Z[j]; for (let i = 0; i < n; i++) r[i] -= z[i] * b[j]; }
  const l1 = lambda * alpha, den = 1 + lambda * (1 - alpha);
  let sweeps = 0;
  const ys = Math.max(1e-12, yc.reduce((s, v) => s + v * v, 0) / n);
  const update = (j: number): number => {
    if (pr.sd[j] === 0) return 0;
    const z = Z[j];
    let g = 0;
    for (let i = 0; i < n; i++) g += z[i] * r[i];
    const old = b[j];
    const nb = soft(g / n + old, l1) / den;
    if (nb === old) return 0;
    const d = nb - old;
    for (let i = 0; i < n; i++) r[i] -= z[i] * d;
    b[j] = nb;
    return d * d;   // change in fit, in units of variance (z has variance 1)
  };
  for (;;) {
    let dmax = 0;
    for (let j = 0; j < p; j++) dmax = Math.max(dmax, update(j));
    sweeps++;
    if (dmax < tol * ys || sweeps >= maxSweeps) break;
    // iterate on the active set
    const act: number[] = [];
    for (let j = 0; j < p; j++) if (b[j] !== 0) act.push(j);
    for (;;) {
      let dm = 0;
      for (const j of act) dm = Math.max(dm, update(j));
      sweeps++;
      if (dm < tol * ys || sweeps >= maxSweeps) break;
    }
  }
  const beta = new Float64Array(p);
  let a0 = pr.ymean, nz = 0;
  for (let j = 0; j < p; j++) {
    if (b[j] === 0 || pr.sd[j] === 0) continue;
    beta[j] = b[j] / pr.sd[j];
    a0 -= beta[j] * pr.mean[j];
    nz++;
  }
  return { lambda, a0, beta, nonzero: nz, sweeps, b };
}

export function predict(fit: Fit, cols: Float64Array[], rows: number[]): Float64Array {
  const out = new Float64Array(rows.length).fill(fit.a0);
  for (let j = 0; j < cols.length; j++) {
    const bj = fit.beta[j];
    if (bj === 0) continue;
    const c = cols[j];
    rows.forEach((i, ii) => { out[ii] += bj * c[i]; });
  }
  return out;
}

export function median(v: ArrayLike<number>): number {
  const s = Array.from(v).sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

/** Median absolute error of predictions p against the truth y[rows]. */
export function mae(p: ArrayLike<number>, y: Float64Array, rows: number[]): number {
  return median(rows.map((i, ii) => Math.abs(p[ii] - y[i])));
}

/** Log-spaced lambda grid from hi down to lo. */
export function lambdaGrid(lo: number, hi: number, m: number): number[] {
  return Array.from({ length: m }, (_, i) => hi * Math.pow(lo / hi, i / (m - 1)));
}

export interface PathPoint { lambda: number; trainMae: number; testMae: number; testR: number; nonzero: number; nAge: number; fit: Fit }

/**
 * Fit the elastic net along a lambda grid (largest first, warm starts) on the selected
 * columns, using only the training rows; evaluate on training and test rows.
 */
export function path(d: Data, cols: number[], alpha: number, grid: number[]): PathPoint[] {
  const C = cols.map(j => d.X[j]);
  const pr = prepare(C, d.y, d.train);
  let b: Float64Array | undefined;
  const out: PathPoint[] = [];
  const ytest = d.test.map(i => d.y[i]);
  for (const lam of grid) {
    const f = enetFit(pr, lam, alpha, b);
    b = Float64Array.from(f.b);
    const pt = predict(f, C, d.train), pe = predict(f, C, d.test);
    let nAge = 0;
    cols.forEach((j, jj) => { if (f.beta[jj] !== 0 && d.isAge[j]) nAge++; });
    out.push({ lambda: lam, trainMae: mae(pt, d.y, d.train), testMae: mae(pe, d.y, d.test),
      testR: pearson(Array.from(pe), ytest), nonzero: f.nonzero, nAge, fit: f });
  }
  return out;
}

export function pearson(a: number[], b: number[]): number {
  const n = a.length;
  const ma = a.reduce((s, v) => s + v, 0) / n, mb = b.reduce((s, v) => s + v, 0) / n;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) ** 2; sbb += (b[i] - mb) ** 2; }
  return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : 0;
}
