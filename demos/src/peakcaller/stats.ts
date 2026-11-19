// Small statistics helpers shared by the peakcaller, crosscorr and fingerprint demos
// (no DOM, testable with node).

const LNFACT: number[] = [0];
for (let i = 1; i <= 256; i++) LNFACT.push(LNFACT[i - 1] + Math.log(i));

/** ln(n!) (table up to 256, Stirling series above). */
export function lnFact(n: number): number {
  if (n <= 256) return LNFACT[n];
  const x = n + 1;
  return (x - 0.5) * Math.log(x) - x + 0.5 * Math.log(2 * Math.PI) + 1 / (12 * x) - 1 / (360 * x ** 3);
}

/** ln P(X = k) for X ~ Poisson(lam). */
export function lnPoisPmf(k: number, lam: number): number {
  if (lam <= 0) return k === 0 ? 0 : -Infinity;
  return k * Math.log(lam) - lam - lnFact(k);
}

/** -log10 P(X >= k) for X ~ Poisson(lam), k an integer >= 0. Accurate far into the tail. */
export function poisUpperScore(k: number, lam: number): number {
  if (k <= 0) return 0;
  if (lam <= 0) return Infinity;
  if (k > lam) {
    // sum_{j>=k} pmf(j) = pmf(k) * (1 + lam/(k+1) + lam^2/((k+1)(k+2)) + ...)
    let t = 1, s = 1;
    for (let j = k + 1; j < k + 10000; j++) {
      t *= lam / j; s += t;
      if (t < 1e-16 * s) break;
    }
    return -(lnPoisPmf(k, lam) + Math.log(s)) / Math.LN10;
  }
  // k <= lam: P(X >= k) = 1 - P(X <= k-1), not small, so plain subtraction is fine
  let t = 1, s = 1;
  for (let j = k - 1; j > 0; j--) {
    t *= j / lam; s += t;
    if (t < 1e-16 * s) break;
  }
  const lower = Math.exp(lnPoisPmf(k - 1, lam) + Math.log(s));
  return -Math.log10(Math.max(1e-300, 1 - lower));
}

/** Poisson random number (Knuth below 30, normal approximation above; used for totals). */
export function rpois(lam: number, unif: () => number, norm: () => number): number {
  if (lam <= 0) return 0;
  if (lam < 30) {
    const L = Math.exp(-lam);
    let k = 0, p = unif();
    while (p > L) { k++; p *= unif(); }
    return k;
  }
  return Math.max(0, Math.round(lam + Math.sqrt(lam) * norm()));
}

/** Gamma(shape, 1) random number (Marsaglia–Tsang). */
export function rgamma(shape: number, unif: () => number, norm: () => number): number {
  if (shape < 1) return rgamma(shape + 1, unif, norm) * Math.pow(unif(), 1 / shape);
  const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
  for (;;) {
    let x: number, v: number;
    do { x = norm(); v = 1 + c * x; } while (v <= 0);
    v = v * v * v;
    const u = unif();
    if (u < 1 - 0.0331 * x ** 4 || Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
  }
}

/**
 * Benjamini–Hochberg on scores s = -log10 p. Returns -log10 q for each entry.
 * q_(i) = min_{j >= i} p_(j) * n / j, with p sorted ascending.
 */
export function bhScores(score: ArrayLike<number>): Float64Array {
  const n = score.length;
  const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => score[b] - score[a]);
  const out = new Float64Array(n);
  let best = 0; // running max of -log10 q from the least significant end
  for (let r = n - 1; r >= 0; r--) {
    const i = idx[r];
    const qs = Math.max(0, score[i] - Math.log10(n / (r + 1)));
    best = r === n - 1 ? qs : Math.max(best, qs);
    out[i] = best;
  }
  return out;
}
