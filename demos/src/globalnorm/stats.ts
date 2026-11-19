// Small statistics helpers for globalnorm (pure, no DOM).

/** log Gamma(x), Lanczos approximation (x > 0). */
export function lgamma(x: number): number {
  const g = 7;
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
  x -= 1;
  let a = c[0];
  const t = x + g + 0.5;
  for (let i = 1; i < g + 2; i++) a += c[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

/** Continued fraction for the regularised incomplete beta (Numerical Recipes betacf). */
function betacf(a: number, b: number, x: number): number {
  const FPMIN = 1e-300;
  let c = 1, d = 1 - (a + b) * x / (a + 1);
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 300; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((a + m2 - 1) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 3e-14) break;
  }
  return h;
}

/** Regularised incomplete beta I_x(a, b). */
export function ibeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x);
  if (x < (a + 1) / (a + b + 2)) return Math.exp(lbt) * betacf(a, b, x) / a;
  return 1 - Math.exp(lbt) * betacf(b, a, 1 - x) / b;
}

/** Two-sided p-value of a t statistic with df degrees of freedom. */
export function tTwoSided(t: number, df: number): number {
  return ibeta(df / (df + t * t), df / 2, 0.5);
}

/** Benjamini–Hochberg adjusted p-values. */
export function bh(p: number[]): number[] {
  const m = p.length;
  const idx = p.map((_, i) => i).sort((a, b) => p[a] - p[b]);
  const out = new Array<number>(m);
  let run = 1;
  for (let r = m - 1; r >= 0; r--) {
    const i = idx[r];
    run = Math.min(run, p[i] * m / (r + 1));
    out[i] = run;
  }
  return out;
}

export function median(v: ArrayLike<number>): number {
  const s = Array.from(v).filter(Number.isFinite).sort((a, b) => a - b);
  if (s.length === 0) return NaN;
  const h = s.length >> 1;
  return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2;
}

/** Poisson random number; normal approximation for large means. */
export function rpois(mu: number, unif: () => number, norm: () => number): number {
  if (mu <= 0) return 0;
  if (mu > 60) return Math.max(0, Math.round(mu + Math.sqrt(mu) * norm()));
  const L = Math.exp(-mu);
  let k = 0, p = 1;
  do { k++; p *= unif(); } while (p > L);
  return k - 1;
}
