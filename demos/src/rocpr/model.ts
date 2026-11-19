// Pure model code for the rocpr demo (no DOM), so it can be tested with node.
// Binormal model of a classifier's scores: closed windows (negatives) ~ N(0, sd),
// open windows (positives) ~ N(d, sd); a fraction `pi` of all windows is positive.
// Theoretical (infinite test set) rates and curves, plus empirical, tie-aware curves
// for a finite simulated test set, computed as in the unit's R code.

export interface Binormal { d: number; sd: number; pi: number }

/** Standard normal CDF, via the Numerical Recipes erfc (fractional error < 1.2e-7). */
export function pnorm(x: number): number {
  const z = Math.abs(x) / Math.SQRT2;
  const t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
    t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? 1 - r / 2 : r / 2;
}

export function dnorm(x: number, mu = 0, sd = 1): number {
  const z = (x - mu) / sd;
  return Math.exp(-0.5 * z * z) / (sd * Math.sqrt(2 * Math.PI));
}

export interface Rates {
  t: number; tpr: number; fpr: number; precision: number; accuracy: number; f1: number;
}

/** Rates when every window with score >= t is called positive. */
export function rates(m: Binormal, t: number): Rates {
  const tpr = 1 - pnorm((t - m.d) / m.sd);
  const fpr = 1 - pnorm(t / m.sd);
  const tp = m.pi * tpr, fp = (1 - m.pi) * fpr;
  const precision = tp + fp > 0 ? tp / (tp + fp) : precisionLimit(m);
  const accuracy = tp + (1 - m.pi) * (1 - fpr);
  const f1 = precision + tpr > 0 ? 2 * precision * tpr / (precision + tpr) : 0;
  return { t, tpr, fpr, precision, accuracy, f1 };
}

/** Precision as the threshold goes to +infinity (equal variances: positives dominate the far tail). */
function precisionLimit(m: Binormal): number { return m.d > 0 ? 1 : m.d < 0 ? 0 : m.pi; }

/** ROC AUC of the binormal model: P(positive score > negative score). */
export function rocAuc(m: Binormal): number { return pnorm(m.d / (m.sd * Math.SQRT2)); }

/** Average precision = area under the PR curve = integral of precision dTPR. */
export function avgPrecision(m: Binormal, steps = 4000): number {
  const lo = Math.min(0, m.d) - 9 * m.sd, hi = Math.max(0, m.d) + 9 * m.sd;
  const h = (hi - lo) / steps;
  let s = 0;
  for (let i = 0; i < steps; i++) {
    const t = lo + (i + 0.5) * h;                       // midpoint rule
    s += rates(m, t).precision * dnorm(t, m.d, m.sd) * h;
  }
  return s;
}

/** Points along the ROC / PR curves from high to low threshold. */
export function curve(m: Binormal, n = 241): Rates[] {
  const lo = Math.min(0, m.d) - 5 * m.sd, hi = Math.max(0, m.d) + 5 * m.sd;
  const out: Rates[] = [];
  for (let i = 0; i < n; i++) out.push(rates(m, hi - (hi - lo) * i / (n - 1)));
  return out;
}

/** Threshold whose (FPR, TPR) or (recall, precision) point is closest to (x, y). */
export function nearestThreshold(m: Binormal, x: number, y: number, kind: "roc" | "pr"): number {
  let best = 0, bd = Infinity;
  const lo = Math.min(0, m.d) - 5 * m.sd, hi = Math.max(0, m.d) + 5 * m.sd;
  for (let i = 0; i <= 600; i++) {
    const t = lo + (hi - lo) * i / 600;
    const r = rates(m, t);
    const px = kind === "roc" ? r.fpr : r.tpr, py = kind === "roc" ? r.tpr : r.precision;
    const dd = (px - x) ** 2 + (py - y) ** 2;
    if (dd < bd) { bd = dd; best = t; }
  }
  return best;
}

// ---------------- empirical (finite test set) ----------------

export interface EmpPoint { threshold: number; tpr: number; fpr: number; precision: number }

/**
 * Walk down the ranked list; tied scores form one threshold (one point), exactly as
 * `curves()` in the unit's R6. Returns the points after each distinct score.
 */
export function empiricalCurve(score: number[], label: number[]): EmpPoint[] {
  const o = score.map((_, i) => i).sort((a, b) => score[b] - score[a]);
  const P = label.reduce((s, l) => s + l, 0), N = label.length - P;
  const out: EmpPoint[] = [];
  let tp = 0, fp = 0;
  for (let j = 0; j < o.length; j++) {
    if (label[o[j]]) tp++; else fp++;
    if (j === o.length - 1 || score[o[j + 1]] !== score[o[j]])
      out.push({ threshold: score[o[j]], tpr: tp / P, fpr: fp / N, precision: tp / (tp + fp) });
  }
  return out;
}

/** ROC AUC by trapezoids over the empirical curve (ties = diagonal steps). */
export function aucTrapezoid(c: EmpPoint[]): number {
  let s = 0, x0 = 0, y0 = 0;
  for (const p of c) { s += (p.fpr - x0) * (p.tpr + y0) / 2; x0 = p.fpr; y0 = p.tpr; }
  return s;
}

/** ROC AUC by counting (positive, negative) pairs; a tie counts 1/2. */
export function aucPairs(score: number[], label: number[]): number {
  let s = 0, n = 0;
  for (let i = 0; i < score.length; i++) {
    if (!label[i]) continue;
    for (let j = 0; j < score.length; j++) {
      if (label[j]) continue;
      n++;
      s += score[i] > score[j] ? 1 : score[i] === score[j] ? 0.5 : 0;
    }
  }
  return s / n;
}

/** Average precision: sum over thresholds of (increase in recall) x precision, as avg_prec() in R6. */
export function averagePrecision(c: EmpPoint[]): number {
  let s = 0, r0 = 0;
  for (const p of c) { s += (p.tpr - r0) * p.precision; r0 = p.tpr; }
  return s;
}

/** A simulated test set from the binormal model: exactly round(n * pi) positives (at least 1). */
export function simulate(m: Binormal, n: number, norm: () => number): { score: number[]; label: number[] } {
  const nPos = Math.max(1, Math.round(n * m.pi));
  const score: number[] = [], label: number[] = [];
  for (let i = 0; i < n; i++) {
    const pos = i < nPos;
    score.push((pos ? m.d : 0) + m.sd * norm());
    label.push(pos ? 1 : 0);
  }
  return { score, label };
}
