// Simulated ChIP-seq library for the fingerprint demo (no DOM, testable with node).
//
// Genome: NBINS bins of BIN bp. A fraction UNMAPPABLE of the bins can receive no reads
// (gaps, repeats); the others get a mild background bias (Gamma(BG_SHAPE), mean 1: GC,
// open chromatin), the same in the ChIP and its input. Enriched bins: the first nE of a
// fixed random order of the mappable bins, each with a strength ~ Exponential(1).
// A fragment comes from the enriched component with probability f (spread over the
// enriched bins in proportion to their strength), else from the background.
//
// Library: M = complexity * N distinct molecules. PCR copies each molecule with an
// efficiency ~ Gamma(PCR_SHAPE)/PCR_SHAPE, and N reads are drawn from the copies, so a
// molecule is read c ~ NegBin(mean N/M, size PCR_SHAPE) times. We only need the molecules
// read at least once: their number is M * P(c >= 1), and their c follows the zero-truncated
// distribution. Each molecule has a position (bin, offset, strand) drawn at random.
//
// Metrics follow ENCODE: positions = (5' end, strand); NRF = distinct positions / reads;
// PBC1 = positions with exactly one read / distinct positions; PBC2 = M1 / M2.
// FRiP = reads in the enriched bins / all reads (the enriched bins play the role of peaks).
import { rng } from "../lib/ui";
import { lnFact } from "../peakcaller/stats";

export const NBINS = 20000, BIN = 1000, UNMAPPABLE = 0.08, BG_SHAPE = 8, PCR_SHAPE = 2;
export const READS = 200000;

export interface Params {
  f: number;          // fraction of fragments from enriched bins, 0..1
  enrichedPct: number; // % of bins that are enriched
  complexity: number; // distinct molecules in the library / reads sequenced
}

export interface Genome { bg: Float64Array; order: Int32Array; strength: Float64Array }

export function makeGenome(seed: number): Genome {
  const r = rng(seed);
  const bg = new Float64Array(NBINS), strength = new Float64Array(NBINS);
  const mappable: number[] = [];
  for (let b = 0; b < NBINS; b++) {
    if (r.unif() < UNMAPPABLE) continue;
    bg[b] = gamma(BG_SHAPE, r) / BG_SHAPE;
    strength[b] = -Math.log(1 - r.unif());
    mappable.push(b);
  }
  for (let i = mappable.length - 1; i > 0; i--) {   // shuffle: the order in which bins become enriched
    const j = Math.floor(r.unif() * (i + 1));
    [mappable[i], mappable[j]] = [mappable[j], mappable[i]];
  }
  return { bg, order: Int32Array.from(mappable), strength };
}

function gamma(shape: number, r: { unif: () => number; norm: () => number }): number {
  if (shape < 1) return gamma(shape + 1, r) * Math.pow(r.unif(), 1 / shape);
  const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
  for (;;) {
    let x: number, v: number;
    do { x = r.norm(); v = 1 + c * x; } while (v <= 0);
    v = v * v * v;
    const u = r.unif();
    if (u < 1 - 0.0331 * x ** 4 || Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
  }
}

/** Probability of each bin and the enriched mask. */
export function binProbs(p: Params, g: Genome): { prob: Float64Array; enriched: Uint8Array } {
  const nE = Math.round(p.enrichedPct / 100 * NBINS);
  const enriched = new Uint8Array(NBINS);
  let sE = 0, sB = 0;
  for (let i = 0; i < Math.min(nE, g.order.length); i++) { enriched[g.order[i]] = 1; sE += g.strength[g.order[i]]; }
  for (let b = 0; b < NBINS; b++) sB += g.bg[b];
  const prob = new Float64Array(NBINS);
  const f = sE > 0 ? p.f : 0;
  for (let b = 0; b < NBINS; b++) prob[b] = (1 - f) * g.bg[b] / sB + (enriched[b] ? f * g.strength[b] / sE : 0);
  return { prob, enriched };
}

/** NegBin(mean mu, size a) pmf for c = 0..: returns the zero-truncated cumulative table. */
function truncNbCdf(mu: number, a: number): { p0: number; cdf: number[] } {
  const lp = (c: number) => lnGammaInt(c + a) - lnGammaInt(a) - lnFact(c) + a * Math.log(a / (a + mu)) + c * Math.log(mu / (a + mu));
  const p0 = Math.exp(lp(0));
  const cdf: number[] = [];
  let acc = 0;
  for (let c = 1; c < 5000; c++) {
    acc += Math.exp(lp(c)) / (1 - p0);
    cdf.push(acc);
    if (1 - acc < 1e-12 && c > mu) break;
  }
  cdf[cdf.length - 1] = 1;
  return { p0, cdf };
}

/** ln Gamma(x) for x = integer + a (Lanczos). */
function lnGammaInt(x: number): number {
  const g = 7, c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lnGammaInt(1 - x);
  x -= 1;
  let s = c[0];
  for (let i = 1; i < g + 2; i++) s += c[i] / (x + i);
  const t = x + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(s);
}

export interface Library {
  reads: number; distinct: number; m1: number; m2: number;
  nrf: number; pbc1: number; pbc2: number;
  countsAll: Float64Array;    // reads per bin, duplicates kept
  countsDedup: Float64Array;  // distinct positions per bin
  perPos: number[];           // number of positions with 1, 2, 3, 4, 5+ reads
  enriched: Uint8Array;
}

export function simulate(p: Params, g: Genome, seed: number, N = READS): Library {
  const r = rng(seed * 104729 + 7);
  const { prob, enriched } = binProbs(p, g);
  const cum = new Float64Array(NBINS + 1);
  for (let b = 0; b < NBINS; b++) cum[b + 1] = cum[b] + prob[b];
  const M = p.complexity * N, mu = N / M;
  const { p0, cdf } = truncNbCdf(mu, PCR_SHAPE);
  // molecules read at least once: Binomial(M, 1 - p0), normal approximation (M is large)
  const mean = M * (1 - p0), sd = Math.sqrt(M * (1 - p0) * p0);
  const D = Math.max(1, Math.round(mean + sd * r.norm()));
  const keys = new Float64Array(D), copies = new Float64Array(D);
  const countsAll = new Float64Array(NBINS);
  let reads = 0;
  for (let i = 0; i < D; i++) {
    const u = r.unif() * cum[NBINS];
    let lo = 0, hi = NBINS;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (cum[m] <= u) lo = m; else hi = m; }
    const off = Math.floor(r.unif() * BIN), strand = r.unif() < 0.5 ? 0 : 1;
    const v = r.unif();
    let c = 0; while (c < cdf.length - 1 && cdf[c] < v) c++;
    c += 1;
    keys[i] = (lo * BIN + off) * 2 + strand;
    copies[i] = c;
    countsAll[lo] += c;
    reads += c;
  }
  // positions: molecules sharing a 5' end and strand count as one position
  const idx = Array.from({ length: D }, (_, i) => i).sort((a, b) => keys[a] - keys[b]);
  const countsDedup = new Float64Array(NBINS);
  const perPos = [0, 0, 0, 0, 0];
  let distinct = 0, m1 = 0, m2 = 0;
  for (let i = 0; i < D;) {
    let j = i, c = 0;
    while (j < D && keys[idx[j]] === keys[idx[i]]) { c += copies[idx[j]]; j++; }
    distinct++;
    if (c === 1) m1++; else if (c === 2) m2++;
    perPos[Math.min(5, c) - 1]++;
    countsDedup[Math.floor(keys[idx[i]] / 2 / BIN)]++;
    i = j;
  }
  return { reads, distinct, m1, m2, nrf: distinct / reads, pbc1: m1 / distinct, pbc2: m2 ? m1 / m2 : Infinity,
    countsAll, countsDedup, perPos, enriched };
}

/** Fingerprint: bins sorted by count; returns cumulative read fraction at npts evenly spaced ranks. */
export function fingerprint(counts: Float64Array, npts = 400): { x: number[]; y: number[] } {
  const s = Float64Array.from(counts).sort();
  const tot = s.reduce((a, b) => a + b, 0) || 1;
  const x: number[] = [0], y: number[] = [0];
  let acc = 0, next = 1;
  for (let i = 0; i < s.length; i++) {
    acc += s[i];
    if ((i + 1) / s.length >= next / npts - 1e-12) { x.push((i + 1) / s.length); y.push(acc / tot); next++; }
  }
  return { x, y };
}

/** Fraction of reads in the top `frac` of bins. */
export function topShare(counts: Float64Array, frac: number): number {
  const s = Float64Array.from(counts).sort();
  const tot = s.reduce((a, b) => a + b, 0) || 1;
  let acc = 0;
  const k = Math.round(frac * s.length);
  for (let i = s.length - k; i < s.length; i++) acc += s[i];
  return acc / tot;
}

export function frip(counts: Float64Array, enriched: Uint8Array): number {
  let a = 0, t = 0;
  for (let b = 0; b < counts.length; b++) { t += counts[b]; if (enriched[b]) a += counts[b]; }
  return t ? a / t : 0;
}

/** ENCODE classes (ChIP-seq). */
export function nrfClass(v: number): string { return v > 0.9 ? "ideal" : v >= 0.8 ? "compliant" : v >= 0.5 ? "acceptable" : "concerning"; }
export function pbc1Class(v: number): string { return v > 0.9 ? "none" : v >= 0.8 ? "mild" : v >= 0.5 ? "moderate" : "severe"; }
export function pbc2Class(v: number): string { return v > 10 ? "none" : v >= 3 ? "mild" : v >= 1 ? "moderate" : "severe"; }
