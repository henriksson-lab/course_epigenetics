// Model for the shufflenull demo (no DOM; testable with node).
//
// A toy chromosome of L kb with regions where no peak can be (an unassembled short arm and
// smaller runs of N, the centromere, blacklisted repeats) and two GC-rich, gene-dense regions.
// Two simulated interval sets:
//   A: 150 intervals (5–15 kb), 75% inside the GC-rich regions, none in gaps or blacklist.
//   B: 100 intervals (10–30 kb), the same GC preference (70% GC-rich), never in gaps or blacklist;
//      in the scenario "real association", 30% of B are in addition placed on top of an A interval.
// Statistic (as `bedtools intersect -u -a A -b B | wc -l`): the number of A intervals that share
// at least one base with some B interval. Shuffle B (same widths) under a chosen null:
//   naive:          anywhere on the chromosome (bedtools shuffle -g genome),
//   keep out:       not overlapping gaps, centromere or blacklist (-excl),
//   match GC:       each interval stays in its own class, GC-rich or not (by its midpoint),
// and compare: enrichment = observed / mean(null); empirical p = (k + 1) / (N + 1), k = number of
// shuffles with an overlap count >= the observed one.
import { rng } from "../lib/ui";

export const L = 10_000;                       // chromosome length, kb
export type Iv = [number, number];             // [start, end) in kb

export const GAPS: Iv[] = [[0, 1000], [5200, 5300], [8600, 8700], [9950, 10000]];
export const CENTRO: Iv = [1000, 1400];
export const BLACK: Iv[] = [[2300, 2330], [3900, 3940], [4700, 4725], [6000, 6040], [7300, 7330], [9100, 9140]];
export const GC: Iv[] = [[1800, 3200], [6500, 7600]];
export const EXCL: Iv[] = [...GAPS, CENTRO, ...BLACK].sort((a, b) => a[0] - b[0]);

export const NA = 150, NB = 100;

export type Scenario = "assoc" | "none";
export interface NullOpts { excl: boolean; gc: boolean }

const overlaps = (a: Iv, b: Iv) => a[0] < b[1] && b[0] < a[1];
const inAny = (x: Iv, set: Iv[]) => set.some(s => overlaps(x, s));
const isGC = (x: Iv) => GC.some(g => (x[0] + x[1]) / 2 >= g[0] && (x[0] + x[1]) / 2 < g[1]);

type U = () => number;

/** Total length of [0, L) minus a set of intervals, per region. */
function complement(set: Iv[], lo = 0, hi = L): Iv[] {
  const out: Iv[] = [];
  let x = lo;
  for (const s of [...set].sort((a, b) => a[0] - b[0])) {
    if (s[1] <= lo || s[0] >= hi) continue;
    if (s[0] > x) out.push([x, Math.min(s[0], hi)]);
    x = Math.max(x, s[1]);
  }
  if (x < hi) out.push([x, hi]);
  return out;
}

export const NONGC: Iv[] = complement(GC);

/** Place an interval of width w uniformly with its start in one of the allowed segments
 *  (the whole interval inside the segment), weighted by the number of possible starts. */
function placeIn(segs: Iv[], w: number, u: U): Iv | null {
  const ok = segs.filter(s => s[1] - s[0] >= w);
  const tot = ok.reduce((a, s) => a + (s[1] - s[0] - w), 0);
  if (tot <= 0) return null;
  let r = u() * tot;
  for (const s of ok) {
    const n = s[1] - s[0] - w;
    if (r < n || s === ok[ok.length - 1]) { const st = s[0] + Math.min(r, n); return [st, st + w]; }
    r -= n;
  }
  return null;
}

/** Segments of `region` that are allowed when gaps/blacklist must be avoided. */
function allowed(region: Iv[], excl: boolean): Iv[] {
  if (!excl) return region;
  const out: Iv[] = [];
  for (const r of region) out.push(...complement(EXCL, r[0], r[1]));
  return out;
}

export interface Data { A: Iv[]; B: Iv[]; bGC: boolean[] }

/** The observed data of a scenario. Fixed seed: the same "experiment" every time. */
export function makeData(sc: Scenario, seed = 3): Data {
  const r = rng(seed), u = r.unif;
  const gcOK = allowed(GC, true), restOK = allowed(NONGC, true);
  const A: Iv[] = [];
  for (let i = 0; i < NA; i++) {
    const w = Math.round(5 + 10 * u());
    A.push(placeIn(u() < 0.75 ? gcOK : restOK, w, u)!);
  }
  A.sort((a, b) => a[0] - b[0]);
  const B: Iv[] = [];
  for (let i = 0; i < NB; i++) {
    const w = Math.round(10 + 20 * u());
    const v = u(), z = u(), y = u();
    if (sc === "assoc" && v < 0.3) {
      // on top of a random A interval: midpoint within the A interval
      const a = A[Math.floor(z * A.length)];
      const mid = a[0] + y * (a[1] - a[0]);
      B.push([Math.round(mid - w / 2), Math.round(mid - w / 2) + w]);
    } else {
      B.push(placeIn(z < 0.7 ? gcOK : restOK, w, u)!);
    }
  }
  B.sort((a, b) => a[0] - b[0]);
  return { A, B, bGC: B.map(isGC) };
}

/** Number of A intervals that overlap at least one B interval. */
export function overlapCount(A: Iv[], B: Iv[]): number {
  const b = [...B].sort((x, y) => x[0] - y[0]);
  let n = 0;
  for (const a of A) {
    for (const x of b) {
      if (x[0] >= a[1]) break;
      if (x[1] > a[0]) { n++; break; }
    }
  }
  return n;
}

/** One shuffle of B (same widths) under the null. */
export function shuffleB(d: Data, opt: NullOpts, u: U): Iv[] {
  const whole = allowed([[0, L]], opt.excl);
  const gcSeg = allowed(GC, opt.excl), restSeg = allowed(NONGC, opt.excl);
  return d.B.map((b, i) => {
    const w = b[1] - b[0];
    const segs = opt.gc ? (d.bGC[i] ? gcSeg : restSeg) : whole;
    return placeIn(segs, w, u)!;
  });
}

export interface Stats { obs: number; n: number; mean: number; k: number; p: number; enrich: number }

export function stats(obs: number, nulls: number[]): Stats {
  const n = nulls.length;
  const mean = n ? nulls.reduce((a, b) => a + b, 0) / n : NaN;
  const k = nulls.filter(x => x >= obs).length;
  return { obs, n, mean, k, p: n ? (k + 1) / (n + 1) : NaN, enrich: n ? obs / mean : NaN };
}

/** A stream of shuffles for one (scenario, null) setting, reproducible from the seed. */
export function shuffler(d: Data, opt: NullOpts, seed: number) {
  const u = rng(seed).unif;
  return () => { const b = shuffleB(d, opt, u); return { b, count: overlapCount(d.A, b) }; };
}

export const fracExcluded = (): number => EXCL.reduce((a, s) => a + s[1] - s[0], 0) / L;
export const fracGC = (): number => GC.reduce((a, s) => a + s[1] - s[0], 0) / L;
export { inAny, isGC };
