// Loop extrusion on a 1D chain of bins, and the contact map it leaves. No DOM, so node
// tests can import it.
//
// The chain has n bins. A cohesin has two legs (bins l < r) and holds the chain between them
// as a loop. Each time step:
//   1. every cohesin is released with probability 1 / lifetime (WAPL);
//   2. new cohesins are loaded (NIPBL) at random free neighbouring bins c, c + 1;
//   3. every active, unstalled leg steps one bin outward (left leg left, right leg right),
//      unless the next bin holds another leg (legs cannot pass); a leg that steps off the
//      end of the chain releases its cohesin (the chain continues beyond the window).
// A CTCF site whose motif points INTO the loop stops a leg that arrives at it, with
// probability pStop (the convergent rule): a "→" site stops left-moving legs, a "←" site
// right-moving legs; a site in the other orientation is passed freely. A stopped leg stays
// until its cohesin is released. One-sided extrusion: only one leg (chosen at loading) moves.
//
// Contacts: in each snapshot the chain is a line graph with a shortcut between the two legs of
// every cohesin (the loop base holds them together). The contact frequency of bins i and j is
// d^-alpha, with d their shortest distance along this folded chain. Without cohesin d = |i - j|
// and the map is the polymer background P(s) = s^-alpha; a loop base itself has d = 1.

export type Dir = 1 | -1;          // 1 = "→" (blocks legs moving left), -1 = "←" (blocks legs moving right)
export interface Site { pos: number; dir: Dir }

export interface Params {
  n: number;           // number of bins
  sites: Site[];       // CTCF sites
  loadRate: number;    // expected loading attempts per time step (whole chain)
  lifetime: number;    // mean residence time in steps (release prob 1/lifetime per step)
  pStop: number;       // probability that a convergent CTCF site stops an arriving leg
  oneSided: boolean;   // one-sided extrusion: only one leg moves
  alpha: number;       // polymer exponent of the contact law d^-alpha
}

export interface Cohesin {
  l: number; r: number;
  lStall: boolean; rStall: boolean;   // stopped at a CTCF site
  lMove: boolean; rMove: boolean;     // which legs extrude (both, or one if one-sided)
}

export type Rand = () => number;

/** Site orientation per bin: 0 none, 1 "→", -1 "←" (last site wins if two share a bin). */
export function siteArray(n: number, sites: Site[]): Int8Array {
  const a = new Int8Array(n);
  for (const s of sites) if (s.pos >= 0 && s.pos < n) a[s.pos] = s.dir;
  return a;
}

/** One cell: the cohesins on its chain and the leg occupancy of every bin. */
export class Cell {
  coh: Cohesin[] = [];
  occ: Uint8Array;
  constructor(public n: number) { this.occ = new Uint8Array(n); }

  step(p: Params, site: Int8Array, rnd: Rand): void {
    const n = this.n;
    // 1. release
    const pOff = 1 / Math.max(1, p.lifetime);
    this.coh = this.coh.filter(c => {
      if (rnd() < pOff) { this.occ[c.l] = 0; this.occ[c.r] = 0; return false; }
      return true;
    });
    // 2. load: floor(rate) attempts plus one more with the fractional probability
    let k = Math.floor(p.loadRate) + (rnd() < p.loadRate - Math.floor(p.loadRate) ? 1 : 0);
    while (k-- > 0) {
      const c = Math.floor(rnd() * (n - 1));
      const side = rnd();
      if (this.occ[c] || this.occ[c + 1]) continue;
      this.occ[c] = 1; this.occ[c + 1] = 1;
      const leftOnly = side < 0.5;
      this.coh.push({ l: c, r: c + 1, lStall: false, rStall: false,
        lMove: !p.oneSided || leftOnly, rMove: !p.oneSided || !leftOnly });
    }
    // 3. extrude; a leg that runs off the end of the chain releases its cohesin
    let fell = false;
    for (const c of this.coh) {
      if (c.lMove && !c.lStall) {
        const t = c.l - 1;
        if (t < 0) { this.occ[c.l] = 0; c.l = -1; fell = true; continue; }
        if (!this.occ[t]) {
          this.occ[c.l] = 0; this.occ[t] = 1; c.l = t;
          if (site[t] === 1 && rnd() < p.pStop) c.lStall = true;
        }
      }
      if (c.rMove && !c.rStall) {
        const t = c.r + 1;
        if (t >= n) { this.occ[c.r] = 0; c.r = -1; fell = true; continue; }
        if (!this.occ[t]) {
          this.occ[c.r] = 0; this.occ[t] = 1; c.r = t;
          if (site[t] === -1 && rnd() < p.pStop) c.rStall = true;
        }
      }
    }
    if (fell) this.coh = this.coh.filter(c => {
      if (c.l >= 0 && c.r >= 0) return true;
      if (c.l >= 0) this.occ[c.l] = 0;
      if (c.r >= 0) this.occ[c.r] = 0;
      return false;
    });
  }
}

/**
 * Shortest distances along the folded chain from every bin (BFS on the line graph plus one
 * shortcut per cohesin). Writes into `dist` (n*n, row-major); `partner`, `queue` are scratch.
 */
export function chainDistances(n: number, coh: Cohesin[], dist: Int16Array,
  partner: Int16Array, queue: Int16Array): void {
  partner.fill(-1);
  for (const c of coh) { partner[c.l] = c.r; partner[c.r] = c.l; }
  for (let s = 0; s < n; s++) {
    const off = s * n;
    for (let j = 0; j < n; j++) dist[off + j] = -1;
    dist[off + s] = 0;
    let qh = 0, qt = 0;
    queue[qt++] = s;
    while (qh < qt) {
      const u = queue[qh++], du = dist[off + u] + 1;
      if (u > 0 && dist[off + u - 1] < 0) { dist[off + u - 1] = du; queue[qt++] = u - 1; }
      if (u < n - 1 && dist[off + u + 1] < 0) { dist[off + u + 1] = du; queue[qt++] = u + 1; }
      const v = partner[u];
      if (v >= 0 && dist[off + v] < 0) { dist[off + v] = du; queue[qt++] = v; }
    }
  }
}

/** Contact law: frequency for folded-chain distance d >= 1. Precomputed table. */
export function powTable(n: number, alpha: number): Float64Array {
  const t = new Float64Array(n + 1);
  for (let d = 1; d <= n; d++) t[d] = Math.pow(d, -alpha);
  return t;
}

/** Add one snapshot's contacts (upper triangle, i < j) to `map` (n*n). */
export function addSnapshot(n: number, coh: Cohesin[], map: Float64Array, pw: Float64Array,
  scratch: { dist: Int16Array; partner: Int16Array; queue: Int16Array }): void {
  if (coh.length === 0) {
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) map[i * n + j] += pw[j - i];
    return;
  }
  chainDistances(n, coh, scratch.dist, scratch.partner, scratch.queue);
  const d = scratch.dist;
  for (let i = 0; i < n; i++) {
    const off = i * n;
    for (let j = i + 1; j < n; j++) map[off + j] += pw[d[off + j]];
  }
}

export function makeScratch(n: number) {
  return { dist: new Int16Array(n * n), partner: new Int16Array(n), queue: new Int16Array(n) };
}

/** P(s): mean of the upper-triangle map along each diagonal, s = 1..n-1 (index s). */
export function contactDecay(n: number, map: Float64Array): Float64Array {
  const ps = new Float64Array(n);
  for (let s = 1; s < n; s++) {
    let sum = 0;
    for (let i = 0; i + s < n; i++) sum += map[i * n + i + s];
    ps[s] = sum / (n - s);
  }
  return ps;
}

/**
 * A population of cells run together. `run` advances all cells by `steps` steps; with
 * `snapEvery` > 0 it adds a snapshot of every cell every `snapEvery` steps to `map`.
 */
export class Population {
  cells: Cell[];
  map: Float64Array;
  nSnap = 0;
  t = 0;
  private site: Int8Array;
  private pw: Float64Array;
  private scratch;
  // statistics over snapshots
  sumCoh = 0; sumLoop = 0; nLoop = 0; sumStalledLegs = 0; sumLegs = 0;

  constructor(public p: Params, nCells: number, public rnd: Rand) {
    this.cells = Array.from({ length: nCells }, () => new Cell(p.n));
    this.map = new Float64Array(p.n * p.n);
    this.site = siteArray(p.n, p.sites);
    this.pw = powTable(p.n, p.alpha);
    this.scratch = makeScratch(p.n);
  }

  run(steps: number, snapEvery: number): void {
    for (let k = 0; k < steps; k++) {
      for (const c of this.cells) c.step(this.p, this.site, this.rnd);
      this.t++;
      if (snapEvery > 0 && this.t % snapEvery === 0) this.snapshot();
    }
  }

  snapshot(): void {
    for (const c of this.cells) {
      addSnapshot(this.p.n, c.coh, this.map, this.pw, this.scratch);
      this.nSnap++;
      this.sumCoh += c.coh.length;
      for (const h of c.coh) {
        this.sumLoop += h.r - h.l; this.nLoop++;
        this.sumLegs += 2;
        if (h.lStall) this.sumStalledLegs++;
        if (h.rStall) this.sumStalledLegs++;
      }
    }
  }

  /** Mean contact map (upper triangle), per snapshot. */
  mean(): Float64Array {
    const m = new Float64Array(this.map.length);
    if (this.nSnap === 0) return m;
    for (let k = 0; k < m.length; k++) m[k] = this.map[k] / this.nSnap;
    return m;
  }
}


export type PairKind = "convergent" | "divergent" | "tandem";
export interface SitePair { a: Site; b: Site; kind: PairKind; oe: number }

/** Orientation of two sites a (left) and b (right): → ← convergent, ← → divergent, else tandem. */
export function pairKind(a: Site, b: Site): PairKind {
  if (a.dir === 1 && b.dir === -1) return "convergent";
  if (a.dir === -1 && b.dir === 1) return "divergent";
  return "tandem";
}

/**
 * Observed/expected at the pixel of every pair of CTCF sites minD..maxD bins apart, sorted
 * from the strongest. `mean` is the upper-triangle map, `ps` its P(s) (contactDecay).
 */
export function sitePairs(n: number, sites: Site[], mean: Float64Array, ps: Float64Array,
  minD: number, maxD: number): SitePair[] {
  const s = [...sites].sort((x, y) => x.pos - y.pos);
  const out: SitePair[] = [];
  for (let i = 0; i < s.length; i++) for (let j = i + 1; j < s.length; j++) {
    const d = s[j].pos - s[i].pos;
    if (d < minD || d > maxD || ps[d] <= 0) continue;
    out.push({ a: s[i], b: s[j], kind: pairKind(s[i], s[j]), oe: mean[s[i].pos * n + s[j].pos] / ps[d] });
  }
  return out.sort((x, y) => y.oe - x.oe);
}
