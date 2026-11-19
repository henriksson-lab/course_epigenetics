// Strand cross-correlation of simulated single-end ChIP-seq (no DOM, testable with node).
//
// Genome: GLEN bp. Mappability: alternating unique and repeat stretches (exponential
// lengths). Each repeat stretch belongs to a repeat whose copies are similar over COPY bp
// (exponential, mean COPY_MEAN): a read of length L whose L-mer starts in that stretch is
// unmappable if L <= COPY (a longer read reaches unique sequence). So longer reads lose
// fewer positions. A + strand read covering [x, x + L) has its 5' end at x; a − strand read
// covering [y − L, y) has its 5' end at y (BED end coordinate), and its L-mer starts at y − L.
// So + strand mappability at x equals − strand mappability at x + L: wherever reads are
// missing on the + strand, they are also missing on the − strand L bp further on. That
// shared pattern correlates the strands at shift L in any library (the "phantom" peak).
//
// Fragments: background fragments uniform along the genome, modulated by a smooth regional
// factor (log-normal per 10-kb block: GC, copy number, open chromatin), which gives the
// positive background correlation at all shifts (the minimum of the curve). Signal
// fragments are centred on NSITES binding sites (SD 15 bp). Fragment lengths ~ Normal(f, 0.12 f).
// A fragment [a, a + len) gives a + read at a or a − read at a + len (random strand).
//
// CC(k) = Pearson correlation of fwd(x) and rev(x + k) over the genome, computed as
// (sum fwd(x) rev(x + k) / GLEN − mean fwd · mean rev) / (sd fwd · sd rev), as in the exercise's R1.
import { rng } from "../lib/ui";

export const GLEN = 1_000_000, NSITES = 300, SITE_SD = 15, BLOCK = 20000, KMAX = 600, REGION_SD = 1.0;
export const MAP = { unique: 12, repeat: 5 };   // mean stretch lengths (bp)
export const COPY_MEAN = 80;
export const LEN_CV = 0.12;

export interface Params {
  reads: number;     // total reads (before the mappability filter)
  signal: number;    // fraction of fragments from binding sites, 0..1
  frag: number;      // mean fragment length
  readLen: number;   // read length L
}

export interface Genome {
  repeatStarts: Int32Array;  // start of each repeat stretch
  repeatId: Int32Array;      // index of the repeat stretch at each position, -1 if unique
  copyLen: Float64Array;     // per repeat stretch: reads up to this length are unmappable
  sites: Int32Array;
  region: Float64Array;      // regional factor per BLOCK
  regionCum: Float64Array;   // cumulative, for sampling
}

export function makeGenome(seed: number): Genome {
  const r = rng(seed);
  const repeatId = new Int32Array(GLEN).fill(-1);
  const rs: number[] = [], cl: number[] = [];
  let x = 0, uniq = true;
  while (x < GLEN) {
    const len = Math.max(1, Math.round(-Math.log(1 - r.unif()) * (uniq ? MAP.unique : MAP.repeat)));
    if (!uniq) {
      repeatId.fill(rs.length, x, Math.min(GLEN, x + len));
      rs.push(x); cl.push(-Math.log(1 - r.unif()) * COPY_MEAN);
    }
    x += len; uniq = !uniq;
  }
  const sites = new Int32Array(NSITES);
  for (let i = 0; i < NSITES; i++) sites[i] = Math.floor(1000 + r.unif() * (GLEN - 2000));
  sites.sort();
  const nb = Math.ceil(GLEN / BLOCK);
  const region = new Float64Array(nb), regionCum = new Float64Array(nb + 1);
  for (let b = 0; b < nb; b++) { region[b] = Math.exp(REGION_SD * r.norm()); regionCum[b + 1] = regionCum[b] + region[b]; }
  return { repeatStarts: Int32Array.from(rs), repeatId, copyLen: Float64Array.from(cl), sites, region, regionCum };
}

/** Can a read of length L whose L-mer starts at x be mapped uniquely? */
export function mappable(g: Genome, x: number, L: number): boolean {
  const id = g.repeatId[x];
  return id < 0 || g.copyLen[id] < L;
}

export interface Reads { fwd: Uint16Array; rev: Uint16Array; nPlus: number; nMinus: number; dropped: number }

/** Simulate p.reads reads (before the mappability filter). */
export function simulate(p: Params, g: Genome, seed: number): Reads {
  const r = rng(seed * 7919 + 13);
  const fwd = new Uint16Array(GLEN + 1), rev = new Uint16Array(GLEN + 1);
  let nPlus = 0, nMinus = 0, dropped = 0;
  const nb = g.region.length, tot = g.regionCum[nb];
  const n = p.reads;
  for (let i = 0; i < n; i++) {
    // the same number of draws for every read, so a slider changes reads, not the stream
    const uSig = r.unif(), zLen = r.norm(), u1 = r.unif(), u2 = r.unif(), zSite = r.norm(), uStr = r.unif();
    const len = Math.max(p.readLen, Math.round(p.frag * (1 + LEN_CV * zLen)));
    let centre: number;
    if (uSig < p.signal) {
      centre = g.sites[Math.floor(u1 * NSITES)] + SITE_SD * zSite;
    } else {
      // pick a block proportional to its regional factor, then a uniform position in it
      const u = u1 * tot;
      let lo = 0, hi = nb;
      while (hi - lo > 1) { const m = (lo + hi) >> 1; if (g.regionCum[m] <= u) lo = m; else hi = m; }
      centre = lo * BLOCK + u2 * BLOCK;
    }
    const a = Math.round(centre - len / 2);
    const plus = uStr < 0.5;
    if (plus) {
      if (a < 0 || a + p.readLen > GLEN) { dropped++; continue; }
      if (!mappable(g, a, p.readLen)) { dropped++; continue; }
      fwd[a]++; nPlus++;
    } else {
      const e = a + len;                      // 5' end of the − read (BED end)
      if (e - p.readLen < 0 || e > GLEN) { dropped++; continue; }
      if (!mappable(g, e - p.readLen, p.readLen)) { dropped++; continue; }
      rev[e]++; nMinus++;
    }
  }
  return { fwd, rev, nPlus, nMinus, dropped };
}

/** CC(k) for k = 0..kmax. */
export function crossCorr(rd: Reads, kmax = KMAX): Float64Array {
  const { fwd, rev } = rd;
  const N = GLEN;
  let sf = 0, sf2 = 0, sr = 0, sr2 = 0;
  let np = 0;
  for (let x = 0; x < N; x++) if (fwd[x]) np++;
  const pos = new Int32Array(np), val = new Float64Array(np);
  np = 0;
  for (let x = 0; x < N; x++) {
    const a = fwd[x], b = rev[x];
    if (a) { sf += a; sf2 += a * a; pos[np] = x; val[np++] = a; }
    if (b) { sr += b; sr2 += b * b; }
  }
  const mf = sf / N, mr = sr / N;
  const sdf = Math.sqrt(sf2 / N - mf * mf), sdr = Math.sqrt(sr2 / N - mr * mr);
  const cc = new Float64Array(kmax + 1);
  for (let k = 0; k <= kmax; k++) {
    let s = 0;
    for (let i = 0; i < np; i++) { const y = pos[i] + k; if (y <= N) s += val[i] * rev[y]; }
    cc[k] = (s / N - mf * mr) / (sdf * sdr);
  }
  return cc;
}

export interface QC { frag: number; ccFrag: number; ccRead: number; ccMin: number; nsc: number; rsc: number }

/**
 * Fragment length = the shift with the highest CC outside 10 .. L + 10 (phantompeakqualtools'
 * default exclusion, -x=10:(readlen+10)). CC(read length) = CC(L).
 * min = the smallest CC over all shifts computed (0 .. KMAX).
 */
export function qc(cc: Float64Array, readLen: number): QC {
  let frag = -1, best = -Infinity, mn = Infinity;
  for (let k = 0; k < cc.length; k++) {
    mn = Math.min(mn, cc[k]);
    if (k < 10 || k <= readLen + 10) continue;
    if (cc[k] > best) { best = cc[k]; frag = k; }
  }
  const ccRead = cc[readLen];
  return { frag, ccFrag: best, ccRead, ccMin: mn, nsc: best / mn, rsc: (best - mn) / (ccRead - mn) };
}

/** Aggregate 5' end counts around anchors, in bins of `bin` bp over [-half, half). */
export function aggregate(counts: Uint16Array, anchors: ArrayLike<number>, half: number, bin: number): Float64Array {
  const nb = Math.ceil(2 * half / bin);
  const out = new Float64Array(nb);
  for (let i = 0; i < anchors.length; i++) {
    const a = anchors[i];
    for (let d = -half; d < half; d++) {
      const x = a + d;
      if (x < 0 || x > GLEN) continue;
      const c = counts[x];
      if (c) out[Math.floor((d + half) / bin)] += c;
    }
  }
  return out;
}
