// Model for the tn5ladder demo (no DOM; testable with node).
//
// A toy chromatin fibre of FIB bp with a nucleosome-depleted region (NDR) in the middle.
// Nucleosome slots (147 bp core) are laid out outwards from the NDR edges, separated by linkers
// whose lengths vary around the mean (SD LINKER_SD bp). Each slot holds a nucleosome with
// probability "occupancy"; an empty slot is free DNA.
// Tn5 insertions are a Poisson process along the DNA with rate
//   rate(x) = TN5_RATE × amount × (1 in free DNA, P_PROT inside a nucleosome core)
//             × optional helical factor (1 + HELIX_AMP cos(2π x / 10.5 + φ))
// (the helical factor is a toy for "Tn5 reaches the backbone on one face of the helix",
// with one random phase φ per fibre).
// The DNA is cut at every insertion, so the fragments are the pieces between consecutive
// insertions (pieces that run off the end of the fibre are not counted). Each fragment end
// carries adapter A or B at random (each insertion joins one adapter to each new end); only
// A–B fragments are amplified and sequenced. Fragments shorter than MIN_LEN are lost in the
// clean-up; fragments longer than MAX_LEN are not shown (they amplify and cluster poorly).
// The 9-bp target duplication and Tn5's sequence bias are ignored.
import { rng } from "../lib/ui";

export const FIB = 5000, CORE = 147, LINKER_SD = 9, P_PROT = 0.1, TN5_RATE = 1 / 70,
  HELIX = 10.5, HELIX_AMP = 0.9, MIN_LEN = 35, MAX_LEN = 1000;

export interface Params { linker: number; ndr: number; occ: number; tn5: number; helix: boolean }
export const DEFAULTS: Params = { linker: 45, ndr: 300, occ: 0.9, tn5: 1, helix: false };

export interface Nuc { start: number; end: number; present: boolean }
export interface Frag { start: number; end: number; a: "A" | "B"; b: "A" | "B"; ok: boolean }
export interface Fibre { nucs: Nuc[]; ndr: [number, number]; ins: number[]; frags: Frag[] }

type U = () => number;

function normal(u: U): number {
  const a = 1 - u(), b = u();
  return Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * b);
}

/** Lay out nucleosome slots outwards from the NDR. */
export function layout(p: Params, u: U): { nucs: Nuc[]; ndr: [number, number] } {
  const c = FIB / 2, ndr: [number, number] = [Math.round(c - p.ndr / 2), Math.round(c + p.ndr / 2)];
  const nucs: Nuc[] = [];
  const link = () => Math.max(3, Math.round(p.linker + LINKER_SD * normal(u)));
  // rightwards: the first nucleosome (+1) sits right at the NDR edge
  for (let x = ndr[1]; x + CORE <= FIB; ) {
    nucs.push({ start: x, end: x + CORE, present: u() < p.occ });
    x += CORE + link();
  }
  for (let x = ndr[0]; x - CORE >= 0; ) {
    nucs.push({ start: x - CORE, end: x, present: u() < p.occ });
    x -= CORE + link();
  }
  nucs.sort((a, b) => a.start - b.start);
  return { nucs, ndr };
}

/** Is position x inside a nucleosome core? (nucs sorted) */
function protectedAt(nucs: Nuc[], x: number): boolean {
  let lo = 0, hi = nucs.length - 1;
  while (lo <= hi) {
    const m = (lo + hi) >> 1, n = nucs[m];
    if (x < n.start) hi = m - 1;
    else if (x >= n.end) lo = m + 1;
    else return n.present;
  }
  return false;
}

/** Tn5 insertion positions by thinning a Poisson process with the maximal rate. */
export function insertions(nucs: Nuc[], p: Params, u: U): number[] {
  const lmax = TN5_RATE * p.tn5 * (p.helix ? 1 + HELIX_AMP : 1);
  const phase = u() * 2 * Math.PI;
  const out: number[] = [];
  for (let x = -Math.log(1 - u()) / lmax; x < FIB; x += -Math.log(1 - u()) / lmax) {
    let r = TN5_RATE * p.tn5 * (protectedAt(nucs, x) ? P_PROT : 1);
    if (p.helix) r *= 1 + HELIX_AMP * Math.cos(2 * Math.PI * x / HELIX + phase);
    if (u() * lmax < r) out.push(Math.floor(x));
  }
  return out;
}

export function fibre(p: Params, u: U): Fibre {
  const { nucs, ndr } = layout(p, u);
  const ins = insertions(nucs, p, u);
  const frags: Frag[] = [];
  // the adapter on the right side of insertion i and on the left side of insertion i
  const right = ins.map(() => (u() < 0.5 ? "A" : "B") as "A" | "B");
  const left = ins.map(() => (u() < 0.5 ? "A" : "B") as "A" | "B");
  for (let i = 0; i + 1 < ins.length; i++) {
    if (ins[i + 1] === ins[i]) continue;
    const a = right[i], b = left[i + 1];
    frags.push({ start: ins[i], end: ins[i + 1], a, b, ok: a !== b });
  }
  return { nucs, ndr, ins, frags };
}

export interface Hist {
  bin: number; counts: number[];      // sequenced fragments (A–B, MIN_LEN..MAX_LEN) per bin
  total: number;                      // all fragments between two insertions
  sameAdapter: number;                // A–A or B–B (not amplified)
  tooShort: number; tooLong: number;  // A–B but < MIN_LEN or > MAX_LEN
  seq: number;                        // sequenced and shown
  nfr: number; mono: number; di: number; // sequenced fragments < 100, 180–250, 315–473 bp
  insPerKb: number;                   // insertions per kb of fibre
  mode: number;                       // the most common length (bp; bin centre)
}

export function histogram(p: Params, nFib: number, seed: number, bin = 2): Hist {
  const u = rng(seed).unif;
  const nb = Math.ceil(MAX_LEN / bin);
  const counts = new Array(nb).fill(0);
  let total = 0, same = 0, short = 0, long = 0, nfr = 0, mono = 0, di = 0, nIns = 0;
  for (let f = 0; f < nFib; f++) {
    const fb = fibre(p, u);
    nIns += fb.ins.length;
    for (const fr of fb.frags) {
      total++;
      if (!fr.ok) { same++; continue; }
      const len = fr.end - fr.start;
      if (len < MIN_LEN) { short++; continue; }
      if (len > MAX_LEN) { long++; continue; }
      counts[Math.min(nb - 1, Math.floor(len / bin))]++;
      if (len < 100) nfr++;
      else if (len >= 180 && len <= 250) mono++;
      else if (len >= 315 && len <= 473) di++;
    }
  }
  let im = 0;
  counts.forEach((c, i) => { if (c > counts[im]) im = i; });
  const seq = counts.reduce((a, b) => a + b, 0);
  return { bin, counts, total, sameAdapter: same, tooShort: short, tooLong: long, seq, nfr, mono, di,
    insPerKb: nIns / (nFib * FIB / 1000), mode: (im + 0.5) * bin };
}

/** Local maxima (within ±40 bp, above 3% of the highest) of a smoothed histogram above `minLen` bp:
 *  where the ladder's peaks are. */
export function peaks(h: Hist, smooth = 15, minLen = MIN_LEN): number[] {
  const k = Math.max(1, Math.round(smooth / h.bin));
  const s = h.counts.map((_, i) => {
    let a = 0, n = 0;
    for (let j = i - k; j <= i + k; j++) if (j >= 0 && j < h.counts.length) { a += h.counts[j]; n++; }
    return a / n;
  });
  const out: number[] = [];
  const w = Math.round(40 / h.bin), floor = 0.03 * Math.max(...s);
  for (let i = 1; i < s.length - 1; i++) {
    const x = (i + 0.5) * h.bin;
    if (x < minLen + 5) continue;
    let isMax = s[i] > floor;
    for (let j = Math.max(0, i - w); j <= Math.min(s.length - 1, i + w); j++) if (s[j] > s[i]) { isMax = false; break; }
    if (isMax && (out.length === 0 || x - out[out.length - 1] > 60)) out.push(Math.round(x));
  }
  return out;
}
