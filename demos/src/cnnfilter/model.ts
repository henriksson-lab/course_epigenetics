// Pure model code for the cnnfilter demo (no DOM), so it can be tested with node.
// One convolutional filter (4 x k weights) scans a one-hot encoded DNA sequence; bias + ReLU,
// global max pooling, then one dense weight + bias and a sigmoid give P(bound).
// Optionally the filter's reverse complement scans too (the other strand), and the pooling
// takes the maximum over both strands. In-silico mutagenesis (ISM): the change in P(bound)
// for every possible single-base substitution.

export const BASES = ["A", "C", "G", "T"] as const;
export type Filter = number[][];          // [4][k], rows A, C, G, T

export interface Net {
  filter: Filter;
  bias: number;        // conv bias, added before ReLU
  w: number;           // dense weight
  b: number;           // dense bias
  both: boolean;       // also scan the reverse complement filter
}

/** Index of a base (A0 C1 G2 T3), or -1. */
export function baseIndex(c: string): number { return "ACGT".indexOf(c.toUpperCase()); }

/** Clean user input to an upper-case ACGT string (other letters dropped). */
export function cleanSeq(s: string): string { return s.toUpperCase().replace(/[^ACGT]/g, ""); }

/** One-hot matrix [4][L]. */
export function oneHot(seq: string): number[][] {
  const m = BASES.map(() => new Array<number>(seq.length).fill(0));
  for (let i = 0; i < seq.length; i++) { const b = baseIndex(seq[i]); if (b >= 0) m[b][i] = 1; }
  return m;
}

export function revComp(seq: string): string {
  const c: Record<string, string> = { A: "T", C: "G", G: "C", T: "A" };
  return seq.split("").reverse().map(x => c[x] ?? "N").join("");
}

/** Reverse-complement filter: reverse the columns, swap rows A<->T and C<->G. */
export function revCompFilter(f: Filter): Filter {
  const k = f[0].length;
  return [3, 2, 1, 0].map(r => Array.from({ length: k }, (_, j) => f[r][k - 1 - j]));
}

/** Convolution (no padding): score at window start i = sum_j W[base(i+j)][j]. Length L - k + 1. */
export function convolve(seq: string, f: Filter): number[] {
  const k = f[0].length, n = seq.length - k + 1;
  const out: number[] = [];
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < k; j++) { const b = baseIndex(seq[i + j]); if (b >= 0) s += f[b][j]; }
    out.push(s);
  }
  return out;
}

export const relu = (x: number): number => (x > 0 ? x : 0);
export const sigmoid = (z: number): number => 1 / (1 + Math.exp(-z));

export interface Forward {
  fwd: number[];          // raw filter scores, + strand filter
  rev: number[] | null;   // raw scores of the reverse-complement filter (null if !both)
  pooled: number;         // global max of ReLU(score + bias) over positions (and strands)
  argmax: number;         // window start of the max (-1 if the sequence is shorter than k)
  strand: "+" | "-";
  z: number;              // dense: w * pooled + b
  p: number;              // sigmoid(z)
}

export function forward(seq: string, net: Net): Forward {
  const fwd = convolve(seq, net.filter);
  const rev = net.both ? convolve(seq, revCompFilter(net.filter)) : null;
  let pooled = 0, argmax = -1, strand: "+" | "-" = "+", best = -Infinity;
  fwd.forEach((s, i) => { if (s > best) { best = s; argmax = i; strand = "+"; } });
  rev?.forEach((s, i) => { if (s > best) { best = s; argmax = i; strand = "-"; } });
  if (argmax >= 0) pooled = relu(best + net.bias);
  const z = net.w * pooled + net.b;
  return { fwd, rev, pooled, argmax, strand, z, p: sigmoid(z) };
}

/** In-silico mutagenesis: delta[b][i] = P(bound | base b at i) - P(bound | reference); 0 at the reference base. */
export function ism(seq: string, net: Net): number[][] {
  const p0 = forward(seq, net).p;
  const out = BASES.map(() => new Array<number>(seq.length).fill(0));
  for (let i = 0; i < seq.length; i++) {
    for (let b = 0; b < 4; b++) {
      if (BASES[b] === seq[i]) continue;
      const mut = seq.slice(0, i) + BASES[b] + seq.slice(i + 1);
      out[b][i] = forward(mut, net).p - p0;
    }
  }
  return out;
}

/** Log-odds (log2) weights from a count matrix, background 0.25, pseudocount `ps` split over the bases. */
export function logOdds(counts: number[][], ps = 1): Filter {
  const k = counts[0].length;
  const tot = Array.from({ length: k }, (_, j) => counts.reduce((s, r) => s + r[j], 0));
  return counts.map(r => r.map((n, j) => Math.log2((n + ps / 4) / (tot[j] + ps) / 0.25)));
}

/** Highest and lowest possible score of a filter (best / worst base in every column). */
export function scoreRange(f: Filter): [number, number] {
  let lo = 0, hi = 0;
  for (let j = 0; j < f[0].length; j++) {
    const col = f.map(r => r[j]);
    lo += Math.min(...col); hi += Math.max(...col);
  }
  return [lo, hi];
}

/** Consensus: the best base in each column. */
export function consensus(f: Filter): string {
  return Array.from({ length: f[0].length }, (_, j) => {
    let b = 0;
    for (let r = 1; r < 4; r++) if (f[r][j] > f[b][j]) b = r;
    return BASES[b];
  }).join("");
}

// JASPAR 2024 MA0139.2 (CTCF), counts, rows A, C, G, T; consensus GCCACCAGGGGGCGC.
export const CTCF_COUNTS: number[][] = [
  [281, 56, 8, 744, 40, 107, 851, 5, 333, 54, 12, 56, 104, 372, 82],
  [49, 800, 903, 13, 528, 433, 11, 0, 3, 12, 0, 8, 733, 13, 482],
  [449, 21, 0, 65, 334, 48, 32, 903, 566, 504, 890, 775, 5, 507, 307],
  [134, 36, 2, 91, 11, 324, 18, 3, 9, 341, 8, 71, 67, 17, 37],
];

// The filter of the exercise's B2 (width 4, consensus CCAG, bias -5).
export const B2_FILTER: Filter = [
  [-1, -1, 2, -1],
  [2, 2, -1, -1],
  [-1, -1, -1, 2],
  [-1, -1, -1, -1],
];

/** A "GC" filter: +1 for C or G, -1 for A or T, in every column. */
export function gcFilter(k: number): Filter {
  return [0, 1, 1, 0].map(gc => new Array<number>(k).fill(gc ? 1 : -1));
}

/** Random filter from a seeded normal generator. */
export function randomFilter(k: number, norm: () => number, sd = 1.5): Filter {
  return BASES.map(() => Array.from({ length: k }, () => Math.round(sd * norm() * 10) / 10));
}

/** A seeded shuffle of the bases of a sequence (same composition). */
export function shuffleSeq(seq: string, unif: () => number): string {
  const a = seq.split("");
  for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(unif() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
  return a.join("");
}
