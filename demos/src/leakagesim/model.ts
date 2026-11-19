// Pure simulation for the leakagesim demo (no DOM), so it can be tested with node.
//
// Simulated genome: nChrom "chromosomes", each with the same number of regions. Every region
// has a label (open/closed, balanced) and a random DNA sequence whose GC content is
//   GC = 0.42 + (gc spread between regions) * N(0,1) + signal/2 * (+1 open, -1 closed),
// so `signal` is the true difference in mean GC between open and closed regions (0 = labels
// carry no information about sequence at all). A fraction `dup` of regions are homologous
// copies (2% of bases mutated) of a region on another chromosome, with the same label
// (a segmental duplication or repeat family). Each region is cut into nWin windows of winLen bp
// that start `step` = winLen * (1 - overlap) apart (overlap 0 = adjacent windows).
//
// Classifier: k-nearest neighbours on 4-mer counts (Euclidean distance), as in the unit's R1.
// Evaluated with three test sets: 20% of windows at random, 20% of regions (all windows of a
// region together), and whole chromosomes (2 of 10). Baseline: window GC as the score (no training).

export interface Params {
  overlap: number;   // 0..0.9, overlap between neighbouring windows of a region
  dup: number;       // 0..0.5, fraction of regions that are homologous copies
  signal: number;    // 0..0.12, true GC difference between open and closed regions
  seed: number;
}

export const NCHROM = 10, PER_CHROM = 30, NWIN = 3, WINLEN = 200, K = 4;
export const GC_SD = 0.035;           // spread of GC between regions of the same class
export const MUT = 0.02;              // divergence of homologous copies
export const TEST_CHROMS = [8, 9];    // "chr9" and "chr10"
export const KNN_K = 15;

export interface Region { chrom: number; label: number; gc: number; dupOf: number; seq: string }
export interface Win { region: number; chrom: number; label: number; offset: number; gc: number; kmers: Uint16Array; counts: Uint8Array; norm2: number }
export interface Data { regions: Region[]; wins: Win[]; step: number }

/** mulberry32: a small seeded generator (same as lib/ui rng, kept here so node tests need no DOM). */
export function mulberry(seed: number): { unif: () => number; norm: () => number } {
  let a = seed >>> 0;
  const unif = (): number => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  const norm = (): number => Math.sqrt(-2 * Math.log(1 - unif())) * Math.cos(2 * Math.PI * unif());
  return { unif, norm };
}

function randomSeq(n: number, gc: number, unif: () => number): string {
  let s = "";
  for (let i = 0; i < n; i++) {
    const u = unif();
    s += u < gc / 2 ? "G" : u < gc ? "C" : u < gc + (1 - gc) / 2 ? "A" : "T";
  }
  return s;
}

/**
 * Build the data set. Every random choice comes from a stream tied to the seed and the region,
 * so moving one slider changes only what that slider controls.
 */
export function makeData(p: Params): Data {
  const n = NCHROM * PER_CHROM;
  const step = Math.max(1, Math.round(WINLEN * (1 - p.overlap)));
  const len = WINLEN + (NWIN - 1) * step;
  const maxLen = WINLEN + (NWIN - 1) * WINLEN;          // sequence drawn at full length, then cut
  const lab = mulberry(p.seed * 1009 + 1);
  const labels = Array.from({ length: n }, (_, i) => (i % 2));
  for (let i = n - 1; i > 0; i--) { const j = Math.floor(lab.unif() * (i + 1)); [labels[i], labels[j]] = [labels[j], labels[i]]; }
  const regions: Region[] = [];
  for (let r = 0; r < n; r++) {
    const g = mulberry(p.seed * 7919 + r * 104729 + 3);
    const z = g.norm(), uDup = g.unif(), uSrc = g.unif();
    const chrom = Math.floor(r / PER_CHROM);
    const isDup = uDup < p.dup;
    let src = -1;
    if (isDup) {
      // a source region on another chromosome that is not itself a copy
      const cand = regions.filter(q => q.chrom !== chrom && q.dupOf < 0);
      if (cand.length) src = regions.indexOf(cand[Math.floor(uSrc * cand.length)]);
    }
    if (src >= 0) {
      const s = regions[src];
      const m = mulberry(p.seed * 31337 + r * 7 + 11);
      let seq = "";
      for (const c of s.seq) seq += m.unif() < MUT ? "ACGT"[Math.floor(m.unif() * 4)] : c;
      regions.push({ chrom, label: s.label, gc: s.gc, dupOf: src, seq });
    } else {
      const y = labels[r];
      const gc = Math.min(0.9, Math.max(0.1, 0.42 + GC_SD * z + (p.signal / 2) * (y ? 1 : -1)));
      const sq = mulberry(p.seed * 4241 + r * 613 + 5);
      regions.push({ chrom, label: y, gc, dupOf: -1, seq: randomSeq(maxLen, gc, sq.unif) });
    }
  }
  const wins: Win[] = [];
  regions.forEach((rg, r) => {
    const seq = rg.seq.slice(0, len);
    for (let w = 0; w < NWIN; w++) {
      const s = seq.slice(w * step, w * step + WINLEN);
      wins.push({ region: r, chrom: rg.chrom, label: rg.label, offset: w * step, ...features(s) });
    }
  });
  return { regions, wins, step };
}

/** 4-mer counts as a sparse vector (sorted indices + counts), its squared norm, and GC fraction. */
export function features(s: string): { gc: number; kmers: Uint16Array; counts: Uint8Array; norm2: number } {
  const code: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };
  const m = new Map<number, number>();
  let gc = 0;
  for (let i = 0; i < s.length; i++) if (s[i] === "G" || s[i] === "C") gc++;
  for (let i = 0; i + K <= s.length; i++) {
    let idx = 0;
    for (let j = 0; j < K; j++) idx = idx * 4 + code[s[i + j]];
    m.set(idx, (m.get(idx) ?? 0) + 1);
  }
  const keys = [...m.keys()].sort((a, b) => a - b);
  const counts = Uint8Array.from(keys.map(k => m.get(k)!));
  let norm2 = 0;
  for (const c of counts) norm2 += c * c;
  return { gc: gc / s.length, kmers: Uint16Array.from(keys), counts, norm2 };
}

export type SplitName = "random" | "region" | "chromosome";
export const SPLITS: SplitName[] = ["random", "region", "chromosome"];

/** Test-set membership of every window for one split. */
export function splitTest(d: Data, split: SplitName, seed: number): boolean[] {
  const r = mulberry(seed * 2027 + (split === "random" ? 1 : 2));
  if (split === "random") return d.wins.map(() => r.unif() < 0.2);
  if (split === "region") {
    const testReg = d.regions.map(() => r.unif() < 0.2);
    return d.wins.map(w => testReg[w.region]);
  }
  return d.wins.map(w => TEST_CHROMS.includes(w.chrom));
}

export interface Result {
  nTest: number; nTrain: number;
  acc1: number;          // 1-NN accuracy
  aucKnn: number;        // ROC AUC of the 15-NN vote (fraction of open neighbours)
  aucGc: number;         // ROC AUC of window GC as the score (no training)
  nnSame: number;        // fraction of test windows whose nearest training window is from the same region
  nnHomolog: number;     // ... is from a homologous copy of its region (or the region it was copied from)
}

/** k-NN evaluation for one split. */
export function evaluate(d: Data, test: boolean[]): Result {
  const tr: number[] = [], te: number[] = [];
  test.forEach((t, i) => (t ? te : tr).push(i));
  const dense = new Float32Array(256);
  const nnLab: number[] = [], vote: number[] = [], lab: number[] = [], gcs: number[] = [];
  let same = 0, homo = 0, right = 0;
  const dist = new Float64Array(tr.length);
  for (const i of te) {
    const a = d.wins[i];
    dense.fill(0);
    for (let j = 0; j < a.kmers.length; j++) dense[a.kmers[j]] = a.counts[j];
    for (let t = 0; t < tr.length; t++) {
      const b = d.wins[tr[t]];
      let dot = 0;
      for (let j = 0; j < b.kmers.length; j++) dot += dense[b.kmers[j]] * b.counts[j];
      dist[t] = a.norm2 + b.norm2 - 2 * dot;
    }
    // the KNN_K smallest distances (ties broken by index, for reproducibility): insertion into a short list
    const best: number[] = [];
    for (let t = 0; t < tr.length; t++) {
      if (best.length === KNN_K && dist[t] >= dist[best[KNN_K - 1]]) continue;
      let j = best.length < KNN_K ? best.length : KNN_K - 1;
      while (j > 0 && dist[best[j - 1]] > dist[t]) { if (j < KNN_K) best[j] = best[j - 1]; j--; }
      best[j] = t;
    }
    const top = best.map(t => d.wins[tr[t]]);
    const nn = top[0];
    if (nn.label === a.label) right++;
    if (nn.region === a.region) same++;
    else {
      const ra = d.regions[a.region], rb = d.regions[nn.region];
      const root = (r: number, q: Region) => (q.dupOf >= 0 ? q.dupOf : r);
      if (root(a.region, ra) === root(nn.region, rb)) homo++;
    }
    nnLab.push(nn.label);
    vote.push(top.reduce((s, w) => s + w.label, 0) / top.length);
    lab.push(a.label);
    gcs.push(a.gc);
  }
  const n = te.length;
  return {
    nTest: n, nTrain: tr.length, acc1: right / n, aucKnn: aucPairs(vote, lab), aucGc: aucPairs(gcs, lab),
    nnSame: same / n, nnHomolog: homo / n,
  };
}

/** ROC AUC by ranks (Mann-Whitney): P(score of a positive > score of a negative), ties count 1/2. */
export function aucPairs(score: number[], label: number[]): number {
  const idx = score.map((_, i) => i).sort((a, b) => score[a] - score[b]);
  const rank = new Float64Array(score.length);
  for (let i = 0; i < idx.length;) {
    let j = i;
    while (j + 1 < idx.length && score[idx[j + 1]] === score[idx[i]]) j++;
    const r = (i + j) / 2 + 1;                       // average rank of the tie group
    for (let t = i; t <= j; t++) rank[idx[t]] = r;
    i = j + 1;
  }
  let n1 = 0, s1 = 0;
  label.forEach((l, i) => { if (l) { n1++; s1 += rank[i]; } });
  const n0 = label.length - n1;
  if (!n1 || !n0) return NaN;
  return (s1 - n1 * (n1 + 1) / 2) / (n1 * n0);
}

export function runAll(p: Params): { data: Data; tests: Record<SplitName, boolean[]>; res: Record<SplitName, Result> } {
  const data = makeData(p);
  const tests = {} as Record<SplitName, boolean[]>, res = {} as Record<SplitName, Result>;
  for (const s of SPLITS) { tests[s] = splitTest(data, s, p.seed); res[s] = evaluate(data, tests[s]); }
  return { data, tests, res };
}
