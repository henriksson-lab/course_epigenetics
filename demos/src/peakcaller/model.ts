// Simulated ChIP-seq and input over a 30-kb stretch of a 2.7-Gb genome, and a MACS-style
// peak caller on the middle 20 kb (no DOM, testable with node).
//
// Reads are fragments of length D = 200 bp (single-end reads extended to d, as MACS does).
// Background fragments fall uniformly with density bg * N / G per bp (bg = local level
// relative to the genome average), plus an optional 2-kb "bump" (open chromatin or a
// copy-number gain) that is in the ChIP and the input alike. True binding sites add
// fragments centred on the site (SD 20 bp) to the ChIP only.
//
// All random draws come from fixed pools at the maximum of every slider, thinned to the
// current setting (keeping a pool read with probability rate / max rate is an exact way to
// draw a Poisson process at the lower rate), so moving a slider adds or removes reads
// instead of redrawing everything.
//
// Caller (MACS2/3): if the input is smaller, the ChIP pileup is scaled by N_input / N_ChIP
// and rounded up. lambda_BG = N_ctrl * d / G; lambda_w = (control 5' ends in a window of
// width w centred on the base) * d / w, for w = d, 1 kb, 10 kb; lambda_local = the max.
// Without an input: max(lambda_BG of the ChIP, lambda_10k of the ChIP). p = P(X >= k) with
// X ~ Poisson(lambda_local), scored every STEP bp; q by Benjamini–Hochberg over those
// positions; peaks = runs above the cutoff, joined across gaps <= GAP, at least MINLEN long.
import { rng } from "../lib/ui";
import { poisUpperScore, rpois, bhScores } from "./stats";

export const G = 2.7e9, D = 200, SIM = 30000, OFF = 5000, SHOW = 20000, STEP = 10;
export const GAP = 50, SITE_SD = 20, MINLEN = 100;
export const SITES = [3000, 8000, 13000, 17000];      // in shown coordinates (0..SHOW)
export const SITE_REL = [1, 0.3, 0.55, 0.15];          // relative strength
export const BUMP = { from: 9500, to: 11500 };          // shown coordinates
export const MAX = { depth: 80, bg: 3, bump: 10, signal: 100 };

export interface Params {
  depth: number;     // ChIP reads, millions
  inputPct: number;  // input depth, % of the ChIP
  signal: number;    // expected fragments at the strongest site per 10 M ChIP reads
  bg: number;        // local background, x genome average
  bump: number;      // background inside the bump, x the local background (1 = no bump)
  control: boolean;  // use the input
  local: boolean;    // use the local windows (else lambda_BG only)
}

export interface Frag { c: number; s: 1 | -1; u: number }   // centre (SIM coords), strand, thinning uniform

export interface Pools { chipBg: Frag[]; inputBg: Frag[]; sites: Frag[][] }

/** Background density at SIM position x relative to the local level (bump = fold inside). */
function bgShape(x: number, bump: number): number {
  const a = BUMP.from + OFF, b = BUMP.to + OFF, ramp = 200;
  let w = 0;
  if (x > a - ramp && x < b + ramp) w = Math.min(1, (x - (a - ramp)) / ramp, ((b + ramp) - x) / ramp);
  return 1 + (bump - 1) * Math.max(0, w);
}

export function makePools(seed: number): Pools {
  const r = rng(seed);
  const maxRate = MAX.depth * 1e6 / G * MAX.bg * MAX.bump;
  const bgPool = (): Frag[] => {
    const n = rpois(maxRate * SIM, r.unif, r.norm), out: Frag[] = [];
    for (let i = 0; i < n; i++) out.push({ c: r.unif() * SIM, s: r.unif() < 0.5 ? 1 : -1, u: r.unif() });
    return out;
  };
  const chipBg = bgPool(), inputBg = bgPool();
  const sites = SITES.map(() => {
    const n = rpois(MAX.signal * MAX.depth / 10, r.unif, r.norm), out: Frag[] = [];
    for (let i = 0; i < n; i++) out.push({ c: SITE_SD * r.norm(), s: r.unif() < 0.5 ? 1 : -1, u: r.unif() });
    return out;
  });
  return { chipBg, inputBg, sites };
}

/** Keep pool fragments for the current parameters. */
export function draw(p: Params, pools: Pools): { chip: Frag[]; input: Frag[] } {
  const maxRate = MAX.depth * MAX.bg * MAX.bump;
  const keepBg = (f: Frag, depth: number) => f.u < depth * p.bg * bgShape(f.c, p.bump) / maxRate;
  const chip = pools.chipBg.filter(f => keepBg(f, p.depth));
  const inDepth = p.depth * p.inputPct / 100;
  const input = pools.inputBg.filter(f => keepBg(f, inDepth));
  const maxSite = MAX.signal * MAX.depth / 10;
  pools.sites.forEach((pool, j) => {
    const want = p.signal * p.depth / 10 * SITE_REL[j];
    for (const f of pool) if (f.u < want / maxSite) chip.push({ c: f.c + SITES[j] + OFF, s: f.s, u: f.u });
  });
  return { chip, input };
}

/** Coverage of extended fragments [c - D/2, c + D/2) at every base of SIM. */
export function pileup(frags: Frag[]): Float64Array {
  const diff = new Float64Array(SIM + 1);
  for (const f of frags) {
    const a = Math.max(0, Math.round(f.c - D / 2)), b = Math.min(SIM, Math.round(f.c + D / 2));
    if (b > a) { diff[a]++; diff[b]--; }
  }
  const cov = new Float64Array(SIM);
  let run = 0;
  for (let i = 0; i < SIM; i++) { run += diff[i]; cov[i] = run; }
  return cov;
}

/** Prefix sums of 5' ends: cum[i] = number of 5' ends at positions < i. */
export function fivePrimeCum(frags: Frag[]): Float64Array {
  const n = new Float64Array(SIM + 1);
  for (const f of frags) {
    const e = Math.floor(f.s === 1 ? f.c - D / 2 : f.c + D / 2 - 1);
    if (e >= 0 && e < SIM) n[e + 1]++;
  }
  for (let i = 1; i <= SIM; i++) n[i] += n[i - 1];
  return n;
}

/** lambda_w at SIM position x from 5' end prefix sums: count in [x - w/2, x + w/2) * d / w. */
export function lambdaWin(cum: Float64Array, x: number, w: number): number {
  const a = Math.max(0, Math.round(x - w / 2)), b = Math.min(SIM, Math.round(x + w / 2));
  return (cum[b] - cum[a]) * D / w;
}

export type Source = "BG" | "d" | "1k" | "10k";

export interface Result {
  pos: Float64Array;       // shown coordinate of each scored position
  rawChip: Float64Array;   // ChIP pileup (unscaled)
  k: Float64Array;         // pileup used in the test (scaled and rounded up if needed)
  inputPile: Float64Array; // input pileup (the input is never the larger library here)
  lam: Float64Array;       // lambda_local
  parts: { BG: number; d: Float64Array; "1k": Float64Array; "10k": Float64Array };
  src: Source[];           // which estimate set lambda_local
  pScore: Float64Array;    // -log10 p
  qScore: Float64Array;    // -log10 q
  scale: number;           // factor applied to the ChIP pileup
  nChip: number; nInput: number;
}

export function call(p: Params, chip: Frag[], input: Frag[]): Result {
  const nPos = SHOW / STEP;
  const Nchip = p.depth * 1e6, Ninput = Nchip * p.inputPct / 100;
  const scale = p.control && Ninput < Nchip ? Ninput / Nchip : 1;
  const cp = pileup(chip), ip = pileup(input);
  const ctrl = p.control ? input : chip;
  const cum = fivePrimeCum(ctrl);
  const lamBG = (p.control ? Math.min(Ninput, Nchip) : Nchip) * D / G;
  const r: Result = {
    pos: new Float64Array(nPos), rawChip: new Float64Array(nPos), k: new Float64Array(nPos),
    inputPile: new Float64Array(nPos), lam: new Float64Array(nPos),
    parts: { BG: lamBG, d: new Float64Array(nPos), "1k": new Float64Array(nPos), "10k": new Float64Array(nPos) },
    src: [], pScore: new Float64Array(nPos), qScore: new Float64Array(nPos), scale,
    nChip: chip.length, nInput: input.length,
  };
  for (let i = 0; i < nPos; i++) {
    const xs = i * STEP + STEP / 2, x = Math.floor(xs + OFF);
    r.pos[i] = xs;
    r.rawChip[i] = cp[x];
    r.k[i] = Math.max(0, Math.ceil(cp[x] * scale - 1e-9));
    r.inputPile[i] = ip[x];
    // windows: with an input all three from the input (already at the smaller depth);
    // without an input only the 10-kb window of the ChIP itself.
    const ld = p.control ? lambdaWin(cum, x, D) : NaN;
    const l1 = p.control ? lambdaWin(cum, x, 1000) : NaN;
    const l10 = lambdaWin(cum, x, 10000);
    r.parts.d[i] = ld; r.parts["1k"][i] = l1; r.parts["10k"][i] = l10;
    let lam = lamBG, src: Source = "BG";
    if (p.local) {
      if (p.control && ld > lam) { lam = ld; src = "d"; }
      if (p.control && l1 > lam) { lam = l1; src = "1k"; }
      if (l10 > lam) { lam = l10; src = "10k"; }
    }
    r.lam[i] = lam; r.src.push(src);
    r.pScore[i] = poisUpperScore(r.k[i], lam);
  }
  r.qScore = bhScores(r.pScore);
  return r;
}

export interface Peak { from: number; to: number; summit: number; best: number; site: number }  // site index or -1

/** Runs of positions with score >= cut, joined across gaps <= GAP, kept if >= MINLEN long. */
export function peaks(r: Result, score: Float64Array, cut: number): Peak[] {
  const runs: [number, number][] = [];
  for (let i = 0; i < score.length; i++) {
    if (score[i] < cut) continue;
    const a = r.pos[i] - STEP / 2, b = r.pos[i] + STEP / 2;
    const last = runs[runs.length - 1];
    if (last && a - last[1] <= GAP) last[1] = b; else runs.push([a, b]);
  }
  const out: Peak[] = [];
  for (const [a, b] of runs) {
    if (b - a < MINLEN) continue;
    let summit = a, best = -1, bestS = 0;
    for (let i = 0; i < score.length; i++) {
      if (r.pos[i] < a || r.pos[i] > b) continue;
      if (r.rawChip[i] > best) { best = r.rawChip[i]; summit = r.pos[i]; }
      bestS = Math.max(bestS, score[i]);
    }
    const site = SITES.findIndex(s => s >= a - D / 2 && s <= b + D / 2);
    out.push({ from: a, to: b, summit, best: bestS, site });
  }
  return out;
}

/** Index of the scored position nearest to shown coordinate x. */
export function indexAt(x: number): number {
  return Math.max(0, Math.min(SHOW / STEP - 1, Math.floor(x / STEP)));
}
