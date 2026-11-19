// polycombmemory model (pure, no DOM). A stochastic read–write model of a Polycomb domain,
// in the spirit of Dodd et al. 2007 (Cell 129:813), reduced to two states per nucleosome.
//
// A row of N nucleosomes, each marked (H3K27me3) or unmarked. In each small time step dt:
//  - a marked nucleosome loses its mark with rate `turnover` (demethylation, histone exchange);
//  - an unmarked nucleosome gains it with rate
//        spont + feedback × f^h + nucleation,   f = fraction of marked nucleosomes within ±range,
//    the read–write loop: PRC2 bound to a marked nucleosome (EED) marks its neighbours.
//    h = 1 is a linear loop; h = 2 is cooperative (two marked neighbours needed, as for PRC2,
//    which is allosterically activated by H3K27me3). Dodd et al. showed that a two-state loop
//    without cooperativity (or a third, active state) is not bistable;
//  - during the nucleation pulse (e.g. cold for FLC), nucleosomes in the central nucleation
//    region are marked with an extra rate `nucRate`;
//  - at each replication (every `cycle` hours, a random phase per cell) each marked
//    nucleosome is kept with probability 1/2 and is otherwise replaced by a new, unmarked one.
// Each simulated cell follows one lineage (one daughter after every division).
import { rng } from "../lib/ui";

export interface Params {
  n: number;          // nucleosomes
  cells: number;
  feedback: number;   // per hour, at 100% marked neighbours
  range: number;      // neighbours on each side that can recruit (>= n: the whole domain)
  coop: boolean;      // recruitment rate ∝ f² (true) or ∝ f (false)
  turnover: number;   // per hour
  spont: number;      // per hour, spontaneous marking anywhere
  cycle: number;      // hours between replications; 0 = no division
  pulse: number;      // hours of nucleation (from t = 0)
  nucRate: number;    // per hour, in the nucleation region during the pulse
  nucWidth: number;   // nucleosomes in the nucleation region (centred)
  startMarked: boolean;
  tEnd: number;       // hours simulated
  dt: number;
  seed: number;
}

export const DEFAULTS: Params = {
  n: 60, cells: 150, feedback: 1.5, range: 4, coop: true, turnover: 0.1, spont: 0, cycle: 22,
  pulse: 48, nucRate: 0.01, nucWidth: 6, startMarked: false, tEnd: 480, dt: 0.5, seed: 1,
};

export interface Result {
  p: Params;
  times: Float64Array;          // record times (every hour)
  frac: Float32Array[];         // [cell][record] fraction marked
  kymo: Uint8Array[];           // example cells: [k] -> records × n states (row-major by time)
  examples: number[];           // which cells are shown as kymographs
  divisions: number[][];        // division times per example cell
  final: Float32Array;          // fraction marked per cell, mean over the last 24 h
}

export function simulate(p: Params, nExamples = 4): Result {
  const r = rng(p.seed * 4801 + 11);
  const N = p.n, nC = p.cells, dt = p.dt;
  const steps = Math.round(p.tEnd / dt);
  const recEvery = Math.max(1, Math.round(1 / dt));
  const nRec = Math.floor(steps / recEvery) + 1;
  const times = new Float64Array(nRec);
  for (let k = 0; k < nRec; k++) times[k] = k * recEvery * dt;
  const frac = Array.from({ length: nC }, () => new Float32Array(nRec));
  const examples = Array.from({ length: Math.min(nExamples, nC) }, (_, k) => k);
  const kymo = examples.map(() => new Uint8Array(nRec * N));
  const divisions: number[][] = examples.map(() => []);
  const final = new Float32Array(nC);
  const pOff = 1 - Math.exp(-p.turnover * dt);
  const nuc0 = Math.floor((N - p.nucWidth) / 2), nuc1 = nuc0 + p.nucWidth;
  const R = Math.min(p.range, N);
  const s = new Uint8Array(N), nxt = new Uint8Array(N), cum = new Int32Array(N + 1);
  // lookup: probability of marking for an unmarked nucleosome, by marked-neighbour count and window size
  const maxW = 2 * R;
  const pOnTab = (nuc: number) => {
    const tab = new Float64Array((maxW + 1) * (maxW + 1));
    for (let w = 1; w <= maxW; w++) for (let c = 0; c <= w; c++)
      tab[w * (maxW + 1) + c] = 1 - Math.exp(-(p.spont + p.feedback * (p.coop ? (c / w) ** 2 : c / w) + nuc) * dt);
    return tab;
  };
  const tabPlain = pOnTab(0), tabNuc = pOnTab(p.nucRate);

  for (let c = 0; c < nC; c++) {
    s.fill(p.startMarked ? 1 : 0);
    const phase = p.cycle > 0 ? r.unif() * p.cycle : Infinity;
    let nextDiv = phase;
    const ex = examples.indexOf(c);
    const record = (k: number) => {
      let m = 0;
      for (let i = 0; i < N; i++) m += s[i];
      frac[c][k] = m / N;
      if (ex >= 0) kymo[ex].set(s, k * N);
    };
    record(0);
    for (let st = 1; st <= steps; st++) {
      const t = st * dt;
      const pulsing = t <= p.pulse;
      cum[0] = 0;
      for (let i = 0; i < N; i++) cum[i + 1] = cum[i] + s[i];
      for (let i = 0; i < N; i++) {
        if (s[i]) { nxt[i] = r.unif() < pOff ? 0 : 1; continue; }
        const lo = Math.max(0, i - R), hi = Math.min(N - 1, i + R);
        const w = hi - lo;                       // neighbours, self excluded
        const cnt = cum[hi + 1] - cum[lo];       // self is unmarked
        const tab = pulsing && i >= nuc0 && i < nuc1 ? tabNuc : tabPlain;
        nxt[i] = r.unif() < tab[w * (maxW + 1) + cnt] ? 1 : 0;
      }
      s.set(nxt);
      if (t >= nextDiv) {
        for (let i = 0; i < N; i++) if (s[i] && r.unif() < 0.5) s[i] = 0;
        if (ex >= 0) divisions[ex].push(t);
        nextDiv += p.cycle;
      }
      if (st % recEvery === 0) record(st / recEvery);
    }
    // "at the end" = mean over the last 24 h, so a division just before tEnd does not count as a loss
    const k0 = Math.max(0, nRec - 1 - Math.round(24 / (recEvery * dt)));
    let m = 0;
    for (let k = k0 + 1; k < nRec; k++) m += frac[c][k];
    final[c] = m / (nRec - 1 - k0);
  }
  return { p, times, frac, kymo, examples, divisions, final };
}

/** Mean over cells of the fraction marked, per record. */
export function meanFrac(res: Result): Float64Array {
  const n = res.times.length, out = new Float64Array(n);
  for (const f of res.frac) for (let k = 0; k < n; k++) out[k] += f[k];
  for (let k = 0; k < n; k++) out[k] /= res.frac.length;
  return out;
}

/** Fraction of cells with at least `thr` of their nucleosomes marked at the end ("silenced"). */
export function fracSilenced(res: Result, thr = 0.5): number {
  let k = 0;
  for (const v of res.final) if (v >= thr) k++;
  return k / res.final.length;
}
