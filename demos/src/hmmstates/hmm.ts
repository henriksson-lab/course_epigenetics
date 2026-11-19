// A ChromHMM-style hidden Markov model: K hidden chromatin states along 200-bp bins, each
// emitting M binary marks, independent given the state:
//   P(o_t | s) = prod_m e[s][m]^o_tm (1 - e[s][m])^(1 - o_tm).
// Simulate, forward–backward (scaled) posteriors, Viterbi (logs), emissions-only decoding.
// No DOM, so node tests can import it.

export interface Model {
  pi: number[];       // initial probabilities, length K
  A: number[][];      // K x K transitions, A[r][s] = P(s_{t+1} = s | s_t = r)
  E: number[][];      // K x M emission probabilities, E[s][m] = P(mark m present | s)
}

/** Transitions with the same self-transition `stay` for every state, the rest split evenly. */
export function stayMatrix(K: number, stay: number): number[][] {
  return Array.from({ length: K }, (_, r) => Array.from({ length: K }, (_, s) => (r === s ? stay : (1 - stay) / (K - 1))));
}

/** P(o | s) for one bin. */
export function emission(m: Model, s: number, o: ArrayLike<number>): number {
  let p = 1;
  const e = m.E[s];
  for (let k = 0; k < e.length; k++) p *= o[k] ? e[k] : 1 - e[k];
  return p;
}

function pick(row: number[], u: number): number {
  let c = 0;
  for (let k = 0; k < row.length - 1; k++) { c += row[k]; if (u < c) return k; }
  return row.length - 1;
}

/**
 * Simulate T bins. `us` holds T*(1+M) uniforms: per bin one for the state, then one per mark.
 * Re-using the same uniforms with other parameters gives a comparable ("same seed") sequence.
 */
export function simulate(m: Model, us: ArrayLike<number>, T: number): { s: number[]; o: number[][] } {
  const M = m.E[0].length, s: number[] = [], o: number[][] = [];
  for (let t = 0; t < T; t++) {
    const base = t * (1 + M);
    s.push(pick(t === 0 ? m.pi : m.A[s[t - 1]], us[base]));
    o.push(m.E[s[t]].map((e, k) => (us[base + 1 + k] < e ? 1 : 0)));
  }
  return { s, o };
}

/** Scaled forward–backward: log P(o) and the posterior P(s_t = k | o), T x K. */
export function forwardBackward(m: Model, o: number[][]): { logLik: number; post: number[][] } {
  const T = o.length, K = m.pi.length;
  if (T === 0) return { logLik: 0, post: [] };
  const al: number[][] = [], c: number[] = [];
  const em = o.map(ot => Array.from({ length: K }, (_, k) => emission(m, k, ot)));
  for (let t = 0; t < T; t++) {
    const a = Array.from({ length: K }, (_, s) => {
      let prior = 0;
      if (t === 0) prior = m.pi[s];
      else for (let r = 0; r < K; r++) prior += al[t - 1][r] * m.A[r][s];
      return prior * em[t][s];
    });
    const sum = a.reduce((x, y) => x + y, 0);
    c.push(sum);
    al.push(a.map(x => x / sum));
  }
  const be: number[][] = new Array(T);
  be[T - 1] = new Array(K).fill(1);
  for (let t = T - 2; t >= 0; t--) {
    be[t] = Array.from({ length: K }, (_, r) => {
      let v = 0;
      for (let s = 0; s < K; s++) v += m.A[r][s] * em[t + 1][s] * be[t + 1][s];
      return v / c[t + 1];
    });
  }
  const post = al.map((a, t) => {
    const g = a.map((x, k) => x * be[t][k]);
    const z = g.reduce((x, y) => x + y, 0);
    return g.map(x => x / z);
  });
  return { logLik: c.reduce((x, v) => x + Math.log(v), 0), post };
}

/** Viterbi in log space: the single most probable state path. Ties go to the lower state. */
export function viterbi(m: Model, o: number[][]): number[] {
  const T = o.length, K = m.pi.length;
  if (T === 0) return [];
  const L = Math.log;
  let d = Array.from({ length: K }, (_, s) => L(m.pi[s]) + L(emission(m, s, o[0])));
  const back: number[][] = [];
  for (let t = 1; t < T; t++) {
    const bp: number[] = [], nd: number[] = [];
    for (let s = 0; s < K; s++) {
      let best = -Infinity, arg = 0;
      for (let r = 0; r < K; r++) {
        const v = d[r] + L(m.A[r][s]);
        if (v > best + 1e-12) { best = v; arg = r; }
      }
      bp.push(arg);
      nd.push(best + L(emission(m, s, o[t])));
    }
    back.push(bp);
    d = nd;
  }
  const z = new Array<number>(T);
  let arg = 0;
  for (let s = 1; s < K; s++) if (d[s] > d[arg] + 1e-12) arg = s;
  z[T - 1] = arg;
  for (let t = T - 1; t > 0; t--) z[t - 1] = back[t - 1][z[t]];
  return z;
}

/** Emissions only: each bin gets the state that makes its own marks most likely (ties: lower state). */
export function emissionOnly(m: Model, o: number[][]): number[] {
  const K = m.pi.length;
  return o.map(ot => {
    let arg = 0, best = -1;
    for (let s = 0; s < K; s++) { const p = emission(m, s, ot); if (p > best + 1e-15) { best = p; arg = s; } }
    return arg;
  });
}

/** Most probable state per bin from the posterior (what ChromHMM reports). */
export function posteriorArgmax(post: number[][]): number[] {
  return post.map(p => p.indexOf(Math.max(...p)));
}

export function segments(path: number[]): number {
  let n = path.length ? 1 : 0;
  for (let t = 1; t < path.length; t++) if (path[t] !== path[t - 1]) n++;
  return n;
}

export function agreement(a: number[], b: number[]): number {
  const n = Math.min(a.length, b.length);
  let k = 0;
  for (let i = 0; i < n; i++) if (a[i] === b[i]) k++;
  return n ? k / n : NaN;
}
