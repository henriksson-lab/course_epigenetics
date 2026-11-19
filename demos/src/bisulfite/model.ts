// Bisulfite conversion, read types, three-letter alignment and methylation calls for a short
// double-stranded sequence (no DOM; testable with node).
//
// Coordinates are those of the top (+) strand, 0-based. A CpG at i has its top C at i and its
// bottom C at i + 1 (opposite the G). Methylation is symmetric (both Cs of a CpG on one molecule
// share the state); non-CpG Cs are unmethylated (as in most mammalian cells).
//
// Read types (Bismark's names): OT, OB = the converted original top / bottom strand;
// CTOT, CTOB = their PCR complements (reverse complements). Directional libraries give only
// OT and OB. "Reference orientation" = how a genome browser shows a read on the + strand:
// OT and CTOT then carry top-strand information (C vs T at reference Cs), OB and CTOB
// bottom-strand information (G vs A at reference Gs).
import { rng } from "../lib/ui";

export type ReadType = "OT" | "OB" | "CTOT" | "CTOB";
export const READ_TYPES: ReadType[] = ["OT", "OB", "CTOT", "CTOB"];
export type Context = "CG" | "CHG" | "CHH";

const COMP: Record<string, string> = { A: "T", C: "G", G: "C", T: "A", N: "N" };
export const comp = (s: string): string => s.split("").map(c => COMP[c] ?? "N").join("");
export const rev = (s: string): string => s.split("").reverse().join("");
export const revcomp = (s: string): string => rev(comp(s));
export const ct = (s: string): string => s.replace(/C/g, "T");
export const ga = (s: string): string => s.replace(/G/g, "A");

/** Clean user input: upper case, ACGT only. */
export function cleanSeq(s: string): string {
  return s.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
}

/** CpG positions (top C) in a top-strand sequence. */
export function cpgs(seq: string): number[] {
  const out: number[] = [];
  for (let i = 0; i + 1 < seq.length; i++) if (seq[i] === "C" && seq[i + 1] === "G") out.push(i);
  return out;
}

/** Context of a C on the top strand at i ("" if not a C, or at the end without context). */
export function topContext(seq: string, i: number): Context | "" {
  if (seq[i] !== "C") return "";
  if (seq[i + 1] === "G") return "CG";
  if (i + 2 < seq.length) return seq[i + 2] === "G" ? "CHG" : "CHH";
  if (i + 1 < seq.length) return "CHH";   // last-but-one base: assume CHH (unknown in Bismark)
  return "";
}

/** Context of the bottom-strand C opposite a top G at i. */
export function bottomContext(seq: string, i: number): Context | "" {
  if (seq[i] !== "G") return "";
  if (seq[i - 1] === "C") return "CG";
  if (i - 2 >= 0) return seq[i - 2] === "C" ? "CHG" : "CHH";
  if (i - 1 >= 0) return "CHH";
  return "";
}

export interface Molecule {
  top: string;            // top strand of this molecule (with its allele)
  meth: Set<number>;      // top positions i of CpGs methylated on this molecule (both strands)
}

/**
 * Bisulfite-convert one strand of a molecule. Returns the converted strand 5'->3'.
 * Unmethylated Cs become T with probability eff (1 = complete conversion); 5mC stays C.
 */
export function convertStrand(m: Molecule, strand: "top" | "bottom", eff: number, unif?: () => number): string {
  const L = m.top.length;
  const out: string[] = [];
  if (strand === "top") {
    for (let i = 0; i < L; i++) {
      const b = m.top[i];
      if (b !== "C") { out.push(b); continue; }
      const isMeth = m.top[i + 1] === "G" && m.meth.has(i);
      out.push(isMeth || (unif ? unif() >= eff : false) ? "C" : "T");
    }
    return out.join("");
  }
  // bottom strand 5'->3' = reverse complement; position k corresponds to top i = L-1-k
  for (let k = 0; k < L; k++) {
    const i = L - 1 - k;
    const b = COMP[m.top[i]];
    if (b !== "C") { out.push(b); continue; }
    const isMeth = m.top[i - 1] === "C" && m.meth.has(i - 1);
    out.push(isMeth || (unif ? unif() >= eff : false) ? "C" : "T");
  }
  return out.join("");
}

/** The read (5'->3', as in the FASTQ) of each type from one molecule. */
export function readsFromMolecule(m: Molecule, eff = 1, unif?: () => number): Record<ReadType, string> {
  const ot = convertStrand(m, "top", eff, unif), ob = convertStrand(m, "bottom", eff, unif);
  return { OT: ot, OB: ob, CTOT: revcomp(ot), CTOB: revcomp(ob) };
}

/** A read of a given type shown in reference (+ strand) orientation. */
export function refOrientation(type: ReadType, read: string): string {
  return type === "OT" || type === "CTOB" ? read : revcomp(read);
}

export interface Candidate { label: string; readConv: "C→T" | "G→A"; genome: "C→T" | "G→A"; strand: "+" | "−"; type: ReadType; mismatches: number; shown: string }

const mism = (a: string, b: string): number => { let k = 0; for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) k++; return k; };

/**
 * Bismark's three-letter alignment of one read against a reference: the read is converted in
 * silico (C→T; for non-directional libraries also G→A) and aligned to the C→T and the G→A
 * converted genome. Each combination stands for one read type. `shown` is the converted read
 * in reference orientation.
 */
export function align(read: string, ref: string, directional: boolean): Candidate[] {
  const gCT = ct(ref), gGA = ga(ref);
  const rc = ct(read), rg = ga(read);
  const c: Candidate[] = [
    { label: "C→T read, C→T genome (+)", readConv: "C→T", genome: "C→T", strand: "+", type: "OT", mismatches: mism(rc, gCT), shown: rc },
    { label: "C→T read, G→A genome (−)", readConv: "C→T", genome: "G→A", strand: "−", type: "OB", mismatches: mism(revcomp(rc), gGA), shown: revcomp(rc) },
  ];
  if (!directional) {
    c.push({ label: "G→A read, C→T genome (−)", readConv: "G→A", genome: "C→T", strand: "−", type: "CTOT", mismatches: mism(revcomp(rg), gCT), shown: revcomp(rg) });
    c.push({ label: "G→A read, G→A genome (+)", readConv: "G→A", genome: "G→A", strand: "+", type: "CTOB", mismatches: mism(rg, gGA), shown: rg });
  }
  return c;
}

/** Mismatches of an unconverted (ordinary) alignment: best of + and − strand. */
export function plainMismatches(read: string, ref: string): number {
  return Math.min(mism(read, ref), mism(revcomp(read), ref));
}

export interface Call { pos: number; strand: "+" | "-"; context: Context; meth: boolean | null; letter: string }

const XM: Record<Context, string> = { CG: "z", CHG: "x", CHH: "h" };

/**
 * Methylation calls of a read of type t (5'->3' sequence) against the reference: one call
 * per reference C on the strand the read came from. meth = null when the base is neither
 * the converted nor the unconverted one (a mismatch, e.g. a SNP). Letters as in Bismark's
 * XM tag (upper case = methylated).
 */
export function callRead(type: ReadType, read: string, ref: string): Call[] {
  const r = refOrientation(type, read);
  const top = type === "OT" || type === "CTOT";
  const out: Call[] = [];
  for (let i = 0; i < ref.length; i++) {
    const cx = top ? topContext(ref, i) : bottomContext(ref, i);
    if (!cx) continue;
    const b = r[i];
    const [mb, ub] = top ? ["C", "T"] : ["G", "A"];
    const meth = b === mb ? true : b === ub ? false : null;
    const l = XM[cx];
    out.push({ pos: i, strand: top ? "+" : "-", context: cx, meth, letter: meth === null ? "." : meth ? l.toUpperCase() : l });
  }
  return out;
}

// ---------------------------------------------------------------- simulation

export type Genotype = "ref" | "het" | "hom";

export interface SimSettings {
  seq: string;
  levels: Map<number, number>;   // CpG top position -> true methylation level (0, 0.5, 1)
  snpPos: number;                // top C of the CpG carrying a C/T SNP (-1: none)
  genotype: Genotype;
  directional: boolean;
  coverage: number;              // reads
  eff: number;                   // conversion efficiency
  seed: number;
}

export interface SimRead { type: ReadType; read: string; ref: string /* ref orientation */; alt: boolean; calls: Call[] }

/** Draw `coverage` reads, each covering the whole sequence (a simplification). */
export function simulateReads(s: SimSettings): SimRead[] {
  const r = rng(s.seed);
  const types: ReadType[] = s.directional ? ["OT", "OB"] : READ_TYPES;
  const out: SimRead[] = [];
  for (let k = 0; k < s.coverage; k++) {
    const alt = s.snpPos >= 0 && (s.genotype === "hom" || (s.genotype === "het" && r.unif() < 0.5));
    const top = alt ? s.seq.slice(0, s.snpPos) + "T" + s.seq.slice(s.snpPos + 1) : s.seq;
    const meth = new Set<number>();
    for (const [p, lev] of s.levels) if (r.unif() < lev) meth.add(p);
    const type = types[Math.floor(r.unif() * types.length)];
    const reads = readsFromMolecule({ top, meth }, s.eff, r.unif);
    const read = reads[type];
    out.push({ type, read, ref: refOrientation(type, read), alt, calls: callRead(type, read, s.seq) });
  }
  return out;
}

export interface CpgSummary { pos: number; level: number; mTop: number; nTop: number; mBot: number; nBot: number; beta: number }

/** Per-CpG counts, strands merged as in coverage2cytosine --merge_CpG. */
export function summarizeCpgs(seq: string, levels: Map<number, number>, reads: SimRead[]): CpgSummary[] {
  return cpgs(seq).map(p => {
    let mT = 0, nT = 0, mB = 0, nB = 0;
    for (const rd of reads) for (const c of rd.calls) {
      if (c.context !== "CG" || c.meth === null) continue;
      if (c.strand === "+" && c.pos === p) { nT++; if (c.meth) mT++; }
      if (c.strand === "-" && c.pos === p + 1) { nB++; if (c.meth) mB++; }
    }
    const n = nT + nB;
    return { pos: p, level: levels.get(p) ?? 0, mTop: mT, nTop: nT, mBot: mB, nBot: nB, beta: n ? (mT + mB) / n : NaN };
  });
}

/** Non-CpG calls: methylated and total (estimates the non-conversion rate). */
export function nonCpg(reads: SimRead[]): { m: number; n: number } {
  let m = 0, n = 0;
  for (const rd of reads) for (const c of rd.calls) {
    if (c.context === "CG" || c.meth === null) continue;
    n++; if (c.meth) m++;
  }
  return { m, n };
}
