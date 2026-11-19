# Real-data figures for hic-processing/lec.qmd ("Hi-C: from reads to contact matrices").
# Run from this folder:  python3 make-hicproc-figs.py
#
# Needs (data/download.sh hic-processing; data/prepare-hic.sh):
#   data/gm12878_hic_R1.fastq.gz, data/gm12878_hic_R2.fastq.gz (read pairs of SRR1658570 whose
#   sides map to chr21), data/gm12878_hic_fullgenome_stats.txt (pairtools stats of the first
#   30 M read pairs aligned to the whole genome), data/hg38_chr21.fa.gz,
#   data/gm12878_hic_chr19_chr21.mcool (ENCODE's map of the same experiment, all replicates).
# On the PATH: bwa, samtools, pairtools 1.1, cooler 0.10. Python: cooltools 0.7, numpy, pandas,
# matplotlib.
#
# The processing (the same commands as in exercise.qmd) runs once into ./hic-work
# (set HIC_WORK to put it elsewhere; ~1 GB, ~10 min on 8 threads). Re-running only redoes
# missing steps.
#
# Writes fig04-fragments.png, fig08-pair-types.png, fig10-orientation.png, fig12-read-funnel.png,
# fig14-resolution-depth.png, fig18-balance-real.png, fig19-pofs.png, fig20-obs-exp.png
# and prints the numbers quoted in the slides.
import os, re, gzip, subprocess, time, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
import cooler, cooltools

warnings.filterwarnings("ignore")
D = os.path.abspath("../data")
W = os.path.abspath(os.environ.get("HIC_WORK", "hic-work"))
T = os.environ.get("THREADS", "8")
os.makedirs(W, exist_ok=True)
plt.rcParams.update({"font.size": 13, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 14})
HIC = LinearSegmentedColormap.from_list("hic", ["#ffffff", "#f4b6a6", "#c0392b", "#5a1209"])
DIV = LinearSegmentedColormap.from_list("div", ["#2a6f97", "#ffffff", "#c0392b"])
GREEN, CRIMSON, INK, GREY, BLUE, ORANGE, RED = "#2e8b57", "#a8324a", "#1b2a3a", "#5b6b7a", "#2a6f97", "#e67e22", "#c0392b"
ENC = f"{D}/gm12878_hic_chr19_chr21.mcool::resolutions/"


def step(out, cmd):
    """run a shell command in the work folder unless its output exists."""
    if os.path.exists(os.path.join(W, out)):
        return
    print("running:", out, flush=True)
    t0 = time.time()
    subprocess.run(f"set -euo pipefail; cd {W}; D={D}; T={T}; {cmd}", shell=True, check=True, executable="/bin/bash")
    print(f"  {(time.time() - t0) / 60:.1f} min", flush=True)


# ---- 1. the pipeline of exercise.qmd ------------------------------------------------
step("chr21.fa.bwt", "zcat $D/hg38_chr21.fa.gz > chr21.fa && bwa index chr21.fa 2> bwa_index.log && "
     "samtools faidx chr21.fa && cut -f1,2 chr21.fa.fai > chr21.sizes")
step("chr21.pairs.gz",
     "bwa mem -SP5M -t $T chr21.fa $D/gm12878_hic_R1.fastq.gz $D/gm12878_hic_R2.fastq.gz 2> bwa.log "
     "| pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --drop-sam --drop-seq "
     "-c chr21.sizes --output-stats parse.stats -o chr21.pairs.gz")
step("chr21.nodups.pairs.gz",
     "pairtools sort --nproc 4 -o chr21.sorted.pairs.gz chr21.pairs.gz && "
     "pairtools dedup --max-mismatch 3 --output-stats dedup.stats --output-dups chr21.dups.pairs.gz "
     "-o chr21.nodups.pairs.gz chr21.sorted.pairs.gz")
step("chr21.valid.pairs.gz",
     "pairtools select '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' "
     "-o chr21.valid.pairs.gz chr21.nodups.pairs.gz")
step("chr21.mcool",
     "cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chr21.sizes:10000 chr21.valid.pairs.gz chr21_10kb.cool && "
     "cooler zoomify -r 10000,50000,100000,250000,500000 --balance --balance-args '--cis-only' "
     "-o chr21.mcool chr21_10kb.cool")


def read_stats(fn):
    s = pd.read_csv(fn, sep="\t", header=None, names=["key", "value"], comment=None)
    return s.set_index("key")["value"]


full = read_stats(f"{D}/gm12878_hic_fullgenome_stats.txt")
dd = read_stats(f"{W}/dedup.stats")
ps = read_stats(f"{W}/parse.stats")
f = lambda k, s=full: float(s[k])
print("full genome (first 30 M read pairs of SRR1658570):")
for k in ["total", "total_unmapped", "total_single_sided_mapped", "total_mapped", "total_dups", "total_nodups",
          "cis", "trans", "cis_1kb+", "cis_20kb+"]:
    print(f"  {k:28s} {int(f(k)):>12,}  ({f(k) / f('total') * 100:5.1f}% of pairs)")
print(f"  duplicates: {f('total_dups') / f('total_mapped') * 100:.1f}% of mapped pairs; "
      f"cis {f('cis') / f('total_nodups') * 100:.1f}% of unique pairs; cis > 20 kb {f('cis_20kb+') / f('total_nodups') * 100:.1f}%;"
      f" trans {f('trans') / f('total_nodups') * 100:.1f}%")
print(f"chr21 subset: {int(ps['total']):,} read pairs; mapped both sides {int(dd['total_mapped']):,}; "
      f"duplicates {int(dd['total_dups']):,} ({float(dd['total_dups']) / float(dd['total_mapped']) * 100:.1f}%); "
      f"unique {int(dd['total_nodups']):,}; cis >1 kb {int(dd['cis_1kb+']):,}")

# ---- fig04: restriction fragments on chr21 ----------------------------------------------
with gzip.open(f"{D}/hg38_chr21.fa.gz", "rt") as fh:
    seq = "".join(l.strip() for l in fh if not l.startswith(">")).upper()
acgt = sum(seq.count(b) for b in "ACGT")
fig, ax = plt.subplots(figsize=(8.5, 4.3))
for site, name, col in [("GATC", "MboI / DpnII (GATC)", RED), ("AAGCTT", "HindIII (AAGCTT)", BLUE)]:
    pos = np.array([m.start() for m in re.finditer(site, seq)])
    frag = np.diff(pos)
    # skip fragments that span assembly gaps (runs of N)
    nn = np.cumsum(np.frombuffer(seq.encode(), dtype=np.uint8) == ord("N"))
    keep = (nn[pos[1:]] - nn[pos[:-1]]) == 0
    frag = frag[keep]
    print(f"fig04: {name}: {len(pos):,} sites on chr21 ({acgt / 1e6:.1f} Mb of ACGT); expected 1 per "
          f"{4 ** len(site):,} bp = {acgt / 4 ** len(site):,.0f}; median fragment {np.median(frag):.0f} bp, mean {frag.mean():.0f} bp")
    ax.hist(np.log10(frag[frag > 0]), bins=80, range=(0, 5.5), color=col, alpha=0.6, label=f"{name}: median {np.median(frag):.0f} bp",
            density=True)
ax.set_xticks([0, 1, 2, 3, 4, 5]); ax.set_xticklabels(["1", "10", "100", "1 k", "10 k", "100 k"])
ax.set_xlabel("restriction fragment length on chr21 (bp)"); ax.set_ylabel("density")
ax.set_title("A 4-cutter makes ~8× shorter fragments than a 6-cutter on chr21")
ax.set_ylim(0, 1.35)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig("fig04-fragments.png", dpi=150); plt.close(fig)

# junction motif in the reads
n = j = 0
pos_j = []
with gzip.open(f"{D}/gm12878_hic_R1.fastq.gz", "rt") as fh:
    for i, line in enumerate(fh):
        if i % 4 == 1:
            n += 1
            k = line.find("GATCGATC")
            if k >= 0:
                j += 1; pos_j.append(k)
        if n >= 200000:
            break
print(f"junction GATCGATC in {j / n * 100:.1f}% of the first {n:,} R1 reads; median position {np.median(pos_j):.0f}")

# ---- fig08: pair types (full genome) -----------------------------------------------------
pt = full[full.index.str.startswith("pair_types/")]
pt.index = pt.index.str.replace("pair_types/", "")
pt = pt.astype(float).drop("DD", errors="ignore").sort_values(ascending=False)   # DD = duplicates (dedup)
fig, ax = plt.subplots(figsize=(8.5, 4.2))
cols = [GREEN if k in ("UU", "UR", "RU") else GREY for k in pt.index]
ax.bar(pt.index, pt.values / pt.sum() * 100, color=cols)
for i, v in enumerate(pt.values / pt.sum() * 100):
    ax.text(i, v + 1, f"{v:.1f}", ha="center", fontsize=10)
ax.set_ylabel("% of read pairs"); ax.set_xlabel("pair type (side 1, side 2)")
ax.set_title("30 M read pairs, whole genome: green = kept as contacts")
fig.tight_layout(); fig.savefig("fig08-pair-types.png", dpi=150); plt.close(fig)
print("pair types (% of all pairs):", (pt / pt.sum() * 100).round(1).to_dict())

# ---- fig10: orientation vs distance (full genome) ------------------------------------------
df = full[full.index.str.startswith("dist_freq/")].astype(float)
rows = []
for k, v in df.items():
    m = re.match(r"dist_freq/(\d+)-?(\d+|\+)?/(..)", k)
    if m:
        rows.append((int(m.group(1)), m.group(3), v))
dfq = pd.DataFrame(rows, columns=["lo", "strand", "n"]).pivot(index="lo", columns="strand", values="n")
dfq = dfq[(dfq.index >= 100) & (dfq.index < 1e6)]
frac = dfq.div(dfq.sum(axis=1), axis=0)
print("fig10: % by orientation:\n", (frac * 100).round(1).to_string())
fig, ax = plt.subplots(figsize=(8.5, 4.5))
names = {"+-": "+ −  inward (dangling ends)", "-+": "− +  outward (self-circles)", "++": "+ +", "--": "− −"}
colors = {"+-": RED, "-+": BLUE, "++": GREY, "--": "#9aa5b1"}
for s in ["+-", "-+", "++", "--"]:
    ax.semilogx(frac.index, frac[s] * 100, "-o", ms=3, color=colors[s], label=names[s], lw=2)
ax.axhline(25, color=INK, ls="--", lw=1)
ax.set_xlabel("distance between the two sides (bp)"); ax.set_ylabel("% of pairs at this distance")
ax.set_title("Orientations equalise above ~1–2 kb")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("fig10-orientation.png", dpi=150); plt.close(fig)
conv = frac[(frac[["+-", "-+", "++", "--"]].sub(0.25).abs() < 0.03).all(axis=1)]
print(f"fig10: orientations within 3 points of 25% from {conv.index.min() if len(conv) else 'n/a'} bp on")

# ---- fig12: read funnel (full genome) --------------------------------------------------------
steps_ = [("read pairs", f("total")),
          ("both sides mapped", f("total_mapped")),
          ("after removing duplicates", f("total_nodups")),
          ("cis (same chromosome)", f("cis")),
          ("cis, > 1 kb apart", f("cis_1kb+")),
          ("cis, > 20 kb apart", f("cis_20kb+"))]
fig, ax = plt.subplots(figsize=(8.8, 4.2))
y = np.arange(len(steps_))[::-1]
ax.barh(y, [v / 1e6 for _, v in steps_], color=[GREY, BLUE, BLUE, GREEN, GREEN, GREEN])
for yy, (lab, v) in zip(y, steps_):
    ax.text(v / 1e6 + 0.3, yy, f"{v / 1e6:.1f} M ({v / f('total') * 100:.0f}%)", va="center", fontsize=12)
ax.set_yticks(y); ax.set_yticklabels([s for s, _ in steps_])
ax.set_xlabel("million read pairs"); ax.set_xlim(0, f("total") / 1e6 * 1.3)
ax.set_title("Where 30 M read pairs go (whole genome)")
fig.tight_layout(); fig.savefig("fig12-read-funnel.png", dpi=150); plt.close(fig)

# ---- our chr21 map vs ENCODE -------------------------------------------------------------
ours = lambda r: cooler.Cooler(f"{W}/chr21.mcool::resolutions/{r}")
enc = lambda r: cooler.Cooler(ENC + str(r))
o100, e100 = ours(100000), enc(100000)
def n_contacts(c):   # each contact once: upper triangle including the diagonal
    m = c.matrix(balance=False, sparse=True).fetch("chr21")
    return int((m.sum() + m.diagonal().sum()) // 2)
n_ours, n_enc = n_contacts(o100), n_contacts(e100)
print(f"contacts on chr21: ours {n_ours:,} (cis, all distances), ENCODE {n_enc:,}: {n_enc / n_ours:.0f}x more")
for r in [10000, 50000, 100000, 250000]:
    for lab, c in (("ours", ours(r)), ("ENCODE", enc(r))):
        m = c.matrix(balance=False, sparse=True).fetch("chr21")
        cov = np.asarray(m.sum(axis=1)).ravel()          # the matrix is symmetric: row sum = contacts of the bin
        cov = cov[cov > 0]
        print(f"  {r // 1000:>4} kb {lab:7s}: {np.mean(cov >= 1000) * 100:5.1f}% of non-empty bins with ≥ 1000 contacts")

# fig14: same region, three resolutions, ours vs ENCODE
reg = "chr21:33,000,000-38,000,000"
fig, axs = plt.subplots(2, 3, figsize=(12, 8))
for j, r in enumerate([250000, 50000, 10000]):
    for i, (lab, c) in enumerate((("our chr21 reads", ours(r)), ("ENCODE, all replicates", enc(r)))):
        m = c.matrix(balance=True).fetch(reg)
        ax = axs[i, j]
        ax.imshow(m, cmap=HIC, norm=LogNorm(vmin=np.nanpercentile(m[m > 0], 10) if np.any(m > 0) else 1e-4,
                                             vmax=np.nanpercentile(m, 99.5)),
                  extent=[33, 38, 38, 33], interpolation="none")
        ax.set_title(f"{lab}\n{r // 1000} kb bins", fontsize=12)
        ax.set_xticks([33, 35, 37]); ax.set_yticks([33, 35, 37])
        for sp in ax.spines.values():
            sp.set_visible(True)
fig.supxlabel("chr21 (Mb)")
fig.tight_layout(); fig.savefig("fig14-resolution-depth.png", dpi=150); plt.close(fig)

# fig18: raw vs balanced (ours, 250 kb), with coverage
c = ours(250000)
raw = c.matrix(balance=False).fetch("chr21").astype(float)
bal = c.matrix(balance=True).fetch("chr21")
L = c.chromsizes["chr21"] / 1e6
cov = raw.sum(axis=0)
w = c.bins().fetch("chr21")["weight"].values
print(f"fig18: 250-kb bins on chr21: {len(w)}, masked (NaN weight) {np.isnan(w).sum()}; coverage of unmasked bins "
      f"from {np.nanmin(cov[~np.isnan(w)]):.0f} to {np.nanmax(cov[~np.isnan(w)]):.0f} contacts")
fig = plt.figure(figsize=(12, 6.2))
gs = fig.add_gridspec(2, 2, height_ratios=[5, 1.1], hspace=0.08, wspace=0.3)
for k, (M, title) in enumerate([(raw, "raw counts"), (bal, "balanced (ICE)")]):
    ax = fig.add_subplot(gs[0, k])
    MM = np.where(M > 0, M, np.nan)
    ax.imshow(MM, cmap=HIC, norm=LogNorm(vmin=np.nanpercentile(MM, 5), vmax=np.nanpercentile(MM, 99.5)),
              extent=[0, L, L, 0], interpolation="none")
    ax.set_title(f"our chr21 map, 250 kb: {title}")
    plt.setp(ax.get_xticklabels(), visible=False)
    a = fig.add_subplot(gs[1, k], sharex=ax)
    x = (np.arange(len(cov)) + 0.5) * 0.25
    far = np.abs(np.subtract.outer(np.arange(len(M)), np.arange(len(M)))) >= 2   # balancing ignores 2 diagonals
    y = np.nansum(np.where(far, M, 0), axis=0)
    y = np.where(np.isnan(w), 0, y) if k == 1 else y
    a.bar(x, y, width=0.25, color=BLUE if k == 0 else GREEN)
    a.set_ylabel("row sum\n(|i−j| ≥ 2)", rotation=0, ha="right", va="center", fontsize=11); a.set_yticks([])
    a.set_xlabel("chr21 (Mb)")
fig.savefig("fig18-balance-real.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# fig19: P(s) ours vs ENCODE (10 kb)
fig, ax = plt.subplots(figsize=(8, 5))
for lab, c, col in (("our chr21 reads", ours(10000), ORANGE), ("ENCODE", enc(10000), BLUE)):
    v = pd.DataFrame({"chrom": ["chr21"], "start": [0], "end": [c.chromsizes["chr21"]], "name": ["chr21"]})
    e = cooltools.expected_cis(c, view_df=v, nproc=2, smooth=True, aggregate_smoothed=False)
    e["s"] = e.dist * 10000
    e = e[(e.dist >= 2) & (e.s <= 2e7)]
    yv = e["balanced.avg.smoothed"] / e["balanced.avg.smoothed"].sum()
    ax.loglog(e.s, yv, color=col, lw=2, label=lab)
    sel = (e.s >= 1e5) & (e.s <= 1e6)
    print(f"fig19: {lab}: P(s) slope 100 kb–1 Mb {np.polyfit(np.log10(e.s[sel]), np.log10(e['balanced.avg'][sel]), 1)[0]:.2f}")
ax.set_xlabel("genomic distance s (bp)"); ax.set_ylabel("P(s) (normalised to sum 1)")
ax.set_title("Distance decay on chr21: few reads, same curve")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("fig19-pofs.png", dpi=150); plt.close(fig)

# fig20: balanced vs observed/expected (ENCODE chr19, 100 kb, one arm region)
c = enc(100000)
v = pd.DataFrame({"chrom": ["chr19"], "start": [0], "end": [c.chromsizes["chr19"]], "name": ["chr19"]})
e = cooltools.expected_cis(c, view_df=v, nproc=2).set_index("dist")["balanced.avg"]
reg, s0, s1 = "chr19:30,000,000-58,617,616", 30.0, 58.617616
m = c.matrix(balance=True).fetch(reg)
nb = m.shape[0]
i0 = int(30_000_000 // 100000)
dist = np.abs(np.subtract.outer(np.arange(nb), np.arange(nb)))
oe = m / e.reindex(np.arange(nb)).values[dist]
fig, (a0, a1) = plt.subplots(1, 2, figsize=(12, 5.6))
a0.imshow(m, cmap=HIC, norm=LogNorm(vmin=np.nanpercentile(m[m > 0], 5), vmax=np.nanpercentile(m, 99.5)),
          extent=[s0, s1, s1, s0], interpolation="none")
a0.set_title("balanced: dominated by distance")
im = a1.imshow(np.log2(oe), cmap=DIV, norm=TwoSlopeNorm(0, -1.5, 1.5), extent=[s0, s1, s1, s0], interpolation="none")
a1.set_title("log2 observed / expected: the checkerboard")
fig.colorbar(im, ax=a1, fraction=0.045)
for a in (a0, a1):
    a.set_xlabel("chr19 (Mb)")
    for sp in a.spines.values():
        sp.set_visible(True)
fig.tight_layout(); fig.savefig("fig20-obs-exp.png", dpi=150); plt.close(fig)
print("done")
