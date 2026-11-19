# Real-data figures for hic-features/lec.qmd ("Calling compartments, TADs and loops").
# Run from this folder:  python3 make-hicfeatures-figs.py
#
# Needs (data/download.sh hic-features; data/prepare-hic.sh matrices):
#   data/gm12878_hic_chr19_chr21.mcool, data/gm12878_hic_loops_hiccups.bedpe.gz,
#   data/gm12878_hic_domains_arrowhead.bedpe.gz, data/gm12878_hic_eigenvector_10kb.bigWig,
#   data/gm12878_ctcf_peaks.bed.gz, data/gm12878_ctcf_fc.bigWig, data/gm12878_h3k27ac_fc.bigWig,
#   data/gm12878_h3k9me3_fc.bigWig, data/hg38.fa.bgz (+ .fai, .gzi);
#   Python: cooler 0.10, cooltools 0.7, bioframe, pyfaidx, pyBigWig, numpy, pandas, matplotlib;
#   samtools on the PATH; R with GenomicRanges, rtracklayer, Rsamtools, JASPAR2024, TFBSTools,
#   motifmatchr (for the CTCF motif orientation at loop anchors, the same code as exercise R1).
#
# Writes fig03-eigenvector-real.png, fig04-eigenvector-marks.png, fig07-saddle-real.png,
# fig09-insulation-real.png, fig10-insulation-windows.png, fig11-boundaries-ctcf.png,
# fig14-dots-real.png, fig15-ctcf-orientation.png, fig17-apa-real.png
# and prints the numbers quoted in the slides. Work files go to ./features-work.
import os, subprocess, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
import cooler, cooltools, bioframe, pyBigWig

warnings.filterwarnings("ignore")
D = os.path.abspath("../data")
W = "features-work"
os.makedirs(W, exist_ok=True)
MC = f"{D}/gm12878_hic_chr19_chr21.mcool::resolutions/"
plt.rcParams.update({"font.size": 13, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 14})
HIC = LinearSegmentedColormap.from_list("hic", ["#ffffff", "#f4b6a6", "#c0392b", "#5a1209"])
DIV = LinearSegmentedColormap.from_list("div", ["#2a6f97", "#ffffff", "#c0392b"])
GREEN, CRIMSON, K27AC, K9, CTCF, INK, GREY, BLUE = ("#2e8b57", "#a8324a", "#e67e22", "#6b4fa0", "#8c564b",
                                                    "#1b2a3a", "#5b6b7a", "#2a6f97")
CHROMS = ["chr19", "chr21"]


def view_of(clr):
    return pd.DataFrame({"chrom": CHROMS, "start": 0, "end": [clr.chromsizes[c] for c in CHROMS], "name": CHROMS})


def bw_mean(path, df):
    b = pyBigWig.open(path)
    v = [b.stats(r.chrom, int(r.start), int(min(r.end, b.chroms(r.chrom))), type="mean")[0] for r in df.itertuples()]
    b.close()
    return np.array([np.nan if x is None else x for x in v], dtype=float)


def read_bedpe(fn):
    t = pd.read_csv(fn, sep="\t", comment="#", header=None).iloc[:, :6]
    t.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    return t


def triangle(ax, mat, start, binsize, depth, norm, cmap=HIC):
    n = mat.shape[0]
    edges = start + np.arange(n + 1) * binsize
    X, Y = np.meshgrid(edges, edges)
    Xr, Yr = (X + Y) / 2 / 1e6, (X - Y) / 2 / 1e6
    m = np.ma.masked_invalid(mat.astype(float))
    m[np.tril_indices(n, -1)] = np.ma.masked
    ax.pcolormesh(Xr, Yr, m, cmap=cmap, norm=norm, shading="flat", rasterized=True)
    ax.set_ylim(0, depth / 2 / 1e6)
    ax.set_xlim(edges[0] / 1e6, edges[-1] / 1e6)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)


# ===================================================================================
# 1. Compartments (100 kb)
# ===================================================================================
c100 = cooler.Cooler(MC + "100000")
view = view_of(c100)
fa = f"{W}/chr19_chr21.fa"
if not os.path.exists(fa):
    subprocess.run(f"samtools faidx {D}/hg38.fa.bgz chr19 chr21 > {fa}", shell=True, check=True)
bins = c100.bins()[:][["chrom", "start", "end"]]
gc = bioframe.frac_gc(bins, bioframe.load_fasta(fa))
ev_vals, ev = cooltools.eigs_cis(c100, view_df=view, phasing_track=gc, n_eigs=3)
ev["GC"] = gc["GC"].values
ev["H3K27ac"] = bw_mean(f"{D}/gm12878_h3k27ac_fc.bigWig", ev)
ev["H3K9me3"] = bw_mean(f"{D}/gm12878_h3k9me3_fc.bigWig", ev)
ev["ENCODE"] = bw_mean(f"{D}/gm12878_hic_eigenvector_10kb.bigWig", ev)
ev.to_csv(f"{W}/eigs_100kb.tsv", sep="\t", index=False)
for ch in CHROMS:
    e = ev[(ev.chrom == ch) & ev.E1.notna()]
    r = {k: e.E1.corr(e[k], method="spearman") for k in ["GC", "H3K27ac", "H3K9me3", "ENCODE"]}
    print(f"compartments {ch}: {len(e)} bins with E1, {np.mean(e.E1 > 0) * 100:.0f}% A; "
          f"Spearman E1 vs GC {r['GC']:.2f}, H3K27ac {r['H3K27ac']:.2f}, H3K9me3 {r['H3K9me3']:.2f}, "
          f"ENCODE eigenvector {r['ENCODE']:.2f}; E2 vs GC {e.E2.corr(e.GC, method='spearman'):.2f}")
print("eigenvalues:", ev_vals[["name", "eigval1", "eigval2", "eigval3"]].round(1).to_string(index=False))

# fig03: chr19 O/E, correlation, E1
ch = "chr19"
exp100 = cooltools.expected_cis(c100, view_df=view, nproc=2)
lo, hi = c100.extent(ch)
A = c100.matrix(balance=True).fetch(ch)
n = A.shape[0]
e = exp100[exp100.region1 == ch].set_index("dist")["balanced.avg"]
dist = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
OE = A / e.reindex(np.arange(n)).values[dist]
OEl = np.log2(OE)
OEl[~np.isfinite(OEl)] = np.nan
good = ~np.all(np.isnan(OEl), axis=1)
C = np.full_like(OEl, np.nan)
Cg = pd.DataFrame(OEl[np.ix_(good, good)]).T.corr().values
C[np.ix_(good, good)] = Cg
L = c100.chromsizes[ch] / 1e6
e1 = ev[ev.chrom == ch].E1.values
fig = plt.figure(figsize=(12, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[5, 1], hspace=0.06, wspace=0.38)
for k, (M, title, norm) in enumerate([(OEl, "log2 observed / expected", TwoSlopeNorm(0, -1.5, 1.5)),
                                      (C, "Pearson correlation of rows", TwoSlopeNorm(0, -0.8, 0.8))]):
    ax = fig.add_subplot(gs[0, k])
    im = ax.imshow(M, cmap=DIV, norm=norm, extent=[0, L, L, 0], interpolation="none")
    ax.set_title(f"{ch}, 100 kb: {title}")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_ylabel("Mb")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    a = fig.add_subplot(gs[1, k], sharex=ax)
    x = (np.arange(n) + 0.5) * 0.1
    a.fill_between(x, 0, np.where(e1 > 0, e1, 0), color=GREEN, step="mid", lw=0)
    a.fill_between(x, 0, np.where(e1 < 0, e1, 0), color=CRIMSON, step="mid", lw=0)
    a.set_yticks([]); a.set_ylabel("E1", rotation=0, ha="right", va="center")
    a.set_xlabel(f"{ch} (Mb)")
    pos = ax.get_position(); p2 = a.get_position()
    a.set_position([pos.x0, p2.y0, pos.width, p2.height])
fig.savefig("fig03-eigenvector-real.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# fig04: E1 against GC, H3K27ac, H3K9me3
fig, axs = plt.subplots(1, 3, figsize=(12, 4.0))
cols = {"chr19": BLUE, "chr21": "#e67e22"}
for ax, (k, lab, logy) in zip(axs, [("GC", "GC content", False), ("H3K27ac", "H3K27ac (fold, log)", True),
                                     ("H3K9me3", "H3K9me3 (fold, log)", True)]):
    for chm in CHROMS:
        d = ev[(ev.chrom == chm) & ev.E1.notna()]
        ax.scatter(d.E1, d[k], s=7, alpha=0.6, color=cols[chm], label=chm, lw=0)
    d = ev[ev.E1.notna()]
    ax.set_title(f"{lab.split(' (')[0]}: Spearman {d.E1.corr(d[k], method='spearman'):.2f}")
    if logy:
        ax.set_yscale("log")
    ax.axvline(0, color=GREY, lw=0.8, ls="--")
    ax.set_xlabel("E1  (B ← → A)"); ax.set_ylabel(lab)
axs[0].legend(frameon=False, markerscale=2.5, loc="upper left")
fig.tight_layout()
fig.savefig("fig04-eigenvector-marks.png", dpi=150)
plt.close(fig)

# fig07: saddle plot (cooltools, 38 quantile bins) + 5-group strength
track = ev[["chrom", "start", "end", "E1"]]
S, Cn = cooltools.saddle(c100, exp100, track, "cis", n_bins=38, qrange=(0.025, 0.975), view_df=view)
sad = (S / Cn)[1:-1, 1:-1]
S5, C5 = cooltools.saddle(c100, exp100, track, "cis", n_bins=5, qrange=(0.0, 1.0), view_df=view)
s5 = (S5 / C5)[1:-1, 1:-1]
strength = (s5[0, 0] + s5[4, 4]) / (s5[0, 4] + s5[4, 0])
print("saddle, 5 E1 groups (rows B -> A), mean O/E:\n", np.round(s5, 2))
print(f"compartment strength (BB + AA) / (BA + AB) = ({s5[0,0]:.2f} + {s5[4,4]:.2f}) / ({s5[0,4]:.2f} + {s5[4,0]:.2f}) = {strength:.2f}")
fig, (a0, a1) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})
im = a0.imshow(np.log2(sad), cmap=DIV, norm=TwoSlopeNorm(0, -1, 1))
a0.set_title("saddle plot, chr19 + chr21 (100 kb)")
a0.set_xlabel("E1 quantile (B → A)"); a0.set_ylabel("E1 quantile (B → A)")
a0.set_xticks([0, 35]); a0.set_xticklabels(["B", "A"]); a0.set_yticks([0, 35]); a0.set_yticklabels(["B", "A"])
fig.colorbar(im, ax=a0, fraction=0.045, label="log2 mean O/E")
im = a1.imshow(np.log2(s5), cmap=DIV, norm=TwoSlopeNorm(0, -1, 1))
for i in range(5):
    for j in range(5):
        a1.text(j, i, f"{s5[i, j]:.2f}", ha="center", va="center", fontsize=11)
a1.set_xticks(range(5)); a1.set_xticklabels(["1\nB", "2", "3", "4", "5\nA"])
a1.set_yticks(range(5)); a1.set_yticklabels(["1 B", "2", "3", "4", "5 A"])
a1.set_title(f"5 groups: strength {strength:.1f}")
fig.tight_layout()
fig.savefig("fig07-saddle-real.png", dpi=150)
plt.close(fig)

# ===================================================================================
# 2. Insulation and boundaries (10 kb)
# ===================================================================================
c10 = cooler.Cooler(MC + "10000")
v10 = view_of(c10)
windows = [100_000, 200_000, 500_000]
ins = cooltools.insulation(c10, windows, view_df=v10, verbose=False)
ins.to_csv(f"{W}/insulation_10kb.tsv", sep="\t", index=False)
ctcf = pd.read_csv(f"{D}/gm12878_ctcf_peaks.bed.gz", sep="\t", header=None).iloc[:, [0, 1, 2, 6, 9]]
ctcf.columns = ["chrom", "start", "end", "signal", "summit"]
ctcf = ctcf[ctcf.chrom.isin(CHROMS)]
summits = pd.DataFrame({"chrom": ctcf.chrom, "start": ctcf.start + ctcf.summit, "end": ctcf.start + ctcf.summit + 1})
dom = read_bedpe(f"{D}/gm12878_hic_domains_arrowhead.bedpe.gz")
dom = dom[dom.chrom1.isin(CHROMS)]
edges = pd.concat([pd.DataFrame({"chrom": dom.chrom1, "start": dom.start1, "end": dom.start1 + 1}),
                   pd.DataFrame({"chrom": dom.chrom1, "start": dom.end1 - 1, "end": dom.end1})]).drop_duplicates().reset_index(drop=True)
valid = ins[~ins.is_bad_bin]
rng = np.random.default_rng(1)
print(f"Arrowhead domains on chr19 + chr21: {len(dom)}; distinct domain edges {len(edges)}")
for w in windows:
    b = ins[ins[f"is_boundary_{w}"]][["chrom", "start", "end"]]
    near_ctcf = bioframe.overlap(b, bioframe.expand(summits, pad=10_000), how="left", return_index=True)
    frac = near_ctcf.groupby("index")["index_"].apply(lambda x: x.notna().any()).mean()
    rb = valid.iloc[rng.choice(len(valid), len(b), replace=False)][["chrom", "start", "end"]]
    near_r = bioframe.overlap(rb, bioframe.expand(summits, pad=10_000), how="left", return_index=True)
    frac_r = near_r.groupby("index")["index_"].apply(lambda x: x.notna().any()).mean()
    near_edge = bioframe.overlap(b, bioframe.expand(edges, pad=20_000), how="left", return_index=True)
    fe = near_edge.groupby("index")["index_"].apply(lambda x: x.notna().any()).mean()
    print(f"insulation window {w // 1000} kb: {len(b)} boundaries; within 10 kb of a CTCF summit {frac * 100:.0f}% "
          f"(random bins {frac_r * 100:.0f}%); within 20 kb of an Arrowhead domain edge {fe * 100:.0f}%")
b200 = ins[ins["is_boundary_200000"]][["chrom", "start", "end"]]
eb = bioframe.overlap(edges, bioframe.expand(b200, pad=20_000), how="left", return_index=True)
print(f"  Arrowhead edges within 20 kb of a 200-kb insulation boundary: "
      f"{eb.groupby('index')['index_'].apply(lambda x: x.notna().any()).mean() * 100:.0f}%")

# fig09: insulation along a chr21 window with boundaries and CTCF
ch, s0, s1 = "chr21", 37_000_000, 42_000_000
m = c10.matrix(balance=True).fetch(f"{ch}:{s0}-{s1}")
fig = plt.figure(figsize=(12, 6.4))
gs = fig.add_gridspec(3, 1, height_ratios=[4, 1.5, 0.8], hspace=0.08)
ax = fig.add_subplot(gs[0])
triangle(ax, m, s0, 10000, 1_500_000, LogNorm(vmin=np.nanpercentile(m[m > 0], 30), vmax=np.nanpercentile(m, 99.5)))
ax.set_title(f"{ch}:{s0 / 1e6:.0f}–{s1 / 1e6:.0f} Mb, 10 kb: insulation score (200-kb window) and boundaries")
plt.setp(ax.get_xticklabels(), visible=False)
a = fig.add_subplot(gs[1], sharex=ax)
iw = ins[(ins.chrom == ch) & (ins.start >= s0) & (ins.end <= s1)]
a.plot((iw.start + 5000) / 1e6, iw["log2_insulation_score_200000"], color=INK, lw=1.2)
bb = iw[iw["is_boundary_200000"]]
a.plot((bb.start + 5000) / 1e6, bb["log2_insulation_score_200000"], "o", color="#c0392b", ms=5)
for x in (bb.start + 5000) / 1e6:
    ax.axvline(x, color="#c0392b", lw=0.6, ls=":")
a.axhline(0, color=GREY, lw=0.6, ls="--")
a.set_ylabel("insulation\n(log2)", rotation=0, ha="right", va="center")
plt.setp(a.get_xticklabels(), visible=False)
a2 = fig.add_subplot(gs[2], sharex=ax)
b = pyBigWig.open(f"{D}/gm12878_ctcf_fc.bigWig")
nb = (s1 - s0) // 5000
cv = np.nan_to_num(np.array(b.stats(ch, s0, s1, nBins=nb, type="max"), dtype=float))
a2.fill_between(s0 / 1e6 + (np.arange(nb) + 0.5) * 0.005, 0, cv, color=CTCF, lw=0)
a2.set_yticks([]); a2.set_ylabel("CTCF", rotation=0, ha="right", va="center")
a2.set_xlim(s0 / 1e6, s1 / 1e6); a2.set_xlabel(f"{ch} position (Mb)")
fig.savefig("fig09-insulation-real.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"fig09: {len(bb)} boundaries (200-kb window) in {ch}:{s0}-{s1}")

# fig10: three windows
fig, ax = plt.subplots(figsize=(12, 3.6))
for w, col in zip(windows, ["#e67e22", INK, BLUE]):
    ax.plot((iw.start + 5000) / 1e6, iw[f"log2_insulation_score_{w}"], color=col, lw=1.2, label=f"{w // 1000} kb")
    bw_ = iw[iw[f"is_boundary_{w}"]]
    ax.plot((bw_.start + 5000) / 1e6, bw_[f"log2_insulation_score_{w}"], "o", color=col, ms=4.5)
    print(f"fig10: window {w // 1000} kb: {len(bw_)} boundaries in the window")
ax.axhline(0, color=GREY, lw=0.6, ls="--")
ax.set_xlim(s0 / 1e6, s1 / 1e6); ax.set_xlabel(f"{ch} position (Mb)"); ax.set_ylabel("log2 insulation")
ax.legend(title="window", frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.25))
fig.tight_layout()
fig.savefig("fig10-insulation-windows.png", dpi=150)
plt.close(fig)

# fig11: CTCF ChIP signal around boundaries (100-kb window: the most boundaries), strong vs weak
bwc = pyBigWig.open(f"{D}/gm12878_ctcf_fc.bigWig")
def profile(bnds, flank=300_000, step=10_000):
    offs = np.arange(-flank, flank + 1, step)
    out = []
    for c, st in zip(bnds.chrom.values, bnds.start.values):
        x = st + 5000
        if x - flank - step < 0 or x + flank + step > bwc.chroms(c):
            continue
        v = bwc.stats(c, int(x - flank - step // 2), int(x + flank + step // 2), nBins=len(offs), type="mean")
        out.append([np.nan if y is None else y for y in v])
    return offs / 1e3, np.nanmean(np.array(out, dtype=float), axis=0)
B = ins[ins["is_boundary_100000"]].copy()
q = B["boundary_strength_100000"].quantile([1 / 3, 2 / 3]).values
fig, ax = plt.subplots(figsize=(7.5, 4.4))
for lab, sel, col in [("strongest third", B["boundary_strength_100000"] > q[1], "#c0392b"),
                      ("weakest third", B["boundary_strength_100000"] <= q[0], "#e8a090")]:
    x, y = profile(B[sel])
    ax.plot(x, y, color=col, lw=2.2, label=f"{lab} (n = {sel.sum()})")
    print(f"fig11: {lab}: mean CTCF fold change at the boundary {y[len(y) // 2]:.2f}, 300 kb away {np.mean([y[0], y[-1]]):.2f}")
x, y = profile(valid.iloc[rng.choice(len(valid), 500, replace=False)])
ax.plot(x, y, color=GREY, lw=1.6, ls="--", label="random bins (n = 500)")
print(f"fig11: random bins: {y[len(y) // 2]:.2f}")
ax.set_xlabel("distance from boundary (kb)"); ax.set_ylabel("mean CTCF ChIP fold change")
ax.set_title("CTCF binding peaks at boundaries (100-kb window)")
ax.set_ylim(top=ax.get_ylim()[1] + 0.1)
ax.legend(frameon=False, loc="upper right", fontsize=11)
fig.tight_layout()
fig.savefig("fig11-boundaries-ctcf.png", dpi=150)
plt.close(fig)

# ===================================================================================
# 3. Dots (10 kb) and comparison with HiCCUPS
# ===================================================================================
exp10 = cooltools.expected_cis(c10, view_df=v10, nproc=4)
dots_fn = f"{W}/dots_10kb.tsv"
if os.path.exists(dots_fn):
    dots = pd.read_csv(dots_fn, sep="\t")
else:
    dots = cooltools.dots(c10, expected=exp10, view_df=v10, max_loci_separation=2_000_000, nproc=4)
    dots.to_csv(dots_fn, sep="\t", index=False)
hic = read_bedpe(f"{D}/gm12878_hic_loops_hiccups.bedpe.gz")
hic = hic[hic.chrom1.isin(CHROMS)]
def match(a, b, tol=20_000):
    """fraction of loops in a with a loop in b whose two anchors are both within tol."""
    ok = 0
    for r in a.itertuples():
        bb_ = b[b.chrom1 == r.chrom1]
        d1 = np.abs((bb_.start1 + bb_.end1) / 2 - (r.start1 + r.end1) / 2)
        d2 = np.abs((bb_.start2 + bb_.end2) / 2 - (r.start2 + r.end2) / 2)
        ok += np.any((d1 <= tol) & (d2 <= tol))
    return ok / len(a)
hic2 = hic[(hic.start2 - hic.start1) <= 2_000_000]
print(f"dots: cooltools found {len(dots)} on chr19 + chr21 (10 kb); HiCCUPS has {len(hic)} "
      f"({len(hic2)} spanning ≤ 2 Mb); {match(dots, hic) * 100:.0f}% of ours match a HiCCUPS loop (±20 kb); "
      f"{match(hic2, dots) * 100:.0f}% of HiCCUPS loops (≤ 2 Mb) are found")
ch, s0, s1 = "chr21", 24_750_000, 26_250_000
m = c10.matrix(balance=True).fetch(f"{ch}:{s0}-{s1}")
fig, ax = plt.subplots(figsize=(7.2, 6.6))
ax.imshow(m, cmap=HIC, norm=LogNorm(vmin=np.nanpercentile(m[m > 0], 30), vmax=np.nanpercentile(m, 99.7)),
          extent=[s0 / 1e6, s1 / 1e6, s1 / 1e6, s0 / 1e6], interpolation="none")
hw = hic[(hic.chrom1 == ch) & (hic.start1 >= s0) & (hic.end2 <= s1)]
dw = dots[(dots.chrom1 == ch) & (dots.start1 >= s0) & (dots.end2 <= s1)]
ax.scatter((hw.start2 + hw.end2) / 2e6, (hw.start1 + hw.end1) / 2e6, s=230, facecolors="none", edgecolors=BLUE,
           lw=1.5, label=f"HiCCUPS (ENCODE, {len(hw)})")
ax.scatter((dw.start1 + dw.end1) / 2e6, (dw.start2 + dw.end2) / 2e6, s=150, facecolors="none", edgecolors=INK,
           marker="s", lw=1.5, label=f"cooltools dots (ours, {len(dw)})")
ax.set_title(f"{ch}:{s0 / 1e6:.2f}–{s1 / 1e6:.2f} Mb, 10 kb")
ax.set_xlabel("Mb"); ax.set_ylabel("Mb")
ax.legend(loc="lower left", frameon=True, fontsize=11)
for sp in ax.spines.values():
    sp.set_visible(True)
fig.tight_layout()
fig.savefig("fig14-dots-real.png", dpi=150)
plt.close(fig)
print(f"fig14: window {ch}:{s0}-{s1}: {len(hw)} HiCCUPS loops (upper triangle), {len(dw)} cooltools dots (lower triangle)")

# ===================================================================================
# 4. CTCF motif orientation at loop anchors (R, motifmatchr) and APA
# ===================================================================================
orient_fn = f"{W}/loop_ctcf_orientation.tsv"
if not os.path.exists(orient_fn):
    R = f'''
suppressPackageStartupMessages({{
  library(GenomicRanges); library(rtracklayer); library(Rsamtools); library(Biostrings)
  library(JASPAR2024); library(TFBSTools); library(RSQLite); library(motifmatchr)
}})
D <- "{D}"
lp <- read.table(file.path(D, "gm12878_hic_loops_hiccups.bedpe.gz"), sep = "\\t", comment.char = "#")[, 1:6]
names(lp) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
a1 <- GRanges(lp$chrom1, IRanges(lp$start1 + 1, lp$end1))
a2 <- GRanges(lp$chrom2, IRanges(lp$start2 + 1, lp$end2))
ctcf <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak")
ctcf <- ctcf[overlapsAny(ctcf, c(a1, a2))]
win <- resize(GRanges(seqnames(ctcf), IRanges(start(ctcf) + ctcf$peak, width = 1)), 100, fix = "center")
ctcf_m <- getMatrixByID(dbConnect(SQLite(), db(JASPAR2024())), "MA0139.2")
pos <- matchMotifs(ctcf_m, getSeq(FaFile(file.path(D, "hg38.fa.bgz")), win), p.cutoff = 1e-4, bg = "even",
                   out = "positions")[[1]]
ctcf$motif <- vapply(pos, function(x) if (length(x) == 0) NA_character_ else
                     as.character(mcols(x)$strand)[which.max(mcols(x)$score)], "")
anchor_motif <- function(a) {{
  h <- findOverlaps(a, ctcf)
  d <- data.frame(q = queryHits(h), s = subjectHits(h), sig = ctcf$signalValue[subjectHits(h)])
  d <- d[!is.na(ctcf$motif[d$s]), ]
  d <- d[order(d$q, -d$sig), ]; d <- d[!duplicated(d$q), ]
  out <- rep(NA_character_, length(a)); out[d$q] <- ctcf$motif[d$s]; out
}}
lp$motif1 <- anchor_motif(a1); lp$motif2 <- anchor_motif(a2)
write.table(lp, "{orient_fn}", sep = "\\t", quote = FALSE, row.names = FALSE)
'''
    subprocess.run(["Rscript", "-e", R], check=True)
lo_ = pd.read_csv(orient_fn, sep="\t")
both = lo_.dropna(subset=["motif1", "motif2"])
cls = (both.motif1 + both.motif2).map({"+-": "convergent → ←", "-+": "divergent ← →", "++": "tandem → →",
                                         "--": "tandem ← ←"})
tab = cls.value_counts()
print(f"orientation: {len(lo_)} loops; {lo_.motif1.notna().mean() * 100:.0f}% / {lo_.motif2.notna().mean() * 100:.0f}% of "
      f"left / right anchors have a CTCF peak with a motif; {len(both)} loops with both:")
print((tab / tab.sum() * 100).round(1).to_string())
fig, ax = plt.subplots(figsize=(7.2, 4.0))
order = ["convergent → ←", "tandem → →", "tandem ← ←", "divergent ← →"]
vals = [tab.get(k, 0) / tab.sum() * 100 for k in order]
ax.barh(order[::-1], vals[::-1], color=[GREY, GREY, GREY, CTCF])
for i, v in enumerate(vals[::-1]):
    ax.text(v + 1, i, f"{v:.0f}%", va="center")
ax.axvline(25, color=INK, ls="--", lw=1)
ax.text(26, -0.45, "25% if random", fontsize=11, color=INK)
ax.set_xlabel(f"% of {len(both):,} loops with a motif at both anchors")
ax.set_title("CTCF motifs at loop anchors point into the loop")
ax.set_xlim(0, 100)
fig.tight_layout()
fig.savefig("fig15-ctcf-orientation.png", dpi=150)
plt.close(fig)

# APA at 10 kb on chr19 + chr21: convergent loops, other loops, and the same loops shifted
lo19 = lo_[lo_.chrom1.isin(CHROMS)].copy()
lo19["cls"] = (lo19.motif1.fillna("0") + lo19.motif2.fillna("0"))
def to_feat(t):
    f = pd.DataFrame({"chrom1": t.chrom1, "start1": (t.start1 + t.end1) // 2 - 5000, "end1": (t.start1 + t.end1) // 2 + 5000,
                      "chrom2": t.chrom2, "start2": (t.start2 + t.end2) // 2 - 5000, "end2": (t.start2 + t.end2) // 2 + 5000})
    return f[(f.start2 - f.start1) >= 150_000]          # keep loops well away from the diagonal
sets = {"convergent CTCF loops": to_feat(lo19[lo19.cls == "+-"]),
        "all other loops": to_feat(lo19[lo19.cls != "+-"])}
shift = to_feat(lo19).copy()
for k in ["start1", "end1", "start2", "end2"]:
    shift[k] = shift[k] + 500_000
shift = shift[shift.end2 < shift.chrom2.map(c10.chromsizes)]
sets["control: loops shifted 500 kb"] = shift
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.2), gridspec_kw={"wspace": 0.3})
for ax, (lab, f) in zip(axs, sets.items()):
    st = cooltools.pileup(c10, f, view_df=v10, expected_df=exp10, flank=100_000)
    P = np.nanmean(st, axis=0)
    apa = P[10, 10] / np.nanmean(P[-6:, :6])
    im = ax.imshow(np.log2(P), cmap=DIV, norm=TwoSlopeNorm(0, -1.2, 1.2))
    ax.set_title(f"{lab}\nn = {len(f)}, APA = {apa:.1f}", fontsize=13)
    ax.set_xticks([0, 10, 20]); ax.set_xticklabels(["−100", "0", "+100"])
    ax.set_yticks([0, 10, 20]); ax.set_yticklabels(["−100", "0", "+100"])
    ax.set_xlabel("anchor 2 (kb)")
    if ax is axs[0]:
        ax.set_ylabel("anchor 1 (kb)")
    print(f"APA {lab}: n = {len(f)}, centre O/E {P[10, 10]:.2f}, lower-left corner {np.nanmean(P[-6:, :6]):.2f}, APA {apa:.2f}")
fig.colorbar(im, ax=axs, fraction=0.02, label="log2 mean O/E")
fig.savefig("fig17-apa-real.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("done")
