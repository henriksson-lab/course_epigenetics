# Real-data figures for genome-3d-intro/lec.qmd.
# Run from this folder:  python3 make-genome3d-figs.py
#
# Needs (data/download.sh genome-3d-intro, data/prepare-hic.sh matrices):
#   data/gm12878_hic_genome_1mb.cool, data/gm12878_hic_chr19_chr21.mcool,
#   data/gm12878_hic_loops_hiccups.bedpe.gz, data/gm12878_hic_domains_arrowhead.bedpe.gz,
#   data/gm12878_h3k27ac_fc.bigWig, data/gm12878_h3k9me3_fc.bigWig, data/gm12878_ctcf_fc.bigWig,
#   data/hg38_chr21.fa.gz is not needed here.
# Python packages: cooler 0.10, cooltools 0.7, pyBigWig, numpy, pandas, matplotlib.
#
# Writes fig03-territories.png, fig05-distance-decay.png, fig06-compartments.png,
# fig09-tads.png, fig10-loops.png and prints the numbers quoted in the slides.
import warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
import cooler, cooltools, pyBigWig

warnings.filterwarnings("ignore")
D = "../data"
MC = f"{D}/gm12878_hic_chr19_chr21.mcool::resolutions/"
plt.rcParams.update({"font.size": 13, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 14})
HIC = LinearSegmentedColormap.from_list("hic", ["#ffffff", "#f4b6a6", "#c0392b", "#5a1209"])
GREEN, CRIMSON, K27AC, K9, CTCF, INK, GREY = "#2e8b57", "#a8324a", "#e67e22", "#6b4fa0", "#8c564b", "#1b2a3a", "#5b6b7a"


def bw_bins(path, chrom, start, end, binsize):
    """mean of a bigWig in fixed bins (NaN-safe)."""
    b = pyBigWig.open(path)
    n = int(np.ceil((end - start) / binsize))
    stop = min(start + n * binsize, b.chroms(chrom))   # the last bin may be shorter
    v = np.array(b.stats(chrom, start, stop, nBins=n, type="mean"), dtype=float)
    b.close()
    return np.nan_to_num(v)


def triangle(ax, mat, start, binsize, depth, norm, cmap=HIC):
    """Draw the upper triangle of a square matrix rotated by 45 degrees (x in Mb, y = distance in Mb)."""
    n = mat.shape[0]
    edges = start + np.arange(n + 1) * binsize
    X, Y = np.meshgrid(edges, edges)
    # rotate: x = (i + j) / 2, y = (j - i) / 2 for row i, column j (positions in bp)
    Xr, Yr = (X + Y) / 2 / 1e6, (X - Y) / 2 / 1e6
    m = np.ma.masked_invalid(np.triu(mat).astype(float))
    m[np.tril_indices(n, -1)] = np.ma.masked
    ax.pcolormesh(Xr, Yr, m, cmap=cmap, norm=norm, shading="flat", rasterized=True)
    ax.set_ylim(0, depth / 2 / 1e6)
    ax.set_xlim(edges[0] / 1e6, edges[-1] / 1e6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


# ---- fig03: territories (genome-wide 1 Mb) -------------------------------------
g = cooler.Cooler(f"{D}/gm12878_hic_genome_1mb.cool")
pix = g.pixels(join=True)[:]          # upper triangle: each contact counted once
chroms = list(g.chromnames)
cis = {}
for c in chroms:
    inv = pix[(pix.chrom1 == c) | (pix.chrom2 == c)]
    cis[c] = inv["count"][inv.chrom1 == inv.chrom2].sum() / inv["count"].sum()
cis = pd.Series(cis)
size = pd.Series(g.chromsizes)[chroms].astype(float)
cis_all = pix["count"][pix.chrom1 == pix.chrom2].sum() / pix["count"].sum()
print("fig03: fraction of the contacts involving each chromosome that stay on it (cis):")
print(cis.round(3).to_string())
print(f"  {cis_all:.3f} of all contacts are within a chromosome; total contacts {pix['count'].sum():,}")
print(f"  expected if partners were picked at random (sum of squared sizes / squared genome): "
      f"{(size**2).sum() / size.sum()**2:.3f}")
bal = g.matrix(balance=True)[:]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.4), gridspec_kw={"width_ratios": [1.05, 1]})
im = a1.imshow(bal, cmap=HIC, norm=LogNorm(vmin=np.nanpercentile(bal[bal > 0], 20), vmax=np.nanpercentile(bal, 99.5)),
               interpolation="none")
ext = [g.extent(c) for c in chroms]
ticks = [(lo + hi) / 2 for lo, hi in ext]
show = ["chr1", "chr3", "chr5", "chr7", "chr9", "chr11", "chr13", "chr15", "chr18", "chr22", "chrX"]
a1.set_xticks([t for t, c in zip(ticks, chroms) if c in show]); a1.set_xticklabels([c[3:] for c in show], fontsize=11)
a1.set_yticks([t for t, c in zip(ticks, chroms) if c in show]); a1.set_yticklabels([c[3:] for c in show], fontsize=11)
for lo, hi in ext:
    a1.axhline(lo - 0.5, color="#d0d6dc", lw=0.4); a1.axvline(lo - 0.5, color="#d0d6dc", lw=0.4)
a1.set_title("GM12878 Hi-C, whole genome, 1 Mb bins")
a1.set_xlabel("chromosome"); a1.set_ylabel("chromosome")
for sp in a1.spines.values():
    sp.set_visible(True)
cb = fig.colorbar(im, ax=a1, fraction=0.04, pad=0.02); cb.set_label("balanced contacts (log)")
order = chroms
a2.bar(range(len(order)), cis[order] * 100, color=GREEN)
a2.set_xticks(range(len(order))); a2.set_xticklabels([c[3:] for c in order], fontsize=9.5)
a2.set_ylabel("% of its contacts within the chromosome")
a2.set_ylim(0, 100)
a2.set_title("Most contacts stay within a chromosome")
a2.set_xlabel("chromosome")
fig.tight_layout()
fig.savefig("fig03-territories.png", dpi=150)
plt.close(fig)

# ---- fig05: distance decay ----------------------------------------------------
c5 = cooler.Cooler(MC + "5000")
view = pd.DataFrame({"chrom": ["chr19", "chr21"], "start": [0, 0],
                     "end": [c5.chromsizes["chr19"], c5.chromsizes["chr21"]], "name": ["chr19", "chr21"]})
exp = cooltools.expected_cis(c5, view_df=view, smooth=False, aggregate_smoothed=False, ignore_diags=2, nproc=4)
exp["s"] = exp["dist"] * 5000
fig, ax = plt.subplots(figsize=(7.6, 5.0))
for chrom, col in (("chr19", "#2a6f97"), ("chr21", "#e67e22")):
    e = exp[(exp.region1 == chrom) & (exp.dist >= 2) & (exp["balanced.avg"] > 0) & (exp.dist * 5000 <= 2e7)]
    # log-binned for a clean curve
    edges = np.unique(np.round(np.logspace(np.log10(10000), np.log10(2e7), 40) / 5000).astype(int))
    idx = np.digitize(e.dist, edges)
    agg = e.groupby(idx).apply(lambda d: pd.Series({"s": np.exp(np.log(d.s).mean()),
                                                     "p": (d["balanced.avg"] * d.n_valid).sum() / d.n_valid.sum()}))
    ax.loglog(agg.s, agg.p, "-", color=col, lw=2.2, label=chrom)
    sel = (e.s >= 1e5) & (e.s <= 1e6)
    slope = np.polyfit(np.log10(e.s[sel]), np.log10(e["balanced.avg"][sel]), 1)[0]
    print(f"fig05: {chrom} P(s) slope between 100 kb and 1 Mb: {slope:.2f}")
    if chrom == "chr21":
        r = e[(e.s >= 9e4) & (e.s <= 1.1e5)]["balanced.avg"].mean() / e[(e.s >= 9e5) & (e.s <= 1.1e6)]["balanced.avg"].mean()
        print(f"fig05: chr21 contact frequency at 100 kb / at 1 Mb: {r:.1f}")
s_ref = np.array([1e5, 1e7])
p0 = exp[(exp.region1 == "chr21") & (exp.s == 100000)]["balanced.avg"].values[0]
ax.loglog(s_ref, p0 * 3 * (s_ref / 1e5) ** -1, "--", color=GREY, lw=1.5)
ax.text(1.6e6, p0 * 3 * 16 ** -1 * 1.25, "slope −1", color=GREY, fontsize=13)
ax.set_xlabel("genomic distance s (bp)")
ax.set_ylabel("contact probability P(s) (balanced)")
ax.set_title("Contacts fall with genomic distance")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("fig05-distance-decay.png", dpi=150)
plt.close(fig)

# ---- fig06: compartments on chr19 with marks ---------------------------------------
res = 250000
c = cooler.Cooler(MC + str(res))
chrom = "chr19"
L = c.chromsizes[chrom]
m = c.matrix(balance=True).fetch(chrom)
v19 = pd.DataFrame({"chrom": [chrom], "start": [0], "end": [L], "name": [chrom]})
c100 = cooler.Cooler(MC + "100000")
ex100 = cooltools.expected_cis(c100, view_df=v19, nproc=1)
# E1 at 100 kb, oriented by H3K27ac (more on choosing the sign in "Calling compartments, TADs and loops")
eig = cooltools.eigs_cis(c100, view_df=v19, n_eigs=1)[1]
e1 = eig["E1"].values
nb = len(e1)
k27 = bw_bins(f"{D}/gm12878_h3k27ac_fc.bigWig", chrom, 0, nb * 100000, 100000)[:nb]
k9 = bw_bins(f"{D}/gm12878_h3k9me3_fc.bigWig", chrom, 0, nb * 100000, 100000)[:nb]
ok = np.isfinite(e1)
if np.corrcoef(e1[ok], k27[ok])[0, 1] < 0:
    e1 = -e1
print(f"fig06: chr19 E1 (100 kb) correlation with H3K27ac {np.corrcoef(e1[ok], k27[ok])[0,1]:.2f}, "
      f"with H3K9me3 {np.corrcoef(e1[ok], k9[ok])[0,1]:.2f}; A bins {np.mean(e1[ok] > 0)*100:.0f}%")
fig = plt.figure(figsize=(9.5, 8.6))
gs = fig.add_gridspec(4, 1, height_ratios=[6, 0.9, 0.9, 0.9], hspace=0.08)
ax = fig.add_subplot(gs[0])
ax.imshow(m, cmap=HIC, norm=LogNorm(vmin=np.nanpercentile(m[m > 0], 5), vmax=np.nanpercentile(m, 99)),
          extent=[0, L / 1e6, L / 1e6, 0], interpolation="none", aspect="auto")
ax.set_title(f"{chrom}, 250 kb bins: the checkerboard of A/B compartments")
ax.set_ylabel("Mb"); plt.setp(ax.get_xticklabels(), visible=False)
x = (np.arange(nb) + 0.5) * 0.1
tracks = [(e1, "E1", None), (k27, "H3K27ac", K27AC), (k9, "H3K9me3", K9)]
for i, (y, lab, col) in enumerate(tracks):
    a = fig.add_subplot(gs[i + 1], sharex=ax)
    if col is None:
        a.fill_between(x, 0, np.where(y > 0, y, 0), color=GREEN, step="mid", lw=0)
        a.fill_between(x, 0, np.where(y < 0, y, 0), color=CRIMSON, step="mid", lw=0)
        a.text(1.01, 0.8, "A", color=GREEN, transform=a.transAxes, ha="left", va="center", fontweight="bold")
        a.text(1.01, 0.2, "B", color=CRIMSON, transform=a.transAxes, ha="left", va="center", fontweight="bold")
    else:
        a.fill_between(x, 0, y, color=col, step="mid", lw=0)
        a.set_ylim(0, np.nanpercentile(y, 99.5))
    a.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=12)
    a.set_yticks([])
    a.set_xlim(0, L / 1e6)
    if i < 2:
        plt.setp(a.get_xticklabels(), visible=False)
a.set_xlabel(f"{chrom} position (Mb)")
fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.07)
fig.savefig("fig06-compartments.png", dpi=150)
plt.close(fig)

# ---- fig09: TADs (chr21, 25 kb) --------------------------------------------------------
res, chrom, s0, s1 = 25000, "chr21", 37_000_000, 42_000_000
c = cooler.Cooler(MC + str(res))
m = c.matrix(balance=True).fetch(f"{chrom}:{s0}-{s1}")
dom = pd.read_csv(f"{D}/gm12878_hic_domains_arrowhead.bedpe.gz", sep="\t", comment="#", header=None)
dom = dom[(dom[0] == chrom) & (dom[1] >= s0) & (dom[2] <= s1)]
alld = pd.read_csv(f"{D}/gm12878_hic_domains_arrowhead.bedpe.gz", sep="\t", comment="#", header=None)
print(f"fig09: {len(dom)} Arrowhead domains in {chrom}:{s0}-{s1}; {len(alld)} genome-wide, median size "
      f"{(alld[2] - alld[1]).median() / 1e3:.0f} kb")
fig, ax = plt.subplots(figsize=(11, 4.6))
depth = 2_000_000
triangle(ax, m, s0, res, depth, LogNorm(vmin=np.nanpercentile(m[m > 0], 30), vmax=np.nanpercentile(m, 99.5)))
for _, d in dom.iterrows():
    a, b = d[1] / 1e6, d[2] / 1e6
    ax.plot([a, (a + b) / 2, b], [0, (b - a) / 2, 0], color=INK, lw=1.1, ls="--")
ax.set_xlabel(f"{chrom} position (Mb)")
ax.set_title("Domains (TADs) in a 5-Mb window, 25 kb bins; dashed: Arrowhead calls")
ax.set_ylabel("distance →")
fig.tight_layout()
fig.savefig("fig09-tads.png", dpi=150)
plt.close(fig)

# ---- fig10: loops (chr21, 5 kb) + CTCF ----------------------------------------------------
res, chrom, s0, s1 = 10000, "chr21", 24_750_000, 26_250_000
c = cooler.Cooler(MC + str(res))
m = c.matrix(balance=True).fetch(f"{chrom}:{s0}-{s1}")
loops = pd.read_csv(f"{D}/gm12878_hic_loops_hiccups.bedpe.gz", sep="\t", comment="#", header=None)
lp = loops[(loops[0] == chrom) & (loops[1] >= s0) & (loops[5] <= s1)]
print(f"fig10: {len(lp)} HiCCUPS loops in {chrom}:{s0}-{s1}; {len(loops)} genome-wide, "
      f"{(loops[0] == 'chr21').sum()} on chr21; median loop size "
      f"{((loops[4] + loops[5]) / 2 - (loops[1] + loops[2]) / 2).median() / 1e3:.0f} kb")
fig = plt.figure(figsize=(11, 5.6))
gs = fig.add_gridspec(2, 1, height_ratios=[5, 1], hspace=0.05)
ax = fig.add_subplot(gs[0])
depth = 1_000_000
triangle(ax, m, s0, res, depth, LogNorm(vmin=np.nanpercentile(m[m > 0], 30), vmax=np.nanpercentile(m, 99.7)))
for _, l in lp.iterrows():
    xa, xb = (l[1] + l[2]) / 2, (l[4] + l[5]) / 2
    ax.plot((xa + xb) / 2 / 1e6, (xb - xa) / 2 / 1e6, "o", mfc="none", mec="#2a6f97", mew=1.3, ms=19)
ax.set_title("Loops: dots at the corners of domains (10 kb bins; circles: HiCCUPS calls)")
ax.set_ylabel("distance →")
plt.setp(ax.get_xticklabels(), visible=False)
a = fig.add_subplot(gs[1], sharex=ax)
nb = (s1 - s0) // 2000
ct = bw_bins(f"{D}/gm12878_ctcf_fc.bigWig", chrom, s0, s0 + nb * 2000, 2000)
a.fill_between(s0 / 1e6 + (np.arange(nb) + 0.5) * 0.002, 0, ct, color=CTCF, lw=0)
a.set_ylabel("CTCF", rotation=0, ha="right", va="center")
a.set_yticks([])
a.set_xlim(s0 / 1e6, s1 / 1e6)
a.set_xlabel(f"{chrom} position (Mb)")
fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.11)
fig.savefig("fig10-loops.png", dpi=150)
plt.close(fig)
print("done")
