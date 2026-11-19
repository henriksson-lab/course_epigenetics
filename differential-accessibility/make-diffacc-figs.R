# Real-data figures for differential-accessibility/lec.qmd.
# Run from this folder:  Rscript make-diffacc-figs.R   (about 15 minutes, 4 GB RAM)
# Needs: data/download.sh differential-accessibility  (Corces 2016 ATAC and RNA
#        counts, hg19ToHg38 chain, GENCODE v46 GTF)
#        R packages data.table, GenomicRanges, rtracklayer, DESeq2, edgeR, apeglm,
#        ggplot2, patchwork, ggrepel.
#
# Writes fig05, fig08, fig10, fig11, fig12, fig15, fig17, fig18, fig19, fig20,
# fig22 (*.png) and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(DESeq2); library(edgeR); library(ggplot2); library(patchwork); library(ggrepel)
})
set.seed(1)
D <- "../data"
colBT <- c(Bcell = "#2a6f97", CD4 = "#e67e22")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")

# ---- load and decode sample names -----------------------------------------------
atac <- fread(file.path(D, "corces2016_atac_counts.txt.gz"))
cn <- names(atac)[-(1:3)]
healthy <- grep("^(Donor)?[0-9]+-", cn, value = TRUE)            # drop AML (SU...) and bulk CD34
re <- "^(?:Donor)?([0-9]+)-(?:[A-Za-z][A-Za-z0-9]*-)?([0-9]+)([A-C]?)(?:-[0-9]+)?$"
m  <- regmatches(healthy, regexec(re, healthy, perl = TRUE))
code2type <- c(`1` = "HSC", `2` = "MPP", `3` = "LMPP", `4` = "CMP", `5` = "GMP", `6` = "MEP",
               `7` = "Mono", `9` = "CD4", `10` = "CD8", `11` = "NK", `13` = "Bcell",
               `14` = "CLP", `15` = "Ery")
samples <- data.frame(sample = healthy, donor = sapply(m, `[`, 2),
                      celltype = unname(code2type[sapply(m, `[`, 3)]), rep = sapply(m, `[`, 4))
stopifnot(!anyNA(samples$celltype))
cat("healthy libraries:", nrow(samples), " donors:", length(unique(samples$donor)),
    " cell types:", length(unique(samples$celltype)), "\n")
lvl <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Mono", "Ery", "CLP", "Bcell", "CD4", "CD8", "NK")
samples$celltype <- factor(samples$celltype, levels = lvl)
ctcol <- c(HSC = "#1b2a3a", MPP = "#4d5d6c", LMPP = "#7f8c99", CMP = "#a8324a", GMP = "#d4a017",
           MEP = "#9c6644", Mono = "#e9c46a", Ery = "#c0392b", CLP = "#8ab17d",
           Bcell = "#2a6f97", CD4 = "#e67e22", CD8 = "#2e8b57", NK = "#6b4fa0")

# ---- fig05: the design: libraries per donor and cell type ------------------------
tab <- as.data.frame(table(celltype = samples$celltype, donor = samples$donor))
tab <- tab[tab$Freq > 0, ]
tab$main <- tab$celltype %in% c("Bcell", "CD4") & tab$donor %in% c("1022", "4983", "5483")
p05 <- ggplot(tab, aes(donor, celltype)) +
  geom_tile(aes(fill = main), colour = "white", linewidth = 1.5) +
  geom_text(aes(label = Freq), size = 5.5, fontface = "bold", colour = ifelse(tab$main, "white", "#1b2a3a")) +
  scale_fill_manual(values = c(`FALSE` = "#dfe6ec", `TRUE` = "#2e8b57"), guide = "none") +
  scale_y_discrete(limits = rev(lvl)) +
  labs(x = "donor", y = NULL, title = "ATAC-seq libraries per donor and cell type (77 healthy)",
       subtitle = "green: the B-cell vs CD4 T-cell comparison of this unit (3 donors, 8 libraries)")
save(p05, "fig05-corces-design.png", w = 9, h = 5.6)

# ---- liftOver hg19 -> hg38 -----------------------------------------------------
gr19 <- GRanges(atac$Chr, IRanges(atac$Start + 1, atac$End))   # file starts are 0-based
chainf <- tempfile(fileext = ".chain")
writeLines(readLines(gzfile(file.path(D, "hg19ToHg38.over.chain.gz"))), chainf)
lo <- liftOver(gr19, import.chain(chainf))
n  <- lengths(lo)
cat("liftOver: 0 pieces", sum(n == 0), " 1 piece", sum(n == 1), " >1 pieces", sum(n > 1), "\n")
keep <- which(n == 1)
gr38 <- unlist(lo[keep])
ok <- width(gr38) == 500 & as.character(seqnames(gr38)) == as.character(seqnames(gr19)[keep])
cat("  of the single pieces, not 500 bp or other chromosome:", sum(!ok), "\n")
keep <- keep[ok]; gr38 <- gr38[ok]
names(gr38) <- sprintf("peak%06d", keep)          # id = row number in the Corces file
cat("peaks kept in hg38:", length(gr38), "\n")

cnt_all <- as.matrix(atac[keep, samples$sample, with = FALSE])
rownames(cnt_all) <- names(gr38)
rm(atac); invisible(gc())

# ---- fig08: reads in peaks per library -------------------------------------------
samples$inpeaks <- colSums(cnt_all) / 1e6
o <- order(samples$celltype, samples$inpeaks)
samples$lab <- factor(samples$sample, levels = samples$sample[o])
p08 <- ggplot(samples, aes(lab, inpeaks, fill = celltype)) + geom_col(width = 0.8) +
  scale_fill_manual(values = ctcol, name = NULL) +
  labs(x = "77 libraries, grouped by cell type", y = "reads in peaks (millions)") +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "right") + guides(fill = guide_legend(ncol = 1))
save(p08, "fig08-reads-in-peaks.png", w = 12, h = 4.2)
cat("reads in peaks per library (M): min", round(min(samples$inpeaks), 2),
    "median", round(median(samples$inpeaks), 2), "max", round(max(samples$inpeaks), 2), "\n")

# ---- fig10: PCA of all 77 libraries -------------------------------------------
cd_all <- data.frame(samples[, c("donor", "celltype")], row.names = samples$sample)
sel_pk <- rowMeans(cnt_all) >= 10
dall <- DESeqDataSetFromMatrix(cnt_all[sel_pk, ], cd_all, ~ 1)
vall <- vst(dall, blind = TRUE)
top  <- order(-rowVars(assay(vall)))[1:10000]
pca  <- prcomp(t(assay(vall)[top, ]))
ve   <- round(100 * summary(pca)$importance[2, 1:2])
pd   <- data.frame(pca$x[, 1:2], cd_all)
cat("PCA on", sum(sel_pk), "peaks (mean >= 10), top 10000 variable; PC1", ve[1], "% PC2", ve[2], "%\n")
cent <- aggregate(cbind(PC1, PC2) ~ celltype, pd, median)
# label each cell type once, at its median; the other points get empty labels so labels avoid them
labd <- rbind(data.frame(PC1 = cent$PC1, PC2 = cent$PC2, celltype = cent$celltype, lab = as.character(cent$celltype)),
              data.frame(pd[, c("PC1", "PC2", "celltype")], lab = ""))
p10 <- ggplot(pd, aes(PC1, PC2, colour = celltype)) +
  geom_point(aes(shape = donor), size = 3, stroke = 1.2) +
  geom_text_repel(data = labd, aes(label = lab), size = 5, fontface = "bold", seed = 1,
                  box.padding = 0.6, point.padding = 0.3, max.overlaps = Inf, min.segment.length = 0.2,
                  show.legend = FALSE) +
  scale_colour_manual(values = ctcol, guide = "none") +
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 1), name = "donor") +
  labs(x = sprintf("PC1 (%d%%)", ve[1]), y = sprintf("PC2 (%d%%)", ve[2]))
save(p10, "fig10-pca-celltypes.png", w = 11, h = 5.4)

# ---- the B vs CD4 comparison: 8 libraries -> 6 columns --------------------------
sel <- samples[samples$celltype %in% c("Bcell", "CD4") & samples$donor %in% c("1022", "4983", "5483"), ]
cat("libraries used:\n"); print(sel[, c("sample", "donor", "celltype", "rep")], row.names = FALSE)
cnt <- cnt_all[, sel$sample]
cd  <- data.frame(donor = factor(sel$donor), celltype = factor(sel$celltype, levels = c("CD4", "Bcell")),
                  row.names = sel$sample)
cd$unit <- paste(cd$donor, cd$celltype, sep = "_")

# fig11: correlation of the 8 libraries (log2 CPM, peaks with mean CPM >= 2)
lcpm <- cpm(cnt, log = TRUE, prior.count = 2)
lcpm <- lcpm[rowMeans(cpm(cnt)) >= 2, ]
cc <- cor(lcpm)
short <- sub("^Donor", "", sub("-1[0-9]{5}$", "", colnames(cc)))
dimnames(cc) <- list(short, short)
ord <- short[order(sel$celltype, sel$donor)]
ccd <- as.data.frame(as.table(cc)); names(ccd) <- c("a", "b", "r")
ccd$a <- factor(ccd$a, levels = ord); ccd$b <- factor(ccd$b, levels = rev(ord))
p11 <- ggplot(ccd, aes(a, b, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4.2, colour = ifelse(ccd$r > 0.6, "white", "#1b2a3a")) +
  scale_fill_gradient(low = "#f4f6f8", high = "#1b2a3a", limits = c(min(cc), 1), name = "Pearson r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save(p11, "fig11-sample-correlation.png", w = 9, h = 6)
cat("peaks for correlation:", nrow(lcpm), " 5483 B A/B r =", round(cc["5483-Bcell-13A", "5483-Bcell-13B"], 3),
    " 5483 CD4 A/B r =", round(cc["5483-CD4-9A", "5483-CD4-9B"], 3), "\n")

dds0 <- DESeqDataSetFromMatrix(cnt, cd, ~ donor + celltype)
rowRanges(dds0) <- gr38
dds <- collapseReplicates(dds0, groupby = dds0$unit)
keep_pk <- rowSums(counts(dds) >= 10) >= 3
cat("peaks with >= 10 reads in >= 3 of 6 samples:", sum(keep_pk), "of", nrow(dds), "\n")
dds <- dds[keep_pk, ]
t0 <- Sys.time(); dds <- DESeq(dds, quiet = TRUE)
cat("DESeq2 run time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- fig12: median of ratios on real data ---------------------------------------
lc  <- log(counts(dds))
lgm <- rowMeans(lc)
fin <- is.finite(lgm)
s1  <- "1022_Bcell"
lr  <- (lc[fin, s1] - lgm[fin]) / log(2)
tot <- colSums(counts(dds))
sf_tot <- tot / exp(mean(log(tot)))
p12 <- ggplot(data.frame(lr = lr), aes(lr)) +
  geom_histogram(bins = 120, fill = "#9fb3c3") +
  geom_vline(xintercept = log2(sizeFactors(dds)[s1]), colour = "#2e8b57", linewidth = 1.3) +
  geom_vline(xintercept = log2(sf_tot[s1]), colour = "#c0392b", linewidth = 1.3, linetype = 2) +
  annotate("label", x = log2(sizeFactors(dds)[s1]) + 0.12, y = Inf, vjust = 1.4, hjust = 0, size = 5,
           label = sprintf("median ratio: size factor %.2f", sizeFactors(dds)[s1]),
           colour = "#2e8b57", fill = "white") +
  annotate("label", x = log2(sf_tot[s1]) - 0.12, y = Inf, vjust = 3.0, hjust = 1, size = 5,
           label = sprintf("total counts: %.2f", sf_tot[s1]), colour = "#c0392b", fill = "white") +
  coord_cartesian(xlim = c(-5, 3)) +
  labs(x = "log2( count in 1022_Bcell / geometric mean over the 6 samples ), one value per peak",
       y = "number of peaks")
save(p12, "fig12-median-ratio.png", w = 12, h = 4.4)

# ---- size factors: DESeq2, TMM, total; and the consequences ---------------------
y <- normLibSizes(DGEList(counts(dds)))
eff <- y$samples$lib.size * y$samples$norm.factors
sf <- rbind(DESeq2 = sizeFactors(dds), TMM = eff / exp(mean(log(eff))), total = sf_tot)
cat("size factors:\n"); print(round(sf, 2))
run_sf <- function(s) {
  d <- dds; sizeFactors(d) <- s; d <- DESeq(d, quiet = TRUE)
  r <- results(d, name = "celltype_Bcell_vs_CD4", alpha = 0.05)
  c(B_up = sum(r$padj < 0.05 & r$log2FoldChange > 1, na.rm = TRUE),
    T_up = sum(r$padj < 0.05 & r$log2FoldChange < -1, na.rm = TRUE),
    median_lfc = median(r$log2FoldChange, na.rm = TRUE))
}
nres <- t(sapply(rownames(sf), function(k) run_sf(sf[k, ])))
cat("DA peaks (padj < 0.05, |MLE LFC| > 1) by normalisation:\n"); print(round(nres, 3))
sfd <- data.frame(sample = rep(colnames(sf), each = 3), method = rep(rownames(sf), 6), sf = as.vector(sf))
sfd$sample <- factor(sfd$sample, levels = colnames(sf))
sfd$method <- factor(sfd$method, levels = c("DESeq2", "TMM", "total"))
mcol <- c(DESeq2 = "#2e8b57", TMM = "#2a6f97", total = "#c0392b")
p15a <- ggplot(sfd, aes(sample, sf, colour = method)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = mcol, name = NULL) +
  labs(x = NULL, y = "size factor", title = "Size factors") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
nd <- data.frame(method = rep(rownames(nres), 2), dir = rep(c("B > T", "T > B"), each = 3),
                 n = c(nres[, "B_up"], nres[, "T_up"]))
nd$method <- factor(nd$method, levels = c("DESeq2", "TMM", "total"))
p15b <- ggplot(nd, aes(method, n, fill = dir)) + geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1fk", n / 1000)), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 4.3) +
  scale_fill_manual(values = c(`B > T` = "#2a6f97", `T > B` = "#e67e22"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "normalisation", y = "differential peaks", title = "Differential peaks") +
  theme(legend.position = "top")
save(p15a + p15b + plot_layout(widths = c(1.3, 1)), "fig15-sizefactors.png", w = 12, h = 4.8)

# ---- results, shrinkage, unpaired design ---------------------------------------
res <- results(dds, name = "celltype_Bcell_vs_CD4", alpha = 0.05)
t0 <- Sys.time()
rs  <- lfcShrink(dds, coef = "celltype_Bcell_vs_CD4", type = "apeglm", quiet = TRUE)
cat("apeglm run time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cls <- ifelse(!is.na(res$padj) & res$padj < 0.05 & rs$log2FoldChange > 1, "B_up",
       ifelse(!is.na(res$padj) & res$padj < 0.05 & rs$log2FoldChange < -1, "T_up", "ns"))
cat("padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "; with |shrunken LFC| > 1: B_up",
    sum(cls == "B_up"), " T_up", sum(cls == "T_up"), "\n")
dd2 <- dds; design(dd2) <- ~ celltype; dd2 <- DESeq(dd2, quiet = TRUE)
r2 <- results(dd2, name = "celltype_Bcell_vs_CD4", alpha = 0.05)
cat("padj < 0.05: paired", sum(res$padj < 0.05, na.rm = TRUE), " unpaired", sum(r2$padj < 0.05, na.rm = TRUE),
    "; both", sum(res$padj < 0.05 & r2$padj < 0.05, na.rm = TRUE), "\n")
cat("median dispersion: paired", round(median(dispersions(dds)), 3), " unpaired", round(median(dispersions(dd2)), 3), "\n")
ddl <- DESeq(dds, test = "LRT", reduced = ~ celltype, quiet = TRUE)
cat("peaks with a donor effect (LRT padj < 0.05):", sum(results(ddl)$padj < 0.05, na.rm = TRUE), "\n")

# fig17: dispersion
dd <- data.frame(mean = mcols(dds)$baseMean, gw = mcols(dds)$dispGeneEst, fit = mcols(dds)$dispFit,
                 fin = dispersions(dds), fit2 = mcols(dd2)$dispFit)
sub <- dd[sample(nrow(dd), 20000), ]
p17 <- ggplot(sub, aes(mean)) +
  geom_point(aes(y = gw), colour = "grey55", size = 0.5, alpha = 0.4) +
  geom_point(aes(y = fin), colour = "#2a6f97", size = 0.5, alpha = 0.4) +
  geom_line(data = dd[order(dd$mean), ], aes(y = fit), colour = "#c0392b", linewidth = 1.2) +
  geom_line(data = dd[order(dd$mean), ], aes(y = fit2), colour = "#c0392b", linewidth = 1.2, linetype = 2) +
  annotate("text", x = 150, y = 5, label = "grey: per-peak estimate\nblue: shrunk (final)\nred: trend, paired (solid)\n       and unpaired (dashed)",
           hjust = 0, vjust = 1, size = 4.6) +
  scale_x_log10() + scale_y_log10() + coord_cartesian(ylim = c(0.01, 10)) +
  labs(x = "mean of normalised counts", y = "dispersion")
save(p17, "fig17-dispersion.png", w = 9, h = 4.8)

# fig18: MA plots, MLE vs apeglm
ma <- data.frame(mean = res$baseMean, mle = res$log2FoldChange, shr = rs$log2FoldChange, cls = cls)
ma <- ma[!is.na(ma$mle), ]
mid <- ma[sample(nrow(ma)), ]
cc3 <- c(B_up = "#2a6f97", T_up = "#e67e22", ns = "grey70")
pm <- function(v, ttl) ggplot(mid, aes(mean, .data[[v]], colour = cls)) +
  geom_point(size = 0.4, alpha = 0.5) + geom_hline(yintercept = 0) +
  scale_colour_manual(values = cc3, guide = "none") + scale_x_log10() +
  coord_cartesian(ylim = c(-9, 9)) + labs(x = "mean of normalised counts", y = "log2 fold change (B / CD4)", title = ttl)
save(pm("mle", "MLE (results)") + pm("shr", "shrunk (lfcShrink, apeglm)"), "fig18-ma-shrink.png", w = 12, h = 4.6)

# fig19: volcano
vo <- data.frame(lfc = rs$log2FoldChange, p = -log10(res$padj), cls = cls)
vo <- vo[!is.na(vo$p), ]
p19 <- ggplot(vo[sample(nrow(vo)), ], aes(lfc, p, colour = cls)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
  scale_colour_manual(values = cc3, guide = "none") +
  annotate("text", x = 10, y = 37, label = sprintf("open in B\n%s peaks", format(sum(cls == "B_up"), big.mark = " ")),
           colour = "#2a6f97", size = 5, fontface = "bold") +
  annotate("text", x = -10.5, y = 37, label = sprintf("open in CD4 T\n%s peaks", format(sum(cls == "T_up"), big.mark = " ")),
           colour = "#e67e22", size = 5, fontface = "bold") +
  coord_cartesian(xlim = c(-12.5, 12.5)) +
  labs(x = "shrunken log2 fold change (B / CD4)", y = "-log10 adjusted p")
save(p19, "fig19-volcano.png", w = 9, h = 4.8)

# ---- annotation: nearest protein-coding TSS -------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), feature.type = c("gene", "transcript"))
tx  <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding"]
tss <- unique(resize(tx, width = 1, fix = "start"))
mcols(tss) <- mcols(tss)["gene_name"]
pk  <- rowRanges(dds)
dn  <- distanceToNearest(pk, tss, ignore.strand = TRUE)
ann <- data.frame(peak = names(pk), gene = tss$gene_name[subjectHits(dn)], dist = mcols(dn)$distance,
                  cls = cls, lfc = rs$log2FoldChange, padj = res$padj)
cat("fraction within 1 kb of a protein-coding TSS: all tested", round(mean(ann$dist <= 1000), 3),
    " B_up", round(mean(ann$dist[cls == "B_up"] <= 1000), 3), " T_up", round(mean(ann$dist[cls == "T_up"] <= 1000), 3), "\n")
mk <- c("PAX5", "CD79A", "MS4A1", "CD19", "EBF1", "BLK", "CD3E", "BCL11B", "TCF7", "LEF1", "GATA3", "CD28")
cat("DA peaks whose nearest TSS is a marker gene:\n")
print(t(sapply(mk, function(g) table(factor(ann$cls[ann$gene == g], levels = c("B_up", "T_up", "ns"))))))

# fig20: browser-like bars at PAX5 and BCL11B
genes <- gtf[gtf$type == "gene"]
nc <- counts(dds, normalized = TRUE)
locus <- function(g, pad) {
  gg <- genes[genes$gene_name == g][1]
  win <- GRanges(seqnames(gg), IRanges(start(gg) - pad, end(gg) + pad))
  h <- which(overlapsAny(pk, win))
  d <- data.frame(pos = (start(pk)[h] + 249) / 1e6, cls = cls[h],
                  Bcell = rowMeans(nc[h, c("1022_Bcell", "4983_Bcell", "5483_Bcell")]),
                  CD4 = rowMeans(nc[h, c("1022_CD4", "4983_CD4", "5483_CD4")]))
  long <- rbind(data.frame(d[, 1:2], track = "B cells", v = d$Bcell), data.frame(d[, 1:2], track = "CD4 T", v = d$CD4))
  long$track <- factor(long$track, levels = c("B cells", "CD4 T"))
  ex <- gtf[gtf$type == "gene" & gtf$gene_name == g]
  gl <- data.frame(x0 = start(ex) / 1e6, x1 = end(ex) / 1e6, str = as.character(strand(ex)))
  ymax <- max(long$v)
  ggplot(long) +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = v, colour = cls), linewidth = 1.6) +
    scale_colour_manual(values = c(B_up = "#2a6f97", T_up = "#e67e22", ns = "grey70"), guide = "none") +
    scale_y_continuous(breaks = function(l) pretty(c(0, l[2]), n = 3)[pretty(c(0, l[2]), n = 3) >= 0],
                       expand = expansion(mult = c(0.12, 0.05))) +
    facet_grid(track ~ ., switch = "y") +
    geom_segment(data = gl, aes(x = x0, xend = x1, y = -0.12 * ymax, yend = -0.12 * ymax), linewidth = 2.5, colour = "#1b2a3a") +
    annotate("text", x = (gl$x0 + gl$x1) / 2, y = -0.27 * ymax, vjust = 0.5, label = paste0(g, " (", gl$str, ")"), fontface = "italic", size = 4.6) +
    labs(x = sprintf("%s position (Mb, hg38)", as.character(seqnames(gg))), y = "mean normalised count",
         title = g) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))
}
save(locus("PAX5", 60000) + locus("BCL11B", 60000), "fig20-marker-loci.png", w = 13, h = 5)

# ---- RNA-seq from the same study: concordance ----------------------------------
rna <- fread(file.path(D, "corces2016_rna_counts.txt.gz"))
rc  <- as.matrix(rna[, -1]); rownames(rc) <- rna[[1]]
rcol <- grep("-(Bcell|CD4Tcell)$", colnames(rc), value = TRUE)
rcd <- data.frame(donor = factor(sub("-.*", "", rcol)),
                  celltype = factor(ifelse(grepl("Bcell", rcol), "Bcell", "CD4"), levels = c("CD4", "Bcell")), row.names = rcol)
cat("RNA-seq samples:", paste(rcol, collapse = ", "), "\n")
ddr <- DESeqDataSetFromMatrix(rc[, rcol], rcd, ~ donor + celltype)
ddr <- DESeq(ddr[rowSums(counts(ddr) >= 10) >= 2, ], quiet = TRUE)
rr  <- results(ddr, name = "celltype_Bcell_vs_CD4")
cat("RNA LFC of markers:\n"); print(round(rr[intersect(mk, rownames(rr)), "log2FoldChange"], 1))
ann$rna <- rr$log2FoldChange[match(ann$gene, rownames(rr))]
ann$where <- ifelse(ann$dist <= 1000, "promoter (<= 1 kb)", "distal (> 1 kb)")
for (k in c("B_up", "T_up")) {
  a <- ann[ann$cls == k & !is.na(ann$rna), ]
  cat(k, ": nearest gene's RNA LFC median", round(median(a$rna), 2), "; same sign", round(mean(sign(a$rna) == sign(a$lfc)), 3),
      "(n =", nrow(a), ")\n")
}
cat("ns: median RNA LFC", round(median(ann$rna[ann$cls == "ns"], na.rm = TRUE), 2), "\n")
for (w in unique(ann$where)) {
  a <- ann[ann$cls != "ns" & ann$where == w & !is.na(ann$rna), ]
  cat(w, ": Spearman ATAC vs RNA LFC over DA peaks", round(cor(a$lfc, a$rna, method = "spearman"), 2), "(n =", nrow(a), ")\n")
}
bx <- ann[!is.na(ann$rna), ]
bx$cls <- factor(bx$cls, levels = c("T_up", "ns", "B_up"), labels = c("open in T", "not DA", "open in B"))
bx$where <- factor(bx$where, levels = c("promoter (<= 1 kb)", "distal (> 1 kb)"))
p22a <- ggplot(bx, aes(cls, rna, fill = cls)) + geom_boxplot(outlier.shape = NA, width = 0.6) +
  facet_wrap(~ where) + coord_cartesian(ylim = c(-8, 8)) + geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(`open in T` = "#e67e22", `not DA` = "grey80", `open in B` = "#2a6f97"), guide = "none") +
  labs(x = "ATAC peak class", y = "RNA log2 FC of nearest gene (B / CD4)")
pr <- ann[ann$cls != "ns" & ann$where == "promoter (<= 1 kb)" & !is.na(ann$rna), ]
lab <- pr[pr$gene %in% mk, ]; lab <- lab[!duplicated(lab$gene), ]
p22b <- ggplot(pr, aes(lfc, rna, colour = cls)) + geom_point(size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey40") + geom_vline(xintercept = 0, colour = "grey40") +
  geom_text_repel(data = lab, aes(label = gene), colour = "#1b2a3a", size = 4.3, fontface = "italic", seed = 1, min.segment.length = 0) +
  scale_colour_manual(values = cc3, guide = "none") +
  labs(x = "ATAC log2 FC (promoter DA peaks)", y = "RNA log2 FC of the gene",
       title = sprintf("Spearman %.2f", cor(pr$lfc, pr$rna, method = "spearman")))
save(p22a + p22b + plot_layout(widths = c(1.4, 1)), "fig22-rna-concordance.png", w = 12.5, h = 4.8)
cat("done\n")
