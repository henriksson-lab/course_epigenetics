# Real-data figures for methylation-differential/lec.qmd.
# Run from this folder:  Rscript make-diffmeth-figs.R   (about 4 minutes, 4 GB RAM)
# Needs: data/download.sh methylation-differential (ENCODE WGBS per-CpG calls for
#        GM12878 and K562, two replicates each, chr21; UCSC CpG islands; GENCODE
#        v46 GTF; ENCODE RNA-seq gene quantifications for both cell lines)
#        R packages data.table, GenomicRanges, rtracklayer, DSS (with bsseq),
#        ggplot2, patchwork.
#
# Writes fig03, fig04, fig05, fig06, fig07, fig09, fig10, fig11, fig12 (*.png)
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(DSS); library(ggplot2); library(patchwork)
})
set.seed(1)
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
col_gm <- "#2a6f97"; col_k <- "#a8324a"

# ---- read the four samples, merge strands ----------------------------------------------
read_cpg <- function(f) {
  x <- fread(file.path(D, f), header = FALSE, select = c(1, 2, 6, 10, 11, 13),
             col.names = c("chr", "start", "strand", "N", "pct", "call"))
  x <- x[call == "CG"][, M := round(N * pct / 100)]
  x[, pos := ifelse(strand == "+", start + 1L, start)]
  x[, .(N = sum(N), X = sum(M)), by = .(chr, pos)]
}
files <- c(gm1 = "gm12878_wgbs_rep1_chr21.bed.gz", gm2 = "gm12878_wgbs_rep2_chr21.bed.gz",
           k1 = "k562_wgbs_rep1_chr21.bed.gz", k2 = "k562_wgbs_rep2_chr21.bed.gz")
s <- lapply(files, read_cpg)
for (n in names(s)) cat(n, "CpGs", nrow(s[[n]]), " median coverage", median(s[[n]]$N),
                        " methylation", round(sum(s[[n]]$X) / sum(s[[n]]$N), 3), "\n")

# ---- fig03: binomial vs beta-binomial -----------------------------------------------------
dbb <- function(k, n, p, phi) {            # beta-binomial with mean p, overdispersion phi = 1/(a+b+1)
  if (phi == 0) return(dbinom(k, n, p))
  ab <- 1 / phi - 1; a <- p * ab; b <- (1 - p) * ab
  exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))
}
bb <- rbindlist(lapply(c(0, 0.05, 0.2), function(ph)
  data.table(k = 0:20, prob = dbb(0:20, 20, 0.5, ph), phi = sprintf("φ = %.2f%s", ph, ifelse(ph == 0, " (binomial)", "")))))
p03 <- ggplot(bb, aes(k, prob, fill = phi)) + geom_col(position = "dodge", width = 0.85) +
  scale_fill_manual(values = c("#1b2a3a", "#2a6f97", "#e67e22"), name = NULL) +
  labs(x = "methylated reads out of 20 (mean level 0.5)", y = "probability",
       title = "Same mean, more spread: the beta-binomial") + theme(legend.position = "top")
save(p03, "fig03-binom-betabinom.png", w = 10, h = 4.4)

# ---- fig04: overdispersion between replicates ------------------------------------------------
zrep <- function(a, b, minN = 10) {
  m <- merge(a, b, by = c("chr", "pos"))[N.x >= minN & N.y >= minN]
  p <- (m$X.x + m$X.y) / (m$N.x + m$N.y); ok <- p > 0 & p < 1
  m <- m[ok]; p <- p[ok]
  m[, z := (X.x / N.x - X.y / N.y) / sqrt(p * (1 - p) * (1 / N.x + 1 / N.y))][]
}
zg <- zrep(s$gm1, s$gm2); zk <- zrep(s$k1, s$k2)
cat("variance of z: GM12878 reps", round(var(zg$z), 2), "(", nrow(zg), "CpGs)  K562 reps", round(var(zk$z), 2),
    "(", nrow(zk), "CpGs)\n")
zz <- rbind(data.table(pair = sprintf("GM12878 rep1 vs rep2 (variance %.2f)", var(zg$z)), z = zg$z),
            data.table(pair = sprintf("K562 rep1 vs rep2 (variance %.2f)", var(zk$z)), z = zk$z))
nd <- data.table(z = seq(-6, 6, 0.05))[, d := dnorm(z)]
p04 <- ggplot(zz, aes(z)) + geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, fill = "#5b6b7a") +
  geom_line(data = nd, aes(z, d), colour = "#c0392b", linewidth = 1) + facet_wrap(~ pair) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(x = "z = (β1 − β2) / binomial standard error", y = "density",
       title = "Replicate differences vs pure sampling noise (red: N(0,1))")
save(p04, "fig04-overdispersion.png", h = 4.4)

# ---- fig05: replicate agreement, CpG level and 100 kb level -----------------------------------
m4 <- Reduce(function(a, b) merge(a, b, by = c("chr", "pos")),
             lapply(names(s), function(n) setnames(copy(s[[n]]), c("N", "X"), paste0(c("N", "X"), "_", n))))
cat("CpGs in all four samples:", nrow(m4), " with >= 10 reads in all:",
    nrow(m4[N_gm1 >= 10 & N_gm2 >= 10 & N_k1 >= 10 & N_k2 >= 10]), "\n")
b <- m4[, .(gm1 = sum(X_gm1) / sum(N_gm1), gm2 = sum(X_gm2) / sum(N_gm2), k1 = sum(X_k1) / sum(N_k1),
            k2 = sum(X_k2) / sum(N_k2), n = .N), by = .(bin = pos %/% 1e5)][n >= 50]
print(round(cor(b[, .(gm1, gm2, k1, k2)]), 3))
cat("100-kb bins: mean gm2 - gm1", round(mean(b$gm2 - b$gm1), 3), " k2 - k1", round(mean(b$k2 - b$k1), 3), "\n")
p05 <- ggplot(b, aes(gm1, gm2)) + geom_abline(colour = "#c0392b", linetype = 2) + geom_point(colour = col_gm, alpha = 0.7) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "GM12878 rep1", y = "GM12878 rep2", title = "GM12878 replicates") |
  ggplot(b, aes(k1, k2)) + geom_abline(colour = "#c0392b", linetype = 2) + geom_point(colour = col_k, alpha = 0.7) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "K562 rep1", y = "K562 rep2", title = "K562 replicates")
save(p05 + plot_annotation(caption = "each point: mean methylation of a 100-kb bin of chr21"), "fig05-replicate-bins.png", w = 10, h = 5)

# ---- DSS -------------------------------------------------------------------------------------------
bs <- makeBSseqData(lapply(s, as.data.frame), names(s))
cat("BSseq object:", nrow(bs), "CpGs x", ncol(bs), "samples\n")
t0 <- Sys.time()
dml <- DMLtest(bs, group1 = c("gm1", "gm2"), group2 = c("k1", "k2"), smoothing = TRUE, ncores = 4)
cat("DMLtest (smoothing):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
dml0 <- DMLtest(bs, group1 = c("gm1", "gm2"), group2 = c("k1", "k2"), smoothing = FALSE)
dml <- as.data.table(dml); dml0 <- as.data.table(dml0)
dmls <- callDML(dml, delta = 0.2, p.threshold = 1e-5)
dmrs <- as.data.table(callDMR(dml, delta = 0.2, p.threshold = 1e-5))
dmrs[, dir := ifelse(diff.Methy > 0, "lower in K562", "higher in K562")]
cat("DMLs (delta 0.2, p < 1e-5):", nrow(dmls), " lower in K562:", sum(dmls$diff > 0), "\n")
cat("DMRs:", nrow(dmrs), "\n"); print(dmrs[, .(n = .N, median_bp = median(length), median_CpG = median(nCG),
                                               total_Mb = sum(length) / 1e6), by = dir])

# ---- fig06: how many CpGs look significant, by test ------------------------------------------------
cnt <- m4[N_gm1 + N_gm2 >= 10 & N_k1 + N_k2 >= 10]
cnt[, fisher := mapply(function(a, b, c, d) fisher.test(matrix(c(a, b - a, c, d - c), 2))$p.value,
                       X_gm1 + X_gm2, N_gm1 + N_gm2, X_k1 + X_k2, N_k1 + N_k2)]
cnt[, fdr := p.adjust(fisher, "BH")]
tab <- data.table(test = c("Fisher, replicates pooled", "DSS, no smoothing", "DSS, smoothing"),
                  n = c(sum(cnt$fdr < 0.05), sum(dml0$fdr < 0.05, na.rm = TRUE), sum(dml$fdr < 0.05, na.rm = TRUE)),
                  tested = c(nrow(cnt), nrow(dml0), nrow(dml)))
print(tab)
# replicate-vs-replicate: the same tests on a null comparison (GM rep1 vs rep2)
nul <- cnt[, .(p = mapply(function(a, b, c, d) fisher.test(matrix(c(a, b - a, c, d - c), 2))$p.value,
                         X_gm1, N_gm1, X_gm2, N_gm2))]
cat("Fisher GM rep1 vs rep2 (no true difference): FDR < 0.05:", sum(p.adjust(nul$p, "BH") < 0.05), "of", nrow(nul), "\n")
tab$test <- factor(tab$test, levels = tab$test)
p06 <- ggplot(tab, aes(test, n / 1e3)) + geom_col(fill = c("#a8324a", "#5b6b7a", "#2a6f97"), width = 0.6) +
  geom_text(aes(label = formatC(n, big.mark = " ", format = "d")), vjust = -0.4, size = 5) +
  scale_y_continuous(limits = c(0, max(tab$n) / 1e3 * 1.15)) +
  labs(x = NULL, y = "thousands of CpGs at FDR < 0.05", title = "GM12878 vs K562, chr21: significant CpGs by test")
save(p06, "fig06-tests.png", w = 9, h = 4.4)

# ---- annotation ------------------------------------------------------------------------------------
cgi <- fread(file.path(D, "hg38_cpgIslandExt.txt.gz"), header = FALSE)[V2 == "chr21"]
cgi <- GRanges(cgi$V2, IRanges(cgi$V3 + 1, cgi$V4))
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding" & seqnames(gtf) == "chr21"]
genes$id <- sub("\\..*", "", genes$gene_id)
prom <- promoters(genes, 1000, 1000)
kbin <- s$k1[, .(m = sum(X) / sum(N), n = .N), by = .(bin = pos %/% 1e5)][n >= 50]
pmd <- reduce(GRanges("chr21", IRanges(kbin[m < 0.1]$bin * 1e5 + 1, width = 1e5)))
dg <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
dmrs[, `:=`(cgi = overlapsAny(dg, cgi), promoter = overlapsAny(dg, prom), k562_low = overlapsAny(dg, pmd))]
ann <- dmrs[, .(DMRs = .N, `CpG island` = mean(cgi), promoter = mean(promoter), `K562 low domain` = mean(k562_low)), by = dir]
print(ann)
# background: same overlaps for all tested CpGs
bg <- GRanges("chr21", IRanges(dml$pos, width = 1))
cat("background (tested CpGs): CGI", round(mean(overlapsAny(bg, cgi)), 3), " promoter", round(mean(overlapsAny(bg, prom)), 3),
    " K562 low", round(mean(overlapsAny(bg, pmd)), 3), "\n")
al <- melt(ann, id.vars = c("dir", "DMRs"), variable.name = "feature", value.name = "frac")
p11 <- ggplot(al, aes(feature, frac, fill = dir)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(`lower in K562` = col_gm, `higher in K562` = col_k), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "fraction of DMRs overlapping", title = "Where the DMRs are (chr21)") + theme(legend.position = "top")

# ---- fig10: DMR sizes ----------------------------------------------------------------------------
p10 <- ggplot(dmrs, aes(length, fill = dir)) + geom_histogram(bins = 40, position = "identity", alpha = 0.75) +
  scale_x_log10() + scale_fill_manual(values = c(`lower in K562` = col_gm, `higher in K562` = col_k), name = NULL) +
  labs(x = "DMR length (bp, log scale)", y = "DMRs", title = sprintf("%d DMRs: %d lower, %d higher in K562",
       nrow(dmrs), sum(dmrs$dir == "lower in K562"), sum(dmrs$dir == "higher in K562"))) +
  theme(legend.position = "top")
save(p10 | p11, "fig10-dmrs.png", w = 13, h = 4.8)

# ---- fig09: volcano of single-CpG tests ---------------------------------------------------------
set.seed(1)
v <- dml[!is.na(pval)][sample(.N, 60000)]
v[, cls := ifelse(fdr < 0.05 & abs(diff) > 0.2, ifelse(diff > 0, "lower in K562", "higher in K562"), "not called")]
p09 <- ggplot(v, aes(-diff, -log10(pmax(pval, 1e-50)), colour = cls)) + geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = c(`lower in K562` = col_gm, `higher in K562` = col_k, `not called` = "#c9d3dc"), name = NULL) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = 2) +
  labs(x = "Δβ (K562 − GM12878), smoothed", y = "−log10 p (capped at 50)", title = "DSS Wald test per CpG (60 000 of 344 000 shown)") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + theme(legend.position = "top")
save(p09, "fig09-volcano.png", w = 9, h = 5)

# ---- RNA ------------------------------------------------------------------------------------------
rna <- function(f) { x <- fread(file.path(D, f))[grepl("^ENSG", gene_id)]; setNames(x$TPM, sub("\\..*", "", x$gene_id)) }
r <- data.table(id = genes$id, gene = genes$gene_name,
                gm = (rna("gm12878_rna_rep1_genes.tsv")[genes$id] + rna("gm12878_rna_rep2_genes.tsv")[genes$id]) / 2,
                k5 = (rna("k562_rna_rep1_genes.tsv")[genes$id] + rna("k562_rna_rep2_genes.tsv")[genes$id]) / 2)
r[, lfc := log2((k5 + 1) / (gm + 1))]
po <- findOverlaps(prom, dg)
pr <- data.table(id = genes$id[queryHits(po)], d = -dmrs$diff.Methy[subjectHits(po)])[, .(d = d[which.max(abs(d))]), by = id]
r <- merge(r, pr, by = "id", all.x = TRUE)
r[, cls := ifelse(is.na(d), "no promoter DMR", ifelse(d > 0, "promoter higher in K562", "promoter lower in K562"))]
print(r[, .(genes = .N, median_lfc = round(median(lfc, na.rm = TRUE), 2), up = sum(lfc > 1, na.rm = TRUE),
            down = sum(lfc < -1, na.rm = TRUE)), by = cls])
r[, cgi := id %in% genes$id[overlapsAny(prom, cgi)]]
print(r[!is.na(d), .(n = .N, spearman = cor(d, lfc, method = "spearman", use = "complete")), by = cgi])
r$cls <- factor(r$cls, levels = c("promoter lower in K562", "no promoter DMR", "promoter higher in K562"))
nlab <- r[!is.na(lfc), .(n = .N), by = cls]
p13 <- ggplot(r[!is.na(lfc)], aes(cls, lfc, colour = cls)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.15, size = 1.6, alpha = 0.7) +
  geom_text(data = nlab, aes(cls, max(r$lfc, na.rm = TRUE) + 1, label = paste0("n = ", n)), colour = "#1b2a3a", size = 4.5) +
  scale_colour_manual(values = c(col_gm, "#5b6b7a", col_k), guide = "none") +
  labs(x = NULL, y = "log2 expression K562 / GM12878", title = "Promoter DMRs and expression (chr21 protein-coding genes)")
save(p13, "fig12-expression.png", w = 10, h = 5)

# ---- fig07 / fig11: loci with raw and smoothed methylation -----------------------------------------
locus_plot <- function(from, to, title, genes_show = TRUE) {
  raw <- rbindlist(lapply(names(s), function(n) s[[n]][pos >= from & pos <= to & N >= 5,
         .(pos, beta = X / N, sample = n, cell = ifelse(grepl("gm", n), "GM12878", "K562"))]))
  sm <- dml[pos >= from & pos <= to, .(pos, GM12878 = mu1, K562 = mu2)]
  sm <- melt(sm, id.vars = "pos", variable.name = "cell", value.name = "beta")
  dd <- dmrs[end >= from & start <= to][, `:=`(start = pmax(start, from), end = pmin(end, to))]
  p <- ggplot() +
    geom_rect(data = dd, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "#d4a017", alpha = 0.2) +
    geom_point(data = raw, aes(pos / 1e6, beta, colour = cell), size = 1, alpha = 0.45) +
    geom_line(data = sm, aes(pos / 1e6, beta, colour = cell), linewidth = 1.3) +
    scale_colour_manual(values = c(GM12878 = col_gm, K562 = col_k), name = NULL) +
    scale_x_continuous(limits = c(from, to) / 1e6, expand = c(0, 0)) +
    labs(x = "chr21 position (Mb)", y = "methylation", title = title) + theme(legend.position = "top")
  if (genes_show) {
    gg <- as.data.table(genes[overlapsAny(genes, GRanges("chr21", IRanges(from, to)))])
    gg[, `:=`(s = pmax(start, from), e = pmin(end, to))]
    p <- p + geom_segment(data = gg, aes(x = s / 1e6, xend = e / 1e6, y = -0.08, yend = -0.08), linewidth = 2, colour = "#1b2a3a") +
      geom_text(data = gg, aes(x = (s + e) / 2e6, y = -0.16, label = gene_name), fontface = "italic", size = 4.3)
  }
  p
}
# example: a CpG-island promoter DMR higher in K562 at a gene expressed in GM12878 (TPM > 2)
# with the largest drop in expression in K562
print(r[cls != "no promoter DMR" & !is.na(lfc)][order(lfc)][, .(gene, gm = round(gm, 1), k5 = round(k5, 1), lfc = round(lfc, 2), d = round(d, 2), cgi)])
ex <- r[cls == "promoter higher in K562" & cgi & gm > 2][order(lfc)][1]
g1 <- genes[genes$id == ex$id]
po1 <- dmrs[overlapsAny(dg, promoters(g1, 1000, 1000))]
cat("example gene:", ex$gene, " TPM GM12878", round(ex$gm, 1), " K562", round(ex$k5, 1), " promoter DMR:",
    po1$start, "-", po1$end, " beta GM", round(po1$meanMethy1, 2), " K562", round(po1$meanMethy2, 2), "\n")
tss <- ifelse(as.character(strand(g1)) == "+", start(g1), end(g1))
p12 <- locus_plot(tss - 5000, tss + 5000, sprintf("A promoter DMR: %s (TPM %.1f in GM12878, %.1f in K562)", ex$gene, ex$gm, ex$k5))
save(p12, "fig11-dmr-locus.png", w = 11, h = 5)
p07 <- locus_plot(36.214e6, 36.226e6, "Raw CpGs (dots) and DSS-smoothed means (lines)")
save(p07, "fig07-smoothing.png", w = 11, h = 4.8)
