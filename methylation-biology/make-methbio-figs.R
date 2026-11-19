# Real-data figures for methylation-biology/lec.qmd.
# Run from this folder:  Rscript make-methbio-figs.R   (about 1 minute, 2 GB RAM)
# Needs: data/download.sh methylation-biology  (ENCODE WGBS per-CpG calls for
#        GM12878 and K562 on chr21, UCSC CpG islands, GENCODE v46 GTF,
#        GM12878 H3K4me3 bigWig, GM12878 RNA-seq gene quantifications)
#        R packages data.table, GenomicRanges, rtracklayer, ggplot2, patchwork.
#
# Writes fig08, fig09, fig11, fig12, fig13, fig14 (*.png) and prints the numbers
# quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(ggplot2); library(patchwork)
})
set.seed(1)
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
col_gm <- "#2a6f97"; col_k <- "#a8324a"

# ---- read per-CpG calls, merge the two strands of each CpG ----------------------
# ENCODE bed9+: col 6 strand, col 10 reads covering the C, col 11 % methylated
# (an integer, rounded down), col 13 the genotype called in this sample.
# The C of a CpG on the - strand is 1 bp to the right
# of the C on the + strand; both get the + strand position (1-based).
read_cpg <- function(f) {
  x <- fread(file.path(D, f), header = FALSE, select = c(1, 2, 6, 10, 11, 13),
             col.names = c("chr", "start", "strand", "N", "pct", "call"))
  x <- x[call == "CG"]                                # drop CpGs where this genome has a SNP
  x[, M := round(N * pct / 100)]                      # exact for N < 50
  x[, pos := ifelse(strand == "+", start + 1L, start)]
  x[, .(N = sum(N), M = sum(M)), by = .(chr, pos)]
}
gm <- read_cpg("gm12878_wgbs_rep1_chr21.bed.gz")
k5 <- read_cpg("k562_wgbs_rep1_chr21.bed.gz")
cat("GM12878 rep1 chr21: CpGs", nrow(gm), " mean coverage", round(mean(gm$N), 1),
    " coverage >= 10:", sum(gm$N >= 10), "\n")
gm10 <- gm[N >= 10][, beta := M / N]
gr <- GRanges(gm10$chr, IRanges(gm10$pos, width = 2)); gr$beta <- gm10$beta

# ---- annotation: CGIs, promoters, gene bodies -----------------------------------
cgi <- fread(file.path(D, "hg38_cpgIslandExt.txt.gz"), header = FALSE)[V2 == "chr21"]
cgi <- GRanges(cgi$V2, IRanges(cgi$V3 + 1, cgi$V4))         # file is 0-based
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding" & seqnames(gtf) == "chr21"]
genes$id <- sub("\\..*", "", genes$gene_id)
prom <- promoters(genes, upstream = 1000, downstream = 1000)
body <- genes[width(genes) > 5000]
body <- resize(body, width(body) - 2000, fix = "end")       # gene body without the first 2 kb
cat("chr21: CGIs", length(cgi), " protein-coding genes", length(genes),
    " promoters overlapping a CGI", sum(overlapsAny(prom, cgi)), "\n")

gr$region <- "intergenic"
gr$region[overlapsAny(gr, genes)] <- "gene body"
gr$region[overlapsAny(gr, prom)] <- "promoter"
gr$cgi <- overlapsAny(gr, cgi)
tab <- as.data.table(mcols(gr))
print(tab[, .(n = .N, mean = round(mean(beta), 3), below20 = round(mean(beta < 0.2), 3),
              above80 = round(mean(beta > 0.8), 3)), by = .(cgi)])
print(tab[, .(n = .N, mean = round(mean(beta), 3), median = round(median(beta), 3)), by = .(region, cgi)])

# ---- fig08: the distribution of per-CpG methylation ------------------------------
k10 <- k5[N >= 10][, beta := M / N]
hd <- rbind(data.table(cell = "GM12878", beta = gm10$beta), data.table(cell = "K562", beta = k10$beta))
cat("fraction < 0.2 / > 0.8: GM12878", round(mean(gm10$beta < 0.2), 2), round(mean(gm10$beta > 0.8), 2),
    "  K562", round(mean(k10$beta < 0.2), 2), round(mean(k10$beta > 0.8), 2), "\n")
cat("mean methylation of CpGs (coverage >= 10): GM12878", round(mean(gm10$beta), 3), " K562", round(mean(k10$beta), 3), "\n")
h <- data.table(beta = gm10$beta, cgi = ifelse(gr$cgi, "in CpG islands", "outside CpG islands"))
p10 <- ggplot(h, aes(beta)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.05), fill = col_gm, colour = "white") +
  facet_wrap(~ cgi, scales = "free_y") +
  labs(x = "methylation level of a CpG (methylated reads / all reads)", y = "density",
       title = "GM12878, chr21: per-CpG methylation (CpGs with >= 10 reads)")
save(p10, "fig08-meth-distribution.png", h = 4.2)

# ---- fig09: methylation by genomic context -------------------------------------
tab[, ctx := factor(paste0(region, ifelse(cgi, "\nCGI", "\nnot CGI")),
                    levels = c("promoter\nCGI", "promoter\nnot CGI", "gene body\nCGI", "gene body\nnot CGI",
                               "intergenic\nCGI", "intergenic\nnot CGI"))]
nlab <- tab[, .(n = .N), by = ctx]
p11 <- ggplot(tab, aes(ctx, beta)) +
  geom_violin(fill = "#eaf1f7", colour = "#5b6b7a", scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  geom_text(data = nlab, aes(ctx, 1.07, label = paste0("n = ", formatC(n, big.mark = " ", format = "d"))),
            size = 4.2, colour = "#5b6b7a") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  labs(x = NULL, y = "CpG methylation", title = "GM12878, chr21: methylation by genomic context")
save(p11, "fig09-meth-context.png", h = 4.6)

# ---- fig11: promoter methylation vs H3K4me3 ------------------------------------
pm <- findOverlaps(prom, gr)
pmeth <- tapply(gr$beta[subjectHits(pm)], queryHits(pm), mean)
pn <- tapply(gr$beta[subjectHits(pm)], queryHits(pm), length)
ok <- as.integer(names(pmeth))[pn >= 5]
pr <- prom[ok]; pr$meth <- as.numeric(pmeth[as.character(ok)])
pr$cgi <- overlapsAny(pr, cgi)
k4 <- rtracklayer::summary(BigWigFile(file.path(D, "gm12878_h3k4me3_fc.bigWig")), pr, type = "mean")
pr$k4 <- sapply(k4, function(z) z$score)
cat("promoters with >= 5 covered CpGs:", length(pr), " Spearman(meth, H3K4me3):",
    round(cor(pr$meth, pr$k4, method = "spearman"), 2), "\n")
cat("  promoters with meth < 0.2:", sum(pr$meth < 0.2), " of which H3K4me3 FC > 5:", sum(pr$meth < 0.2 & pr$k4 > 5),
    "; meth > 0.6:", sum(pr$meth > 0.6), " of which H3K4me3 FC > 5:", sum(pr$meth > 0.6 & pr$k4 > 5), "\n")
d12 <- data.table(meth = pr$meth, k4 = pr$k4, cgi = ifelse(pr$cgi, "CGI promoter", "non-CGI promoter"))
p12 <- ggplot(d12, aes(meth, k4 + 1, colour = cgi)) + geom_point(size = 2.4, alpha = 0.8) +
  scale_y_log10() + scale_colour_manual(values = c(`CGI promoter` = "#c0392b", `non-CGI promoter` = "#5b6b7a"), name = NULL) +
  labs(x = "mean CpG methylation, TSS ± 1 kb", y = "H3K4me3 fold enrichment + 1",
       title = "GM12878, chr21: protein-coding promoters") +
  theme(legend.position = "top")
save(p12, "fig11-promoter-k4me3.png", w = 7.5, h = 5.2)

# ---- fig12: methylation along genes by expression (metagene) ---------------------
rna <- fread(file.path(D, "gm12878_rna_rep1_genes.tsv"))[grepl("^ENSG", gene_id)]
rna[, id := sub("\\..*", "", gene_id)]
g <- genes[width(genes) > 5000]
g$tpm <- rna$TPM[match(g$id, rna$id)]
g <- g[!is.na(g$tpm)]
g$expr <- cut(g$tpm, c(-Inf, 0.5, 10, Inf), labels = c("not expressed (TPM < 0.5)", "low (0.5-10)", "high (> 10)"))
print(table(g$expr))
# 20 bins upstream (5 kb), 40 bins in the body (scaled), 20 bins downstream (5 kb):
# bins -20..-1 upstream, 0..39 gene body, 40..59 downstream
meta <- rbindlist(lapply(seq_along(g), function(i) {
  gi <- g[i]; s <- start(gi); e <- end(gi); st <- as.character(strand(gi))
  ov <- subjectHits(findOverlaps(gi + 5000, gr))
  if (!length(ov)) return(NULL)
  p <- start(gr)[ov]
  rel <- if (st == "+") p - s else e - p                  # distance from TSS along the gene
  len <- e - s
  x <- ifelse(rel < 0, rel / 5000 * 20, ifelse(rel > len, 40 + (rel - len) / 5000 * 20, rel / len * 40))
  data.table(bin = pmin(floor(x), 59), beta = gr$beta[ov], expr = gi$expr)
}))
mm <- meta[bin >= -20, .(beta = mean(beta), n = .N), by = .(bin, expr)]
p13 <- ggplot(mm, aes(bin + 0.5, beta, colour = expr)) + geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("#5b6b7a", "#d4a017", "#2e8b57"), name = NULL) +
  scale_x_continuous(breaks = c(-20, 0, 40, 60), labels = c("-5 kb", "TSS", "TES", "+5 kb")) +
  geom_vline(xintercept = c(0, 40), linetype = 2, colour = "grey60") +
  labs(x = "gene body scaled to the same length", y = "mean CpG methylation",
       title = "GM12878, chr21: methylation along protein-coding genes (> 5 kb)") +
  theme(legend.position = "top")
save(p13, "fig12-metagene-expression.png", w = 10, h = 4.8)
cat("gene body (bins 8-40 of 40) mean methylation by expression:\n")
print(meta[bin >= 8 & bin < 40, .(beta = round(mean(beta), 3)), by = expr])
cat("TSS bin (0-1) mean methylation by expression:\n")
print(meta[bin %in% 0:1, .(beta = round(mean(beta), 3)), by = expr])

# ---- fig14: partially methylated domains along chr21 ------------------------------
bins <- rbind(gm[, .(meth = sum(M) / sum(N), n = .N), by = .(bin = pos %/% 1e5)][, cell := "GM12878"],
              k5[, .(meth = sum(M) / sum(N), n = .N), by = .(bin = pos %/% 1e5)][, cell := "K562"])
bins <- bins[n >= 50]
full <- CJ(bin = seq(min(bins$bin), max(bins$bin)), cell = c("GM12878", "K562"))
bins <- merge(full, bins, by = c("bin", "cell"), all.x = TRUE)   # NA breaks the line at gaps
p14 <- ggplot(bins, aes((bin + 0.5) / 10, meth, colour = cell)) + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c(GM12878 = col_gm, K562 = col_k), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "chr21 position (Mb)", y = "mean methylation, 100-kb bins",
       title = "Partially methylated domains: chr21 in GM12878 and K562") +
  theme(legend.position = "top")
save(p14, "fig14-pmd-chr21.png", h = 4.4)
w <- dcast(bins[!is.na(meth)], bin ~ cell, value.var = "meth")[!is.na(GM12878) & !is.na(K562)]
cat("100-kb bins:", nrow(w), " K562 < 0.1:", sum(w$K562 < 0.1), " their mean in GM12878:",
    round(mean(w$GM12878[w$K562 < 0.1]), 2), " other bins GM12878:", round(mean(w$GM12878[w$K562 >= 0.1]), 2), "\n")
cat("chr21 global methylation (all CpGs, read-weighted): GM12878", round(sum(gm$M) / sum(gm$N), 3),
    " K562", round(sum(k5$M) / sum(k5$N), 3), "\n")

# ---- fig13: one locus: DNMT3L, AIRE, PFKL ---------------------------------------
reg <- GRanges("chr21", IRanges(44245000, 44330000))
x <- gm[pos >= start(reg) & pos <= end(reg) & N >= 5][, beta := M / N]
cat("locus CpGs with >= 5 reads:", nrow(x), "\n")
k4 <- import(file.path(D, "gm12878_h3k4me3_fc.bigWig"), which = reg)
k4 <- data.table(start = start(k4), end = end(k4), score = k4$score)
cg <- as.data.table(subsetByOverlaps(cgi, reg))
ex <- gtf[gtf$type == "exon" & overlapsAny(gtf, reg) & gtf$gene_type == "protein_coding"]
ex <- unlist(reduce(split(ex, ex$gene_name)))            # exons of all transcripts, merged
ex$gene_name <- names(ex)
gl <- as.data.table(genes[overlapsAny(genes, reg)])[, .(gene_name, start, end, strand)]
gl[, `:=`(start = pmax(start, start(reg)), end = pmin(end, end(reg)))]
gl[, y := c(1, 2, 1)[seq_len(.N)]]                       # stack genes in two rows
exd <- as.data.table(ex)[, .(gene_name, start, end)][gl[, .(gene_name, y)], on = "gene_name"]
exd[, `:=`(start = pmax(start, start(reg)), end = pmin(end, end(reg)))]
lim <- c(start(reg), end(reg)) / 1e6
pg <- ggplot() +
  geom_segment(data = gl, aes(start / 1e6, y, xend = end / 1e6, yend = y), colour = "#1b2a3a") +
  geom_rect(data = exd, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = y - 0.25, ymax = y + 0.25), fill = "#1b2a3a") +
  geom_text(data = gl, aes((start + end) / 2e6, y + 0.62, label = paste0(gene_name, ifelse(strand == "+", " →", " ←"))),
            size = 5, fontface = "italic") +
  geom_rect(data = cg, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -0.3, ymax = 0.1), fill = "#2e8b57") +
  annotate("text", x = lim[1], y = -0.1, label = "CpG islands", hjust = 0, size = 4.2, colour = "#2e8b57", vjust = 0.5) +
  scale_x_continuous(limits = lim, expand = c(0, 0)) + scale_y_continuous(limits = c(-0.4, 2.9)) +
  theme_void() + labs(title = "GM12878, chr21:44.245-44.330 Mb")+
  theme(plot.title = element_text(face = "bold", size = 15))
pm2 <- ggplot(x, aes(pos / 1e6, beta)) + geom_point(size = 1.3, colour = "#1b2a3a", alpha = 0.6) +
  scale_x_continuous(limits = lim, expand = c(0, 0)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "CpG\nmethylation") + theme(axis.text.x = element_blank())
pk <- ggplot(k4) + geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = score), fill = "#c0392b", colour = NA) +
  scale_x_continuous(limits = lim, expand = c(0, 0)) +
  labs(x = "chr21 position (Mb)", y = "H3K4me3\nfold enrichment")
p15 <- pg / pm2 / pk + plot_layout(heights = c(1.1, 1.3, 0.9))
save(p15, "fig13-locus-aire-pfkl.png", w = 12, h = 6)
