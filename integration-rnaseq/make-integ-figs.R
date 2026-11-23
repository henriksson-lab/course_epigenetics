# Real-data figures for integration-rnaseq/lec.qmd ("Linking regulatory elements to genes").
# Run from this folder:  Rscript make-integ-figs.R   (about 4 minutes, 6 GB RAM)
# Needs (data/download.sh integration-rnaseq; data/prepare-integ.sh for the K562 table):
#   crispr_benchmark_k562.tsv.gz, crispr_benchmark_heldout.tsv.gz, k562_abc_elements.tsv.gz,
#   k562_rna_rep1_genes.tsv, gm12878_rna_rep1_genes.tsv, gm12878_atac_peaks.bed.gz,
#   gm12878_atac_fc.bigWig, gm12878_h3k27ac_fc.bigWig, gencode.v46.basic.annotation.gtf.gz,
#   hg38_blacklist.bed.gz, corces2016_atac_counts.txt.gz, corces2016_rna_counts.txt.gz,
#   hg19ToHg38.over.chain.gz
# R packages: data.table, GenomicRanges, rtracklayer, ggplot2, patchwork, ggrepel.
#
# Writes fig03, fig04, fig07, fig08, fig10, fig12, fig13, fig15, fig17, fig18 (*.png)
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(ggplot2); library(patchwork); library(ggrepel)
})
set.seed(1)
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
colPN <- c(`TRUE` = "#c0392b", `FALSE` = "#5b6b7a")
chroms <- paste0("chr", c(1:22, "X"))

pr_curve <- function(score, truth) {
  o <- order(score, decreasing = TRUE)
  tp <- cumsum(truth[o])
  # A threshold includes every pair with the same score: never split a tie.
  last <- !duplicated(score[o], fromLast = TRUE)
  k <- which(last)
  data.frame(threshold = score[o][k], recall = tp[k] / sum(truth), precision = tp[k] / k)
}
auprc <- function(score, truth) { p <- pr_curve(score, truth); sum(diff(c(0, p$recall)) * p$precision) }
prec_at_recall <- function(score, truth, r = 0.7) { p <- pr_curve(score, truth); p$precision[which(p$recall >= r)[1]] }

# ---- the K562 CRISPRi benchmark ----------------------------------------------------
cr <- fread(file.path(D, "crispr_benchmark_k562.tsv.gz"))
cat("pairs", nrow(cr), " positives", sum(cr$Regulated), " elements", uniqueN(paste(cr$chrom, cr$chromStart)),
    " genes", uniqueN(cr$measuredGeneSymbol), " fraction positive", round(mean(cr$Regulated), 4), "\n")
print(cr[, .(pairs = .N, pos = sum(Regulated)), by = Dataset])
print(cr[, .(median_kb = median(distanceToTSS) / 1000, q90_kb = quantile(distanceToTSS, 0.9) / 1000), by = Regulated])
cat("median effect of positives", median(cr$EffectSize[cr$Regulated]), "\n")

# fig03: distance of tested pairs, and the fraction that is positive by distance
bins <- c(0, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 3e6)
cr[, dbin := cut(distanceToTSS, bins, labels = c("<10", "10-20", "20-50", "50-100", "100-200", "200-500", "500-1000", ">1000"))]
fb <- cr[, .(n = .N, pos = sum(Regulated), frac = mean(Regulated)), by = dbin][order(dbin)]
print(fb)
p03a <- ggplot(cr, aes(distanceToTSS / 1000, fill = Regulated)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2) + scale_x_log10(breaks = c(1, 10, 100, 1000)) +
  facet_wrap(~ ifelse(Regulated, "positive pairs (471)", "negative pairs (9885)"), ncol = 1, scales = "free_y") +
  scale_fill_manual(values = colPN, guide = "none") +
  labs(x = "distance from element to the gene's TSS (kb)", y = "pairs", title = "Tested element–gene pairs")
p03b <- ggplot(fb, aes(dbin, frac)) + geom_col(fill = "#c0392b") +
  geom_text(aes(label = paste0(pos, "/", n)), vjust = -0.3, size = 4) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "distance (kb)", y = "pairs that are positive", title = "Positives are close") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save(p03a + p03b + plot_layout(widths = c(1, 1.1)), "fig03-crispr-distance.png", w = 12, h = 4.8)

# ---- nearest gene --------------------------------------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), feature.type = "gene")
pc  <- gtf[gtf$gene_type == "protein_coding" & seqnames(gtf) %in% chroms]
tss <- resize(pc, width = 1, fix = "start"); tss$ensg <- sub("\\..*", "", tss$gene_id)
rna <- fread(file.path(D, "k562_rna_rep1_genes.tsv"))
tss$tpm <- rna$TPM[match(tss$ensg, sub("\\..*", "", rna$gene_id))]
expr <- tss[!is.na(tss$tpm) & tss$tpm >= 1]
cat("protein-coding genes", length(tss), " expressed in K562 (TPM >= 1)", length(expr), "\n")
tested <- GRanges(cr$chrom, IRanges(cr$chromStart + 1, cr$chromEnd))
mid <- resize(tested, width = 1, fix = "center")
cr[, nearest_pc   := tss$gene_name[nearest(mid, tss)]  == measuredGeneSymbol]
cr[, nearest_expr := expr$gene_name[nearest(mid, expr)] == measuredGeneSymbol]
pr_point <- function(pred, truth) c(TP = sum(pred & truth), FP = sum(pred & !truth),
  precision = sum(pred & truth) / sum(pred), recall = sum(pred & truth) / sum(truth))
print(round(rbind(nearest_pc = pr_point(cr$nearest_pc, cr$Regulated),
                  nearest_expr = pr_point(cr$nearest_expr, cr$Regulated)), 3))

# fig04: rank of the true target among expressed genes, by TSS distance from the element
pos <- cr[Regulated == TRUE]
pmid <- resize(GRanges(pos$chrom, IRanges(pos$chromStart + 1, pos$chromEnd)), 1, fix = "center")
rk <- sapply(seq_len(nrow(pos)), function(i) {
  e <- expr[as.character(seqnames(expr)) == pos$chrom[i]]
  d <- abs(start(e) - start(pmid[i])); t <- e$gene_name == pos$measuredGeneSymbol[i]
  if (!any(t)) return(NA)                                   # target not expressed (TPM < 1)
  sum(d < min(d[t])) + 1 })                                 # rank of the target by TSS distance
cat("positives whose target has TPM < 1:", sum(is.na(rk)), "\n"); rk <- rk[!is.na(rk)]
rk_tab <- table(factor(pmin(rk, 6), 1:6, labels = c("1st", "2nd", "3rd", "4th", "5th", "6th or\nfurther")))
print(rk_tab); print(round(prop.table(rk_tab), 3))
p04 <- ggplot(data.frame(rk_tab), aes(Var1, Freq)) + geom_col(fill = c("#2e8b57", rep("#c0392b", 5))) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", Freq, 100 * Freq / sum(Freq))), vjust = -0.3, size = 4.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "the true target is the ... closest expressed gene to the enhancer", y = "positive pairs",
       title = sprintf("CRISPRi-validated pairs in K562 (%d with an expressed target)", length(rk)))
save(p04, "fig04-nearest-rank.png", w = 10, h = 4.6)

# ---- K562 candidate elements and ABC -------------------------------------------------
el  <- fread(file.path(D, "k562_abc_elements.tsv.gz"))
el[, `:=`(activity = sqrt(dnase * h3k27ac), activity_atac = sqrt(atac * h3k27ac), mid = (start + end) %/% 2)]
elg <- GRanges(el$chr, IRanges(el$start + 1, el$end))
cat("K562 elements", nrow(el), " promoters", sum(el$promoter), "\n")
print(round(cor(el[, .(atac, dnase, h3k27ac)], method = "spearman"), 2))
h  <- findOverlaps(tested, elg)
cat("tested elements overlapping a candidate", round(mean(seq_along(tested) %in% queryHits(h)), 3),
    " positives", round(mean(which(cr$Regulated) %in% queryHits(h)), 3), "\n")
best <- function(v) { s <- numeric(length(tested)); m <- tapply(v[subjectHits(h)], queryHits(h), max); s[as.integer(names(m))] <- m; s }
cr[, A := best(el$activity)]
genes <- cr[, .(chr = chrTSS[1], tss = startTSS[1] + 1), by = .(gene = measuredGeneSymbol)]
win <- GRanges(genes$chr, IRanges(pmax(genes$tss - 5e6, 1), genes$tss + 5e6))
hw  <- findOverlaps(win, elg)
abc <- data.table(gene = genes$gene[queryHits(hw)], e = subjectHits(hw))
abc[, d := abs(el$mid[e] - genes$tss[queryHits(hw)])]
abc[, C := (pmax(d, 5000) + 1)^-1.024]
abc[, ABC := el$activity[e] * C / sum(el$activity[e] * C), by = gene]
abc[, ABC_atac := el$activity_atac[e] * C / sum(el$activity_atac[e] * C), by = gene]
hm <- merge(data.table(pair = queryHits(h), e = subjectHits(h), gene = cr$measuredGeneSymbol[queryHits(h)]),
            abc[, .(gene, e, ABC, ABC_atac)], by = c("gene", "e"))
s <- hm[, .(ABC = max(ABC), ABC_atac = max(ABC_atac)), by = pair]
cr[, `:=`(ABC = 0, ABC_atac = 0)]; cr[s$pair, `:=`(ABC = s$ABC, ABC_atac = s$ABC_atac)]
cr[, inv_dist := 1 / distanceToTSS]

# ---- correlation across 13 blood cell types (Corces 2016) -----------------------------
atac <- fread(file.path(D, "corces2016_atac_counts.txt.gz"))
code2type <- c(`1` = "HSC", `2` = "MPP", `3` = "LMPP", `4` = "CMP", `5` = "GMP", `6` = "MEP", `7` = "Mono",
               `9` = "CD4", `10` = "CD8", `11` = "NK", `13` = "Bcell", `14` = "CLP", `15` = "Ery")
healthy <- grep("^(Donor)?[0-9]+-", names(atac), value = TRUE)
re <- "^(?:Donor)?([0-9]+)-(?:[A-Za-z][A-Za-z0-9]*-)?([0-9]+)([A-C]?)(?:-[0-9]+)?$"
ct <- code2type[sapply(regmatches(healthy, regexec(re, healthy, perl = TRUE)), `[`, 3)]
cpm_by_type <- function(m, type) sapply(split(seq_along(type), type), function(i)
  log2(rowMeans(sweep(m[, i, drop = FALSE], 2, colSums(m[, i, drop = FALSE]), "/")) * 1e6 + 1))
Acor <- cpm_by_type(as.matrix(atac[, ..healthy]), ct)
chainf <- tempfile(fileext = ".chain")
writeLines(readLines(gzfile(file.path(D, "hg19ToHg38.over.chain.gz"))), chainf)
lo <- liftOver(GRanges(atac$Chr, IRanges(atac$Start + 1, atac$End)), import.chain(chainf))
one <- lengths(lo) == 1; peaks38 <- unlist(lo[one]); Acor <- Acor[one, ]
rnac <- fread(file.path(D, "corces2016_rna_counts.txt.gz"))
rcol <- grep("^[0-9]+-", names(rnac), value = TRUE)
rtype <- sub("^NKcell$", "NK", sub("Tcell$", "", sub("^[0-9]+-", "", rcol)))
Rcor <- cpm_by_type(as.matrix(rnac[, ..rcol]), rtype); rownames(Rcor) <- rnac[[1]]; Rcor <- Rcor[, colnames(Acor)]
h2 <- findOverlaps(tested, peaks38)
g  <- match(cr$measuredGeneSymbol, rownames(Rcor))
rr <- suppressWarnings(mapply(function(q, s) if (is.na(g[q])) NA else cor(Acor[s, ], Rcor[g[q], ]), queryHits(h2), subjectHits(h2)))
best_r <- suppressWarnings(tapply(rr, queryHits(h2), max, na.rm = TRUE))
cr[, r := NA_real_]; cr[as.integer(names(best_r)), r := best_r]; cr[!is.finite(r), r := NA]
cr[, r_score := fifelse(is.na(r), -1, r)]
print(cr[, .(with_r = sum(!is.na(r)), median_r = median(r, na.rm = TRUE), above_0.8 = mean(r > 0.8, na.rm = TRUE)), by = Regulated])
nullr <- suppressWarnings(replicate(5000, cor(Acor[sample(nrow(Acor), 1), ], Rcor[sample(nrow(Rcor), 1), ])))
cat("random peak-gene pairs: quantiles of r\n"); print(round(quantile(nullr, c(0.5, 0.9, 0.95, 0.99), na.rm = TRUE), 2))

# fig07: an example pair and a random pair, across the 13 cell types
ex <- cr[Regulated == TRUE & !is.na(r)][order(-r)]
print(head(ex[, .(measuredGeneSymbol, chrom, chromStart, distanceToTSS, EffectSize, r)], 8))
pick <- ex[measuredGeneSymbol %in% c("GATA1", "HBG1", "HBG2", "TAL1", "LMO2", "MYB", "KLF1", "GFI1B", "NFE2", "HBE1")][1]
if (is.na(pick$measuredGeneSymbol)) pick <- ex[1]
pe <- subjectHits(h2)[queryHits(h2) == which(cr$chrom == pick$chrom & cr$chromStart == pick$chromStart & cr$measuredGeneSymbol == pick$measuredGeneSymbol)[1]]
pe <- pe[which.max(sapply(pe, function(s) cor(Acor[s, ], Rcor[pick$measuredGeneSymbol, ])))]
d1 <- data.frame(celltype = colnames(Acor), atac = Acor[pe, ], rna = Rcor[pick$measuredGeneSymbol, ])
# a random gene on another chromosome that correlates just as well with the same peak
other <- rownames(Rcor)[!rownames(Rcor) %in% tss$gene_name[as.character(seqnames(tss)) == pick$chrom] & rownames(Rcor) %in% tss$gene_name]
rc_other <- suppressWarnings(cor(Acor[pe, ], t(Rcor[other, ])))[1, ]
g2 <- names(which.max(rc_other))
cat("example:", pick$measuredGeneSymbol, " r =", round(cor(d1$atac, d1$rna), 2),
    "; best trans gene", g2, "r =", round(max(rc_other, na.rm = TRUE), 2),
    "; trans genes with r > 0.8:", sum(rc_other > 0.8, na.rm = TRUE), "of", sum(!is.na(rc_other)), "\n")
d2 <- data.frame(celltype = colnames(Acor), atac = Acor[pe, ], rna = Rcor[g2, ])
sc <- function(d, gname, ttl) ggplot(d, aes(atac, rna)) + geom_point(size = 3, colour = "#6b4fa0") +
  geom_text_repel(aes(label = celltype), size = 4, seed = 1) +
  labs(x = "ATAC at the element (log2 CPM)", y = paste0(gname, " RNA (log2 CPM)"),
       title = sprintf("%s  (r = %.2f)", ttl, cor(d$atac, d$rna)))
save(sc(d1, pick$measuredGeneSymbol, paste0("CRISPRi target: ", pick$measuredGeneSymbol, ", ", round(pick$distanceToTSS / 1000), " kb")) +
     sc(d2, g2, paste0("another chromosome: ", g2)), "fig07-correlation-example.png", w = 12, h = 4.8)

# fig08: correlation for positive and negative pairs vs random pairs
dd <- rbind(data.frame(r = cr$r[cr$Regulated & !is.na(cr$r)], set = "CRISPRi positive pairs"),
            data.frame(r = cr$r[!cr$Regulated & !is.na(cr$r)], set = "CRISPRi negative pairs"),
            data.frame(r = nullr[!is.na(nullr)], set = "random peak + random gene"))
dd$set <- factor(dd$set, levels = unique(dd$set))
p08 <- ggplot(dd, aes(r, colour = set)) + geom_density(linewidth = 1.2, adjust = 0.8) +
  scale_colour_manual(values = c("#c0392b", "#5b6b7a", "#d4a017"), name = NULL) +
  labs(x = "Pearson r across 13 cell types (accessibility vs expression)", y = "density",
       title = "Correlation of accessibility and expression") + theme(legend.position = "right")
save(p08, "fig08-correlation-null.png", w = 12, h = 4.4)

# ---- fig10: the power-law contact used by ABC ------------------------------------------
dd <- data.frame(d = 10^seq(3, 6.7, length.out = 200))
dd$C <- (pmax(dd$d, 5000) + 1)^-1.024
dd$C <- dd$C / dd$C[which.min(abs(dd$d - 5000))]
p10 <- ggplot(dd, aes(d / 1000, C)) + geom_line(linewidth = 1.3, colour = "#2a6f97") +
  scale_x_log10(breaks = c(1, 5, 10, 100, 1000, 5000), labels = c("1", "5", "10", "100", "1000", "5000")) +
  scale_y_log10(breaks = c(1, 0.1, 0.01, 0.001), labels = c("1", "0.1", "0.01", "0.001")) +
  geom_vline(xintercept = c(5, 5000), linetype = 2, colour = "#5b6b7a") +
  annotate("text", x = 5.4, y = 0.002, label = "flat below 5 kb", hjust = 0, size = 4.5, colour = "#5b6b7a") +
  annotate("text", x = 4800, y = 0.5, label = "ABC stops at 5 Mb", hjust = 1, size = 4.5, colour = "#5b6b7a") +
  annotate("text", x = 120, y = 0.12, label = "slope ≈ −1: 10× farther,\n10× fewer contacts", hjust = 0, size = 4.5, colour = "#2a6f97") +
  labs(x = "distance from the TSS (kb)", y = "relative contact", title = "ABC's power law: contact ≈ distance to the power −1.02")
save(p10, "fig10-powerlaw.png", w = 9, h = 4.6)

# ---- GM12878: candidate elements, activity, ABC at the PPIF locus ------------------------
prom <- resize(resize(pc, 1, fix = "start"), 500, fix = "center")
p  <- fread(file.path(D, "gm12878_atac_peaks.bed.gz"))[V1 %in% chroms]
pk <- GRanges(p$V1, IRanges(p$V2 + p$V10 + 1 - 250, width = 500))
gm <- reduce(c(granges(pk), granges(prom)), ignore.strand = TRUE)
gm <- gm[!overlapsAny(gm, import(file.path(D, "hg38_blacklist.bed.gz")))]
gm$promoter <- overlapsAny(gm, prom)
signal <- function(bw, gr) {
  cv <- import(bw, which = reduce(gr), as = "RleList"); s <- numeric(length(gr))
  for (ch in unique(as.character(seqnames(gr)))) { i <- which(as.character(seqnames(gr)) == ch)
    s[i] <- viewMeans(Views(cv[[ch]], ranges(gr)[i])) }
  s }
gm$atac    <- signal(file.path(D, "gm12878_atac_fc.bigWig"), gm) * width(gm) / 500
gm$h3k27ac <- signal(file.path(D, "gm12878_h3k27ac_fc.bigWig"), resize(gm, width(gm) + 500, fix = "center")) * (width(gm) + 500) / 1000
gm$activity <- sqrt(gm$atac * gm$h3k27ac)
cat("GM12878 elements", length(gm), " promoters", sum(gm$promoter), "\n")
cat("Spearman ATAC vs H3K27ac (GM12878 elements):", round(cor(gm$atac, gm$h3k27ac, method = "spearman"), 2), "\n")

# fig12: ATAC vs H3K27ac in GM12878 elements
dg <- data.frame(atac = gm$atac, k27 = gm$h3k27ac, a = gm$activity, type = ifelse(gm$promoter, "promoter", "distal"))
cat("elements with zero activity:", sum(dg$a == 0), " (promoters:", sum(dg$a == 0 & dg$type == "promoter"), ")\n")
cat("activity quantiles (non-zero):\n"); print(round(quantile(dg$a[dg$a > 0], c(0.05, 0.5, 0.95)), 2))
dg <- dg[dg$atac > 0.05 & dg$k27 > 0.05, ]
p12 <- ggplot(dg[sample(nrow(dg), 30000), ], aes(atac, k27)) + geom_point(aes(colour = type), size = 0.4, alpha = 0.35) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = c(distal = "#2e8b57", promoter = "#c0392b"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "ATAC (fold change, summed over the element)", y = "H3K27ac (fold change, summed)",
       title = sprintf("GM12878: %s candidate elements", format(length(gm), big.mark = " ")))
p12b <- ggplot(data.frame(a = gm$activity, type = ifelse(gm$promoter, "promoter", "distal"))[gm$activity > 0, ], aes(a, fill = type)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.6) + scale_x_log10() +
  scale_fill_manual(values = c(distal = "#2e8b57", promoter = "#c0392b"), name = NULL) +
  labs(x = "activity A = √(ATAC × H3K27ac)", y = "elements", title = "Promoters: on or off")
save(p12 + p12b + plot_layout(widths = c(1.15, 1)), "fig12-activity.png", w = 12, h = 4.8)

ho <- fread(file.path(D, "crispr_benchmark_heldout.tsv.gz"))[CellType == "GM12878"]
gh <- ho[, .(chr = chrTSS[1], tss = startTSS[1] + 1), by = .(gene = measuredGeneSymbol)]
gmid <- (start(gm) + end(gm)) %/% 2
abc_gene <- function(chr, tss) {
  j <- which(as.character(seqnames(gm)) == chr & abs(gmid - tss) <= 5e6)
  ac <- gm$activity[j] * (pmax(abs(gmid[j] - tss), 5000) + 1)^-1.024
  data.table(e = j, dist = abs(gmid[j] - tss), ABC = ac / sum(ac)) }
gabc <- gh[, abc_gene(chr, tss), by = gene]
hh <- findOverlaps(GRanges(ho$chrom, IRanges(ho$chromStart + 1, ho$chromEnd)), gm)
hm2 <- merge(data.table(pair = queryHits(hh), e = subjectHits(hh), gene = ho$measuredGeneSymbol[queryHits(hh)]), gabc, by = c("gene", "e"))
ho[, ABC := 0]; s2 <- hm2[, .(ABC = max(ABC)), by = pair]; ho[s2$pair, ABC := s2$ABC]
cat("GM12878 held-out pairs:", nrow(ho), " positives", sum(ho$Regulated), "\n")
print(table(predicted = ho$ABC >= 0.02, CRISPRi = ho$Regulated))

# fig13: the PPIF locus
reg <- GRanges("chr10", IRanges(78.93e6, 79.40e6))
trk <- function(bw) { x <- import(file.path(D, bw), which = reg); x <- x[order(start(x))]
  data.frame(pos = (start(x) + end(x)) / 2, score = x$score) }
binned <- function(d, w = 1000) { d$b <- floor(d$pos / w) * w + w / 2; aggregate(score ~ b, d, max) }
ta <- binned(trk("gm12878_atac_fc.bigWig")); tk <- binned(trk("gm12878_h3k27ac_fc.bigWig"))
ppif <- gh[gene == "PPIF"]
pab <- gabc[gene == "PPIF"]; pab$start <- start(gm)[pab$e]; pab$end <- end(gm)[pab$e]
pab <- pab[start >= start(reg) & end <= end(reg)]
hp <- ho[measuredGeneSymbol == "PPIF"]
gl <- pc[overlapsAny(pc, reg)]
gdf <- data.frame(s = pmax(start(gl), start(reg)), e = pmin(end(gl), end(reg)), name = gl$gene_name, strand = as.character(strand(gl)))
gdf$y <- rep(c(0, -1), length.out = nrow(gdf))
xl <- c(start(reg), end(reg)) / 1e6
p13a <- ggplot(ta, aes(b / 1e6, score)) + geom_area(fill = "#2e8b57") + labs(y = "ATAC", x = NULL) +
  coord_cartesian(xlim = xl) + theme(axis.text.x = element_blank())
p13b <- ggplot(tk, aes(b / 1e6, score)) + geom_area(fill = "#e67e22") + labs(y = "H3K27ac", x = NULL) +
  coord_cartesian(xlim = xl) + theme(axis.text.x = element_blank())
p13c <- ggplot(pab, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = ABC)) + geom_rect(fill = "#2a6f97", colour = "#2a6f97") +
  geom_hline(yintercept = 0.02, linetype = 2) + scale_y_sqrt(breaks = c(0, 0.02, 0.1, 0.5)) +
  labs(y = "ABC", x = NULL) + coord_cartesian(xlim = xl) + theme(axis.text.x = element_blank())
p13d <- ggplot(hp, aes(x = (chromStart + chromEnd) / 2e6, y = EffectSize, colour = Regulated)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_point(size = 2.6) + scale_colour_manual(values = colPN, guide = "none") +
  labs(y = "CRISPRi", x = NULL) + coord_cartesian(xlim = xl) + theme(axis.text.x = element_blank())
p13e <- ggplot(gdf) + geom_segment(aes(x = s / 1e6, xend = e / 1e6, y = y, yend = y), linewidth = 3, colour = "#1b2a3a") +
  geom_text(aes(x = (s + e) / 2e6, y = y - 0.45, label = name), size = 4, fontface = "italic") +
  coord_cartesian(xlim = xl, ylim = c(-1.7, 0.4)) + labs(x = "chr10 position (Mb, GRCh38)", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
save(p13a / p13b / p13c / p13d / p13e + plot_layout(heights = c(1, 1, 1.1, 1.1, 0.9)), "fig13-abc-locus.png", w = 12, h = 6.6)
print(hp[, .(kb = round(distanceToTSS / 1000), effect = round(EffectSize, 3), Regulated, ABC = signif(ABC, 2))][order(kb)])

# ---- fig15: the benchmark composition ---------------------------------------------------
cmp <- cr[, .(pairs = .N, positives = sum(Regulated)), by = Dataset]
cmp[, label := c(Nasser2021 = "CRISPRi-FlowFISH\n(Fulco 2019 and others)", Gasperini2019 = "single-cell screen\n(Gasperini 2019)",
                 Schraivogel2020 = "targeted single-cell\n(Schraivogel 2020)")[Dataset]]
p15a <- ggplot(melt(cmp, id.vars = c("Dataset", "label")), aes(label, value, fill = variable)) +
  geom_col(position = "dodge") + geom_text(aes(label = value), position = position_dodge(0.9), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c(pairs = "#5b6b7a", positives = "#c0392b"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + labs(x = NULL, y = "element–gene pairs",
  title = "K562 benchmark: 10 356 pairs, 471 positive") + theme(legend.position = "top")
p15b <- ggplot(cr[Regulated == TRUE], aes(-EffectSize)) + geom_histogram(bins = 30, fill = "#c0392b") +
  scale_x_continuous(labels = scales::percent) + labs(x = "decrease in expression when the element is silenced",
  y = "positive pairs", title = "Most effects are modest")
save(p15a + p15b + plot_layout(widths = c(1.2, 1)), "fig15-benchmark.png", w = 12, h = 4.6)
cat("positives with < 20% decrease:", round(mean(cr$EffectSize[cr$Regulated] > -0.2), 2), "\n")

# ---- fig17: precision-recall curves -----------------------------------------------------
scores <- list(`1 / distance` = cr$inv_dist, `activity only` = cr$A, `correlation (13 cell types)` = cr$r_score,
               `ABC (ATAC × H3K27ac)` = cr$ABC_atac, `ABC (DNase × H3K27ac)` = cr$ABC)
res <- t(sapply(scores, function(v) c(AUPRC = auprc(v, cr$Regulated), precision_at_70pct_recall = prec_at_recall(v, cr$Regulated))))
print(round(res, 3))
thr <- cr[, .(called = sum(ABC >= 0.02), TP = sum(ABC >= 0.02 & Regulated))]
cat("ABC >= 0.02: called", thr$called, " TP", thr$TP, " precision", round(thr$TP / thr$called, 3),
    " recall", round(thr$TP / sum(cr$Regulated), 3), "\n")
pc_all <- rbindlist(lapply(names(scores), function(n) data.table(pr_curve(scores[[n]], cr$Regulated), score = n)))
pc_all[, score := factor(score, levels = names(scores))]
lab <- setNames(sprintf("%s  (AUPRC %.2f)", names(scores), res[, 1]), names(scores))
np <- pr_point(cr$nearest_pc, cr$Regulated)
p17 <- ggplot(pc_all, aes(recall, precision, colour = score)) + geom_path(linewidth = 1.1) +
  geom_hline(yintercept = mean(cr$Regulated), linetype = 2, colour = "grey50") +
  annotate("text", x = 0.3, y = mean(cr$Regulated) - 0.025, label = "random guess: 4.5% precision", hjust = 0, size = 4, colour = "grey40") +
  annotate("point", x = np["recall"], y = np["precision"], size = 4, shape = 18) +
  annotate("text", x = np["recall"] - 0.02, y = np["precision"] - 0.04, label = "nearest gene", hjust = 1, size = 4.5) +
  annotate("point", x = thr$TP / sum(cr$Regulated), y = thr$TP / thr$called, size = 4, shape = 21, fill = "#2e8b57") +
  annotate("text", x = thr$TP / sum(cr$Regulated) + 0.02, y = thr$TP / thr$called + 0.04, label = "ABC ≥ 0.02", hjust = 0, size = 4.5, colour = "#2e8b57") +
  scale_colour_manual(values = c("#5b6b7a", "#d4a017", "#6b4fa0", "#2a6f97", "#2e8b57"), labels = lab, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) + labs(title = "Which score finds the CRISPRi-validated pairs? (K562)")
save(p17, "fig17-pr-curves.png", w = 12, h = 5.2)

# ---- fig18: ABC is distance, modulated by activity -----------------------------------------
d18 <- cr[ABC > 0]
p18 <- ggplot(d18[order(Regulated)], aes(distanceToTSS / 1000, ABC, colour = Regulated)) +
  geom_point(aes(size = Regulated), alpha = 0.7) + scale_size_manual(values = c(`TRUE` = 1.8, `FALSE` = 0.6), guide = "none") +
  scale_x_log10() + scale_y_log10() + geom_hline(yintercept = 0.02, linetype = 2) +
  scale_colour_manual(values = colPN, labels = c(`TRUE` = "positive", `FALSE` = "negative"), name = NULL) +
  labs(x = "distance to TSS (kb)", y = "ABC score", title = "ABC falls with distance; activity moves pairs up or down")
cat("pairs with ABC > 0:", nrow(d18), "; positives with ABC = 0:", sum(cr$Regulated & cr$ABC == 0), "\n")
fn <- cr[Regulated == TRUE & ABC < 0.02]
cat("missed positives (ABC < 0.02):", nrow(fn), "; median distance", round(median(fn$distanceToTSS) / 1000), "kb",
    "; found positives median distance", round(median(cr[Regulated == TRUE & ABC >= 0.02]$distanceToTSS) / 1000), "kb\n")
save(p18, "fig18-abc-vs-distance.png", w = 10, h = 4.8)
