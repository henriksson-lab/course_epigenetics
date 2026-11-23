# Real-data figures for motifs-footprinting/lec.qmd.
# Run from this folder:  Rscript make-motifs-figs.R
# Needs: data/download.sh motifs-footprinting differential-accessibility
#        (Corces 2016 ATAC counts, hg19ToHg38 chain, hg38.fa.bgz + index,
#        GM12878 ATAC reads rep1 + bowtie2 index of chr21/chrM, GM12878 ATAC and
#        CTCF peaks, hg38 blacklist); plus differential_peaks.tsv saved by
#        Differential accessibility R11 (DIFFERENTIAL_PEAKS overrides its path).
#        bowtie2, samtools and fastp on the PATH;
#        R packages SummarizedExperiment, rtracklayer, Rsamtools, GenomicAlignments,
#        JASPAR2024, TFBSTools, motifmatchr, chromVAR, universalmotif,
#        ggplot2, patchwork, ggrepel, pheatmap, BiocParallel.
#
# Writes fig05-family-logos.png, fig06-background-real.png, fig19-motif-vs-binding.png,
# fig07-enrichment-volcano.png, fig08-or-scatter.png, fig11-chromvar-heatmap.png,
# fig12-chromvar-variability.png, fig15-ctcf-footprint.png, fig16-tn5-bias.png
# and prints the numbers quoted in the slides. About 15 minutes on 8 cores.
# Set MOTIF_CACHE=/some/dir to keep intermediate results between runs, and
# ATAC_BAM=/path/to/gm12878_atac_rep1.bam to reuse the BAM made in the exercise.

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(rtracklayer); library(GenomicRanges); library(Rsamtools)
  library(GenomicAlignments); library(Biostrings); library(Matrix)
  library(JASPAR2024); library(TFBSTools); library(motifmatchr); library(chromVAR)
  library(universalmotif); library(BiocParallel)
  library(ggplot2); library(patchwork); library(ggrepel); library(pheatmap)
})
register(MulticoreParam(8))
D <- "../data"
set.seed(1)
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
colB <- "#2a6f97"; colT <- "#e67e22"
cache_dir <- Sys.getenv("MOTIF_CACHE", "")
if (nzchar(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cached <- function(name, expr) {
  if (cache_dir == "") return(expr)
  f <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr; saveRDS(x, f); x
}
genome <- FaFile(file.path(D, "hg38.fa.bgz"))

# ---- A. Reuse the saved Differential accessibility result --------------------
diff_file <- Sys.getenv("DIFFERENTIAL_PEAKS", file.path(D, "differential_peaks.tsv"))
if (!file.exists(diff_file)) stop("Run differential-accessibility/make-diffacc-figs.R first, or set DIFFERENTIAL_PEAKS to the table from exercise R11")
result <- read.delim(diff_file, check.names = FALSE)
stopifnot(all(c("peak_id", "chr", "start", "end", "log2FoldChange", "padj", "class") %in% names(result)),
          !anyDuplicated(result$peak_id), !anyNA(result$class),
          all(result$class %in% c("B_up", "T_up", "ns")))
tested <- GRanges(result$chr, IRanges(result$start, result$end))  # 1-based closed
names(tested) <- result$peak_id
tested$lfc <- result$log2FoldChange; tested$padj <- result$padj
stopifnot(all(width(tested) == 500))
b_up <- result$class == "B_up"
t_up <- result$class == "T_up"
ns <- result$class == "ns"
# Include the input checksum in caches so a changed peak table cannot reuse old matches.
input_key <- unname(tools::md5sum(diff_file))
cat("tested", length(tested), "B-up", sum(b_up), "T-up", sum(t_up), "ns", sum(ns), "\n")

# ---- B. JASPAR motifs, matches, enrichment -------------------------------------
jdb <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024()))
pfms <- getMatrixSet(jdb, opts = list(collection = "CORE", tax_group = "vertebrates",
                                      all_versions = FALSE))
names(pfms) <- paste(name(pfms), ID(pfms), sep = "_")
fam <- sapply(pfms, function(p) paste(tags(p)$family, collapse = "; "))
cat("motifs:", length(pfms), "\n")
hit <- cached(paste0("hit-", input_key), motifMatches(matchMotifs(pfms, tested, genome = genome, p.cutoff = 5e-5)))
cat("median fraction of peaks with a match:", median(colMeans(hit)), "\n")

enrich <- function(target, hit, bg) {
  stopifnot(!any(target & bg))
  a <- colSums(hit[target, , drop = FALSE]); b <- sum(target) - a
  c <- colSums(hit[bg, , drop = FALSE]);     d <- sum(bg) - c
  # one-sided Fisher = hypergeometric upper tail; log10 so that tiny p-values survive
  lp <- phyper(a - 1, a + c, b + d, a + b, lower.tail = FALSE, log.p = TRUE) / log(10)
  data.frame(motif = colnames(hit), tf = sub("_MA.*", "", colnames(hit)),
             family = fam[colnames(hit)], a, b, c, d,
             OR = (a * d) / (b * c), log10p = lp,
             padj = p.adjust(10^lp, "BH"))
}
eB <- enrich(b_up, hit, bg = ns); eT <- enrich(t_up, hit, bg = ns)
for (e in list(B = eB, T = eT)) {
  print(head(e[order(e$log10p), c("tf", "a", "OR", "log10p")], 12), row.names = FALSE, digits = 3)
  cat("significant (BH < 0.05, OR > 1):", sum(e$padj < 0.05 & e$OR > 1),
      " Bonferroni:", sum(10^e$log10p < 0.05 / nrow(e) & e$OR > 1), "\n")
}

# ---- fig05: motif families (JASPAR logos, aligned) -----------------------------
fams <- list(ETS = c("Spi1_MA0080.7", "SPIB_MA0081.3", "ETS1_MA0098.4", "GABPA_MA0062.4"),
             RUNX = c("Runx1_MA0002.3", "RUNX2_MA0511.2", "RUNX3_MA0684.3"),
             PAX = c("PAX5_MA0014.4", "PAX2_MA0067.3", "PAX8_MA2094.1"),
             "TCF/LEF" = c("TCF7_MA0769.3", "Lef1_MA0768.3", "TCF7L2_MA0523.2"))
stopifnot(all(unlist(fams) %in% names(pfms)))
logo <- lapply(names(fams), function(f) {
  m <- convert_motifs(pfms[fams[[f]]])
  for (i in seq_along(m)) m[[i]]@name <- sub("_MA.*", "", fams[[f]][i])
  view_motifs(m, show.positions = FALSE, names.pos = "top", text.size = 13) +
    ggtitle(paste(f, "family")) + theme(plot.title = element_text(face = "bold", size = 16))
})
ggsave("fig05-family-logos.png", wrap_plots(logo, nrow = 1), width = 15, height = 5.6, dpi = 110,
       bg = "white")

# ---- fig06: the same motifs against three backgrounds ---------------------------
gc_of <- function(gr) {
  s <- getSeq(genome, gr)
  as.numeric(letterFrequency(s, "GC") / width(s))
}
gc <- cached(paste0("gc-", input_key), gc_of(tested))
# random 500-bp windows on chr1-22, X without N and blacklisted regions
rnd <- cached(paste0("rnd-", input_key), {
  set.seed(1)
  sz <- read.table(file.path(D, "hg38.chrom.sizes"), col.names = c("chr", "len"))
  sz <- sz[sz$chr %in% paste0("chr", c(1:22, "X")), ]
  bl <- import(file.path(D, "hg38_blacklist.bed.gz"))
  n <- 3 * sum(b_up)
  chr <- sample(sz$chr, n, replace = TRUE, prob = sz$len)
  st <- sapply(chr, function(ch) sample.int(sz$len[sz$chr == ch] - 500, 1))
  r <- GRanges(chr, IRanges(st, width = 500))
  r <- r[!overlapsAny(r, bl)]
  s <- getSeq(genome, r)
  r <- r[letterFrequency(s, "N")[, 1] == 0]
  r[seq_len(sum(b_up))]
})
gc_rnd <- gc_of(rnd)
hit_rnd <- motifMatches(matchMotifs(pfms, rnd, genome = genome, p.cutoff = 5e-5))
# GC-matched: sample ns peaks without replacement within each 2% GC bin
set.seed(1)
bin <- cut(gc, seq(0, 1, 0.02), include.lowest = TRUE)
pool <- split(which(ns), bin[ns])
gcm <- unlist(lapply(split(which(b_up), bin[b_up]), function(i) {
  p <- pool[[as.character(bin[i[1]])]]
  if (length(p) == 0) integer(0) else p[sample.int(length(p), min(length(p), length(i)))]
}))
gcm_set <- seq_along(tested) %in% gcm
cat("GC: B-up", round(mean(gc[b_up]), 3), " ns", round(mean(gc[ns]), 3),
    " GC-matched", round(mean(gc[gcm]), 3), " random", round(mean(gc_rnd), 3), "\n")
or_vs <- function(hb) {          # OR of B-up vs a background match matrix hb
  a <- colSums(hit[b_up, ]); b <- sum(b_up) - a
  c <- colSums(hb); d <- nrow(hb) - c
  (a * d) / (b * c)
}
show <- c("SPIB_MA0081.3", "POU2F2_MA0507.3", "IRF4_MA1419.2", "EBF1_MA0154.5",
          "CTCF_MA0139.2", "KLF4_MA0039.5", "SP1_MA0079.5", "Nrf1_MA0506.3")
stopifnot(all(show %in% names(pfms)))
bgs <- list("ns peaks" = or_vs(hit[ns, ]),
            "GC-matched peaks" = or_vs(hit[gcm, ]),
            "random genome windows" = or_vs(hit_rnd))
d6 <- do.call(rbind, lapply(names(bgs), function(n)
  data.frame(bg = n, motif = sub("_MA.*", "", show), OR = bgs[[n]][show])))
d6$bg <- factor(d6$bg, levels = names(bgs))
d6$motif <- factor(d6$motif, levels = rev(sub("_MA.*", "", show)))
print(reshape(d6, idvar = "motif", timevar = "bg", direction = "wide"), digits = 3)
gcd <- rbind(data.frame(set = "B-up peaks", gc = gc[b_up]),
             data.frame(set = "ns peaks", gc = gc[ns]),
             data.frame(set = "random genome windows", gc = gc_rnd))
bgcol <- c("B-up peaks" = colB, "ns peaks" = "#5b6b7a", "GC-matched peaks" = "#2e8b57",
           "random genome windows" = "#a8324a")
p6a <- ggplot(gcd, aes(gc, colour = set)) + geom_density(linewidth = 1.2) +
  scale_colour_manual(values = bgcol, name = NULL) +
  coord_cartesian(xlim = c(0.15, 1)) + scale_x_continuous(breaks = seq(0.2, 1, 0.2)) +
  labs(x = "GC content of the 500-bp peak", y = "density", title = "GC content") +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1), legend.key.size = unit(0.9, "lines"))
p6b <- ggplot(d6, aes(OR, motif, colour = bg)) +
  geom_vline(xintercept = 1, colour = "grey40", linetype = 2) +
  geom_point(size = 4) + scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
  scale_colour_manual(values = bgcol, name = "B-up peaks compared with") +
  labs(x = "odds ratio (log scale)", y = NULL, title = "Enrichment depends on the background")
ggsave("fig06-background-real.png", p6a + p6b + plot_layout(widths = c(1.25, 1.5)),
       width = 14, height = 5.2, dpi = 110, bg = "white")

# ---- fig07: enrichment volcano, B-up and T-up ----------------------------------
famcol <- c("ETS" = "#2a6f97", "RUNX" = "#e67e22", "POU (OCT)" = "#6b4fa0", "IRF" = "#2e8b57",
            "TCF/LEF" = "#c0392b", "EBF" = "#d4a017", "CTCF" = "#8c564b", "other" = "#b8c2cc")
fam_short <- function(f, tf) {
  out <- rep("other", length(f))
  out[grepl("Ets-related", f) & !grepl("::", tf)] <- "ETS"
  out[grepl("Runt", f)] <- "RUNX"
  out[grepl("POU domain", f)] <- "POU (OCT)"
  out[grepl("Interferon", f)] <- "IRF"
  out[grepl("TCF-7", f)] <- "TCF/LEF"
  out[grepl("Early B-Cell", f)] <- "EBF"
  out[grepl("^CTCF", tf)] <- "CTCF"
  factor(out, levels = names(famcol))
}
vol <- function(e, title, col) {
  e$fs <- fam_short(e$family, e$tf)
  e$y <- -e$log10p
  o <- e[order(e$log10p), ]
  lab <- rbind(head(o[o$fs != "other", ], 8), head(o[o$fs == "other", ], 3),
               e[e$motif %in% c("EBF1_MA0154.5", "PAX5_MA0014.4", "CTCF_MA0139.2"), ])
  lab <- lab[!duplicated(lab$motif), ]
  ggplot(e, aes(log2(OR), y, colour = fs)) +
    geom_point(data = e[e$fs == "other", ], size = 1.6, alpha = 0.6) +
    geom_point(data = e[e$fs != "other", ], size = 2.6) +
    geom_hline(yintercept = -log10(0.05 / nrow(e)), linetype = 2, colour = "grey40") +
    geom_text_repel(data = lab, aes(label = tf), colour = "#1b2a3a", size = 4.2, max.overlaps = 50,
                    show.legend = FALSE, min.segment.length = 0, seed = 1, box.padding = 0.4) +
    scale_colour_manual(values = famcol, name = "family", drop = FALSE) +
    scale_y_sqrt(breaks = c(5, 50, 100, 200, 400, 600, 800)) +
    labs(x = "log2 odds ratio", y = "−log10 p (square-root scale)", title = title) +
    theme(plot.title = element_text(colour = col))
}
p7 <- vol(eB, sprintf("Motifs in B-up peaks (n = %s)", format(sum(b_up), big.mark = " ")), colB) +
  vol(eT, sprintf("Motifs in T-up peaks (n = %s)", format(sum(t_up), big.mark = " ")), colT) +
  plot_layout(guides = "collect")
ggsave("fig07-enrichment-volcano.png", p7, width = 14, height = 6, dpi = 110, bg = "white")

# ---- fig08: B-up vs T-up odds ratios ---------------------------------------------
d8 <- data.frame(tf = eB$tf, B = log2(eB$OR), T = log2(eT$OR[match(eB$motif, eT$motif)]),
                 fs = fam_short(eB$family, eB$tf))
d8$motif <- eB$motif
lab8 <- d8[d8$motif %in% c("Spi1_MA0080.7", "SPIB_MA0081.3", "ETS1_MA0098.4", "GABPA_MA0062.4",
                           "Runx1_MA0002.3", "TCF7_MA0769.3", "POU2F2_MA0507.3", "IRF4_MA1419.2",
                           "EBF1_MA0154.5", "CTCF_MA0139.2", "Bcl11B_MA1989.2", "PAX5_MA0014.4",
                           "KLF4_MA0039.5"), ]
p8 <- ggplot(d8, aes(B, T, colour = fs)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(data = d8[d8$fs == "other", ], size = 1.6, alpha = 0.6) +
  geom_point(data = d8[d8$fs != "other", ], size = 2.8) +
  geom_text_repel(data = lab8, aes(label = tf), colour = "#1b2a3a", size = 4.5, show.legend = FALSE,
                  min.segment.length = 0, seed = 2, max.overlaps = 50, box.padding = 0.5) +
  scale_colour_manual(values = famcol, name = "family") +
  labs(x = "log2 OR in B-up peaks", y = "log2 OR in T-up peaks")
ggsave("fig08-or-scatter.png", p8, width = 9.5, height = 6.2, dpi = 110, bg = "white")

# Set ENRICHMENT_ONLY=1 to update only the figures that depend on the saved classes.
if (Sys.getenv("ENRICHMENT_ONLY") == "1") quit(save = "no")

# ---- D. chromVAR across 13 cell types --------------------------------------------
A <- cached("chromvar-input", {
  cnt <- read.delim(file.path(D, "corces2016_atac_counts.txt.gz"), check.names = FALSE)
  peaks19 <- GRanges(cnt$Chr, IRanges(cnt$Start + 1, cnt$End))
  cnt <- as.matrix(cnt[, -(1:3)])
  code <- c("1" = "HSC", "2" = "MPP", "3" = "LMPP", "4" = "CMP", "5" = "GMP", "6" = "MEP",
            "7" = "Mono", "9" = "CD4", "10" = "CD8", "11" = "NK", "13" = "B", "14" = "CLP",
            "15" = "Ery")
  nm <- colnames(cnt); ok <- grepl("^(Donor)?[0-9]{4}-", nm)
  info <- data.frame(sample = nm[ok],
                     donor = sub("^(Donor)?([0-9]{4})-.*", "\\2", nm[ok]),
                     celltype = code[sub(".*-([0-9]{1,2})[A-C]?(-[0-9]{6})?$", "\\1", nm[ok])])
  chain <- tempfile(fileext = ".chain")
  writeLines(readLines(file.path(D, "hg19ToHg38.over.chain.gz")), chain)
  lo <- liftOver(peaks19, import.chain(chain))
  keep <- which(lengths(lo) == 1)
  peaks <- unlist(lo[keep])
  ok <- width(peaks) == 500 & as.character(seqnames(peaks)) == as.character(seqnames(peaks19)[keep])
  keep <- keep[ok]; peaks <- peaks[ok]
  cnt <- cnt[keep, info$sample]
  list(info = info, peaks = peaks, cnt = cnt)
})
dev <- cached("dev", {
  info <- A$info
  grp <- paste(info$celltype, info$donor, sep = "_")
  cm <- sapply(split(info$sample, grp), function(s) rowSums(A$cnt[, s, drop = FALSE]))
  top <- order(rowSums(cm), decreasing = TRUE)[1:100000]
  se <- SummarizedExperiment(assays = list(counts = cm[top, ]), rowRanges = A$peaks[top],
                             colData = DataFrame(celltype = sub("_.*", "", colnames(cm))))
  se <- addGCBias(se, genome = genome)
  mm <- matchMotifs(pfms, se, genome = genome, p.cutoff = 5e-5)
  bg <- getBackgroundPeaks(se, niterations = 50)
  computeDeviations(se, annotations = mm, background_peaks = bg)
})
z <- deviationScores(dev)
v <- computeVariability(dev); v <- v[order(-v$variability), ]
cat("chromVAR samples:", ncol(dev), "\n"); print(table(dev$celltype))
print(head(v[, c("name", "variability")], 25), digits = 3)
# top variable motifs, skipping near-duplicates (|r| > 0.9 with a motif already chosen)
pick <- character(0)
for (m in rownames(v)) {
  if (length(pick) == 0 || all(abs(cor(z[m, ], t(z[pick, , drop = FALSE]))) < 0.9)) pick <- c(pick, m)
  if (length(pick) == 15) break
}
print(v[pick, c("name", "variability")], digits = 3)
# heatmap: motifs of well-known haematopoietic regulators (plus CTCF)
known <- c("HLF_MA0043.4", "HOXA9_MA0594.3", "MEIS1_MA0498.3", "GATA1::TAL1_MA0140.3",
           "KLF1_MA0493.3", "NFE2_MA0841.2", "CEBPA_MA0102.5", "CEBPB_MA0466.4", "Spi1_MA0080.7",
           "SPIB_MA0081.3", "IRF8_MA0652.2", "EBF1_MA0154.5", "TCF3_MA0522.4", "POU2F2_MA0507.3",
           "PAX5_MA0014.4", "MEF2C_MA0497.2", "NFKB1_MA0105.4", "Runx1_MA0002.3", "ETS1_MA0098.4",
           "TCF7_MA0769.3", "Lef1_MA0768.3", "Gata3_MA0037.5", "TBX21_MA0690.3", "EOMES_MA0800.2",
           "FOS::JUN_MA0099.4", "CTCF_MA0139.2")
stopifnot(all(known %in% rownames(z)))
ord_ct <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Ery", "Mono", "CLP", "B", "CD4", "CD8", "NK")
o <- order(match(dev$celltype, ord_ct), colnames(dev))
zz <- z[known, o]; rownames(zz) <- sub("_MA.*", "", rownames(zz))
zz <- pmax(pmin(zz, 15), -15)
ann <- data.frame(celltype = factor(dev$celltype[o], levels = ord_ct), row.names = colnames(zz))
ctcol <- setNames(c("#5b6b7a", "#7f8c99", "#a3adb7", "#6b4fa0", "#9b7fc4", "#a8324a", "#d4626f",
                    "#2e8b57", "#d4a017", colB, colT, "#f0b27a", "#8c564b"), ord_ct)
png("fig11-chromvar-heatmap.png", width = 1500, height = 760, res = 110)
pheatmap(zz, cluster_cols = FALSE, cluster_rows = TRUE, show_colnames = FALSE, annotation_col = ann,
         annotation_colors = list(celltype = ctcol), gaps_col = cumsum(table(ann$celltype)),
         color = colorRampPalette(c("#2a6f97", "white", "#c0392b"))(101),
         breaks = seq(-15, 15, length.out = 102), fontsize = 13, fontsize_row = 12,
         treeheight_row = 25, main = "chromVAR deviation z-scores (capped at ±15)")
dev.off()

# what the motif looks like (JASPAR names the TF the matrix was measured for, not the blood TF)
like <- c("GATA1::TAL1" = "GATA", SPIC = "PU.1-type ETS", "ERF::FIGLA" = "E-box", FIGLA = "E-box",
          BHLHE22 = "E-box", ASCL1 = "E-box", ZBTB7A = "GC-rich zinc finger", IRF9 = "IRF",
          "STAT1::STAT2" = "IRF / ISRE", Atf3 = "AP-1", "MAFG::NFE2L1" = "NFE2 / MAF",
          ETV1 = "ETS", ETV2 = "ETS", ETV7 = "ETS", TBR1 = "T-box")
d12 <- data.frame(tf = sub("_MA.*", "", pick), var = v[pick, "variability"])
d12$tf <- ifelse(d12$tf %in% names(like), paste0(d12$tf, "  (", like[d12$tf], ")"), d12$tf)
d12$tf <- factor(d12$tf, levels = rev(d12$tf))
p12a <- ggplot(d12, aes(var, tf)) + geom_col(fill = "#5b6b7a", width = 0.7) +
  labs(x = "variability (sd of z across samples)", y = NULL,
       title = "Most variable motifs (one per cluster)")
zb <- data.frame(celltype = factor(dev$celltype, levels = ord_ct),
                 GATA1 = z["GATA1::TAL1_MA0140.3", ], SPI1 = z["Spi1_MA0080.7", ],
                 EBF1 = z["EBF1_MA0154.5", ], TCF7 = z["TCF7_MA0769.3", ], CEBPA = z["CEBPA_MA0102.5", ])
zl <- reshape(zb, direction = "long", varying = names(zb)[-1], v.names = "z", timevar = "motif",
              times = names(zb)[-1])
zl$motif <- factor(zl$motif, levels = names(zb)[-1])
p12b <- ggplot(zl, aes(celltype, z, colour = celltype)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 2.2, position = position_jitter(width = 0.15, seed = 1)) +
  facet_wrap(~ motif, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_colour_manual(values = ctcol, guide = "none") +
  labs(x = NULL, y = "deviation z-score", title = "One dot per donor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
ggsave("fig12-chromvar-variability.png", p12a + p12b + plot_layout(widths = c(1, 1.3)),
       width = 14, height = 6.6, dpi = 110, bg = "white")

# ---- E. CTCF footprint in GM12878 ATAC-seq (chr21) ----------------------------------
# The BAM is made with the commands of the exercise (C1) unless ATAC_BAM points to one.
bam <- Sys.getenv("ATAC_BAM", file.path(tempdir(), "gm12878_atac_rep1.bam"))
if (!file.exists(bam)) {
  tmp <- tempfile("atac"); dir.create(tmp)
  cmd <- sprintf(paste(
    "fastp -i %1$s/gm12878_atac_rep1_R1.fastq.gz -I %1$s/gm12878_atac_rep1_R2.fastq.gz",
    "-o %2$s/R1.fq.gz -O %2$s/R2.fq.gz --detect_adapter_for_pe -w 8 -j %2$s/fp.json -h %2$s/fp.html &&",
    "bowtie2 -p 8 -X 2000 --very-sensitive -x %1$s/bowtie2/hg38_chr21_chrM",
    "-1 %2$s/R1.fq.gz -2 %2$s/R2.fq.gz 2> %2$s/bt2.log",
    "| samtools fixmate -m -u - - | samtools view -u -f 2 -q 30 - | samtools sort -u -",
    "| samtools markdup -r -s - %3$s 2> %2$s/markdup.log && samtools index %3$s"), D, tmp, bam)
  stopifnot(system(cmd) == 0)
}
chr21 <- GRanges("chr21", IRanges(1, 46709983))
ga <- readGAlignments(bam, param = ScanBamParam(which = chr21))
ins <- GRanges("chr21", IRanges(ifelse(strand(ga) == "+", start(ga) + 4, end(ga) - 4), width = 1),
               seqinfo = Seqinfo("chr21", 46709983))
cat("Tn5 insertions on chr21:", length(ins), "\n")
cov_ins <- coverage(ins)$chr21

ctcf_pk <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak")
atac_pk <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
ctcf_pfm <- pfms[["CTCF_MA0139.2"]]
sites <- matchMotifs(ctcf_pfm, chr21, genome = genome, out = "positions", p.cutoff = 1e-4)[[1]]
sites <- sites[!overlapsAny(sites, import(file.path(D, "hg38_blacklist.bed.gz")))]
grp <- ifelse(overlapsAny(sites, ctcf_pk), "in CTCF ChIP peak",
              ifelse(overlapsAny(sites, atac_pk), "open, no ChIP peak", "closed, no ChIP peak"))
lev <- c("in CTCF ChIP peak", "open, no ChIP peak", "closed, no ChIP peak")
ns <- table(factor(grp, levels = lev)); print(ns)
pk21 <- ctcf_pk[seqnames(ctcf_pk) == "chr21"]
with_motif <- overlapsAny(pk21, sites)
cat("CTCF motif sites on chr21:", length(sites), " fraction in a ChIP peak:",
    round(mean(grp == lev[1]), 3), "\n")
cat("chr21 ChIP peaks:", length(pk21), " with a motif site:", sum(with_motif),
    round(mean(with_motif), 3), "\n")
prof <- function(s) {                # mean insertions per site at -100..+100, motif-oriented
  w <- promoters(resize(s, 1, fix = "center"), upstream = 100, downstream = 101)
  m <- as.matrix(Views(cov_ins, ranges(w)))
  minus <- as.logical(strand(w) == "-")
  m[minus, ] <- m[minus, ncol(m):1]
  colMeans(m)
}
fp <- do.call(rbind, lapply(lev, function(g)
  data.frame(grp = g, pos = -100:100, ins = prof(sites[grp == g]))))
fp$smooth <- unlist(lapply(split(fp$ins, factor(fp$grp, levels = lev)),
                           function(x) as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2))))
nlab <- setNames(sprintf("%s (n = %s)", lev, formatC(as.numeric(ns[lev]), big.mark = " ", format = "d")), lev)
fp$lab <- factor(nlab[fp$grp], levels = nlab)
for (g in lev) {
  b <- fp[fp$grp == g, ]
  cat(sprintf("%-22s flank(30-60) %.3f  centre(+-10) %.3f  far(>=90) %.3f  centre/flank %.2f\n", g,
              mean(b$ins[abs(b$pos) >= 30 & abs(b$pos) <= 60]), mean(b$ins[abs(b$pos) <= 10]),
              mean(b$ins[abs(b$pos) >= 90]),
              mean(b$ins[abs(b$pos) <= 10]) / mean(b$ins[abs(b$pos) >= 30 & abs(b$pos) <= 60])))
}
mlen <- ncol(TFBSTools::Matrix(ctcf_pfm))
fpcol <- c("#8c564b", "#d4a017", "#b8c2cc")
p15 <- ggplot(fp, aes(pos, colour = lab)) +
  annotate("rect", xmin = -mlen / 2, xmax = mlen / 2, ymin = -Inf, ymax = Inf, fill = "#8c564b",
           alpha = 0.12) +
  annotate("text", x = 0, y = max(fp$ins) * 1.05, label = sprintf("CTCF motif (%d bp)", mlen), colour = "#8c564b",
           size = 4.6) +
  geom_line(aes(y = ins), linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = smooth), linewidth = 1.3, na.rm = TRUE) +
  scale_colour_manual(values = fpcol, name = NULL) +
  labs(x = "position relative to the motif centre (bp; motif read 5' to 3')",
       y = "Tn5 insertions per site", title = "GM12878 ATAC-seq, chr21: CTCF motif matches") +
  theme(legend.position = "bottom")
ggsave("fig15-ctcf-footprint.png", p15, width = 11, height = 5.4, dpi = 110, bg = "white")

# fig19: motif sites vs ChIP peaks on chr21
d19 <- rbind(data.frame(bar = "CTCF motif matches", what = lev, n = as.numeric(ns[lev])),
             data.frame(bar = "CTCF ChIP peaks",
                        what = c("contain a motif match", "no motif match"),
                        n = c(sum(with_motif), sum(!with_motif))))
d19$bar <- factor(d19$bar, levels = c("CTCF ChIP peaks", "CTCF motif matches"))
d19$what <- factor(d19$what, levels = c(lev, "contain a motif match", "no motif match"))
d19 <- do.call(rbind, lapply(split(d19, d19$bar), function(x) {
  x$frac <- x$n / sum(x$n); x$mid <- 1 - (cumsum(x$frac) - x$frac / 2); x }))
fill19 <- c(fpcol[1], fpcol[2], "#8a97a5", "#8c564b", "#c8a99f")
small <- d19[d19$frac <= 0.06, ]
small$x <- c(0, 0.36)[seq_len(nrow(small))]
small$col <- fill19[as.integer(small$what)]
p19 <- ggplot(d19, aes(frac, bar, fill = what)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE), colour = "white") +
  geom_text(data = d19[d19$frac > 0.06, ],
            aes(x = 1 - mid, label = sprintf("%s\n%s (%.0f%%)", what, format(n, big.mark = " "), 100 * frac)),
            colour = "white", size = 4.3, lineheight = 0.9) +
  geom_text(data = small, aes(x = x, y = as.numeric(bar) + 0.42,
                              label = sprintf("%s: %s (%.1f%%)", what, n, 100 * frac)),
            colour = small$col, size = 4.3, hjust = 0, fontface = "bold") +
  scale_fill_manual(values = fill19, guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = NULL, y = NULL, title = "chr21, GM12878: CTCF motif matches (p < 1e-4) and CTCF ChIP-seq peaks")
ggsave("fig19-motif-vs-binding.png", p19, width = 12, height = 4.2, dpi = 110, bg = "white")

# fig16: Tn5 sequence preference around insertion sites, oriented by read strand
set.seed(2)
k <- sample(length(ga), 200000)
st <- as.character(strand(ga))[k]
pos <- ifelse(st == "+", start(ga)[k] + 4, end(ga)[k] - 4)
win <- GRanges("chr21", IRanges(pos - 30, pos + 30), strand = st)
cm <- consensusMatrix(getSeq(genome, win), as.prob = TRUE)[c("A", "C", "G", "T"), ]
colnames(cm) <- -30:30
# compare with the same reads 20-30 bp away: removes the regional GC richness of open chromatin
bgf <- rowMeans(cm[, abs(-30:30) >= 20])
d16 <- data.frame(pos = rep(-15:15, each = 4), base = rep(c("A", "C", "G", "T"), 31),
                  enr = log2(as.vector(cm[, as.character(-15:15)]) / rep(bgf, 31)))
p16 <- ggplot(d16, aes(pos, enr, colour = base)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 1.1) + geom_point(size = 2) +
  scale_colour_manual(values = c(A = "#2e8b57", C = "#2a6f97", G = "#d4a017", T = "#c0392b")) +
  scale_x_continuous(breaks = seq(-15, 15, 3)) +
  labs(x = "position relative to the Tn5 insertion (bp)", y = "log2(frequency / 20-30 bp away)",
       title = "Bases around 200 000 Tn5 insertions (GM12878, chr21)", colour = NULL)
ggsave("fig16-tn5-bias.png", p16, width = 10, height = 4.6, dpi = 110, bg = "white")
print(round(log2(cm[, as.character(-6:6)] / bgf), 2))
cat("largest |log2 enrichment| within 10 bp:", round(max(abs(d16$enr[abs(d16$pos) <= 10])), 2), "\n")
