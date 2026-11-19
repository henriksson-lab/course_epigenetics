# Real-data figures for chromatin-states/lec.qmd.
# Run from this folder:  Rscript make-chromstates-figs.R
#
# Needs: data/prepare-cuttag.sh (GM12878 ChIP-seq BAMs gm12878_*_chr21.bam), data/download.sh
#        chromatin-states (Roadmap E116 segmentation and 18-state emissions, GTF, GM12878
#        fold-change bigWigs, ATAC and CTCF peaks); Java and ChromHMM (CHROMHMM=/path/to/ChromHMM.jar,
#        default ChromHMM/ChromHMM.jar; its COORDS folder is used for CpG islands);
#        R packages GenomicRanges, rtracklayer, ggplot2, patchwork.
#
# ChromHMM runs once into ./states-work (set STATES_WORK to change; ~50 MB, a few minutes for all K).
#
# Writes fig09-tracks-states.png, fig10-combinations.png, fig11-emissions.png, fig12-transitions.png,
# fig13-enrichment.png, fig14-tss-neighbourhood.png, fig15-nstates.png, fig16-vs-roadmap.png,
# fig17-coverage.png and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(ggplot2); library(patchwork)
})
D <- normalizePath("../data")
W <- Sys.getenv("STATES_WORK", "states-work")
dir.create(W, showWarnings = FALSE, recursive = TRUE); W <- normalizePath(W)
CHMM <- normalizePath(Sys.getenv("CHROMHMM", "ChromHMM/ChromHMM.jar"))
COORDS <- file.path(dirname(CHMM), "COORDS", "hg38")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
marks <- c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3", "H3K27me3", "H3K9me3")
mcol <- c(H3K4me3 = "#c0392b", H3K27ac = "#e67e22", H3K4me1 = "#d4a017", H3K36me3 = "#2e8b57",
          H3K27me3 = "#2a6f97", H3K9me3 = "#6b4fa0")
L21 <- 46709983
KS <- c(4, 6, 8, 10, 12, 15, 18)

# ---- 1. ChromHMM (exercise C1-C2) --------------------------------------------
sh <- function(cmd) { rc <- system2("bash", c("-c", shQuote(paste0("set -euo pipefail; cd ", W, "; ", cmd))))
                      if (rc != 0) stop("failed: ", cmd) }
if (!file.exists(file.path(W, "binary/GM12878_chr21_binary.txt"))) {
  writeLines("chr21\t46709983", file.path(W, "chr21.sizes"))
  writeLines(sprintf("GM12878\t%s\tgm12878_%s_chip_chr21.bam\tgm12878_input_chr21.bam", marks, tolower(marks)),
             file.path(W, "cellmarkfiletable.txt"))
  sh(sprintf("java -mx4G -jar %s BinarizeBam -b 200 chr21.sizes %s cellmarkfiletable.txt binary > binarize.log", CHMM, D))
}
for (K in KS) if (!file.exists(file.path(W, sprintf("model_%d/emissions_%d.txt", K, K))))
  sh(sprintf("java -mx4G -jar %s LearnModel -b 200 -l chr21.sizes -noautoopen binary model_%d %d hg38 > learn_%d.log",
             CHMM, K, K, K))

# ---- 2. read -------------------------------------------------------------------
b <- read.table(file.path(W, "binary/GM12878_chr21_binary.txt"), skip = 1, header = TRUE)[, marks]
cat("== binarised: bins", nrow(b), "; % of bins called 1:\n"); print(round(100 * colMeans(b), 2))
cat("no mark:", round(100 * mean(rowSums(b) == 0), 1), "%\n")
reads <- sapply(c(tolower(marks), "input"), function(m) {
  f <- file.path(D, if (m == "input") "gm12878_input_chr21.bam" else sprintf("gm12878_%s_chip_chr21.bam", m))
  as.numeric(system2("samtools", c("view", "-c", f), stdout = TRUE)) })
cat("reads on chr21:\n"); print(reads)
rd <- function(f) read.table(file.path(W, f), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
emis <- function(K) { e <- rd(sprintf("model_%d/emissions_%d.txt", K, K)); e[, marks] }
tran <- function(K) as.matrix(rd(sprintf("model_%d/transitions_%d.txt", K, K)))
loglik <- function(K) { l <- readLines(file.path(W, sprintf("learn_%d.log", K)))
  l <- grep("^\\s+[0-9]+\\s+-[0-9.]+", l, value = TRUE); v <- strsplit(trimws(tail(l, 1)), "\\s+")[[1]]
  c(iterations = as.numeric(v[1]), loglik = as.numeric(v[2])) }

# Name states from their emissions (rules in the spirit of the Roadmap mnemonics).
name_states <- function(e) {
  nm <- apply(e, 1, function(x) {
    if (x["H3K4me3"] >= 0.5 && x["H3K27me3"] >= 0.25) "TssBiv"
    else if (x["H3K4me3"] >= 0.5 && x["H3K27ac"] >= 0.5 && x["H3K4me1"] >= 0.5) "TssFlnkU"
    else if (x["H3K4me3"] >= 0.5 && x["H3K27ac"] >= 0.5) "TssA"
    else if (x["H3K4me3"] >= 0.4) "TssFlnk"
    else if (x["H3K27ac"] >= 0.5 && x["H3K4me1"] >= 0.5) "EnhA"
    else if (x["H3K4me1"] >= 0.3 || x["H3K27ac"] >= 0.3) "EnhWk"
    else if (max(x) < 0.02) "Quies"
    else {
      top <- names(which.max(x[c("H3K36me3", "H3K27me3", "H3K9me3")]))
      if (top == "H3K36me3") { if (x[top] >= 0.2) "Tx" else "TxWk" }
      else if (top == "H3K27me3") { if (x[top] >= 0.2) "ReprPC" else "ReprPCWk" }
      else { if (x[top] >= 0.25) "ZNF/Rpts" else "Het" }
    }
  })
  nm <- as.character(nm)
  for (d in unique(nm[duplicated(nm)])) { i <- which(nm == d); nm[i] <- paste0(d, c("", seq_along(i)[-1])) }
  paste0("E", seq_along(nm), " ", nm)
}
E8 <- emis(8); N8 <- name_states(E8); names(N8) <- paste0("E", 1:8)
cat("\n== 8-state model: emissions and names\n"); print(cbind(round(E8, 2), name = N8))
for (K in KS) cat("K =", K, " iterations / final log-likelihood:", loglik(K), "\n")

seg <- import(file.path(W, "model_8/GM12878_8_segments.bed"), format = "bed")
bins <- GRanges("chr21", IRanges(seq(1, by = 200, length.out = nrow(b)), width = 200))   # ChromHMM's bins
bins$ours <- seg$name[findOverlaps(bins, seg, select = "first")]
bins$oursN <- factor(N8[bins$ours], levels = N8)
rm <- import(file.path(D, "roadmap_E116_18state_hg38.bed.gz"), format = "bed"); rm <- rm[seqnames(rm) == "chr21"]
rmlev <- c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2", "9_EnhA1",
           "10_EnhA2", "11_EnhWk", "12_ZNF/Rpts", "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC", "17_ReprPCWk", "18_Quies")
bins$rm <- factor(rm$name[findOverlaps(resize(bins, 1, fix = "center"), rm, select = "first")], levels = rmlev)
cat("bins without a Roadmap state (not lifted over):", sum(is.na(bins$rm)), "\n")
scol <- c(TssA = "#ff0000", TssFlnk = "#ff4500", TssFlnkU = "#ff4500", TssBiv = "#cd5c5c", EnhA = "#ffc34d",
          EnhWk = "#ffff00", Tx = "#008000", TxWk = "#006400", ReprPC = "#808080", ReprPCWk = "#c0c0c0",
          "ZNF/Rpts" = "#66cdaa", Het = "#8a91d0", Quies = "#e8e8e8")   # Roadmap colours
col_of <- function(n) { base <- sub("^E[0-9]+ ", "", n); base <- sub("[0-9]+$", "", base); unname(scol[base]) }
ncol8 <- setNames(col_of(N8), N8)
rmcol <- c("1_TssA" = "#ff0000", "2_TssFlnk" = "#ff4500", "3_TssFlnkU" = "#ff4500", "4_TssFlnkD" = "#ff4500",
           "5_Tx" = "#008000", "6_TxWk" = "#006400", "7_EnhG1" = "#c2e105", "8_EnhG2" = "#c2e105",
           "9_EnhA1" = "#ffc34d", "10_EnhA2" = "#ffc34d", "11_EnhWk" = "#ffff00", "12_ZNF/Rpts" = "#66cdaa",
           "13_Het" = "#8a91d0", "14_TssBiv" = "#cd5c5c", "15_EnhBiv" = "#bdb76b", "16_ReprPC" = "#808080",
           "17_ReprPCWk" = "#c0c0c0", "18_Quies" = "#e8e8e8")   # Roadmap colours

# ---- fig10: combinations ----------------------------------------------------------
combo <- apply(b, 1, paste, collapse = "")
tab <- sort(table(combo), decreasing = TRUE)
nz <- tab[names(tab) != "000000"]
cat("\n== combinations: distinct", length(tab), "of 64; no mark", round(100 * tab[["000000"]] / sum(tab), 1),
    "%; top 10 of non-empty cover", round(100 * sum(head(nz, 10)) / sum(nz), 1), "%\n")
top <- head(nz, 12)
cd <- data.frame(combo = factor(names(top), levels = names(top)), pct = 100 * as.numeric(top) / sum(nz))
print(cd)
dots <- do.call(rbind, lapply(seq_along(top), function(i) data.frame(combo = cd$combo[i], mark = factor(marks, levels = rev(marks)),
                                                                   on = strsplit(names(top)[i], "")[[1]] == "1")))
p10a <- ggplot(cd, aes(combo, pct)) + geom_col(fill = "#5b6b7a", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.3, size = 3.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "% of bins with ≥ 1 mark", title = sprintf("Top 12 of %d combinations (%.0f%% of bins have no mark)",
                                                             length(tab), 100 * tab[["000000"]] / sum(tab))) +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
p10b <- ggplot(dots, aes(combo, mark)) + geom_point(aes(colour = ifelse(on, as.character(mark), "off")), size = 4.5) +
  scale_colour_manual(values = c(mcol, off = "#e3e8ee"), guide = "none") + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())
ggsave("fig10-combinations.png", p10a / p10b + plot_layout(heights = c(2.2, 1.4)), width = 11, height = 5.6, dpi = 150, bg = "white")

# ---- fig11: emissions (ours and Roadmap) -----------------------------------------
heat <- function(e, title, rlab) {
  d <- data.frame(state = factor(rep(rlab, ncol(e)), levels = rev(rlab)), mark = factor(rep(colnames(e), each = nrow(e)), levels = marks),
                  p = as.vector(as.matrix(e)))
  ggplot(d, aes(mark, state, fill = p)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", p), colour = p > 0.5), size = 3.4) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "#1b2a3a"), guide = "none") +
    scale_fill_gradient(low = "white", high = "#1b3a6b", limits = c(0, 1), name = "P(mark)") +
    labs(x = NULL, y = NULL, title = title) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"), panel.grid = element_blank())
}
rmE <- read.table(file.path(D, "roadmap_18state_emissions.txt"), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)[, marks]
p11 <- heat(E8, "Our 8 states (GM12878, chr21)", N8) + theme(legend.position = "none") +
  heat(rmE, "Roadmap 18-state model", rownames(rmE)) + plot_layout(widths = c(1, 1))
ggsave("fig11-emissions.png", p11, width = 12.5, height = 6.2, dpi = 150, bg = "white")

# ---- fig12: transitions ---------------------------------------------------------
T8 <- tran(8)
td <- data.frame(from = factor(rep(N8, 8), levels = rev(N8)), to = factor(rep(N8, each = 8), levels = N8), p = as.vector(T8))
p12 <- ggplot(td, aes(to, from, fill = p)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(p >= 0.005, sprintf("%.2f", p), "")), size = 3.6, colour = ifelse(td$p > 0.5, "white", "#1b2a3a")) +
  scale_fill_gradient(low = "white", high = "#1b3a6b", limits = c(0, 1), name = "P") +
  labs(x = "to (next bin)", y = "from") + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
ggsave("fig12-transitions.png", p12, width = 9.5, height = 6, dpi = 150, bg = "white")
cat("\n== transitions (diagonal):", round(diag(T8), 3), "\nexpected run length (kb):", round(0.2 / (1 - diag(T8)), 1), "\n")
off <- T8; diag(off) <- 0
for (i in 1:8) cat(N8[i], "-> most likely next other state:", N8[which.max(off[i, ])], round(max(off[i, ]), 3), "\n")

# ---- fig13: enrichments ----------------------------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz")); gtf <- gtf[seqnames(gtf) == "chr21"]
txp <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding"]
tss <- unique(resize(txp, 1, fix = "start"))
cgi <- import(file.path(COORDS, "CpGIsland.hg38.bed.gz"), format = "bed"); cgi <- cgi[seqnames(cgi) == "chr21"]
pk <- function(f) { x <- import(file.path(D, f), format = "narrowPeak"); x[seqnames(x) == "chr21"] }
ann <- list("CpG island" = cgi, "TSS" = tss, "TSS ± 2 kb" = resize(tss, 4000, fix = "center"),
            "exon" = gtf[gtf$type == "exon" & gtf$transcript_type == "protein_coding"],
            "gene body" = reduce(txp), "ATAC peak" = pk("gm12878_atac_peaks.bed.gz"), "CTCF peak" = pk("gm12878_ctcf_peaks.bed.gz"))
frac_state <- table(bins$oursN) / length(bins)
fe <- sapply(ann, function(a) { sel <- overlapsAny(bins, a); (table(bins$oursN[sel]) / sum(sel)) / frac_state })
fed <- data.frame(state = factor(rep(N8, ncol(fe)), levels = rev(N8)), ann = factor(rep(colnames(fe), each = 8), levels = colnames(fe)),
                  fe = as.vector(fe))
gp <- data.frame(state = factor(N8, levels = rev(N8)), pct = 100 * as.numeric(frac_state))
p13a <- ggplot(gp, aes(1, state)) + geom_tile(fill = "#f4f6f8", colour = "white") + geom_text(aes(label = sprintf("%.1f", pct)), size = 3.8) +
  scale_x_continuous(breaks = 1, labels = "% of chr21") + labs(x = NULL, y = NULL) + theme(panel.grid = element_blank())
p13b <- ggplot(fed, aes(ann, state, fill = log2(pmax(fe, 0.01)))) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(fe >= 10, sprintf("%.0f", fe), sprintf("%.1f", fe))), size = 3.8) +
  scale_fill_gradient2(low = "#2a6f97", mid = "white", high = "#c0392b", midpoint = 0, limits = c(-5, 7), oob = scales::squish, name = "log2 FE") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
ggsave("fig13-enrichment.png", p13a + p13b + plot_layout(widths = c(1, 7)), width = 12, height = 5.4, dpi = 150, bg = "white")
cat("\n== fold enrichments\n"); print(round(fe, 2)); cat("% of chr21 per state:\n"); print(round(100 * frac_state, 2))

# ---- fig14: around the TSS (our own neighbourhood enrichment, strand-aware) --------
st <- as.character(strand(tss))
pos <- seq(-5000, 5000, 200)
nb <- do.call(rbind, lapply(pos, function(p) {
  q <- GRanges("chr21", IRanges(start(tss) + ifelse(st == "-", -p, p), width = 1))
  s <- bins$oursN[findOverlaps(q, bins, select = "first")]
  data.frame(pos = p, state = levels(bins$oursN), fe = as.numeric(table(s) / length(s) / frac_state))
}))
nb$state <- factor(nb$state, levels = N8)
p14 <- ggplot(nb, aes(pos / 1000, fe, colour = state)) + geom_hline(yintercept = 1, colour = "grey60") + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = setNames(c("#2e8b57", "#c9a400", "#e67e22", "#c0392b", "#8c564b", "#9aa5b1", "#3cb4a0", "#2a6f97")[seq_along(N8)], N8)) +
  scale_y_log10() + labs(x = "position relative to protein-coding TSS (kb, transcription →)", y = "fold enrichment (log scale)", colour = NULL)
ggsave("fig14-tss-neighbourhood.png", p14, width = 11, height = 4.8, dpi = 150, bg = "white")
cat("\n== TSS neighbourhood: max fold enrichment and its position (kb)\n")
print(do.call(rbind, lapply(split(nb, nb$state), function(d) data.frame(max = round(max(d$fe), 1), at = d$pos[which.max(d$fe)] / 1000))))

# ---- fig15: number of states ---------------------------------------------------------
ll <- t(sapply(KS, loglik)); lld <- data.frame(K = KS, ll = ll[, "loglik"])
p15a <- ggplot(lld, aes(K, ll)) + geom_line(colour = "#1b2a3a") + geom_point(size = 3, colour = "#1b2a3a") +
  scale_x_continuous(breaks = KS) + labs(x = "number of states K", y = "log-likelihood (ChromHMM)", title = "A. Likelihood always increases")
E4 <- emis(4); E12 <- emis(12)
p15b <- heat(E4, "B. K = 4", name_states(E4)) + theme(legend.position = "none", axis.text.x = element_text(size = 9, angle = 30, hjust = 1))
p15c <- heat(E12, "C. K = 12", name_states(E12)) + theme(legend.position = "none", axis.text.x = element_text(size = 9, angle = 30, hjust = 1))
ggsave("fig15-nstates.png", p15a + p15b + p15c + plot_layout(widths = c(1.1, 1, 1)), width = 14, height = 5.6, dpi = 150, bg = "white")
cat("\n== K = 4 names:", name_states(E4), "\n== K = 12 names:", name_states(E12), "\n")

# ---- fig16: ours vs Roadmap ------------------------------------------------------
ok <- !is.na(bins$rm)
ct <- table(ours = bins$oursN[ok], roadmap = bins$rm[ok])
pr <- 100 * prop.table(ct, 1)
cd <- as.data.frame(pr); names(cd)[3] <- "pct"
cd$ours <- factor(cd$ours, levels = rev(N8))
p16 <- ggplot(cd, aes(roadmap, ours, fill = pct)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct >= 5, sprintf("%.0f", pct), "")), size = 3.5, colour = ifelse(cd$pct > 50, "white", "#1b2a3a")) +
  scale_fill_gradient(low = "white", high = "#1b3a6b", limits = c(0, 100), name = "% of our\nstate's bins") +
  labs(x = "Roadmap E116 (GM12878), 18-state model", y = "our state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave("fig16-vs-roadmap.png", p16, width = 12, height = 5.6, dpi = 150, bg = "white")
cat("\n== ours vs Roadmap (% of our state's bins; top 3 per state)\n")
for (s in rownames(pr)) { x <- sort(pr[s, ], decreasing = TRUE)[1:3]; cat(s, ":", paste(sprintf("%s %.0f%%", names(x), x), collapse = ", "), "\n") }
cat("Roadmap state shares on chr21 (%):\n"); print(round(100 * prop.table(table(bins$rm[ok])), 2))
# a coarse agreement: map both to five classes
cls <- function(x) { x <- as.character(x)
  ifelse(grepl("Tss", x), "promoter", ifelse(grepl("Enh", x), "enhancer", ifelse(grepl("Tx", x), "transcribed",
  ifelse(grepl("ReprPC", x), "Polycomb", ifelse(grepl("Het|ZNF", x), "heterochromatin", "quiescent"))))) }
agree <- mean(cls(bins$oursN[ok]) == cls(bins$rm[ok]))
cat("agreement in 6 coarse classes (all bins):", round(100 * agree, 1), "%; excluding bins quiescent in both:",
    round(100 * mean((cls(bins$oursN[ok]) == cls(bins$rm[ok]))[!(cls(bins$oursN[ok]) == "quiescent" & cls(bins$rm[ok]) == "quiescent")]), 1), "%\n")

# ---- fig17: coverage ---------------------------------------------------------------
cov1 <- data.frame(state = factor(N8, levels = N8), pct = 100 * as.numeric(frac_state), model = "ours (8 states)")
cov2 <- data.frame(state = factor(rmlev, levels = rmlev), pct = 100 * as.numeric(prop.table(table(bins$rm[ok]))), model = "Roadmap (18 states)")
p17a <- ggplot(cov1, aes(state, pct, fill = state)) + geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.3, size = 3.6) +
  scale_fill_manual(values = ncol8, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "% of chr21", title = "A. Our 8 states") + theme(axis.text.x = element_text(angle = 35, hjust = 1))
p17b <- ggplot(cov2, aes(state, pct, fill = state)) + geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = rmcol, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = NULL, title = "B. Roadmap E116, 18 states") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
ggsave("fig17-coverage.png", p17a + p17b + plot_layout(widths = c(1, 1.8)), width = 13, height = 4.8, dpi = 150, bg = "white")

# ---- fig09: tracks, binary, states at one locus ----------------------------------------
# the interferon-receptor genes (also the locus of the CUT&Tag tracks in "CUT&RUN and CUT&Tag")
reg <- GRanges("chr21", IRanges(33250001, 33550000))
cat("\n== fig09 locus:", as.character(reg), "\n")
bw <- function(m) { x <- import(file.path(D, sprintf("gm12878_%s_fc.bigWig", tolower(m))), which = reg, as = "RleList")[["chr21"]]
  v <- as.numeric(x[start(reg):end(reg)]); bb <- (seq_along(v) - 1) %/% 1000
  data.frame(mark = m, x = start(reg) + tapply(seq_along(v), bb, mean), y = tapply(v, bb, mean)) }
tdf <- do.call(rbind, lapply(marks, bw)); tdf$mark <- factor(tdf$mark, levels = marks)
xl <- c(start(reg), end(reg))
pT <- ggplot(tdf, aes(x, y, fill = mark)) + geom_area() + facet_grid(mark ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = mcol, guide = "none") + scale_x_continuous(limits = xl, expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "ENCODE fold change over input") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 11), strip.placement = "outside",
        axis.text.x = element_blank(), axis.text.y = element_text(size = 8), panel.grid = element_blank())
ib <- which(overlapsAny(bins, reg))
bd <- do.call(rbind, lapply(marks, function(m) data.frame(mark = factor(m, levels = rev(marks)), x = start(bins)[ib], on = b[ib, m] == 1)))
pB <- ggplot(bd[bd$on, ], aes(x + 100, mark, fill = mark)) + geom_tile(width = 200, height = 0.8) +
  scale_fill_manual(values = mcol, guide = "none") + scale_x_continuous(limits = xl, expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "binarised (200-bp bins)") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.text.y = element_text(size = 9))
sd <- data.frame(x = start(bins)[ib] + 100, ours = bins$oursN[ib], rm = bins$rm[ib])
pS <- ggplot(sd) + geom_tile(aes(x, 2, fill = ours), width = 200, height = 0.8) +
  geom_tile(aes(x, 1, fill = rm), width = 200, height = 0.8) +
  scale_fill_manual(values = c(ncol8, rmcol), breaks = N8, labels = N8, name = "our states") +
  scale_y_continuous(breaks = 1:2, labels = c("Roadmap 18", "ours 8")) + scale_x_continuous(limits = xl, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "right",
                                   legend.text = element_text(size = 10), legend.key.size = unit(0.4, "cm"))
gx <- subsetByOverlaps(gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"], reg)
gd <- data.frame(s = pmax(start(gx), xl[1]), e = pmin(end(gx), xl[2]), g = paste0(gx$gene_name, ifelse(strand(gx) == "+", " →", " ←")),
                 row = rep(1:2, length.out = length(gx)))
pG <- ggplot(gd) + geom_segment(aes(x = s, xend = e, y = -row, yend = -row), linewidth = 2.5, colour = "#1b2a3a") +
  geom_text(aes(x = (s + e) / 2, y = -row - 0.6, label = g), size = 3.5) +
  scale_x_continuous(limits = xl, expand = c(0, 0), labels = function(v) sprintf("%.2f Mb", v / 1e6)) +
  scale_y_continuous(limits = c(-3, -0.5), breaks = NULL) + labs(x = "chr21", y = NULL) +
  coord_cartesian(clip = "off") + theme(panel.grid = element_blank())
ggsave("fig09-tracks-states.png", pT / pB / pS / pG + plot_layout(heights = c(6, 2.2, 1, 1.2), guides = "collect"),
       width = 12.5, height = 7, dpi = 150, bg = "white")
cat("genes:", unique(gx$gene_name), "\n\nDone.\n")
