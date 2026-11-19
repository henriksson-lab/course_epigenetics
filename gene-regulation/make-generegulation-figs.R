# Real-data figures for gene-regulation/lec.qmd.
# Run from this folder:  Rscript make-generegulation-figs.R   (about 12 minutes)
# Needs: data/download.sh gene-regulation  (peaks, bigWigs, GTF, blacklist and the
#        bgzipped, indexed hg38 genome hg38.fa.bgz + .fai + .gzi)
#        R packages GenomicRanges, rtracklayer, Rsamtools, Biostrings, ggplot2,
#        patchwork, ggrepel, scales, JASPAR2024, TFBSTools, RSQLite, motifmatchr,
#        universalmotif.
#
# Writes fig04-cpg-promoters.png, fig09-bcell-enhancers.png, fig10-super-enhancers.png,
# fig16-ctcf-logo.png, fig17-chance-matches.png, fig18-ctcf-motif-in-peaks.png,
# fig19-ctcf-centrality.png, fig20-unbound-matches.png, fig24-bcell-motifs.png
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Rsamtools); library(Biostrings)
  library(ggplot2); library(patchwork)
  library(JASPAR2024); library(TFBSTools); library(RSQLite); library(motifmatchr)
  library(universalmotif)
})
D <- "../data"
set.seed(2026)
col <- c(H3K4me3 = "#c0392b", H3K27ac = "#e67e22", H3K4me1 = "#d4a017", CTCF = "#8c564b",
         ATAC = "#2e8b57")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
main_chr <- paste0("chr", c(1:22, "X"))

fa  <- FaFile(file.path(D, "hg38.fa.bgz"))           # bgzipped + faidx: random access
jdb <- dbConnect(SQLite(), db(JASPAR2024()))
ctcf_m <- getMatrixByID(jdb, "MA0139.2")               # CTCF, JASPAR 2024 CORE

gtf   <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
genes <- genes[seqnames(genes) %in% main_chr]
tss   <- resize(genes, 1, fix = "start")

np <- function(f) {
  p <- import(file.path(D, f), format = "narrowPeak")
  p[seqnames(p) %in% main_chr]
}
summits <- function(p) GRanges(seqnames(p), IRanges(start(p) + p$peak, width = 1))

# ---- helpers ----------------------------------------------------------------
# Random genomic windows (main chromosomes, outside the blacklist, no N).
blacklist <- import(file.path(D, "hg38_blacklist.bed.gz"))
random_windows <- function(n, width) {
  len <- seqlengths(fa)[main_chr]
  chr <- sample(main_chr, n, replace = TRUE, prob = len)
  st  <- floor(runif(n, 1, len[chr] - width))
  gr  <- GRanges(chr, IRanges(st, width = width))
  gr  <- gr[!overlapsAny(gr, blacklist)]
  s   <- getSeq(fa, gr)
  ok  <- letterFrequency(s, "N")[, 1] == 0
  list(gr = gr[ok], seq = s[ok])
}
gc_of <- function(s) letterFrequency(s, "GC", as.prob = TRUE)[, 1]
# Pick background sequences with the same GC distribution as the target (2% bins).
gc_match <- function(target, pool) {
  tb <- cut(gc_of(target), seq(0, 1, 0.02), include.lowest = TRUE)
  pb <- cut(gc_of(pool), seq(0, 1, 0.02), include.lowest = TRUE)
  idx <- unlist(lapply(levels(tb), function(b) {
    need <- sum(tb == b); have <- which(pb == b)
    if (need == 0 || length(have) == 0) return(integer(0))
    have[sample.int(length(have), need, replace = need > length(have))]
  }))
  pool[idx]
}
# Fraction of sequences with at least one match (p < 1e-4, both strands).
has_match <- function(m, s, p = 1e-4)
  as.matrix(motifMatches(matchMotifs(m, s, p.cutoff = p, bg = "even")))
# Positions of all matches in many equal-width sequences. Scanning one long
# sequence is much faster than many short ones, so join them with runs of N.
scan_positions <- function(m, s, p = 1e-4, gap = 50) {
  w <- unique(width(s)); stopifnot(length(w) == 1)
  big <- unlist(xscat(s, DNAString(strrep("N", gap))))
  h <- matchMotifs(m, big, p.cutoff = p, bg = "even", out = "positions")[[1]][[1]]
  i <- (start(h) - 1) %/% (w + gap) + 1
  rel <- start(h) - (i - 1) * (w + gap)
  keep <- rel + width(h) - 1 <= w
  data.frame(seq = i[keep], start = rel[keep], end = rel[keep] + width(h)[keep] - 1,
             strand = mcols(h)$strand[keep], score = mcols(h)$score[keep])
}

# ---- fig04: CpG-island and non-CpG promoters --------------------------------
prom <- promoters(genes, upstream = 500, downstream = 500)
ps <- getSeq(fa, prom)
nC <- letterFrequency(ps, "C")[, 1]; nG <- letterFrequency(ps, "G")[, 1]
nCG <- dinucleotideFrequency(ps)[, "CG"]
oe <- nCG * width(ps) / pmax(nC * nG, 1)
gc <- (nC + nG) / width(ps)
k4 <- np("gm12878_h3k4me3_peaks.bed.gz")
marked <- overlapsAny(prom, k4)
cgi <- oe > 0.6 & gc > 0.55
cat(sprintf("fig04: %d protein-coding promoters; CpG-island-like %.0f%%; H3K4me3 peak in GM12878: CGI %.0f%%, non-CGI %.0f%%\n",
            length(prom), 100 * mean(cgi), 100 * mean(marked[cgi]), 100 * mean(marked[!cgi])))
d4 <- data.frame(oe = pmin(oe, 1.3),
                 status = factor(ifelse(marked, "H3K4me3 peak in GM12878", "no H3K4me3 peak"),
                                 levels = c("no H3K4me3 peak", "H3K4me3 peak in GM12878")))
p4 <- ggplot(d4, aes(oe, fill = status)) +
  geom_histogram(binwidth = 0.025, boundary = 0, colour = NA) +
  geom_vline(xintercept = 0.6, linetype = 2, colour = "grey30") +
  annotate("text", x = 0, y = Inf, vjust = 1.3, hjust = 0, size = 5,
           label = sprintf("non-CpG promoters: %.0f%%\n%.0f%% with H3K4me3", 100 * mean(!cgi), 100 * mean(marked[!cgi]))) +
  annotate("text", x = 1.3, y = Inf, vjust = 1.3, hjust = 1, size = 5,
           label = sprintf("CpG-island promoters: %.0f%%\n(o/e > 0.6, GC > 55%%)\n%.0f%% with H3K4me3", 100 * mean(cgi), 100 * mean(marked[cgi]))) +
  scale_fill_manual(values = c("grey75", unname(col["H3K4me3"])), name = NULL) +
  labs(x = "CpG observed / expected, TSS ± 500 bp", y = "protein-coding genes") +
  theme(legend.position = "bottom")
ggsave("fig04-cpg-promoters.png", p4, width = 11, height = 4.6, dpi = 150, bg = "white")

# ---- fig09: an enhancer cluster at PAX5 / ZCCHC7 ------------------------------
reg <- GRanges("chr9", IRanges(36880000, 37380000))
atac <- np("gm12878_atac_peaks.bed.gz")
bw_track <- function(f, name) {
  x <- import(file.path(D, f), which = reg, as = "NumericList")[[1]]
  pos <- seq(start(reg), end(reg))
  bins <- split(as.numeric(x), ceiling(seq_along(x) / 500))   # 500-bp bins
  data.frame(track = name, x = pos[seq(1, length(pos), 500)][seq_along(bins)] + 250,
             y = vapply(bins, max, 0))
}
tr <- rbind(bw_track("gm12878_h3k27ac_fc.bigWig", "H3K27ac"),
            bw_track("gm12878_h3k4me1_fc.bigWig", "H3K4me1"),
            bw_track("gm12878_h3k4me3_fc.bigWig", "H3K4me3"),
            bw_track("gm12878_ctcf_fc.bigWig", "CTCF"))
tr$track <- factor(tr$track, levels = c("H3K27ac", "H3K4me1", "H3K4me3", "CTCF"))
a9 <- subsetByOverlaps(atac, reg)
g9 <- subsetByOverlaps(genes, reg)
ex9 <- subsetByOverlaps(gtf[gtf$type == "exon" & gtf$gene_name %in% g9$gene_name], reg)
ex9 <- unlist(reduce(split(ex9, ex9$gene_name)))           # union of exons per gene
ex9$gene_name <- names(ex9); names(ex9) <- NULL
cat(sprintf("fig09: %d ATAC peaks and %d H3K27ac peaks in %s\n", length(a9),
            length(subsetByOverlaps(np("gm12878_h3k27ac_peaks.bed.gz"), reg)), as.character(reg)))
xl <- c(start(reg), end(reg)) / 1e6
pt <- ggplot(tr, aes(x / 1e6, y, fill = track)) +
  geom_area(colour = NA) +
  facet_grid(track ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = col, guide = "none") +
  scale_x_continuous(limits = xl, expand = c(0, 0)) +
  labs(x = NULL, y = "fold change over input") +
  theme(strip.text.y.left = element_text(angle = 0, face = "bold", size = 14),
        strip.placement = "outside", axis.text.x = element_blank(),
        panel.grid.major.x = element_blank())
pa <- ggplot(as.data.frame(a9)) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 1), fill = col["ATAC"]) +
  scale_x_continuous(limits = xl, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + ggtitle(NULL) +
  annotate("text", x = xl[1] + 0.002, y = 1.5, label = "ATAC peaks", hjust = 0, size = 4.6,
           fontface = "bold", colour = col["ATAC"]) +
  coord_cartesian(ylim = c(0, 2), clip = "off") + theme_void()
gdf <- as.data.frame(g9); gdf$mid <- (pmax(gdf$start, start(reg)) + pmin(gdf$end, end(reg))) / 2
gdf$row <- seq_len(nrow(gdf)) %% 2
edf <- as.data.frame(ex9); edf$row <- gdf$row[match(edf$gene_name, gdf$gene_name)]
pg <- ggplot() +
  geom_segment(data = gdf, aes(x = pmax(start, start(reg)) / 1e6, xend = pmin(end, end(reg)) / 1e6,
                               y = row, yend = row), linewidth = 0.6) +
  geom_rect(data = edf, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = row - 0.18, ymax = row + 0.18),
            fill = "#1b2a3a") +
  geom_text(data = gdf, aes(x = mid / 1e6, y = row - 0.5,
                            label = paste0(gene_name, ifelse(strand == "+", " →", " ←"))),
            size = 4.6, fontface = "italic") +
  scale_x_continuous(limits = xl, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.8, 1.3)) +
  labs(x = "chr9 position (Mb, hg38)", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
p9 <- pt / pa / pg + plot_layout(heights = c(4, 0.35, 1))
ggsave("fig09-bcell-enhancers.png", p9, width = 12, height = 5.8, dpi = 150, bg = "white")

# ---- fig10: "super-enhancers" from stitched H3K27ac (ROSE-like) --------------
k27 <- np("gm12878_h3k27ac_peaks.bed.gz")
k27d <- k27[!overlapsAny(k27, tss + 2500)]                 # distal peaks
st <- reduce(k27d, min.gapwidth = 12500)                   # stitch within 12.5 kb
bw <- BigWigFile(file.path(D, "gm12878_h3k27ac_fc.bigWig"))
st$signal <- unlist(summary(bw, st, type = "mean"))$score * width(st)
st <- st[order(st$signal)]                                  # increasing, as in ROSE
st$gene <- tss$gene_name[nearest(st, tss)]
ys <- st$signal / max(st$signal); xs <- seq_along(ys) / length(ys)
cut_i <- which.min(ys - xs)                                 # tangent of slope 1
se <- st$signal >= st$signal[cut_i]
cat(sprintf("fig10: %d distal H3K27ac peaks -> %d stitched regions; %d above the cutoff (%.1f%%) hold %.0f%% of the signal\n",
            length(k27d), length(st), sum(se), 100 * mean(se), 100 * sum(st$signal[se]) / sum(st$signal)))
top <- rev(tail(seq_along(st), 12))
print(data.frame(region = as.character(st[top]), kb = round(width(st[top]) / 1e3), gene = st$gene[top]))
d10 <- data.frame(rank = seq_along(st), signal = st$signal / 1e3, se = se)
lab_genes <- c("ETV6", "PLCG2", "BCL2", "INPP5D", "TRAF3", "PAX5")   # B-cell genes near SEs
lab <- d10[st$gene %in% lab_genes & se, ]
lab$gene <- st$gene[lab$rank]
lab <- lab[!duplicated(lab$gene), ]
p10 <- ggplot(d10, aes(rank, signal, colour = se)) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = st$signal[cut_i] / 1e3, linetype = 2, colour = "grey40") +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 4.6, colour = "#1b2a3a",
                           fontface = "italic", nudge_x = -2500, direction = "y", hjust = 1,
                           segment.colour = "grey50", min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = c("grey60", unname(col["H3K27ac"])), guide = "none") +
  annotate("text", x = 1, y = st$signal[cut_i] / 1e3, vjust = -0.6, hjust = 0, size = 4.8,
           label = sprintf("%d \"super-enhancers\" (%.0f%% of regions, %.0f%% of signal)",
                           sum(se), 100 * mean(se), 100 * sum(st$signal[se]) / sum(st$signal))) +
  labs(x = "stitched distal H3K27ac regions, ranked by signal",
       y = "H3K27ac signal (fold change × kb)")
ggsave("fig10-super-enhancers.png", p10, width = 11, height = 4.8, dpi = 150, bg = "white")

# ---- fig16: the CTCF motif as a sequence logo --------------------------------
um <- convert_motifs(ctcf_m, "universalmotif")
cnt <- Matrix(ctcf_m)
fr <- (cnt + 0.25) / colSums(cnt + 0.25)                    # pseudocount 1 per column, 0.25 per base
ic <- 2 + colSums(fr * log2(fr))
cat(sprintf("fig16: CTCF MA0139.2, %d sites, %d bp, total IC %.1f bits (columns %s)\n",
            sum(cnt[, 1]), ncol(cnt), sum(ic), paste(round(ic, 2), collapse = " ")))
p16 <- view_motifs(um, use.type = "ICM") +
  labs(title = sprintf("CTCF (JASPAR MA0139.2): %d sites, %.1f bits in total", sum(cnt[, 1]), sum(ic))) +
  theme(plot.title = element_text(face = "bold", size = 16), text = element_text(size = 15))
ggsave("fig16-ctcf-logo.png", p16, width = 11, height = 3.6, dpi = 150, bg = "white")

# ---- CTCF peaks: sequences around the summit and three backgrounds ----------
ctcf <- np("gm12878_ctcf_peaks.bed.gz")
s200 <- getSeq(fa, resize(summits(ctcf), 200, fix = "center"))
pool <- random_windows(600000, 200)
bg_gc  <- gc_match(s200, pool$seq)
bg_rand <- pool$seq[sample.int(length(pool$seq), length(s200))]
bg_shuf <- shuffle_sequences(s200, k = 2, rng.seed = 1)   # same dinucleotides
fr18 <- c(`CTCF peaks` = mean(has_match(ctcf_m, s200)),
          `GC-matched\ngenome` = mean(has_match(ctcf_m, bg_gc)),
          `random\ngenome` = mean(has_match(ctcf_m, bg_rand)),
          `shuffled\npeaks` = mean(has_match(ctcf_m, bg_shuf)))
cat(sprintf("fig18: %d CTCF peaks; GC peaks %.3f, GC-matched %.3f, random %.3f\n", length(ctcf),
            mean(gc_of(s200)), mean(gc_of(bg_gc)), mean(gc_of(bg_rand))))
print(round(fr18, 3))
hit <- has_match(ctcf_m, s200)[, 1]
q <- cut(rank(ctcf$signalValue, ties.method = "first"), 4,
         labels = c("weakest 25%", "2nd", "3rd", "strongest 25%"))
fq <- tapply(hit, q, mean); print(round(fq, 3))
d18a <- data.frame(set = factor(names(fr18), levels = names(fr18)), f = fr18)
p18a <- ggplot(d18a, aes(set, f, fill = set == "CTCF peaks")) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * f)), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c("grey65", unname(col["CTCF"])), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "sequences with a CTCF match", title = "200 bp around the summit") +
  theme(axis.text.x = element_text(size = 13))
d18b <- data.frame(q = factor(names(fq), levels = names(fq)), f = as.numeric(fq))
p18b <- ggplot(d18b, aes(q, f)) + geom_col(width = 0.7, fill = col["CTCF"]) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * f)), vjust = -0.4, size = 5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "peaks ranked by signalValue", y = NULL, title = "by peak strength") +
  theme(axis.text.x = element_text(size = 13))
ggsave("fig18-ctcf-motif-in-peaks.png", p18a + p18b + plot_layout(widths = c(1.15, 1)),
       width = 12, height = 4.6, dpi = 150, bg = "white")

# ---- fig19: where in the peak is the motif? ----------------------------------
s500 <- getSeq(fa, resize(summits(ctcf), 500, fix = "center"))
pos <- scan_positions(ctcf_m, s500)
pos <- pos[order(-pos$score), ]; best <- pos[!duplicated(pos$seq), ]
best$centre <- (best$start + best$end) / 2 - 250.5
bgp <- scan_positions(ctcf_m, gc_match(s500, random_windows(400000, 500)$seq))
bgp <- bgp[order(-bgp$score), ]; bgb <- bgp[!duplicated(bgp$seq), ]
bgb$centre <- (bgb$start + bgb$end) / 2 - 250.5
cat(sprintf("fig19: %d of %d peaks (%.1f%%) have a match within ±250 bp; best match within 50 bp of the summit: %.1f%%; median |distance| %.0f bp; GC-matched background with a match: %.1f%%\n",
            nrow(best), length(s500), 100 * nrow(best) / length(s500), 100 * mean(abs(best$centre) <= 50),
            median(abs(best$centre)), 100 * nrow(bgb) / length(s500)))
d19 <- rbind(data.frame(set = "CTCF peaks", x = best$centre),
             data.frame(set = "GC-matched background", x = bgb$centre))
p19 <- ggplot(d19, aes(x, fill = set)) +
  geom_histogram(binwidth = 10, boundary = 0, colour = NA, position = "identity", alpha = 0.85) +
  scale_fill_manual(values = c(unname(col["CTCF"]), "grey60"), name = NULL) +
  labs(x = "position of the best motif match relative to the peak summit (bp)",
       y = "sequences") + theme(legend.position = c(0.82, 0.8))
ggsave("fig19-ctcf-centrality.png", p19, width = 11, height = 4.4, dpi = 150, bg = "white")

# ---- fig17 + fig20: every CTCF match on chr21, bound or not -----------------
chr21 <- getSeq(fa, GRanges("chr21", IRanges(1, seqlengths(fa)["chr21"])))
acgt <- sum(letterFrequency(chr21, "ACGT"))
cnt21 <- sapply(c(1e-4, 1e-5, 1e-6), function(p)
  length(matchMotifs(ctcf_m, chr21, p.cutoff = p, bg = "even", out = "positions")[[1]][[1]]))
exp21 <- 2 * acgt * c(1e-4, 1e-5, 1e-6)
cat(sprintf("fig17: chr21 %.1f Mb of ACGT; matches at p<1e-4/1e-5/1e-6: %s; expected under uniform background: %s\n",
            acgt / 1e6, paste(cnt21, collapse = "/"), paste(round(exp21), collapse = "/")))
d17 <- data.frame(p = factor(rep(c("p < 1e-4", "p < 1e-5", "p < 1e-6"), 2), levels = c("p < 1e-4", "p < 1e-5", "p < 1e-6")),
                  what = rep(c("found on chr21", "expected by chance: 2 × L × p"), each = 3),
                  n = c(cnt21, exp21))
p17 <- ggplot(d17, aes(p, n, fill = what)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = formatC(round(n), big.mark = " ", format = "d")), position = position_dodge(0.75),
            vjust = -0.4, size = 4.6) +
  scale_fill_manual(values = c("grey65", unname(col["CTCF"])), name = NULL) +
  scale_y_log10(labels = function(x) formatC(x, big.mark = " ", format = "d"),
                expand = expansion(mult = c(0, 0.12))) +
  labs(x = "match threshold", y = "CTCF motif matches, chr21 (log scale)") +
  theme(legend.position = c(0.75, 0.8))
ggsave("fig17-chance-matches.png", p17, width = 11, height = 4.4, dpi = 150, bg = "white")

m21 <- matchMotifs(ctcf_m, chr21, p.cutoff = 1e-4, bg = "even", out = "positions")[[1]][[1]]
g21 <- GRanges("chr21", IRanges(start(m21), end(m21)), score = mcols(m21)$score)
c21 <- ctcf[seqnames(ctcf) == "chr21"]; a21 <- atac[seqnames(atac) == "chr21"]
inpk <- overlapsAny(g21, c21); inat <- overlapsAny(g21, a21)
sq <- cut(rank(g21$score, ties.method = "first"), 4, labels = c("weakest 25%", "2nd", "3rd", "best 25%"))
cat(sprintf("fig20: %d matches on chr21 (p<1e-4); %d CTCF peaks on chr21; matches in a CTCF peak %d (%.1f%%); in ATAC %.1f%%; P(peak | ATAC) %.1f%%, P(peak | not ATAC) %.2f%%; chr21 CTCF peaks with a match %.1f%%\n",
            length(g21), length(c21), sum(inpk), 100 * mean(inpk), 100 * mean(inat),
            100 * mean(inpk[inat]), 100 * mean(inpk[!inat]), 100 * mean(overlapsAny(c21, g21))))
print(round(tapply(inpk, sq, mean), 4))
d20a <- data.frame(g = factor(c("all matches", "in an ATAC peak", "not in an ATAC peak"),
                              levels = c("all matches", "in an ATAC peak", "not in an ATAC peak")),
                   f = c(mean(inpk), mean(inpk[inat]), mean(inpk[!inat])),
                   n = c(length(inpk), sum(inat), sum(!inat)))
p20a <- ggplot(d20a, aes(g, f)) + geom_col(width = 0.7, fill = col["CTCF"]) +
  geom_text(aes(label = sprintf("%.1f%%\nof %s", 100 * f, formatC(n, big.mark = " ", format = "d"))),
            vjust = -0.2, size = 4.6) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.5)) +
  labs(x = NULL, y = "matches inside a CTCF peak", title = "accessibility") +
  theme(axis.text.x = element_text(size = 13))
d20b <- data.frame(q = factor(levels(sq), levels = levels(sq)), f = as.numeric(tapply(inpk, sq, mean)))
p20b <- ggplot(d20b, aes(q, f)) + geom_col(width = 0.7, fill = col["CTCF"]) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * f)), vjust = -0.4, size = 4.6) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3)) +
  labs(x = "matches ranked by motif score", y = NULL, title = "motif score") +
  theme(axis.text.x = element_text(size = 13))
ggsave("fig20-unbound-matches.png", p20a + p20b, width = 12, height = 4.6, dpi = 150, bg = "white")

# ---- fig24: which motifs are enriched in B-cell open chromatin? -------------
ids <- c(`PU.1 (SPI1)` = "MA0080.7", EBF1 = "MA0154.5", PAX5 = "MA0014.4", IRF4 = "MA1419.2",
         `OCT2 (POU2F2)` = "MA0507.3", MEF2B = "MA0660.1", RUNX3 = "MA0684.3", CTCF = "MA0139.2",
         `HNF4A (liver)` = "MA0114.5", `FOXA1 (liver)` = "MA0148.5", `GATA1 (erythroid)` = "MA0035.5")
ms <- do.call(PFMatrixList, lapply(ids, function(i) getMatrixByID(jdb, i)))
names(ms) <- names(ids)
distal <- atac[!overlapsAny(summits(atac), tss + 2000)]
sa <- getSeq(fa, resize(summits(distal), 200, fix = "center"))
sb <- gc_match(sa, pool$seq)
fa_pk <- colMeans(has_match(ms, sa)); fa_bg <- colMeans(has_match(ms, sb))
d24 <- data.frame(tf = names(ids), peaks = fa_pk, bg = fa_bg, enr = fa_pk / fa_bg,
                  group = ifelse(grepl("liver|erythroid", names(ids)), "not expressed in B cells",
                                 ifelse(names(ids) == "CTCF", "insulator", "B-cell TF")))
cat(sprintf("fig24: %d ATAC peaks, %d distal (summit > 2 kb from a protein-coding TSS); GC %.3f vs background %.3f\n",
            length(atac), length(distal), mean(gc_of(sa)), mean(gc_of(sb))))
print(d24[order(-d24$enr), c("tf", "peaks", "bg", "enr")], digits = 3)
d24$tf <- factor(d24$tf, levels = d24$tf[order(d24$enr)])
p24 <- ggplot(d24, aes(enr, tf, fill = group)) + geom_col(width = 0.7) +
  geom_vline(xintercept = 1, colour = "grey30") +
  geom_text(aes(x = pmax(enr, 1), label = sprintf("%.0f%% vs %.0f%%", 100 * peaks, 100 * bg)),
            hjust = -0.08, size = 4.4) +
  scale_x_continuous(trans = "log2", breaks = c(0.5, 1, 2, 4), expand = expansion(mult = c(0.05, 0.25))) +
  scale_fill_manual(values = c(`B-cell TF` = "#2e8b57", insulator = unname(col["CTCF"]),
                               `not expressed in B cells` = "grey65"), name = NULL) +
  labs(x = "enrichment: fraction of ATAC peaks with a match / fraction of GC-matched background (log scale)",
       y = NULL) + theme(legend.position = "bottom")
ggsave("fig24-bcell-motifs.png", p24, width = 11, height = 5.4, dpi = 150, bg = "white")
