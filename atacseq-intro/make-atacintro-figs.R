# Real-data figures for atacseq-intro/lec.qmd: GM12878 ATAC-seq (ENCODE ENCSR095QNB,
# replicate 1), reads from chr21 and chrM.
#
# Run from this folder:  Rscript make-atacintro-figs.R
# Needs: data/download.sh atacseq-intro (reads, reference, bowtie2 index, GTF,
#        H3K4me3/H3K27ac bigWigs, ATAC peaks), and the outputs of the command-line
#        part of the exercise (C1-C5) in the folder ATAC_WORK (default ~/epi/atac):
#        rep1.idxstats.txt, rep1.fraglen.txt, rep1.fragments.bed
#
# Writes fig08-fraglen.png, fig10-periodicity.png, fig11-tn5-bias.png,
# fig13-chrm-fraction.png, fig15-tss-profile.png, fig16-nfr-mono-tss.png,
# fig17-browser.png, and prints the numbers quoted on the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Biostrings)
  library(ggplot2); library(patchwork)
})
set.seed(1)
D <- "../data"
W <- path.expand(Sys.getenv("ATAC_WORK", "~/epi/atac"))
INK <- "#1b2a3a"; GREY <- "#5b6b7a"; GREEN <- "#2e8b57"; BLUE <- "#2a6f97"; RED <- "#c0392b"
theme_set(theme_minimal(base_size = 15) +
          theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15)))

# ---- inputs -------------------------------------------------------------------
fl <- read.table(file.path(W, "rep1.fraglen.txt"), col.names = c("len", "n"))
fl <- merge(data.frame(len = 1:2000), fl, all.x = TRUE); fl$n[is.na(fl$n)] <- 0
frag <- read.table(file.path(W, "rep1.fragments.bed"), col.names = c("chr", "s", "e", "len"))
fr <- GRanges(frag$chr, IRanges(frag$s + 1, frag$e))          # BED -> 1-based
cat("fragments on chr21:", length(fr), "\n")

# ---- fig08: fragment-length histogram with log-scale inset --------------------
mode_in <- function(lo, hi) { x <- fl[fl$len >= lo & fl$len <= hi, ]; x$len[which.max(x$n)] }
smooth <- function(v, k) as.numeric(stats::filter(v, rep(1 / k, k), sides = 2))
fl$sm <- smooth(fl$n, 11)
mode_sm <- function(lo, hi) { x <- fl[fl$len >= lo & fl$len <= hi, ]; x$len[which.max(x$sm)] }
fl$sm21 <- smooth(fl$n, 21)
mode_sm21 <- function(lo, hi) { x <- fl[fl$len >= lo & fl$len <= hi, ]; x$len[which.max(x$sm21)] }
modes <- c(nfr = mode_in(30, 120), mono = mode_sm(150, 280), di = mode_sm(300, 480), tri = mode_sm21(520, 700))
print(modes)
tot <- sum(fl$n)
cat(sprintf("fraction < 100 bp: %.1f%%, 180-250 bp: %.1f%%, > 500 bp: %.1f%%, median %d\n",
            100 * sum(fl$n[fl$len < 100]) / tot, 100 * sum(fl$n[fl$len >= 180 & fl$len <= 250]) / tot,
            100 * sum(fl$n[fl$len > 500]) / tot, fl$len[which(cumsum(fl$n) >= tot / 2)[1]]))
lab <- data.frame(x = modes[1:3], m = modes[1:3], y = fl$sm[modes[1:3]] + c(0, 120, 120),
                  t = c("nucleosome-free", "mono", "di"))     # the tri peak is too weak to label
lab$y[1] <- fl$n[modes[1]]; lab$x[1] <- modes[1] + 40
p_main <- ggplot(fl[fl$len <= 1000, ], aes(len, n / 1000)) +
  geom_area(fill = GREEN, alpha = 0.85) +
  geom_text(data = lab, aes(x, y / 1000, label = paste0(t, "\n", m, " bp")), vjust = -0.15, size = 4.6,
            lineheight = 0.9, colour = INK) +
  scale_x_continuous(breaks = seq(0, 1000, 100), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(c(0, 0.18))) +
  labs(x = "fragment length (bp)", y = "fragments (thousands)")
p_log <- ggplot(fl[fl$len <= 1000 & fl$n > 0, ], aes(len, n)) +
  geom_line(colour = GREEN, linewidth = 0.6) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
  labs(x = NULL, y = NULL, title = "log scale") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = "grey80"),
        plot.title = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank())
p08 <- p_main + inset_element(p_log, left = 0.45, bottom = 0.35, right = 0.99, top = 0.99)
ggsave("fig08-fraglen.png", p08, width = 11, height = 5, dpi = 150, bg = "white")

# ---- fig10: the 10.5-bp periodicity --------------------------------------------
fl$trend <- smooth(fl$n, 31)
per <- fl[fl$len >= 40 & fl$len <= 200, ]
per$ratio <- per$n / per$trend
r <- per$ratio
pk <- per$len[which(diff(sign(diff(r))) == -2) + 1]
pk <- pk[sapply(pk, function(p) r[per$len == p] == max(r[abs(per$len - p) <= 4]) && r[per$len == p] > 1)]   # highest within ±4 bp, above trend
pk_short <- pk[pk <= 140]
cat("periodicity maxima:", pk, "\n")
cat(sprintf("maxima 40-140 bp: %s; mean spacing %.1f bp\n", paste(pk_short, collapse = " "), (max(pk_short) - min(pk_short)) / (length(pk_short) - 1)))
cat(sprintf("all maxima 40-200: mean spacing %.1f bp\n", (max(pk) - min(pk)) / (length(pk) - 1)))
sp <- spec.pgram(per$ratio[per$len <= 150] - 1, taper = 0.1, plot = FALSE, detrend = TRUE)
cat(sprintf("dominant period (periodogram, 40-150 bp): %.1f bp\n", 1 / sp$freq[which.max(sp$spec)]))
p10 <- ggplot(per, aes(len, ratio)) +
  geom_hline(yintercept = 1, colour = "grey60") +
  geom_vline(xintercept = pk, colour = GREY, linetype = 3) +
  geom_line(colour = GREEN, linewidth = 0.9) +
  scale_x_continuous(breaks = seq(40, 200, 10), expand = c(0.01, 0)) +
  labs(x = "fragment length (bp)", y = "counts / trend")
ggsave("fig10-periodicity.png", p10, width = 11, height = 2.6, dpi = 150, bg = "white")

# ---- fig11: Tn5 sequence preference around insertion sites ---------------------
fa <- readDNAStringSet(file.path(D, "hg38_chr21_chrM.fa"))
names(fa) <- sub(" .*", "", names(fa))
chr21 <- fa[["chr21"]]
ins_left <- sample(start(fr), 200000)                          # first base of + strand reads
v <- Views(chr21, start = ins_left - 10, end = ins_left + 19)
m <- consensusMatrix(DNAStringSet(v), as.prob = TRUE)[c("A", "C", "G", "T"), ]
bg <- alphabetFrequency(chr21)[c("A", "C", "G", "T")]; bg <- bg / sum(bg)
bias <- data.frame(pos = rep(-10:19, each = 4), base = rownames(m), f = as.vector(m))
cat("genome base composition:", round(bg, 3), "\n")
cat(sprintf("max deviation from genome average: %.3f\n", max(abs(m - bg))))
basecol <- c(A = "#2e8b57", C = "#2a6f97", G = "#d4a017", T = "#c0392b")
p11 <- ggplot(bias, aes(pos, f, colour = base)) +
  annotate("rect", xmin = -0.5, xmax = 8.5, ymin = -Inf, ymax = Inf, fill = "#e9c46a", alpha = 0.25) +
  annotate("text", x = 4, y = max(bias$f) + 0.03, label = "9-bp duplicated stretch", size = 4.5, colour = INK) +
  geom_hline(data = data.frame(base = names(bg), f = bg), aes(yintercept = f, colour = base),
             linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  geom_vline(xintercept = 4, colour = GREY, linetype = 3) +
  scale_colour_manual(values = basecol) +
  scale_x_continuous(breaks = seq(-10, 19, 2)) +
  labs(x = "position relative to the 5' end of the + strand read (0 = first base)",
       y = "base frequency", colour = NULL) +
  theme(legend.position = "right")
ggsave("fig11-tn5-bias.png", p11, width = 11, height = 4.2, dpi = 150, bg = "white")

# ---- fig13: chrM fraction ------------------------------------------------------
idx <- read.table(file.path(W, "rep1.idxstats.txt"), col.names = c("chr", "len", "mapped", "unmapped"))
n21 <- idx$mapped[idx$chr == "chr21"]; nM <- idx$mapped[idx$chr == "chrM"]
f21 <- 46709983 / 3.1e9                     # chr21 share of the (haploid, nuclear) genome length
sub <- nM / (n21 + nM); est <- nM / (nM + n21 / f21)
cat(sprintf("chr21 %d, chrM %d mapped reads: chrM %.1f%% of subset, %.1f%% estimated genome-wide (f21 = %.4f)\n",
            n21, nM, 100 * sub, 100 * est, f21))
cm <- data.frame(what = factor(c("our files\n(chr21 + chrM only)", "scaled to the\nwhole genome",
                                 "ENCODE QC\n(all reads)"), levels = c("our files\n(chr21 + chrM only)",
                                 "scaled to the\nwhole genome", "ENCODE QC\n(all reads)")),
                 pct = c(100 * sub, 100 * est, 13.9), col = c("a", "b", "b"))
p13 <- ggplot(cm, aes(what, pct, fill = col)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 5.5) +
  scale_fill_manual(values = c(a = "#b9c6d3", b = RED), guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = NULL, y = "reads on chrM (%)") +
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 14, colour = INK))
ggsave("fig13-chrm-fraction.png", p13, width = 9, height = 4.4, dpi = 150, bg = "white")

# ---- TSSs and profiles -----------------------------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
tx <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding" & seqnames(gtf) == "chr21"]
tss <- unique(granges(resize(tx, 1, fix = "start")))
cat("protein-coding TSSs on chr21:", length(tss), "\n")

# insertions per bp at -flank..+flank around oriented TSSs
profile <- function(pos, tss, flank = 2000) {
  pts <- GRanges("chr21", IRanges(pos, width = 1))
  win <- resize(tss, 2 * flank + 1, fix = "center")
  h <- findOverlaps(pts, win, ignore.strand = TRUE)
  d <- pos[queryHits(h)] - start(tss)[subjectHits(h)]
  d <- ifelse(as.character(strand(tss))[subjectHits(h)] == "-", -d, d)
  tabulate(d + flank + 1, 2 * flank + 1)
}
ins <- c(start(fr), end(fr))                    # both ends of each fragment (unshifted)
cat("insertion sites:", length(ins), "\n")
prof <- profile(ins, tss)
bgv <- mean(prof[c(1:100, 3902:4001)])
norm <- prof / bgv
bin10 <- tapply(norm, (seq_along(norm) - 1) %/% 10, mean)
tsse_centre <- mean(norm[2001 + (-5:5)])
tsse_max <- max(bin10)
cat(sprintf("TSS enrichment: background %.2f insertions/bp; centre (±5 bp) %.1f; max of 10-bp bins %.1f at %d bp\n",
            bgv, tsse_centre, tsse_max, (which.max(bin10) - 1) * 10 - 2000 + 5))
tp <- data.frame(d = -2000:2000, v = norm, s = smooth(norm, 11))
p15 <- ggplot(tp, aes(d, s)) +
  annotate("rect", xmin = c(-2000, 1900), xmax = c(-1900, 2000), ymin = -Inf, ymax = Inf, fill = BLUE, alpha = 0.2) +
  geom_hline(yintercept = c(1, 7), colour = c("grey50", RED), linetype = c(1, 2)) +
  annotate("text", x = 1950, y = 7, label = "ENCODE 'ideal' > 7", hjust = 1, vjust = -0.5, colour = RED, size = 4.5) +
  geom_line(colour = GREEN, linewidth = 1) +
  annotate("text", x = 250, y = tsse_max, hjust = 0, size = 5,
           label = sprintf("TSS enrichment = %.1f", tsse_max)) +
  annotate("text", x = -1880, y = 5, label = "flanks: 100 bp each", hjust = 0, colour = BLUE, size = 4.3) +
  scale_x_continuous(breaks = seq(-2000, 2000, 500)) +
  labs(x = "distance to TSS (bp, oriented by strand)", y = "insertions / flank mean")
ggsave("fig15-tss-profile.png", p15, width = 10, height = 4.6, dpi = 150, bg = "white")

# ---- fig16: nucleosome-free vs mono-nucleosome fragment centres ------------------
mid <- (start(fr) + end(fr)) %/% 2
w <- width(fr)
cls <- list("nucleosome-free (< 100 bp)" = w < 100, "mono-nucleosome (180-250 bp)" = w >= 180 & w <= 250)
pf <- do.call(rbind, lapply(names(cls), function(k) {
  p <- profile(mid[cls[[k]]], tss, flank = 1000)
  p <- p / mean(p[c(1:100, 1902:2001)])
  data.frame(cls = k, d = -1000:1000, v = smooth(p, 21))
}))
pf$cls <- factor(pf$cls, levels = names(cls))
mono <- pf[pf$cls == names(cls)[2], ]
cat(sprintf("mono-nucleosome centres: max downstream at %d bp, upstream at %d bp; at TSS %.2f\n",
            mono$d[mono$d > 0][which.max(mono$v[mono$d > 0])], mono$d[mono$d < 0][which.max(mono$v[mono$d < 0])],
            mono$v[mono$d == 0]))
nfr <- pf[pf$cls == names(cls)[1], ]
cat(sprintf("NFR centres: max at %d bp (%.1f x flanks)\n", nfr$d[which.max(nfr$v)], max(nfr$v, na.rm = TRUE)))
p16 <- ggplot(pf, aes(d, v, colour = cls)) +
  geom_hline(yintercept = 1, colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 3) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~cls, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(GREEN, BLUE), guide = "none") +
  scale_x_continuous(breaks = seq(-1000, 1000, 500)) +
  labs(x = "fragment centre, distance to TSS (bp)", y = "fragments / flank mean") +
  theme(strip.text = element_text(face = "bold", size = 14))
ggsave("fig16-nfr-mono-tss.png", p16, width = 11, height = 4.2, dpi = 150, bg = "white")

# ---- fig17: browser view, ATAC vs H3K4me3 and H3K27ac at ITGB2 --------------------
region <- GRanges("chr21", IRanges(44870001, 44950000))
xl <- c(start(region), end(region))
cov <- coverage(fr)[["chr21"]] / length(fr) * 1e6           # fragment coverage, CPM
bins <- function(v, nbin = 800) {
  b <- cut(seq_along(v), nbin, labels = FALSE)
  data.frame(x = xl[1] + tapply(seq_along(v), b, mean) - 1, y = tapply(v, b, mean))
}
atac <- bins(as.numeric(cov[xl[1]:xl[2]])); atac$track <- "ATAC-seq (CPM)"
bw <- function(m, lab) {
  x <- import(file.path(D, sprintf("gm12878_%s_fc.bigWig", m)), which = region, as = "RleList")[["chr21"]]
  d <- bins(as.numeric(x[xl[1]:xl[2]])); d$track <- lab; d
}
tr <- rbind(atac, bw("h3k4me3", "H3K4me3 (FC)"), bw("h3k27ac", "H3K27ac (FC)"))
tr$track <- factor(tr$track, levels = c("ATAC-seq (CPM)", "H3K4me3 (FC)", "H3K27ac (FC)"))
cols <- c("ATAC-seq (CPM)" = GREEN, "H3K4me3 (FC)" = "#c0392b", "H3K27ac (FC)" = "#e67e22")
ptr <- ggplot(tr, aes(x, y, fill = track)) + geom_area() +
  facet_grid(track ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = xl, expand = c(0, 0)) + scale_y_continuous(n.breaks = 3) +
  labs(x = NULL, y = NULL) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
        strip.placement = "outside", axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 9))
pk <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
pk <- subsetByOverlaps(pk, region)
cat("ENCODE ATAC IDR peaks in the browser window:", length(pk), "\n")
pg <- gtf[gtf$type == "transcript" & gtf$gene_type == "protein_coding"]
pg <- subsetByOverlaps(pg, region); pg <- pg[order(-width(pg))]; pg <- pg[!duplicated(pg$gene_name)]
pg$row <- seq_along(pg) %% 2 + 1
ex <- gtf[gtf$type == "exon" & gtf$transcript_id %in% pg$transcript_id]
exd <- data.frame(s = start(ex), e = end(ex), row = pg$row[match(ex$transcript_id, pg$transcript_id)])
txd <- data.frame(s = pmax(start(pg), xl[1]), e = pmin(end(pg), xl[2]), row = pg$row,
                  g = paste0(pg$gene_name, ifelse(strand(pg) == "+", " →", " ←")))
txd$lx <- pmin(pmax((txd$s + txd$e) / 2, xl[1] + 4000), xl[2] - 4000)   # keep labels inside the window
pgene <- ggplot() +
  geom_rect(data = data.frame(s = start(pk), e = end(pk)), aes(xmin = s, xmax = e, ymin = 0.2, ymax = 0.6), fill = GREEN) +
  annotate("text", x = xl[1] + 300, y = 0.75, label = "ATAC peaks (ENCODE IDR)", hjust = 0, vjust = 0, size = 3.8, colour = GREEN) +
  geom_segment(data = txd, aes(x = s, xend = e, y = -row, yend = -row), colour = INK) +
  geom_rect(data = exd, aes(xmin = s, xmax = e, ymin = -row - 0.22, ymax = -row + 0.22), fill = INK) +
  geom_text(data = txd, aes(x = lx, y = -row - 0.5, label = g), size = 4.2, vjust = 1) +
  scale_x_continuous(limits = xl, expand = c(0, 0), breaks = seq(44.88e6, 44.94e6, 0.02e6),
                     labels = function(b) sprintf("%.2f Mb", b / 1e6)) +
  scale_y_continuous(breaks = NULL, limits = c(-3.1, 1.3)) +
  labs(x = "chr21", y = NULL) + theme(panel.grid = element_blank())
p17 <- ptr / pgene + plot_layout(heights = c(4, 1.5))
ggsave("fig17-browser.png", p17, width = 12, height = 5.6, dpi = 150, bg = "white")
cat("genes in window:", paste(txd$g, collapse = ", "), "\n")
