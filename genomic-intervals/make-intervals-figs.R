# Real-data figures for genomic-intervals/lec.qmd.
# Run from this folder:  Rscript make-intervals-figs.R
# Needs: data/download.sh genomic-intervals chromatin-intro  (peaks, GTF, CTCF bigWig)
#        bedtools on the PATH (for the shuffles), R packages GenomicRanges,
#        rtracklayer, ggplot2, patchwork.
#
# Writes fig10-browser-ctcf.png, fig16-distance-to-tss.png, fig18-shuffle-real.png
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(ggplot2); library(patchwork)
})
D <- "../data"
col <- c(H3K4me3 = "#c0392b", H3K27ac = "#e67e22", CTCF = "#8c564b")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))

peaks <- list(
  H3K4me3 = import(file.path(D, "gm12878_h3k4me3_peaks.bed.gz"), format = "narrowPeak"),
  H3K27ac = import(file.path(D, "gm12878_h3k27ac_peaks.bed.gz"), format = "narrowPeak"),
  CTCF    = import(file.path(D, "gm12878_ctcf_peaks.bed.gz"),    format = "narrowPeak"))

gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
tx  <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding"]
tss <- resize(tx, width = 1, fix = "start")
mcols(tss) <- mcols(tss)["gene_name"]
tss <- tss[!duplicated(paste(as.character(tss), tss$gene_name))]

# ---- fig16: distance from each peak to the nearest protein-coding TSS -------
dist <- do.call(rbind, lapply(names(peaks), function(m) {
  d <- mcols(distanceToNearest(peaks[[m]], tss, ignore.strand = TRUE))$distance
  data.frame(mark = m, d = d)
}))
dist$mark <- factor(dist$mark, levels = names(col))
frac <- aggregate(d ~ mark, dist, function(x) mean(x <= 1000))
print(frac)
lab <- data.frame(mark = frac$mark, txt = sprintf("%.0f%% within 1 kb", 100 * frac$d))
p16 <- ggplot(dist, aes(d + 1, fill = mark)) +
  geom_histogram(bins = 60, colour = NA) +
  geom_vline(xintercept = 1001, linetype = 2, colour = "grey40") +
  geom_label(data = lab, aes(x = 3e6, y = Inf, label = txt), vjust = 1.3, hjust = 1,
             size = 5, label.size = 0, fill = "white", inherit.aes = FALSE) +
  facet_wrap(~mark, ncol = 3) +
  scale_x_log10(breaks = c(1, 1e2, 1e4, 1e6), labels = c("0", "100 bp", "10 kb", "1 Mb")) +
  scale_fill_manual(values = col, guide = "none") +
  labs(x = "distance to nearest protein-coding TSS (log scale)", y = "number of peaks")
ggsave("fig16-distance-to-tss.png", p16, width = 12, height = 4.2, dpi = 150, bg = "white")

# ---- fig18: CTCF peaks overlapping H3K27ac, observed vs shuffled -------------
tmp <- tempfile("gi"); dir.create(tmp)
export(peaks$CTCF, file.path(tmp, "ctcf.bed"))
export(peaks$H3K27ac, file.path(tmp, "k27.bed"))
sizes <- read.table(file.path(D, "hg38.chrom.sizes"))
sizes <- sizes[grepl("^chr([0-9]+|X|Y)$", sizes$V1), ]
write.table(sizes, file.path(tmp, "genome.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
sh <- function(cmd) as.integer(system(cmd, intern = TRUE))
sort_ctcf <- sprintf("sort -k1,1 -k2,2n %s/ctcf.bed > %s/ctcf.s.bed", tmp, tmp); system(sort_ctcf)
obs <- sh(sprintf("bedtools intersect -u -a %1$s/ctcf.s.bed -b %1$s/k27.bed | wc -l", tmp))
system(sprintf("zcat %s/hg38_blacklist.bed.gz > %s/bl.bed", normalizePath(D), tmp))
shuf <- sapply(1:100, function(i) sh(sprintf(paste(
  "bedtools shuffle -i %1$s/k27.bed -g %1$s/genome.txt -chrom -excl %1$s/bl.bed -seed %2$d",
  "| sort -k1,1 -k2,2n | bedtools intersect -u -a %1$s/ctcf.s.bed -b - | wc -l"), tmp, i)))
cat(sprintf("observed %d, shuffled mean %.0f (range %d-%d), enrichment %.1f, CTCF peaks %d (%.0f%% overlap)\n",
            obs, mean(shuf), min(shuf), max(shuf), obs / mean(shuf), length(peaks$CTCF),
            100 * obs / length(peaks$CTCF)))
p18 <- ggplot(data.frame(n = shuf), aes(n)) +
  geom_histogram(binwidth = 20, fill = "grey45", colour = NA) +
  geom_vline(xintercept = obs, colour = "#c0392b", linewidth = 1.3) +
  annotate("text", x = obs, y = Inf, vjust = 1.5, hjust = 1.05, colour = "#c0392b", size = 5.2,
           label = sprintf("observed: %d CTCF peaks\noverlap an H3K27ac peak", obs)) +
  annotate("text", x = mean(shuf), y = Inf, vjust = 1.5, hjust = -0.15, size = 5.2, colour = "grey30",
           label = sprintf("100 shuffles of H3K27ac\nmean %.0f  (%.1f-fold enrichment)", mean(shuf), obs / mean(shuf))) +
  scale_x_continuous(limits = c(0, obs * 1.08), expand = c(0, 0)) +
  labs(x = "number of CTCF peaks overlapping an H3K27ac peak", y = "number of shuffles")
ggsave("fig18-shuffle-real.png", p18, width = 11, height = 4.2, dpi = 150, bg = "white")

# ---- fig10: browser view around ICOSLG / DNMT3L / AIRE on chr21 --------------
region <- GRanges("chr21", IRanges(44240001, 44320000))
bw  <- import(file.path(D, "gm12878_ctcf_fc.bigWig"), which = region, as = "RleList")$chr21
v   <- as.numeric(bw[start(region):end(region)])        # zero where the bigWig has no data
b   <- cut(seq_along(v), 1600, labels = FALSE)
sig <- data.frame(x = start(region) - 1 + tapply(seq_along(v), b, mean), y = tapply(v, b, max))
pk  <- subsetByOverlaps(peaks$CTCF, region)
pkdf <- data.frame(s = start(pk), e = end(pk), summit = start(pk) + pk$peak)
ex  <- subsetByOverlaps(gtf[gtf$type == "exon" & gtf$transcript_type == "protein_coding"], region)
txr <- subsetByOverlaps(tx, region)
# one transcript per gene (the longest) to keep the panel readable
txr <- txr[order(-width(txr))]; txr <- txr[!duplicated(txr$gene_name)]
txr$row <- seq_along(txr)
exr <- ex[ex$transcript_id %in% txr$transcript_id]
exdf <- data.frame(s = start(exr), e = end(exr), row = txr$row[match(exr$transcript_id, txr$transcript_id)])
txdf <- data.frame(s = pmax(start(txr), start(region)), e = pmin(end(txr), end(region)), row = txr$row,
                   g = txr$gene_name, str = as.character(strand(txr)))
xl <- c(start(region), end(region))
mb <- scale_x_continuous(limits = xl, expand = c(0, 0),
                         labels = function(b) sprintf("%.2f Mb", b / 1e6))
noax <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              panel.grid.major.x = element_blank())
pa <- ggplot(sig, aes(x, y)) + geom_area(fill = col["CTCF"], colour = col["CTCF"], linewidth = 0.2) +
  mb + labs(y = "CTCF\nfold change", title = "bigWig: signal") + noax
pb <- ggplot(pkdf) + geom_rect(aes(xmin = s, xmax = e, ymin = 0, ymax = 1), fill = col["CTCF"]) +
  geom_segment(aes(x = summit, xend = summit, y = -0.25, yend = 1.25), linewidth = 0.6) +
  mb + scale_y_continuous(breaks = NULL) + labs(y = NULL, title = "narrowPeak: peaks (line = summit)") + noax
pc <- ggplot() +
  geom_segment(data = txdf, aes(x = s, xend = e, y = -row, yend = -row), colour = "#1b2a3a") +
  geom_rect(data = exdf, aes(xmin = s, xmax = e, ymin = -row - 0.25, ymax = -row + 0.25), fill = "#1b2a3a") +
  geom_text(data = txdf, aes(x = pmin(pmax((s + e) / 2, xl[1] + 4000), xl[2] - 4000), y = -row - 0.32, vjust = 1,
                             label = paste0(g, ifelse(str == "+", "  →", "  ←"))), size = 4.3) +
  mb + scale_y_continuous(breaks = NULL, expand = expansion(add = 0.5)) +
  labs(x = "chr21", y = NULL, title = "GTF: genes (GENCODE)") +
  theme(plot.margin = margin(5, 30, 5, 5))
p10 <- pa / pb / pc + plot_layout(heights = c(3, 0.8, 2.8))
ggsave("fig10-browser-ctcf.png", p10, width = 12, height = 6, dpi = 150, bg = "white")
cat("CTCF peaks in browser region:", length(pk), "\n")
