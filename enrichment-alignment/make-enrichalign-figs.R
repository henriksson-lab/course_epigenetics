# Real-data figures for enrichment-alignment/lec.qmd.
# Run from this folder:  Rscript make-enrichalign-figs.R
# Needs: data/download.sh enrichment-alignment (CTCF ChIP and input FASTQ, chr21+chrM
#        reference and bowtie2 index, ENCODE CTCF fold change and peaks, blacklist);
#        on the PATH: perl, fastp, bowtie2, samtools, bedtools, bamCoverage, bamCompare;
#        R packages GenomicAlignments, rtracklayer, ggplot2, patchwork.
#
# The script first runs the exercise's command-line chain (the same commands) into a
# cache folder (default: a temporary folder; set ENRICHALIGN_CACHE to keep it), then
# reads the results. Writes fig04, fig11, fig14, fig15, fig17, fig20, fig21 (*.png)
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicAlignments); library(rtracklayer); library(ggplot2); library(patchwork)
})
set.seed(1)
D <- normalizePath("../data")
cache <- Sys.getenv("ENRICHALIGN_CACHE", file.path(tempdir(), "enrichalign"))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
here <- getwd()
col <- c(ChIP = "#8c564b", input = "#5b6b7a")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))

# ---- run the command-line chain (as in the exercise) ------------------------
sh <- function(cmd) {
  t0 <- Sys.time()
  status <- system2("bash", c("-c", shQuote(paste("set -euo pipefail;", cmd))))
  if (status != 0) stop("failed: ", cmd)
  cat(sprintf("[%5.1f s] %s\n", as.numeric(Sys.time() - t0, units = "secs"), substr(cmd, 1, 90)))
}
if (!file.exists(file.path(cache, "counts.txt"))) {
  setwd(cache)
  for (s in c("chip", "input")) {
    sh(sprintf("zcat %s/gm12878_ctcf_%s.fastq.gz | perl -pe 'tr/@-~/!-_/ if $. %% 4 == 0' | gzip > %s.fastq.gz", D, s, s))
    sh(sprintf("fastp -w 8 -i %1$s.fastq.gz -o %1$s.trim.fastq.gz -h %1$s.fastp.html -j %1$s.fastp.json 2> %1$s.fastp.log", s))
    sh(sprintf("bowtie2 -p 8 -x %1$s/bowtie2/hg38_chr21_chrM -U %2$s.trim.fastq.gz 2> %2$s.bowtie2.log | samtools sort -@ 4 -o %2$s.sorted.bam -", D, s))
    sh(sprintf("samtools index %s.sorted.bam", s))
  }
  sh(sprintf("zcat %s/hg38_blacklist.bed.gz > blacklist.bed", D))
  for (s in c("chip", "input")) {
    sh(sprintf("samtools view -b -q 30 -F 1804 -o %1$s.q30.bam %1$s.sorted.bam", s))
    sh(sprintf("samtools markdup -s %1$s.q30.bam %1$s.markdup.bam 2> %1$s.markdup.txt && samtools index %1$s.markdup.bam", s))
    sh(sprintf("samtools view -b -F 1024 -o %1$s.dedup.bam %1$s.markdup.bam chr21", s))
    sh(sprintf("bedtools intersect -v -abam %1$s.dedup.bam -b blacklist.bed > %1$s.final.bam && samtools index %1$s.final.bam", s))
    sh(sprintf("bamCoverage -b %1$s.final.bam -o %1$s.cpm.bw --normalizeUsing CPM --extendReads 100 --binSize 10 -p 8 > /dev/null 2>&1", s))
  }
  sh(paste("bamCompare -b1 chip.final.bam -b2 input.final.bam -o chip_vs_input.ratio.bw --operation ratio",
           "--scaleFactorsMethod readCount --pseudocount 1 --extendReads 100 --binSize 50 -p 8 > /dev/null 2>&1"))
  sh(paste("bamCompare -b1 chip.final.bam -b2 input.final.bam -o chip_vs_input.log2.bw --operation log2",
           "--scaleFactorsMethod readCount --pseudocount 1 --extendReads 100 --binSize 50 -p 8 > /dev/null 2>&1"))
  sh(paste("for s in chip input; do echo $s",
           "$(zcat $s.fastq.gz | awk 'END { print NR / 4 }')",
           "$(zcat $s.trim.fastq.gz | awk 'END { print NR / 4 }')",
           "$(samtools view -c -F 4 $s.sorted.bam) $(samtools view -c $s.q30.bam)",
           "$(samtools view -c -F 1024 $s.markdup.bam) $(samtools view -c $s.dedup.bam)",
           "$(samtools view -c $s.final.bam); done > counts.txt"))
  setwd(here)
}
f <- function(x) file.path(cache, x)

# ---- fig15: reads left after each filtering step -----------------------------
cnt <- read.table(f("counts.txt"), col.names = c("sample", "raw", "fastp", "mapped", "q30", "dedup", "chr21", "final"))
print(cnt)
steps <- c(raw = "raw reads", fastp = "after fastp", mapped = "aligned", q30 = "MAPQ ≥ 30, -F 1804",
           dedup = "duplicates removed", chr21 = "chr21 only (no chrM)", final = "blacklist removed")
long <- do.call(rbind, lapply(seq_len(nrow(cnt)), function(i)
  data.frame(sample = ifelse(cnt$sample[i] == "chip", "ChIP", "input"), step = names(steps),
             n = unlist(cnt[i, names(steps)]))))
long$step <- factor(steps[long$step], levels = rev(steps))
long$lab <- sprintf("%s (%.0f%%)", format(long$n, big.mark = " "), 100 * long$n / ave(long$n, long$sample, FUN = max))
p15 <- ggplot(long, aes(n / 1000, step, fill = sample)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = lab), hjust = -0.05, size = 4.6) +
  facet_wrap(~sample, ncol = 2) +
  scale_fill_manual(values = col, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.62))) +
  labs(x = "reads (thousands)", y = NULL)
ggsave("fig15-filter-steps.png", p15, width = 12, height = 4.4, dpi = 150, bg = "white")

# ---- fig04: base quality along the read ---------------------------------------
qual <- function(file, n = 1e5) {
  l <- readLines(gzfile(file), n = 4 * n)
  q <- l[seq(4, length(l), 4)]
  m <- do.call(rbind, lapply(q, function(s) utf8ToInt(s) - 64L))   # Phred+64 in the raw file
  data.frame(cycle = seq_len(ncol(m)), q10 = apply(m, 2, quantile, 0.1), q90 = apply(m, 2, quantile, 0.9), q25 = apply(m, 2, quantile, 0.25),
             med = apply(m, 2, median), q75 = apply(m, 2, quantile, 0.75), mean = colMeans(m),
             below20 = colMeans(m < 20))
}
qq <- rbind(cbind(sample = "ChIP (36 bp, 2011)", qual(file.path(D, "gm12878_ctcf_chip.fastq.gz"))),
            cbind(sample = "input (28 bp, 2009)", qual(file.path(D, "gm12878_ctcf_input.fastq.gz"))))
print(aggregate(cbind(mean, below20) ~ sample, qq, function(x) round(c(first = x[1], last = x[length(x)]), 3)))
p04 <- ggplot(qq, aes(cycle)) +
  geom_hline(yintercept = c(20, 30), linetype = 2, colour = "grey60") +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#2a6f97", alpha = 0.18) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2a6f97", alpha = 0.3) +
  geom_line(aes(y = med), colour = "#2a6f97", linewidth = 1.2) +
  geom_line(aes(y = mean), colour = "#c0392b", linewidth = 1, linetype = 1) +
  annotate("text", x = 1, y = 21.5, label = "Q20", hjust = 0, colour = "grey40", size = 4.2) +
  annotate("text", x = 1, y = 31.5, label = "Q30", hjust = 0, colour = "grey40", size = 4.2) +
  facet_wrap(~sample, scales = "free_x") +
  coord_cartesian(ylim = c(0, 41)) +
  labs(x = "position in read (cycle)", y = "base quality Q",
       caption = "blue line: median; dark band: 25-75%; light band: 10-90%; red: mean. First 100 000 reads.")
ggsave("fig04-quality.png", p04, width = 12, height = 4.2, dpi = 150, bg = "white")

# ---- fig11: MAPQ of the reads on chr21 ----------------------------------------
mq <- do.call(rbind, lapply(c("chip", "input"), function(s) {
  a <- readGAlignments(f(paste0(s, ".sorted.bam")), param = ScanBamParam(what = "mapq"))
  data.frame(sample = ifelse(s == "chip", "ChIP (36 bp)", "input (28 bp)"),
             chr = as.character(seqnames(a)), mapq = mcols(a)$mapq)
}))
mq21 <- mq[mq$chr == "chr21", ]
tab <- aggregate(mapq ~ sample, mq21, function(x) c(n = length(x), le1 = mean(x <= 1), ge30 = mean(x >= 30)))
print(tab)
lab11 <- data.frame(sample = tab$sample,
                    txt = sprintf("MAPQ 0-1: %.0f%%\nMAPQ ≥ 30: %.0f%%", 100 * tab$mapq[, "le1"], 100 * tab$mapq[, "ge30"]))
p11 <- ggplot(mq21, aes(mapq, fill = sample)) +
  geom_bar(width = 0.9) +
  geom_vline(xintercept = 29.5, linetype = 2, colour = "#c0392b") +
  geom_label(data = lab11, aes(x = 15, y = Inf, label = txt), vjust = 1.2, size = 4.8,
             linewidth = 0, fill = "white", inherit.aes = FALSE) +
  facet_wrap(~sample) +
  scale_fill_manual(values = c("ChIP (36 bp)" = col[["ChIP"]], "input (28 bp)" = col[["input"]]), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x / 1000, "k")) +
  labs(x = "MAPQ (bowtie2)", y = "reads on chr21")
ggsave("fig11-mapq.png", p11, width = 12, height = 4.0, dpi = 150, bg = "white")

# ---- fig14: a blacklisted region, and chrM ------------------------------------
bl <- import(file.path(D, "hg38_blacklist.bed.gz"))
win <- GRanges("chr21", IRanges(9995001, 10035000))
blw <- subsetByOverlaps(bl, win); seqlevels(blw, pruning.mode = "coarse") <- "chr21"
covtrack <- function(s) {
  a <- readGAlignments(f(paste0(s, ".sorted.bam")), param = ScanBamParam(which = win))
  cv <- coverage(resize(granges(a), 100))$chr21[start(win):end(win)]
  data.frame(sample = ifelse(s == "chip", "ChIP", "input"), pos = start(win) + seq_along(cv) - 1, cov = as.numeric(cv))
}
ct <- rbind(covtrack("chip"), covtrack("input"))
ct <- ct[ct$pos %% 20 == 0, ]
inbl <- sapply(c("chip", "input"), function(s)
  length(readGAlignments(f(paste0(s, ".sorted.bam")), param = ScanBamParam(which = blw))))
nall <- sapply(c("chip", "input"), function(s) sum(mq$sample == ifelse(s == "chip", "ChIP (36 bp)", "input (28 bp)") & mq$chr == "chr21"))
expct <- nall * width(blw) / 40088619
cat(sprintf("blacklist %s: reads ChIP %d input %d; expected %.1f %.1f\n", as.character(blw),
            inbl[1], inbl[2], expct[1], expct[2]))
p14a <- ggplot(ct, aes(pos / 1e6, cov, fill = sample)) +
  annotate("rect", xmin = start(blw) / 1e6, xmax = end(blw) / 1e6, ymin = -Inf, ymax = Inf, fill = "#c0392b", alpha = 0.15) +
  geom_area(colour = NA) +
  facet_wrap(~sample, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = col, guide = "none") +
  labs(x = "chr21 position (Mb)", y = "reads covering (all MAPQ)",
       title = sprintf("A blacklisted region: %d input reads where ~%.0f are expected", inbl[2], expct[2])) +
  theme(plot.title = element_text(size = 14), strip.placement = "outside")
idx <- do.call(rbind, lapply(c("chip", "input"), function(s) {
  x <- table(mq$chr[mq$sample == ifelse(s == "chip", "ChIP (36 bp)", "input (28 bp)")])
  data.frame(sample = ifelse(s == "chip", "ChIP", "input"), chr = names(x), n = as.numeric(x))
}))
idx$frac <- idx$n / ave(idx$n, idx$sample, FUN = sum)
print(idx)
p14b <- ggplot(idx, aes(sample, frac, fill = chr)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%s\n%.0f%%", chr, 100 * frac)), position = position_stack(vjust = 0.5),
            colour = "white", size = 4.3, lineheight = 0.9) +
  scale_fill_manual(values = c(chr21 = "#2a6f97", chrM = "#a8324a"), guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "aligned reads", title = "chrM: 16.6 kb") +
  theme(plot.title = element_text(size = 14))
p14 <- p14a + p14b + plot_layout(widths = c(3.2, 1))
ggsave("fig14-blacklist-chrm.png", p14, width = 12, height = 4.4, dpi = 150, bg = "white")

# ---- fig17: + and - strand reads around CTCF summits -----------------------------
pk <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak")
pk <- pk[seqnames(pk) == "chr21"]
summ <- GRanges("chr21", IRanges(start(pk) + pk$peak, width = 1))   # 1-based summit base
a <- granges(readGAlignments(f("chip.final.bam")))
prof <- function(gr, label) {
  w <- promoters(summ, 300, 301)
  cv <- coverage(gr)$chr21
  v <- colSums(do.call(rbind, lapply(seq_along(w), function(i) as.numeric(cv[start(w)[i]:end(w)[i]]))))
  data.frame(pos = -300:300, n = v, what = label)
}
five <- resize(a, 1, fix = "start")
pr <- rbind(prof(five[strand(five) == "+"], "+ strand read 5' ends"),
            prof(five[strand(five) == "-"], "− strand read 5' ends"),
            prof(resize(a, 100, fix = "start"), "reads extended to 100 bp"))
pk5 <- sapply(split(pr, pr$what), function(d) d$pos[which.max(d$n)])
print(pk5)
cat("chr21 summits", length(summ), "\n")
pr$panel <- ifelse(grepl("extended", pr$what), "extended to the fragment", "5' ends of the reads")
pr$panel <- factor(pr$panel, levels = c("5' ends of the reads", "extended to the fragment"))
p17 <- ggplot(pr, aes(pos, n, colour = what)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = c("+ strand read 5' ends" = "#2a6f97", "− strand read 5' ends" = "#e67e22",
                                 "reads extended to 100 bp" = "#8c564b"), name = NULL) +
  labs(x = sprintf("distance from CTCF peak summit (bp), %d summits on chr21", length(summ)),
       y = "reads, summed over summits") +
  theme(legend.position = "bottom")
ggsave("fig17-strands.png", p17, width = 12, height = 4.4, dpi = 150, bg = "white")

# ---- fig20: tracks at a locus -------------------------------------------------
loc <- GRanges("chr21", IRanges(43775001, 43825000))
rd <- function(file, name) {
  x <- import(file, which = loc, as = "RleList")$chr21[start(loc):end(loc)]
  data.frame(track = name, pos = start(loc) + seq_along(x) - 1, val = as.numeric(x))
}
tr <- rbind(rd(f("input.cpm.bw"), "input, CPM"), rd(f("chip.cpm.bw"), "ChIP, CPM"),
            rd(f("chip_vs_input.log2.bw"), "log2(ChIP / input)"),
            rd(file.path(D, "gm12878_ctcf_fc.bigWig"), "ENCODE fold change"))
tr <- tr[tr$pos %% 10 == 0, ]
tr$track <- factor(tr$track, levels = c("input, CPM", "ChIP, CPM", "log2(ChIP / input)", "ENCODE fold change"))
pkl <- as.data.frame(subsetByOverlaps(pk, loc))
p20 <- ggplot(tr, aes(pos / 1e6, val)) +
  geom_area(aes(fill = track), colour = NA) +
  facet_wrap(~track, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = c("input, CPM" = col[["input"]], "ChIP, CPM" = col[["ChIP"]],
                               "log2(ChIP / input)" = "#1b2a3a", "ENCODE fold change" = "#2a6f97"), guide = "none") +
  geom_segment(data = pkl, aes(x = start / 1e6, xend = end / 1e6, y = -Inf, yend = -Inf),
               inherit.aes = FALSE, linewidth = 4, colour = "#c0392b") +
  labs(x = "chr21 position (Mb)  (red bars: ENCODE peaks)", y = NULL) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13), strip.placement = "outside",
        panel.spacing = unit(0.6, "lines"))
ggsave("fig20-locus.png", p20, width = 12, height = 4.8, dpi = 150, bg = "white")

# ---- fig21: agreement with ENCODE over 1-kb bins --------------------------------
bins <- tileGenome(c(chr21 = 46709983), tilewidth = 1000, cut.last.tile.in.chrom = TRUE)
enc  <- import(file.path(D, "gm12878_ctcf_fc.bigWig"), which = GRanges("chr21", IRanges(1, 46709983)), as = "RleList")
mine <- import(f("chip_vs_input.ratio.bw"), as = "RleList")
e <- binnedAverage(bins, enc["chr21"], "v")$v
m <- binnedAverage(bins, mine["chr21"], "v")$v
ok <- !overlapsAny(bins, bl) & e > 0 & m > 0
inpk <- overlapsAny(bins, pk)
r_all <- cor(log2(e[ok]), log2(m[ok])); r_pk <- cor(log2(e[ok & inpk]), log2(m[ok & inpk]))
cat(sprintf("1-kb bins: %d used of %d; Pearson r (log2) %.2f; Spearman %.2f; bins with a peak: %d, r %.2f\n",
            sum(ok), length(bins), r_all, cor(e[ok], m[ok], method = "spearman"), sum(ok & inpk), r_pk))
bd <- data.frame(e = e[ok], m = m[ok], peak = ifelse(inpk[ok], "contains an ENCODE peak", "other bins"))
bd <- bd[order(bd$peak == "contains an ENCODE peak"), ]
p21 <- ggplot(bd, aes(e, m, colour = peak)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
  geom_point(size = 0.8, alpha = 0.35) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("other bins" = "grey55", "contains an ENCODE peak" = "#c0392b"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  annotate("text", x = min(bd$e), y = max(bd$m), hjust = 0, vjust = 1, size = 5,
           label = sprintf("%s bins of 1 kb on chr21\nPearson r (log scale) = %.2f\nbins with a peak: r = %.2f",
                           format(sum(ok), big.mark = " "), r_all, r_pk)) +
  labs(x = "ENCODE fold change (pooled replicates, MACS2)", y = "our ChIP / input (rep 1, bamCompare)") +
  theme(legend.position = "bottom")
ggsave("fig21-bins.png", p21, width = 8, height = 5.6, dpi = 150, bg = "white")

cat("cache:", cache, "\n")
