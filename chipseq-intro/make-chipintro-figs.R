# Real-data figures for chipseq-intro/lec.qmd.
# Run from this folder:  Rscript make-chipintro-figs.R
#
# Needs: data/gm12878_ctcf_{chip,input}.fastq.gz (data/prepare-reads.sh),
#        data/gm12878_h3k36me3_{chip,input}.fastq.gz (data/prepare-chip.sh),
#        data/hg38_chr21_chrM.fa with its bowtie2 index data/bowtie2/hg38_chr21_chrM,
#        the ENCODE fold-change bigWigs, blacklist and GENCODE annotation (data/download.sh);
#        on the PATH: fastp, bowtie2, samtools (>= 1.10), bedtools;
#        R packages GenomicRanges, GenomicAlignments, rtracklayer, ggplot2, patchwork.
#
# The processing (the same commands as in exercise.qmd) runs once into a work folder,
# by default ./chip-work (set CHIP_WORK to put it elsewhere; about 300 MB, 2-3 min on 8 cores).
#
# Writes fig06-input-tss.png, fig13-tracks.png, fig14-autocorrelation.png,
# fig15-tss-profiles.png, fig18-crosscor.png, fig19-h3k36me3-gene.png
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicAlignments); library(rtracklayer)
  library(ggplot2); library(patchwork)
})
set.seed(1)
D <- normalizePath("../data")
W <- Sys.getenv("CHIP_WORK", "chip-work")
dir.create(W, showWarnings = FALSE, recursive = TRUE)
W <- normalizePath(W)
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15)))
col_mark <- c(CTCF = "#8c564b", H3K4me3 = "#c0392b", H3K27ac = "#e67e22", H3K4me1 = "#d4a017",
              H3K36me3 = "#2e8b57", H3K27me3 = "#2a6f97", H3K9me3 = "#6b4fa0")
col_input <- "#5b6b7a"

# ---- 1. processing: the alignment chain of "From reads to coverage" -----------
chain <- paste(
  "chain() { local s=$1;",
  "if [ \"$3\" = 64 ]; then zcat $2 | perl -pe 'tr/@-~/!-_/ if $. % 4 == 0' | gzip > $s.fastq.gz; else cp $2 $s.fastq.gz; fi;",
  "fastp -w $T -i $s.fastq.gz -o $s.trim.fastq.gz -j $s.fastp.json -h $s.fastp.html 2> $s.fastp.log;",
  "bowtie2 -p $T -x $D/bowtie2/hg38_chr21_chrM -U $s.trim.fastq.gz 2> $s.bowtie2.log | samtools sort -@ 4 -o $s.sorted.bam -;",
  "samtools index $s.sorted.bam;",
  "samtools view -b -q 30 -F 1804 -o $s.q30.bam $s.sorted.bam;",
  "samtools markdup -s $s.q30.bam $s.markdup.bam 2> $s.markdup.txt;",
  "samtools index $s.markdup.bam;",
  "samtools view -b -F 1024 -o $s.dedup.bam $s.markdup.bam chr21;",
  "bedtools intersect -v -abam $s.dedup.bam -b blacklist.bed > $s.final.bam;",
  "samtools index $s.final.bam; };")
step <- function(out, cmd) {
  if (all(file.exists(file.path(W, out)))) return(invisible())
  cat("running:", out, "\n")
  rc <- system2("bash", c("-c", shQuote(paste0("set -euo pipefail; cd ", W, "; D=", D, "; T=", THREADS,
                                                "; zcat $D/hg38_blacklist.bed.gz > blacklist.bed; ", chain, cmd))))
  if (rc != 0) stop("step failed: ", out)
}
samples <- data.frame(name = c("ctcf", "ctcf_input", "h3k36me3", "h3k36me3_input"),
                      fq = c("gm12878_ctcf_chip", "gm12878_ctcf_input", "gm12878_h3k36me3_chip", "gm12878_h3k36me3_input"),
                      phred = c(64, 64, 33, 33))
for (i in seq_len(nrow(samples)))
  step(paste0(samples$name[i], ".final.bam.bai"),
       sprintf("chain %s $D/%s.fastq.gz %d", samples$name[i], samples$fq[i], samples$phred[i]))

rd <- function(s) granges(readGAlignments(file.path(W, paste0(s, ".final.bam"))))
reads <- lapply(setNames(samples$name, samples$name), rd)
print(sapply(reads, length))

# ---- annotation and ENCODE fold-change tracks --------------------------------
chr21 <- GRanges("chr21", IRanges(1, 46709983))
G21 <- 40088619                                   # non-N length of chr21
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), which = chr21)
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
seqlevels(genes, pruning.mode = "coarse") <- "chr21"
tss <- resize(genes, 1, fix = "start")
marks <- c(CTCF = "ctcf", H3K4me3 = "h3k4me3", H3K27ac = "h3k27ac", H3K4me1 = "h3k4me1",
           H3K36me3 = "h3k36me3", H3K27me3 = "h3k27me3", H3K9me3 = "h3k9me3")
fc <- lapply(marks, function(m)
  import(file.path(D, sprintf("gm12878_%s_fc.bigWig", m)), which = chr21, as = "RleList")["chr21"])
k4 <- viewMeans(Views(fc$H3K4me3$chr21, ranges(promoters(tss, 1000, 1001))))
genes$active <- k4 > 5; tss$active <- k4 > 5
genes$silent <- k4 < 1; tss$silent <- k4 < 1
cat("protein-coding genes on chr21:", length(genes), " active (H3K4me3 FC > 5 at TSS):", sum(genes$active),
    " silent (FC < 1):", sum(genes$silent), "\n")

# ---- fig06: input reads at active and silent TSSs -----------------------------
prof_cov <- function(gr, pos, str, w, ext = 100) {
  cv <- coverage(resize(gr, ext, fix = "start"))$chr21
  avg <- length(gr) * ext / G21
  m <- sapply(seq_along(pos), function(i) { v <- as.numeric(cv[(pos[i] - w):(pos[i] + w)]); if (str[i] == "-") rev(v) else v })
  list(mean = rowMeans(m) / avg, median = apply(m, 1, median) / avg)
}
w <- 3000
d06 <- do.call(rbind, lapply(c("ctcf_input", "h3k36me3_input"), function(s) {
  do.call(rbind, lapply(c("active", "silent"), function(g) {
    sel <- mcols(tss)[[g]]
    p <- prof_cov(reads[[s]], start(tss)[sel], as.character(strand(tss))[sel], w)
    data.frame(input = s, group = g, x = -w:w, mean = p$mean, median = p$median)
  }))
}))
d06$input <- factor(d06$input, c("ctcf_input", "h3k36me3_input"),
                    c("input of the CTCF ChIP (Yale, 28 bp)", "input of the H3K36me3 ChIP (Broad, 51 bp)"))
d06$group <- factor(d06$group, c("active", "silent"),
                    c(sprintf("active TSSs (n = %d)", sum(tss$active)), sprintf("silent TSSs (n = %d)", sum(tss$silent))))
print(aggregate(cbind(mean, median) ~ input + group, data = d06[d06$x == 0, ], FUN = identity))
print(aggregate(cbind(mean, median) ~ input + group, data = d06[abs(d06$x) == 2000, ], FUN = mean))
d06s <- d06[d06$x %% 20 == 0, ]
p06 <- ggplot(d06s, aes(x / 1000, mean, colour = group)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 1.1) + facet_wrap(~ input) +
  scale_colour_manual(values = c("#2e8b57", "#a8324a"), name = NULL) +
  labs(x = "distance from TSS (kb)", y = "input coverage / chr21 average") +
  theme(legend.position = "top", strip.text = element_text(size = 14))
ggsave("fig06-input-tss.png", p06, width = 11, height = 4.6, dpi = 150, bg = "white")

# ---- fig13: tracks at one locus ----------------------------------------------
win <- GRanges("chr21", IRanges(33000000, 33620000))
d13 <- do.call(rbind, lapply(names(marks), function(m) {
  b <- tile(win, width = 500)[[1]]
  v <- binnedAverage(b, fc[[m]], "v")$v
  data.frame(mark = m, x = (start(b) + 250) / 1e6, v = v)
}))
d13$mark <- factor(d13$mark, names(marks))
gw <- subsetByOverlaps(genes, win)
gw <- gw[width(gw) > 2000 & !grepl("^ENSG|-", gw$gene_name)]
gw <- sort(gw)
gdf <- data.frame(s = pmax(start(gw), start(win)) / 1e6, e = pmin(end(gw), end(win)) / 1e6,
                  name = gw$gene_name, strand = as.character(strand(gw)),
                  active = ifelse(gw$active, "active", ifelse(gw$silent, "silent", "other")))
gdf$y <- rep(c(0.3, 0.65, 1.0), length.out = nrow(gdf))
p13a <- ggplot(d13, aes(x, v, fill = mark)) + geom_area(colour = NA) +
  facet_grid(mark ~ ., scales = "free_y", switch = "y") + scale_fill_manual(values = col_mark, guide = "none") +
  scale_y_continuous(breaks = function(l) c(0, floor(l[2] * 0.9)), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = NULL) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13, face = "bold"),
        strip.placement = "outside", axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 9))
p13b <- ggplot(gdf) +
  geom_segment(aes(x = s, xend = e, y = y, yend = y, colour = active), linewidth = 3) +
  geom_text(aes(x = (s + e) / 2, y = y + 0.17, label = name), size = 3.6) +
  scale_colour_manual(values = c(active = "#2e8b57", silent = "#a8324a", other = "#9aa5b1"), guide = "none") +
  scale_y_continuous(limits = c(0.1, 1.35)) +
  labs(x = "chr21 position (Mb)", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
ggsave("fig13-tracks.png", p13a / p13b + plot_layout(heights = c(7, 1.1)), width = 12, height = 6.6, dpi = 150, bg = "white")
cat("fig13 genes:", paste0(gdf$name, " (", gdf$active, ")", collapse = ", "), "\n")

# ---- fig14: autocorrelation of each mark ---------------------------------------
b100 <- tileGenome(c(chr21 = 46709983), tilewidth = 100, cut.last.tile.in.chrom = TRUE)
b100 <- b100[start(b100) > 13e6]                     # leave out the short arm and centromere
lags <- unique(round(10^seq(0, log10(1000), length.out = 30)))
d14 <- do.call(rbind, lapply(names(marks), function(m) {
  v <- log2(binnedAverage(b100, fc[[m]], "v")$v + 0.5)
  data.frame(mark = m, lag = lags * 100,
             r = sapply(lags, function(l) cor(v[1:(length(v) - l)], v[(1 + l):length(v)])))
}))
d14$mark <- factor(d14$mark, names(marks))
for (L in c(200, 1000, 10000, 50000)) { z <- d14[d14$lag == d14$lag[which.min(abs(d14$lag - L))], ]; cat("autocorrelation at", z$lag[1], "bp:", paste(z$mark, round(z$r, 2)), "\n") }
d14p <- d14[d14$mark %in% c("CTCF", "H3K4me3", "H3K36me3", "H3K27me3"), ]
p14 <- ggplot(d14p, aes(lag, r, colour = mark)) + geom_line(linewidth = 1.3) +
  scale_colour_manual(values = col_mark, name = NULL) +
  scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000, 30000, 1e5),
                labels = c("100 bp", "300 bp", "1 kb", "3 kb", "10 kb", "30 kb", "100 kb")) +
  labs(x = "distance between two 100-bp bins", y = "correlation of the signal") +
  theme(legend.position = c(0.85, 0.8), legend.text = element_text(size = 14))
ggsave("fig14-autocorrelation.png", p14, width = 9, height = 5, dpi = 150, bg = "white")

# ---- fig15: profiles at active TSSs, and a scaled metagene ---------------------
prof_fc <- function(rle, pos, str, w) {
  m <- sapply(seq_along(pos), function(i) { v <- as.numeric(rle[(pos[i] - w):(pos[i] + w)]); if (str[i] == "-") rev(v) else v })
  rowMeans(m)
}
d15a <- do.call(rbind, lapply(c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3"), function(m)
  data.frame(mark = m, x = -w:w, v = prof_fc(fc[[m]]$chr21, start(tss)[tss$active], as.character(strand(tss))[tss$active], w))))
d15a <- d15a[d15a$x %% 10 == 0, ]
d15a$mark <- factor(d15a$mark, names(marks))
for (m in levels(droplevels(d15a$mark))) {
  z <- d15a[d15a$mark == m, ]
  cat(m, "at TSS:", round(z$v[z$x == 0], 1), " max:", round(max(z$v), 1), "at", z$x[which.max(z$v)], "bp\n")
}
p15a <- ggplot(d15a, aes(x / 1000, v, colour = mark)) + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = col_mark, name = NULL) +
  labs(title = sprintf("Around %d active TSSs", sum(tss$active)), x = "distance from TSS (kb)", y = "mean fold change over input") +
  theme(legend.position = c(0.84, 0.8))
# metagene: 5 kb upstream, gene body scaled to 20 bins, 5 kb downstream
meta <- function(rle, g) {
  t(sapply(seq_along(g), function(i) {
    s <- start(g)[i]; e <- end(g)[i]
    up <- binnedAverage(GRanges("chr21", IRanges(seq(s - 5000, s - 500, 500), width = 500)), RleList(chr21 = rle), "v")$v
    body <- binnedAverage(GRanges("chr21", IRanges(round(seq(s, e, length.out = 21))[-21], width = diff(round(seq(s, e, length.out = 21))))), RleList(chr21 = rle), "v")$v
    dn <- binnedAverage(GRanges("chr21", IRanges(seq(e + 1, e + 4501, 500), width = 500)), RleList(chr21 = rle), "v")$v
    v <- c(up, body, dn); if (as.character(strand(g))[i] == "-") rev(v) else v
  }))
}
gl <- genes[width(genes) > 5000 & start(genes) > 5001 & (genes$active | genes$silent)]
d15b <- do.call(rbind, lapply(c("H3K36me3"), function(m) do.call(rbind, lapply(c("active", "silent"), function(gr) {
  mm <- meta(fc[[m]]$chr21, gl[mcols(gl)[[gr]]])
  data.frame(mark = m, group = gr, pos = seq_len(ncol(mm)), v = colMeans(mm))
}))))
d15b$group <- factor(d15b$group, c("active", "silent"),
                     c(sprintf("active genes (n = %d)", sum(gl$active)), sprintf("silent genes (n = %d)", sum(gl$silent))))
print(aggregate(v ~ mark + group, d15b[d15b$pos %in% 11:30, ], mean))
p15b <- ggplot(d15b, aes(pos, v, colour = mark, linetype = group)) + geom_line(linewidth = 1.1) +
  annotate("rect", xmin = 10.5, xmax = 30.5, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = col_mark, guide = "none") + scale_linetype_manual(values = c(1, 2), name = NULL) +
  scale_x_continuous(breaks = c(1, 10.5, 30.5, 40), labels = c("-5 kb", "TSS", "TES", "+5 kb")) +
  labs(title = "H3K36me3 over whole genes (body scaled)", x = NULL, y = "mean fold change over input") +
  theme(legend.position = c(0.22, 0.88))
ggsave("fig15-tss-profiles.png", p15a + p15b, width = 12.5, height = 4.9, dpi = 150, bg = "white")

# ---- fig18: strand cross-correlation, sharp vs broad --------------------------
crosscor <- function(gr, maxk = 600, from = 13e6, to = 46709983) {
  # Pearson correlation of + strand 5' end counts at x with - strand counts at x + k,
  # over chr21 without the short arm; computed from the non-zero positions only (fast)
  five <- resize(gr, 1, fix = "start")
  five <- five[start(five) >= from & start(five) <= to - maxk]
  n <- to - from + 1 - maxk
  p <- tabulate(start(five)[as.character(strand(five)) == "+"] - from + 1, n + maxk)
  m <- tabulate(start(five)[as.character(strand(five)) == "-"] - from + 1, n + maxk)
  ip <- which(p[1:n] > 0)
  mp <- mean(p[1:n]); sp <- sd(p[1:n]); mm <- mean(m); sm <- sd(m)   # a shift changes these negligibly
  sapply(0:maxk, function(k) (sum(p[ip] * m[ip + k]) / n - mp * mm) / (sp * sm))
}
d18 <- do.call(rbind, lapply(samples$name, function(s) data.frame(sample = s, k = 0:600, cc = crosscor(reads[[s]]))))
for (s in samples$name) { z <- d18[d18$sample == s & d18$k > 40, ]; cat(s, "max CC", signif(max(z$cc), 2), "at k =", z$k[which.max(z$cc)], "; CC at read length:", signif(d18$cc[d18$sample == s & d18$k == ifelse(grepl("ctcf", s), 36, 51)], 2), "\n") }
d18 <- do.call(rbind, lapply(split(d18, d18$sample), function(z) {   # 25-bp running mean (the broad mark is noisy per bp)
  z$cc <- as.numeric(stats::filter(z$cc, rep(1 / 25, 25), sides = 2)); z[!is.na(z$cc), ] }))
d18$assay <- ifelse(grepl("ctcf", d18$sample), "CTCF (TF), 36 bp reads", "H3K36me3 (broad mark), 51 bp reads")
d18$what <- ifelse(grepl("input", d18$sample), "input", "ChIP")
p18 <- ggplot(d18, aes(k, cc, colour = what)) + geom_line(linewidth = 1.1) +
  facet_wrap(~ assay, scales = "free_y") +
  scale_colour_manual(values = c(ChIP = "#1b2a3a", input = col_input), name = NULL) +
  labs(x = "shift of the − strand (bp)", y = "strand cross-correlation (25-bp running mean)") +
  theme(legend.position = "top", strip.text = element_text(size = 14))
ggsave("fig18-crosscor.png", p18, width = 11, height = 4.6, dpi = 150, bg = "white")

# ---- fig19: our own H3K36me3 reads over transcribed genes ---------------------
win19 <- GRanges("chr21", IRanges(33220000, 33600000))
cpm <- function(gr) coverage(resize(gr, 200, fix = "start"))$chr21 * 1e6 / length(gr)
b19 <- tile(win19, width = 1000)[[1]]
d19 <- do.call(rbind, lapply(c("h3k36me3", "h3k36me3_input", "ctcf"), function(s)
  data.frame(sample = s, x = (start(b19) + 500) / 1e6, v = binnedAverage(b19, RleList(chr21 = cpm(reads[[s]])), "v")$v)))
d19$sample <- factor(d19$sample, c("h3k36me3", "h3k36me3_input", "ctcf"),
                     c("H3K36me3 ChIP", "its input", "CTCF ChIP"))
g19 <- subsetByOverlaps(genes, win19); g19 <- sort(g19[width(g19) > 2000 & !grepl("^ENSG|-", g19$gene_name)])
g19df <- data.frame(s = pmax(start(g19), start(win19)) / 1e6, e = pmin(end(g19), end(win19)) / 1e6, name = g19$gene_name,
                    str = as.character(strand(g19)), active = ifelse(g19$active, "active", ifelse(g19$silent, "silent", "other")))
g19df$y <- rep(c(0.3, 0.65, 1.0), length.out = nrow(g19df))
p19a <- ggplot(d19, aes(x, v, fill = sample)) + geom_area() +
  facet_grid(sample ~ ., switch = "y", scales = "free_y") +
  scale_fill_manual(values = c("#2e8b57", col_input, "#8c564b"), guide = "none") +
  labs(x = NULL, y = "CPM, 1-kb bins") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13, face = "bold"), strip.placement = "outside",
        axis.text.x = element_blank())
p19b <- ggplot(g19df) + geom_segment(aes(x = s, xend = e, y = y, yend = y, colour = active), linewidth = 3) +
  geom_text(aes(x = (s + e) / 2, y = y + 0.17, label = paste0(name, ifelse(str == "+", " →", " ←"))), size = 3.6) +
  scale_colour_manual(values = c(active = "#2e8b57", silent = "#a8324a", other = "#9aa5b1"), guide = "none") +
  scale_y_continuous(limits = c(0.1, 1.35)) + labs(x = "chr21 position (Mb)", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
ggsave("fig19-h3k36me3-gene.png", p19a / p19b + plot_layout(heights = c(4, 1.1)), width = 12, height = 5.2, dpi = 150, bg = "white")
# reads in active vs silent gene bodies, per kb, relative to the input
gb <- genes[width(genes) > 5000]
dens <- function(gr, g) countOverlaps(g, resize(gr, 1, fix = "start"), ignore.strand = TRUE) / width(g) * 1000 / length(gr) * 1e6
for (grp in c("active", "silent")) {
  g <- gb[mcols(gb)[[grp]]]
  cat(grp, "gene bodies: H3K36me3", round(sum(dens(reads$h3k36me3, g) * width(g)) / sum(width(g)), 2),
      " input", round(sum(dens(reads$h3k36me3_input, g) * width(g)) / sum(width(g)), 2), "reads per kb per million\n")
}
cat("chr21 average:", round(1000 / G21 * 1e6, 2), "reads per kb per million\n")
