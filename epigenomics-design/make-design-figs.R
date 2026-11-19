# Real-data figures for epigenomics-design/lec.qmd.
# Run from this folder:  Rscript make-design-figs.R
#
# Needs: data/download.sh epigenomics-design (reads, bowtie2 index, peaks,
#        blacklist, GENCODE), and on the PATH: fastp, bowtie2, samtools (>= 1.12),
#        bedtools, macs3. R packages GenomicRanges, GenomicAlignments, Rsamtools,
#        rtracklayer, ggplot2, patchwork.
#
# The aligned, filtered BAM files are made with exactly the commands of the
# exercise (alignment block) in a work folder: $DESIGN_WORK if set, otherwise a
# temporary folder. Aligning the ATAC-seq reads takes about 10-20 minutes with
# 8 threads; set DESIGN_WORK to keep the BAMs between runs.
#
# Writes fig15-chip-vs-input.png, fig16-depth-saturation.png,
# fig17-complexity-real.png, fig18-frip-real.png, fig19-fingerprint-real.png,
# fig20-crosscor-real.png, fig21-replicates-real.png, fig22-atac-qc-real.png
# and prints every number quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicAlignments); library(Rsamtools)
  library(rtracklayer); library(ggplot2); library(patchwork)
})
set.seed(1)
D <- normalizePath("../data")
W <- Sys.getenv("DESIGN_WORK", file.path(tempdir(), "design"))
dir.create(W, showWarnings = FALSE, recursive = TRUE)
threads <- 8
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
col <- c(CTCF = "#8c564b", input = "#7f8c99", ATAC1 = "#2e8b57", ATAC2 = "#7fbf9b")
chr21 <- GRanges("chr21", IRanges(1, 46709983))
sh <- function(cmd) { message("+ ", substr(cmd, 1, 120)); stopifnot(system(cmd) == 0) }
bam <- function(s, kind = "final") file.path(W, sprintf("%s.%s.bam", s, kind))

# ---- 0. alignment and filtering (the exercise's alignment block) -------------
se_block <- '
for s in ctcf_chip ctcf_input; do
  fastp --phred64 -w $T -i $D/gm12878_$s.fastq.gz -o $s.trim.fq.gz \\
        -j $s.fastp.json -h $s.fastp.html 2> $s.fastp.log
  bowtie2 -p $T -x $D/bowtie2/hg38_chr21_chrM -U $s.trim.fq.gz 2> $s.bowtie2.log \\
    | samtools view -u -q 30 -F 0x904 - \\
    | samtools sort -@ $T -o $s.filt.bam -
  samtools index $s.filt.bam
  samtools markdup -r -s $s.filt.bam $s.dedup.bam 2> $s.markdup.txt
  samtools index $s.dedup.bam
  samtools view -b $s.dedup.bam chr21 \\
    | bedtools intersect -v -abam - -b $D/hg38_blacklist.bed.gz > $s.final.bam
  samtools index $s.final.bam
done'
pe_block <- '
for r in rep1 rep2; do
  s=atac_$r
  fastp -w $T -i $D/gm12878_atac_${r}_R1.fastq.gz -I $D/gm12878_atac_${r}_R2.fastq.gz \\
        -o $s.R1.fq.gz -O $s.R2.fq.gz -j $s.fastp.json -h $s.fastp.html 2> $s.fastp.log
  bowtie2 -p $T -X 2000 -x $D/bowtie2/hg38_chr21_chrM \\
          -1 $s.R1.fq.gz -2 $s.R2.fq.gz 2> $s.bowtie2.log \\
    | samtools view -u -q 30 -f 2 -F 0x904 - \\
    | samtools sort -n -u -@ $T - \\
    | samtools fixmate -m -u - - \\
    | samtools view -u -f 2 - \\
    | samtools sort -@ $T -o $s.filt.bam -
  samtools index $s.filt.bam
  samtools markdup -r -s $s.filt.bam $s.dedup.bam 2> $s.markdup.txt
  samtools index $s.dedup.bam
  samtools view -b $s.dedup.bam chr21 \\
    | bedtools intersect -v -abam - -b $D/hg38_blacklist.bed.gz > $s.final.bam
  samtools index $s.final.bam
done'
run_block <- function(block) {
  f <- file.path(W, "block.sh")
  writeLines(c("set -euo pipefail", sprintf("D=%s; T=%d; cd %s", D, threads, W), block), f)
  sh(paste("bash", f))
}
if (!file.exists(bam("ctcf_input"))) run_block(se_block)
if (!file.exists(bam("atac_rep2"))) run_block(pe_block)

# ---- 1. numbers: alignment, duplicates, complexity ---------------------------
samples <- c("ctcf_chip", "ctcf_input", "atac_rep1", "atac_rep2")
lab <- c(ctcf_chip = "CTCF ChIP", ctcf_input = "input", atac_rep1 = "ATAC rep1", atac_rep2 = "ATAC rep2")
count <- function(f, which = NULL) {
  p <- if (is.null(which)) ScanBamParam() else ScanBamParam(which = which)
  sum(countBam(f, param = p)$records)
}
# reads per start position (SE: 5' end and strand; PE: fragment start and end),
# on chr21, filtered reads with duplicates still in
position_counts <- function(s) {
  pe <- grepl("atac", s)
  p <- ScanBamParam(which = chr21)
  if (!pe) {
    r <- readGAlignments(bam(s, "filt"), param = p)
    key <- paste(ifelse(as.character(strand(r)) == "+", start(r), end(r)), strand(r))
  } else {
    r <- readGAlignmentPairs(bam(s, "filt"), param = p)
    g <- granges(r, on.discordant.seqnames = "drop")
    key <- paste(start(g), end(g))
  }
  as.integer(table(key))
}
qc <- do.call(rbind, lapply(samples, function(s) {
  n <- position_counts(s)
  filt <- count(bam(s, "filt")); mito <- count(bam(s, "filt"), GRanges("chrM", IRanges(1, 16569)))
  data.frame(sample = s, filt_reads = filt, chrM_frac = mito / filt,
             chr21_units = sum(n), distinct = length(n), M1 = sum(n == 1), M2 = sum(n == 2),
             NRF = length(n) / sum(n), PBC1 = sum(n == 1) / length(n), PBC2 = sum(n == 1) / sum(n == 2),
             dup_frac = 1 - length(n) / sum(n), final_reads = count(bam(s)))
}))
print(qc, digits = 3)
pos_tab <- do.call(rbind, lapply(samples, function(s) {
  n <- pmin(position_counts(s), 10)
  data.frame(sample = lab[[s]], k = as.integer(names(table(n))), frac = as.numeric(table(n)) / length(n))
}))

# ---- fig17: reads per position ------------------------------------------------
pos_tab$sample <- factor(pos_tab$sample, levels = lab)
p17 <- ggplot(pos_tab, aes(factor(k), frac, fill = sample)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  scale_fill_manual(values = setNames(col, lab), name = NULL) +
  scale_x_discrete(labels = c(1:9, "10+")) +
  scale_y_continuous(trans = "sqrt", breaks = c(0.01, 0.05, 0.2, 0.5, 0.9)) +
  labs(x = "reads (SE) or fragments (PE) at the same position", y = "fraction of positions (sqrt scale)") +
  theme(legend.position = "top")
ggsave("fig17-complexity-real.png", p17, width = 11, height = 3.3, dpi = 150, bg = "white")

# ---- 2. FRiP ------------------------------------------------------------------
ctcf_pk <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak")
atac_pk <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
frip <- function(s, pk) {
  r <- readGAlignments(bam(s))
  c(reads = length(r), in_peaks = sum(suppressWarnings(overlapsAny(r, pk, ignore.strand = TRUE))))
}
fr <- rbind(
  data.frame(sample = "CTCF ChIP", peaks = "CTCF peaks", t(frip("ctcf_chip", ctcf_pk))),
  data.frame(sample = "input", peaks = "CTCF peaks", t(frip("ctcf_input", ctcf_pk))),
  data.frame(sample = "ATAC rep1", peaks = "ATAC peaks", t(frip("atac_rep1", atac_pk))),
  data.frame(sample = "ATAC rep2", peaks = "ATAC peaks", t(frip("atac_rep2", atac_pk))),
  data.frame(sample = "input", peaks = "ATAC peaks", t(frip("ctcf_input", atac_pk))))
fr$FRiP <- fr$in_peaks / fr$reads
cov21 <- function(pk) sum(width(reduce(pk[seqnames(pk) == "chr21"]))) / 40.1e6   # assembled chr21, ~40.1 Mb
fr$genome_frac <- ifelse(fr$peaks == "CTCF peaks", cov21(ctcf_pk), cov21(atac_pk))
print(fr, digits = 3)
cat(sprintf("CTCF peaks on chr21: %d; ATAC peaks on chr21: %d\n",
            sum(seqnames(ctcf_pk) == "chr21"), sum(seqnames(atac_pk) == "chr21")))
fr$label <- paste0(fr$sample, "\n(", fr$peaks, ")")
fr$label <- factor(fr$label, levels = fr$label)
fr$fill <- c(col["CTCF"], col["input"], col["ATAC1"], col["ATAC2"], col["input"])
p18 <- ggplot(fr, aes(label, FRiP)) +
  geom_col(aes(fill = I(fill)), width = 0.65) +
  geom_point(aes(y = genome_frac), shape = 95, size = 14, colour = "#c0392b") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * FRiP)), vjust = -0.4, size = 5) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "fraction of reads in peaks",
       subtitle = "red mark: fraction of chr21 covered by the peaks (expected FRiP if reads were uniform)") +
  theme(plot.subtitle = element_text(colour = "#c0392b"))
ggsave("fig18-frip-real.png", p18, width = 10, height = 4.2, dpi = 150, bg = "white")

# ---- 3. fingerprint (1 kb bins over chr21, read 5' ends) ---------------------
bins <- tileGenome(c(chr21 = 46709983), tilewidth = 1000, cut.last.tile.in.chrom = TRUE)
bin_counts <- function(s) {
  r <- granges(readGAlignments(bam(s)))
  countOverlaps(bins, resize(r, 1, fix = "start"))
}
bc <- sapply(c("ctcf_chip", "ctcf_input", "atac_rep1", "atac_rep2"), bin_counts)
fp <- do.call(rbind, lapply(c("ctcf_chip", "ctcf_input", "atac_rep1"), function(s) {
  v <- sort(bc[, s]); data.frame(sample = lab[[s]], x = seq_along(v) / length(v), y = cumsum(v) / sum(v))
}))
fp$sample <- factor(fp$sample, levels = lab)
for (s in c("ctcf_chip", "ctcf_input", "atac_rep1")) {
  v <- sort(bc[, s]); y <- cumsum(v) / sum(v)
  cat(sprintf("fingerprint %s: zero bins %.2f; top 1%% of bins hold %.1f%% of reads; top 10%% hold %.1f%%\n",
              s, mean(v == 0), 100 * (1 - y[ceiling(0.99 * length(v))]), 100 * (1 - y[ceiling(0.9 * length(v))])))
}
fpd <- fp[seq(1, nrow(fp), by = 10), ]
p19 <- ggplot(fpd, aes(x, y, colour = sample)) +
  geom_abline(linetype = 2, colour = "grey60") +
  geom_line(linewidth = 1.3) +
  scale_colour_manual(values = setNames(col[1:3], lab[1:3]), name = NULL) +
  coord_equal() +
  labs(x = "rank of 1 kb bin (fraction, fewest reads first)", y = "cumulative fraction of reads") +
  theme(legend.position = c(0.25, 0.8), legend.background = element_rect(fill = "white", colour = NA))
ggsave("fig19-fingerprint-real.png", p19, width = 6.2, height = 5.6, dpi = 150, bg = "white")

# ---- fig15: ChIP vs input, a window and the whole chromosome ------------------
cpm_bins <- function(s, b) {
  r <- resize(granges(readGAlignments(bam(s))), 1, fix = "start")
  countOverlaps(b, r) / length(r) * 1e6
}
win <- GRanges("chr21", IRanges(33800001, 34200000))
wb <- tile(win, width = 1000)[[1]]
wdf <- rbind(data.frame(pos = start(wb) / 1e6, v = cpm_bins("ctcf_chip", wb), s = "CTCF ChIP"),
             data.frame(pos = start(wb) / 1e6, v = cpm_bins("ctcf_input", wb), s = "input"))
pk_win <- ctcf_pk[overlapsAny(ctcf_pk, win)]
big <- tileGenome(c(chr21 = 46709983), tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)
bdf <- rbind(data.frame(pos = start(big) / 1e6, v = cpm_bins("ctcf_chip", big), s = "CTCF ChIP"),
             data.frame(pos = start(big) / 1e6, v = cpm_bins("ctcf_input", big), s = "input"))
pa <- ggplot(wdf, aes(pos, v, fill = s)) +
  geom_col(width = 0.001, colour = NA) +
  geom_point(data = data.frame(pos = (start(pk_win) + pk_win$peak) / 1e6, s = "CTCF ChIP"),
             aes(pos, -max(wdf$v) * 0.06), shape = 17, size = 2.6, colour = "#c0392b", inherit.aes = FALSE) +
  facet_wrap(~s, ncol = 1) +
  scale_fill_manual(values = c("CTCF ChIP" = col[["CTCF"]], input = col[["input"]]), guide = "none") +
  labs(x = "chr21 position (Mb)", y = "reads per million\nper 1 kb", title = "400 kb window (triangles: ENCODE peaks)")
pb <- ggplot(bdf, aes(pos, v, colour = s)) +
  geom_step(linewidth = 0.7) +
  scale_colour_manual(values = c("CTCF ChIP" = col[["CTCF"]], input = col[["input"]]), name = NULL) +
  labs(x = "chr21 position (Mb)", y = "reads per million\nper 100 kb", title = "all of chr21, 100 kb bins") +
  theme(legend.position = c(0.12, 0.85))
ggsave("fig15-chip-vs-input.png", (pa | pb) & theme(plot.title = element_text(size = 14)), width = 13, height = 4.4, dpi = 150, bg = "white")
cat(sprintf("fig15: ChIP and input CPM in 100 kb bins, correlation %.2f; bins with zero input reads: %d of %d\n",
            cor(bdf$v[bdf$s == "CTCF ChIP"], bdf$v[bdf$s == "input"]), sum(bdf$v[bdf$s == "input"] == 0), length(big)))

# ---- 4. strand cross-correlation ---------------------------------------------
strand_cc <- function(f, shifts = 0:1000, len = 46709983) {
  r <- readGAlignments(f)
  plus <- as.character(strand(r)) == "+"
  fwd <- tabulate(start(r)[plus], len)         # 5' ends of + strand reads
  rev <- tabulate(end(r)[!plus], len)          # 5' ends of - strand reads
  p <- which(fwd > 0)
  m <- mean(fwd) * mean(rev); s <- sd(fwd) * sd(rev)
  sapply(shifts, function(k) {                 # Pearson correlation of fwd(x) and rev(x + k)
    q <- p + k; ok <- q <= len
    (sum(fwd[p[ok]] * rev[q[ok]]) / len - m) / s
  })
}
k <- 0:1000
cc <- list(chip = strand_cc(bam("ctcf_chip")), input = strand_cc(bam("ctcf_input")))
read_len <- c(chip = 36, input = 28)
ccs <- do.call(rbind, lapply(names(cc), function(n) {
  v <- cc[[n]]; frag <- k[k > 50][which.max(v[k > 50])]
  data.frame(sample = n, frag = frag, cc_frag = max(v[k > 50]), cc_read = v[k == read_len[n]],
             min = min(v), NSC = max(v[k > 50]) / min(v),
             RSC = (max(v[k > 50]) - min(v)) / (v[k == read_len[n]] - min(v)))
}))
print(ccs, digits = 3)
ccdf <- rbind(data.frame(k = k, v = cc$chip, s = "CTCF ChIP"), data.frame(k = k, v = cc$input, s = "input"))
ann <- data.frame(s = c("CTCF ChIP", "input"),
                  txt = c(sprintf("fragment peak %d bp\nNSC %.1f   RSC %.2f", ccs$frag[1], ccs$NSC[1], ccs$RSC[1]),
                          sprintf("no fragment peak\nRSC %.2f", ccs$RSC[2])))
vl <- data.frame(s = c("CTCF ChIP", "input"), rl = read_len, fr = c(ccs$frag[1], NA))
p20 <- ggplot(ccdf[ccdf$k <= 500, ], aes(k, v, colour = s)) +
  geom_vline(data = vl, aes(xintercept = rl), linetype = 3, colour = "grey40") +
  geom_vline(data = vl, aes(xintercept = fr), linetype = 2, colour = "#c0392b") +
  geom_line(linewidth = 1.2) +
  geom_text(data = ann, aes(x = 490, y = Inf, label = txt), hjust = 1, vjust = 1.3, size = 4.8,
            colour = "grey20", inherit.aes = FALSE) +
  facet_wrap(~s, scales = "free_y") +
  scale_colour_manual(values = c("CTCF ChIP" = col[["CTCF"]], input = col[["input"]]), guide = "none") +
  labs(x = "strand shift (bp)", y = "cross-correlation")
ggsave("fig20-crosscor-real.png", p20, width = 12, height = 4.2, dpi = 150, bg = "white")

# ---- 5. replicate concordance (ATAC) -----------------------------------------
b10 <- tileGenome(c(chr21 = 46709983), tilewidth = 1e4, cut.last.tile.in.chrom = TRUE)
c10 <- sapply(c("atac_rep1", "atac_rep2", "ctcf_chip"), function(s)
  countOverlaps(b10, resize(granges(readGAlignments(bam(s))), 1, fix = "start")))
keep <- rowSums(c10[, 1:2]) > 0
pk21 <- atac_pk[seqnames(atac_pk) == "chr21"]
cpk <- sapply(c("atac_rep1", "atac_rep2"), function(s)
  countOverlaps(pk21, resize(granges(readGAlignments(bam(s))), 1, fix = "start")))
r10 <- c(pearson_log = cor(log1p(c10[keep, 1]), log1p(c10[keep, 2])),
         spearman = cor(c10[keep, 1], c10[keep, 2], method = "spearman"),
         pearson_log_atac_vs_ctcf = cor(log1p(c10[keep, 1]), log1p(c10[keep, 3])))
r1k <- c(pearson_log_1kb = cor(log1p(bc[, "atac_rep1"]), log1p(bc[, "atac_rep2"])),
         spearman_1kb = cor(bc[, "atac_rep1"], bc[, "atac_rep2"], method = "spearman"))
rpk <- c(pearson_log_peaks = cor(log1p(cpk[, 1]), log1p(cpk[, 2])),
         spearman_peaks = cor(cpk[, 1], cpk[, 2], method = "spearman"))
print(round(c(r10, r1k, rpk), 3)); cat("10 kb bins with reads:", sum(keep), " ATAC peaks on chr21:", length(pk21), "\n")
d10 <- data.frame(x = c10[keep, 1], y = c10[keep, 2])
dpk <- data.frame(x = cpk[, 1], y = cpk[, 2])
p21a <- ggplot(d10, aes(x + 1, y + 1)) + geom_point(alpha = 0.35, size = 1.2, colour = col[["ATAC1"]]) +
  scale_x_log10() + scale_y_log10() + coord_equal() +
  annotate("text", x = 1.3, y = max(d10$y) * 0.8, hjust = 0, size = 5,
           label = sprintf("r = %.2f (log counts)\nSpearman %.2f", r10[1], r10[2])) +
  labs(x = "ATAC rep1 (reads + 1)", y = "ATAC rep2 (reads + 1)", title = "10 kb bins on chr21")
p21b <- ggplot(dpk, aes(x + 1, y + 1)) + geom_point(alpha = 0.35, size = 1.2, colour = col[["ATAC1"]]) +
  scale_x_log10() + scale_y_log10() + coord_equal() +
  annotate("text", x = 1.3, y = max(dpk$y) * 0.8, hjust = 0, size = 5,
           label = sprintf("r = %.2f (log counts)\nSpearman %.2f", rpk[1], rpk[2])) +
  labs(x = "ATAC rep1 (reads + 1)", y = "ATAC rep2 (reads + 1)", title = "ENCODE ATAC peaks on chr21")
ggsave("fig21-replicates-real.png", p21a | p21b, width = 11.5, height = 5, dpi = 150, bg = "white")

# ---- 6. ATAC: fragment sizes and TSS enrichment ------------------------------
frag_sizes <- function(s) {
  r <- readGAlignmentPairs(bam(s), param = ScanBamParam(which = chr21))
  width(granges(r, on.discordant.seqnames = "drop"))
}
fs <- lapply(c(atac_rep1 = "atac_rep1", atac_rep2 = "atac_rep2"), frag_sizes)
fsd <- do.call(rbind, lapply(names(fs), function(s) {
  h <- tabulate(fs[[s]], 1000); data.frame(sample = lab[[s]], len = 1:1000, frac = h / sum(h))
}))
for (s in names(fs)) cat(sprintf("%s fragments %d: < 100 bp %.2f, 180-247 bp %.2f, median %d\n", s, length(fs[[s]]),
                                 mean(fs[[s]] < 100), mean(fs[[s]] >= 180 & fs[[s]] <= 247), median(fs[[s]])))
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
tx <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding" & seqnames(gtf) == "chr21"]
tss <- unique(resize(tx, 1, fix = "start"))
tss_profile <- function(s, flank = 2000) {
  r <- readGAlignments(bam(s))
  plus <- as.character(strand(r)) == "+"
  cuts <- c(start(r)[plus] + 4, end(r)[!plus] - 5)          # Tn5 insertion sites
  cv <- tabulate(cuts, 46709983)
  w <- -flank:flank
  m <- sapply(seq_along(tss), function(i) {
    p <- start(tss)[i]; if (as.character(strand(tss)[i]) == "+") cv[p + w] else rev(cv[p + w])
  })
  prof <- rowSums(m)
  bg <- mean(prof[c(1:100, (length(w) - 99):length(w))])  # the outer 100 bp on each side
  data.frame(sample = lab[[s]], pos = w, enr = as.numeric(stats::filter(prof / bg, rep(1 / 51, 51))))
}
tp <- rbind(tss_profile("atac_rep1"), tss_profile("atac_rep2"))
tse <- tapply(tp$enr, tp$sample, max, na.rm = TRUE)
print(round(tse, 1)); cat("protein-coding TSSs on chr21:", length(tss), "\n")
fsd$sample <- factor(fsd$sample, levels = lab[3:4]); tp$sample <- factor(tp$sample, levels = lab[3:4])
p22a <- ggplot(fsd[fsd$len <= 800, ], aes(len, frac, colour = sample)) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = setNames(col[3:4], lab[3:4]), name = NULL) +
  annotate("text", x = c(75, 200, 400), y = c(0.0072, 0.0038, 0.0021),
           label = c("nucleosome-free", "mono-nucleosome", "di-"), size = 4.5, hjust = 0) +
  labs(x = "fragment length (bp)", y = "fraction of fragments", title = "Fragment sizes") +
  theme(legend.position = c(0.8, 0.8))
p22b <- ggplot(tp, aes(pos, enr, colour = sample)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = setNames(col[3:4], lab[3:4]), guide = "none") +
  annotate("text", x = 1950, y = max(tp$enr, na.rm = TRUE) * 0.95, hjust = 1, size = 5,
           label = sprintf("TSS enrichment\nrep1 %.1f, rep2 %.1f", tse[1], tse[2])) +
  labs(x = "distance to TSS (bp)", y = "Tn5 insertions\n(relative to flanks)", title = "TSS enrichment")
ggsave("fig22-atac-qc-real.png", p22a | p22b, width = 12.5, height = 4.4, dpi = 150, bg = "white")

# ---- 7. depth and saturation: peaks called on subsamples (MACS3) -------------
macs_peaks <- function(s, frac, ctrl = NULL, pe = FALSE) {
  f <- file.path(W, "sub.bam"); out <- file.path(W, "sat")
  if (frac < 1) sh(sprintf("samtools view -b -s %d%s %s > %s", 7, sub("^0", "", sprintf("%.2f", frac)), bam(s), f))
  else file.copy(bam(s), f, overwrite = TRUE)
  sh(sprintf("macs3 callpeak -t %s %s -f %s -g 40e6 -q 0.01 -n sat --outdir %s 2> /dev/null",
             f, if (is.null(ctrl)) "" else paste("-c", bam(ctrl)), if (pe) "BAMPE" else "BAM", out))
  pk <- import(file.path(out, "sat_peaks.narrowPeak"), format = "narrowPeak")
  c(reads = count(f), peaks = length(pk))
}
fracs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
sat <- rbind(
  data.frame(sample = "CTCF ChIP", frac = fracs, t(sapply(fracs, function(f) macs_peaks("ctcf_chip", f, ctrl = "ctcf_input")))),
  data.frame(sample = "ATAC rep1", frac = fracs, t(sapply(fracs, function(f) macs_peaks("atac_rep1", f, pe = TRUE)))))
print(sat)
sat$sample <- factor(sat$sample, levels = c("CTCF ChIP", "ATAC rep1"))
p16 <- ggplot(sat, aes(reads / 1e3, peaks, colour = sample)) + geom_line(linewidth = 1.1) + geom_point(size = 2.4) +
  facet_wrap(~sample, scales = "free") +
  scale_colour_manual(values = c("CTCF ChIP" = col[["CTCF"]], "ATAC rep1" = col[["ATAC1"]]), guide = "none") +
  expand_limits(y = 0, x = 0) +
  labs(x = "reads used, chr21 (thousands)", y = "MACS3 peaks on chr21 (q < 0.01)")
ggsave("fig16-depth-saturation.png", p16, width = 11.5, height = 4.2, dpi = 150, bg = "white")
