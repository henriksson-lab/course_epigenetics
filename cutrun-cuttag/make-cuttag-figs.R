# Real-data figures for cutrun-cuttag/lec.qmd.
# Run from this folder:  Rscript make-cuttag-figs.R
#
# Needs: data/prepare-cuttag.sh (CUT&Tag reads k562_cuttag_*_R{1,2}.fastq.gz, ecoli_U00096.3.fa,
#        ENCODE K562 ChIP alignments k562_*_chr21.bam), data/download.sh cutrun-cuttag
#        (ENCODE K562 peaks, GENCODE GTF, chromosome sizes) and the chr21 + chrM Bowtie2 index;
#        on the PATH: bowtie2, samtools (>= 1.16), bedtools, macs3, R; SEACR_1.3.sh
#        (github.com/FredHutch/SEACR; set SEACR=/path/to/SEACR_1.3.sh);
#        R packages GenomicRanges, rtracklayer, GenomicAlignments, ggplot2, patchwork.
#
# The processing (the same commands as in exercise.qmd) runs once into ./cuttag-work
# (set CUTTAG_WORK to change; ~100 MB, under a minute). Re-running only redoes missing steps.
#
# Writes fig09-tracks.png, fig10-fragment-sizes.png, fig11-duplicates-real.png, fig12-igg-real.png,
# fig13-peaks-compare.png, fig14-fingerprint.png, fig15-depth.png, fig16-bias-real.png
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(GenomicAlignments)
  library(ggplot2); library(patchwork)
})
D <- normalizePath("../data")
W <- Sys.getenv("CUTTAG_WORK", "cuttag-work")
dir.create(W, showWarnings = FALSE, recursive = TRUE); W <- normalizePath(W)
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
SEACR <- Sys.getenv("SEACR", "SEACR/SEACR_1.3.sh")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
col <- c("CUT&Tag H3K4me3" = "#c0392b", "ChIP H3K4me3" = "#e8998f",
         "CUT&Tag H3K27me3" = "#2a6f97", "ChIP H3K27me3" = "#8fb8d4",
         "CUT&Tag IgG" = "#5b6b7a", "ChIP input" = "#b8c2cc")
L21 <- 46709983

# ---- 1. processing (exercise C1-C4) -----------------------------------------
step <- function(out, cmd) {
  if (all(file.exists(file.path(W, out)))) return(invisible())
  cat("running:", out, "\n")
  rc <- system2("bash", c("-c", shQuote(paste0("set -euo pipefail; cd ", W, "; D=", D, "; T=", THREADS,
                                                "; SEACR=", SEACR, "; ", cmd))))
  if (rc != 0) stop("step failed: ", out)
}
S3 <- c("h3k27me3", "h3k4me3", "igg")
step("ecoli.1.bt2", "bowtie2-build --threads $T $D/ecoli_U00096.3.fa ecoli > /dev/null")
step("chrom.sizes", "grep -w chr21 $D/hg38.chrom.sizes > chrom.sizes")
for (s in S3) {
  step(paste0(s, ".raw.bam.bai"), sprintf(paste(
    "bowtie2 -p $T --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700",
    "-x $D/bowtie2/hg38_chr21_chrM -1 $D/k562_cuttag_%1$s_R1.fastq.gz -2 $D/k562_cuttag_%1$s_R2.fastq.gz",
    "2> %1$s.bowtie2.log | samtools fixmate -m -u - - | samtools sort -@ 4 -o %1$s.raw.bam;",
    "samtools index %1$s.raw.bam;",
    "bowtie2 -p $T --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant",
    "--phred33 -I 10 -X 700 -x ecoli -1 $D/k562_cuttag_%1$s_R1.fastq.gz -2 $D/k562_cuttag_%1$s_R2.fastq.gz",
    "2> %1$s.ecoli.log > /dev/null"), s))
  step(paste0(s, ".bam.bai"), sprintf(paste(
    "samtools view -u -q 30 -f 2 -F 1804 %1$s.raw.bam chr21 | samtools collate -u -O - |",
    "samtools fixmate -m -u - - | samtools view -u -f 2 - | samtools sort -u - |",
    "samtools markdup -s - %1$s.bam 2> %1$s.markdup.txt; samtools index %1$s.bam"), s))
  step(paste0(s, ".bedgraph"), sprintf(paste(
    "F=0; [ %1$s = igg ] && F=1024;",                      # drop duplicates from the IgG only
    "samtools view -u -F $F %1$s.bam | samtools sort -n -u - | bedtools bamtobed -bedpe -i - 2> /dev/null",
    "| awk -v OFS='\\t' '$1 == $4 { s = ($2 < $5 ? $2 : $5); e = ($3 > $6 ? $3 : $6); print $1, s, e }'",
    "| sort -k1,1 -k2,2n > %1$s.frag.bed;",
    "bedtools genomecov -bg -i %1$s.frag.bed -g chrom.sizes > %1$s.bedgraph"), s))
}
for (s in c("h3k27me3", "h3k4me3")) {
  step(sprintf("peaks/%s_seacr_igg.stringent.bed", s), sprintf(
    "mkdir -p peaks; bash $SEACR %1$s.bedgraph igg.bedgraph norm stringent peaks/%1$s_seacr_igg > /dev/null", s))
  step(sprintf("peaks/%s_seacr_top1.stringent.bed", s), sprintf(
    "bash $SEACR %1$s.bedgraph 0.01 non stringent peaks/%1$s_seacr_top1 > /dev/null", s))
}
step("peaks/h3k4me3_macs_peaks.narrowPeak", paste(
  "macs3 callpeak -t h3k4me3.bam -c igg.bam -f BAMPE -g 4.0e7 --keep-dup all -q 0.01",
  "-n h3k4me3_macs --outdir peaks 2> peaks/h3k4me3_macs.log"))
step("peaks/h3k27me3_macs_peaks.broadPeak", paste(
  "macs3 callpeak -t h3k27me3.bam -c igg.bam -f BAMPE -g 4.0e7 --keep-dup all --broad --broad-cutoff 0.1",
  "-n h3k27me3_macs --outdir peaks 2> peaks/h3k27me3_macs.log"))

# ---- 2. read everything --------------------------------------------------------
fr <- lapply(setNames(S3, S3), function(s) import(file.path(W, paste0(s, ".frag.bed"))))
rd <- function(f) granges(readGAlignments(file.path(D, paste0(f, "_chr21.bam"))))
chip <- list(h3k4me3 = rd("k562_h3k4me3_chip"), h3k27me3 = rd("k562_h3k27me3_chip"), input = rd("k562_input"))
mid <- function(g) resize(g, 1, fix = "center")      # count a fragment / read once, at its centre
enc <- lapply(c(h3k4me3 = "h3k4me3", h3k27me3 = "h3k27me3"), function(m) {
  x <- import(file.path(D, sprintf("k562_%s_chip_peaks.bed.gz", m)), format = "narrowPeak")
  x[seqnames(x) == "chr21"]
})
encr <- lapply(enc, reduce)
sea <- function(f) { x <- read.table(file.path(W, f)); GRanges(x$V1, IRanges(x$V2 + 1, x$V3), auc = x$V4) }
lib <- c(h3k27me3 = 2984630, h3k4me3 = 1581710, igg = 2127635)    # read pairs per library (GEO/SRA)
num <- function(f, pat) as.numeric(sub("^\\s*(\\d+).*", "\\1", grep(pat, readLines(file.path(W, f)), value = TRUE)[1]))

cat("\n== alignment, duplicates, E. coli\n")
tab <- t(sapply(S3, function(s) {
  pairs <- num(paste0(s, ".bowtie2.log"), "reads; of these")
  eco <- num(paste0(s, ".ecoli.log"), "concordantly exactly 1 time") + num(paste0(s, ".ecoli.log"), "concordantly >1 times")
  idx <- read.table(text = system2("samtools", c("idxstats", file.path(W, paste0(s, ".raw.bam"))), stdout = TRUE))
  frag <- as.numeric(system2("samtools", c("view", "-c", "-f", "64", file.path(W, paste0(s, ".bam"))), stdout = TRUE))
  dup <- as.numeric(system2("samtools", c("view", "-c", "-f", "64", "-f", "1024", file.path(W, paste0(s, ".bam"))), stdout = TRUE))
  c(library = lib[[s]], subset_pairs = pairs, chr21_fragments = frag, duplicates = dup, dup_pct = 100 * dup / frag,
    ecoli = eco, ecoli_pct = 100 * eco / lib[[s]], scale = 10000 / eco,
    chrM_pct = 100 * idx$V3[idx$V1 == "chrM"] / 2 / lib[[s]])
}))
print(round(tab, 3))
cat("ChIP reads on chr21:", sapply(chip, length), "\n")

# ---- fig09: tracks at a locus --------------------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
gtf <- gtf[seqnames(gtf) == "chr21"]
region <- GRanges("chr21", IRanges(33000001, 33500000))
bins <- tile(region, width = 2000)[[1]]
trk <- function(g, name) {
  n <- countOverlaps(bins, mid(g)); data.frame(track = name, x = start(bins), cpm = n / length(g) * 1e6)
}
td <- rbind(trk(fr$h3k4me3, "CUT&Tag H3K4me3"), trk(chip$h3k4me3, "ChIP H3K4me3"),
            trk(fr$h3k27me3, "CUT&Tag H3K27me3"), trk(chip$h3k27me3, "ChIP H3K27me3"),
            trk(fr$igg, "CUT&Tag IgG"), trk(chip$input, "ChIP input"))
td$track <- factor(td$track, levels = names(col))
tx <- subsetByOverlaps(gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"], region)
gd <- data.frame(s = pmax(start(tx), start(region)), e = pmin(end(tx), end(region)), g = tx$gene_name,
                 row = rep(1:2, length.out = length(tx)))
p_tr <- ggplot(td, aes(x, cpm, fill = track)) + geom_col(width = 2000, colour = NA) +
  facet_grid(track ~ ., scales = "free_y", switch = "y") + scale_fill_manual(values = col, guide = "none") +
  scale_x_continuous(limits = c(start(region), end(region)), expand = c(0, 0), labels = NULL) +
  scale_y_continuous(n.breaks = 3) + labs(x = NULL, y = NULL) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
        strip.placement = "outside", panel.grid.major.x = element_blank(), axis.text.y = element_text(size = 9),
        plot.margin = margin(5, 40, 5, 5))
p_g <- ggplot(gd) + geom_segment(aes(x = s, xend = e, y = -row, yend = -row), linewidth = 2.5, colour = "#1b2a3a") +
  geom_text(aes(x = (s + e) / 2, y = -row - 0.6, label = g), size = 3.6) +
  scale_x_continuous(limits = c(start(region), end(region)), expand = c(0, 0), labels = function(b) sprintf("%.2f Mb", b / 1e6)) +
  scale_y_continuous(limits = c(-3, -0.5), breaks = NULL) + labs(x = "chr21 (2-kb bins, counts per million)", y = NULL) +
  coord_cartesian(clip = "off") + theme(panel.grid = element_blank(), plot.margin = margin(5, 40, 5, 5))
ggsave("fig09-tracks.png", p_tr / p_g + plot_layout(heights = c(7, 1.4)), width = 12, height = 6.4, dpi = 150, bg = "white")
cat("\n== fig09 locus", as.character(region), "genes:", paste(unique(tx$gene_name), collapse = " "), "\n")

# ---- fig10: fragment sizes -----------------------------------------------------
fs <- do.call(rbind, lapply(S3, function(s) data.frame(sample = s, len = width(fr[[s]]))))
fs$sample <- factor(fs$sample, levels = S3, labels = c("H3K27me3", "H3K4me3", "IgG"))
h <- do.call(rbind, lapply(split(fs, fs$sample), function(d) {
  t <- table(factor(pmin(d$len, 800), levels = 1:800)); data.frame(sample = d$sample[1], len = 1:800, f = as.numeric(t) / nrow(d))
}))
h$sm <- ave(h$f, h$sample, FUN = function(v) stats::filter(v, rep(1 / 5, 5), sides = 2))
fcol <- c(H3K27me3 = "#2a6f97", H3K4me3 = "#c0392b", IgG = "#5b6b7a")
p10a <- ggplot(h[h$len < 700, ], aes(len, sm * 100, colour = sample)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = fcol) + labs(x = "fragment length (bp)", y = "% of fragments (5-bp smoothed)", colour = NULL,
                                            title = "A. Nucleosome ladder") + theme(legend.position = c(0.8, 0.8))
# 10-bp periodicity: detrend the 1-bp histogram of 140-300 bp fragments, take the strongest period
per <- function(v) { x <- tabulate(v, 1000)[140:300]; tr <- stats::filter(x, rep(1 / 21, 21), sides = 2)
  r <- (x - tr)[11:150]; sp <- spec.pgram(r, plot = FALSE, taper = 0.1, detrend = TRUE); 1 / sp$freq[which.max(sp$spec)] }
period <- per(width(fr$h3k27me3))
h1 <- h[h$sample == "H3K27me3" & h$len >= 140 & h$len <= 300, ]
h1$s3 <- stats::filter(h1$f, rep(1 / 3, 3), sides = 2)
p10b <- ggplot(h1, aes(len, s3 * 100)) + geom_line(colour = "#2a6f97", linewidth = 0.9) +
  scale_x_continuous(breaks = seq(140, 300, 20)) +
  labs(x = "fragment length (bp)", y = "% of fragments (3-bp smoothed)",
       title = sprintf("B. H3K27me3: a %.1f-bp ripple", period))
ggsave("fig10-fragment-sizes.png", p10a + p10b + plot_layout(widths = c(1.3, 1)), width = 12, height = 4.8, dpi = 150, bg = "white")
cat("\n== fragment lengths (median, % < 120, % 120-300)\n")
print(sapply(split(fs$len, fs$sample), function(v) round(c(median = median(v), lt120 = 100 * mean(v < 120),
                                                         mono = 100 * mean(v >= 120 & v < 300), di = 100 * mean(v >= 300)), 1)))
pk <- sapply(split(h, h$sample), function(d) d$len[d$len > 100][which.max(d$sm[d$len > 100])])
cat("mode > 100 bp:", pk, "\n")
cat("periodicity 140-300 bp: H3K27me3", round(period, 1), "bp; H3K4me3", round(per(width(fr$h3k4me3)), 1), "bp\n")

# ---- fig11: duplicates and mitochondria ---------------------------------------
dd <- data.frame(sample = factor(c("H3K27me3", "H3K4me3", "IgG"), levels = c("H3K27me3", "H3K4me3", "IgG")),
                 dup = tab[, "dup_pct"], chrM = tab[, "chrM_pct"], eco = tab[, "ecoli_pct"])
bar <- function(y, lab, ttl) ggplot(dd, aes(sample, .data[[y]], fill = sample)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf(lab, .data[[y]])), vjust = -0.4, size = 4.6) +
  scale_fill_manual(values = fcol, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL, title = ttl)
p11 <- bar("dup", "%.0f%%", "A. Duplicates, % (chr21)") + bar("chrM", "%.2f%%", "B. chrM, % of pairs") +
  bar("eco", "%.2f%%", "C. E. coli, % of pairs")
ggsave("fig11-duplicates-real.png", p11, width = 12, height = 4.4, dpi = 150, bg = "white")

# ---- TSS classes for fig12 and fig16 -------------------------------------------
txp <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding"]
tss <- unique(resize(txp, 1, fix = "start"))
act <- tss[overlapsAny(tss, encr$h3k4me3, maxgap = 500) & !overlapsAny(tss, encr$h3k27me3, maxgap = 2000)]
rep <- tss[overlapsAny(tss, encr$h3k27me3, maxgap = 500)]
cat("\n== TSS:", length(tss), "protein-coding TSS positions;", length(act), "active (H3K4me3, no H3K27me3);",
    length(rep), "in H3K27me3 peaks\n")
S6 <- list("CUT&Tag H3K4me3" = fr$h3k4me3, "ChIP H3K4me3" = chip$h3k4me3, "CUT&Tag H3K27me3" = fr$h3k27me3,
           "ChIP H3K27me3" = chip$h3k27me3, "CUT&Tag IgG" = fr$igg, "ChIP input" = chip$input)
profile <- function(g, t, flank = 5000, bw = 500) {
  c_ <- mid(g); st <- as.character(strand(t))
  br <- seq(-flank, flank, bw)
  ov <- findOverlaps(resize(t, 2 * flank, fix = "center"), c_)
  d <- start(c_)[subjectHits(ov)] - start(t)[queryHits(ov)]
  d <- ifelse(st[queryHits(ov)] == "-", -d, d)
  n <- as.numeric(table(cut(d, br, include.lowest = TRUE)))
  data.frame(pos = head(br, -1) + bw / 2, enr = (n / length(t) / bw) / (length(g) / L21))
}
pr <- do.call(rbind, lapply(c("CUT&Tag IgG", "ChIP input"), function(n)
  cbind(track = n, profile(S6[[n]], act))))
pr$track <- factor(pr$track, levels = names(col))
p12 <- ggplot(pr, aes(pos / 1000, enr, colour = track)) + geom_hline(yintercept = 1, colour = "grey60") +
  geom_line(linewidth = 1.1) + scale_colour_manual(values = col) + expand_limits(y = 0) +
  labs(x = "distance to TSS (kb)", y = "enrichment over a uniform genome", colour = NULL,
       title = sprintf("%d active TSSs on chr21 (H3K4me3, no H3K27me3)", length(act))) +
  theme(legend.position = "right")
ggsave("fig12-igg-real.png", p12, width = 11, height = 4.8, dpi = 150, bg = "white")

enr <- function(g, t) { r <- reduce(resize(t, 1000, fix = "center")); mean(overlapsAny(mid(g), r)) / (sum(width(r)) / L21) }
E <- sapply(S6, function(g) c(active = enr(g, act), repressed = enr(g, rep)))
cat("\n== enrichment in 1 kb around TSS classes\n"); print(round(E, 2))
cat("IgG / input at active TSS, peak of profile:", round(tapply(pr$enr, pr$track, max)[c("CUT&Tag IgG", "ChIP input")], 1), "\n")

# ---- fig13: SEACR vs MACS3 -----------------------------------------------------
calls <- list("H3K4me3 SEACR (IgG)" = sea("peaks/h3k4me3_seacr_igg.stringent.bed"),
              "H3K4me3 SEACR (top 1%)" = sea("peaks/h3k4me3_seacr_top1.stringent.bed"),
              "H3K4me3 MACS3" = import(file.path(W, "peaks/h3k4me3_macs_peaks.narrowPeak"), format = "narrowPeak"),
              "H3K27me3 SEACR (IgG)" = sea("peaks/h3k27me3_seacr_igg.stringent.bed"),
              "H3K27me3 SEACR (top 1%)" = sea("peaks/h3k27me3_seacr_top1.stringent.bed"),
              "H3K27me3 MACS3 broad" = import(file.path(W, "peaks/h3k27me3_macs_peaks.broadPeak"), format = "broadPeak"))
pc <- do.call(rbind, lapply(names(calls), function(n) {
  p <- calls[[n]]; m <- if (grepl("H3K4me3", n)) "h3k4me3" else "h3k27me3"
  data.frame(call = n, mark = sub(" .*", "", n), caller = sub("^\\S+ ", "", n), peaks = length(p),
             median_width = median(width(p)), in_encode = 100 * mean(overlapsAny(p, encr[[m]])),
             encode_found = 100 * mean(overlapsAny(encr[[m]], p)), frip = 100 * mean(overlapsAny(mid(fr[[m]]), p)),
             igg_frip = 100 * mean(overlapsAny(mid(fr$igg), p)))
}))
cat("\n== peak calls\n"); print(pc[, -1], digits = 3)
pl <- reshape(pc[, c("call", "mark", "caller", "peaks", "median_width", "in_encode", "encode_found", "frip")],
              direction = "long", varying = list(c("peaks", "median_width", "in_encode", "encode_found", "frip")),
              v.names = "value", timevar = "metric", times = c("number of peaks", "median width (bp)",
                                                              "% of peaks in ENCODE peaks", "% of ENCODE peaks found",
                                                              "FRiP (%)"))
pl$metric <- factor(pl$metric, levels = unique(pl$metric))
pl$mark <- factor(pl$mark, levels = c("H3K4me3", "H3K27me3"))
pl$caller <- factor(pl$caller, levels = c("SEACR (IgG)", "SEACR (top 1%)", "MACS3", "MACS3 broad"),
                    labels = c("SEACR\n(IgG)", "SEACR\n(top 1%)", "MACS3", "MACS3"))
p13 <- ggplot(pl, aes(caller, value, fill = mark)) + geom_col(width = 0.7) +
  geom_text(aes(label = round(value)), vjust = -0.3, size = 3.4) +
  facet_grid(metric ~ mark, scales = "free", switch = "y") +
  scale_fill_manual(values = c(H3K4me3 = "#c0392b", H3K27me3 = "#2a6f97"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + labs(x = NULL, y = NULL) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11), strip.placement = "outside",
        strip.text.x = element_text(face = "bold", size = 13), axis.text.y = element_blank(), panel.grid = element_blank())
ggsave("fig13-peaks-compare.png", p13, width = 11, height = 6.2, dpi = 150, bg = "white")

# ---- fig14: fingerprint and FRiP -------------------------------------------------
b5 <- tile(GRanges("chr21", IRanges(1, L21)), width = 5000)[[1]]
fp <- do.call(rbind, lapply(names(S6), function(n) {
  cnt <- sort(countOverlaps(b5, mid(S6[[n]])))
  data.frame(track = n, x = seq_along(cnt) / length(cnt), y = cumsum(cnt) / sum(cnt))
}))
fp$track <- factor(fp$track, levels = names(col))
fp <- fp[fp$x > 0.3 | round(fp$x * 1e4) %% 10 == 0, ]
p14a <- ggplot(fp, aes(x, y, colour = track)) + geom_abline(linetype = 2, colour = "grey70") + geom_line(linewidth = 1) +
  scale_colour_manual(values = col) + labs(x = "fraction of 5-kb bins (sorted by count)", y = "cumulative fraction of reads",
                                           colour = NULL, title = "A. Fingerprint, chr21") +
  theme(legend.position = c(0.3, 0.68), legend.text = element_text(size = 11))
top1 <- sapply(names(S6), function(n) { cnt <- sort(countOverlaps(b5, mid(S6[[n]])), decreasing = TRUE)
  sum(cnt[seq_len(ceiling(0.01 * length(cnt)))]) / sum(cnt) })
frip_enc <- c(sapply(names(S6)[1:2], function(n) mean(overlapsAny(mid(S6[[n]]), encr$h3k4me3))),
              sapply(names(S6)[3:4], function(n) mean(overlapsAny(mid(S6[[n]]), encr$h3k27me3))),
              "CUT&Tag IgG" = mean(overlapsAny(mid(fr$igg), encr$h3k4me3)), "ChIP input" = mean(overlapsAny(mid(chip$input), encr$h3k4me3)))
fd <- data.frame(track = factor(names(frip_enc), levels = names(col)), frip = 100 * frip_enc)
p14b <- ggplot(fd, aes(track, frip, fill = track)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", frip)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = col, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "% of reads in ENCODE peaks", title = "B. FRiP in the same (ENCODE) peaks") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11))
ggsave("fig14-fingerprint.png", p14a + p14b + plot_layout(widths = c(1.1, 1)), width = 12, height = 5, dpi = 150, bg = "white")
cat("\n== FRiP in ENCODE peaks (IgG, input: H3K4me3 peaks)\n"); print(round(100 * frip_enc, 1))
cat("ENCODE chr21 peaks:", sapply(enc, length), "regions after merging:", sapply(encr, length),
    "% of chr21:", round(100 * sapply(encr, function(x) sum(width(x))) / L21, 2), "\n")
cat("H3K27me3 FRiP in ENCODE H3K27me3 peaks: IgG", round(100 * mean(overlapsAny(mid(fr$igg), encr$h3k27me3)), 1),
    "input", round(100 * mean(overlapsAny(mid(chip$input), encr$h3k27me3)), 1), "\n")
cat("fraction of reads in the top 1% of 5-kb bins:\n"); print(round(top1, 3))

# ---- fig15: reads needed ------------------------------------------------------------
detect <- function(g, p, n, seed) {
  set.seed(seed); g <- mid(g[sample(length(g), n)])
  k <- countOverlaps(p, g); bg <- sum(!overlapsAny(g, p)) / (L21 - sum(width(p)))
  mean(ppois(k - 1, bg * width(p), lower.tail = FALSE) < 1e-3)
}
grid <- c(1000, 2000, 5000, 10000, 20000, 45000, 90000, 140000)
dp <- do.call(rbind, lapply(c("h3k4me3", "h3k27me3"), function(m) do.call(rbind, lapply(c("CUT&Tag", "ChIP"), function(a) {
  g <- if (a == "CUT&Tag") fr[[m]] else chip[[m]]
  do.call(rbind, lapply(grid[grid <= length(g)], function(n)
    data.frame(mark = ifelse(m == "h3k4me3", "H3K4me3", "H3K27me3"), assay = a, n = n,
               recall = mean(sapply(1:5, function(sd) detect(g, encr[[m]], n, sd))))))
}))))
dp$mark <- factor(dp$mark, levels = c("H3K4me3", "H3K27me3"))
p15 <- ggplot(dp, aes(n, 100 * recall, colour = assay, shape = assay)) + geom_line(linewidth = 1) + geom_point(size = 2.5) +
  facet_wrap(~ mark) + scale_x_log10(labels = function(b) format(b, big.mark = " ", scientific = FALSE)) +
  scale_colour_manual(values = c("CUT&Tag" = "#1b2a3a", ChIP = "#e67e22")) +
  labs(x = "fragments / reads on chr21 (subsampled)", y = "% of ENCODE peaks detected", colour = NULL, shape = NULL) +
  theme(legend.position = "top", strip.text = element_text(face = "bold", size = 14))
ggsave("fig15-depth.png", p15, width = 11, height = 4.8, dpi = 150, bg = "white")
cat("\n== recall of ENCODE peaks vs depth (mean of 5 subsamples)\n"); print(reshape(dp, idvar = c("mark", "assay"), timevar = "n", direction = "wide"), digits = 2)
q <- enc$h3k4me3; q$ct <- countOverlaps(q, mid(fr$h3k4me3)); q$chip <- countOverlaps(q, mid(chip$h3k4me3))
q$quart <- cut(rank(-q$signalValue, ties.method = "first"), 4, labels = c("strongest", "2nd", "3rd", "weakest"))
cat("H3K4me3 ENCODE peaks by ChIP strength: % with >= 3 CUT&Tag fragments\n"); print(round(100 * tapply(q$ct >= 3, q$quart, mean)))
cat("median ChIP reads per peak:\n"); print(tapply(q$chip, q$quart, median))
cat("% within 1 kb of a protein-coding TSS:\n"); print(round(100 * tapply(overlapsAny(q, tss, maxgap = 1000), q$quart, mean)))

# ---- fig16: H3K27me3 at active and repressed promoters ------------------------------
ed <- data.frame(track = factor(rep(colnames(E), each = 2), levels = names(col)),
                 class = factor(rep(c("active TSS", "H3K27me3 TSS"), ncol(E)), levels = c("active TSS", "H3K27me3 TSS")),
                 enr = as.vector(E))
ed <- ed[!grepl("H3K4me3", ed$track), ]
p16 <- ggplot(ed, aes(track, enr, fill = track)) + geom_hline(yintercept = 1, colour = "grey50") + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2g", enr)), vjust = ifelse(ed$enr >= 1, -0.3, 1.3), size = 4.2) +
  facet_wrap(~ class) + scale_y_log10(expand = expansion(mult = 0.15)) + scale_fill_manual(values = col, guide = "none") +
  labs(x = NULL, y = "enrichment in TSS ± 500 bp (log scale)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), strip.text = element_text(face = "bold", size = 14))
ggsave("fig16-bias-real.png", p16, width = 11, height = 4.8, dpi = 150, bg = "white")
cat("\nDone.\n")
