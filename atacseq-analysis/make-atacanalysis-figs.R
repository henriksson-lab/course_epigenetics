# Real-data figures for atacseq-analysis/lec.qmd.
# Run from this folder:  Rscript make-atacanalysis-figs.R
#
# Needs: the ATAC-seq reads data/gm12878_atac_rep{1,2}_R{1,2}.fastq.gz (data/prepare-reads.sh),
#        data/hg38_chr21_chrM.fa with its bowtie2 index data/bowtie2/hg38_chr21_chrM,
#        data/download.sh atacseq-analysis genomic-intervals (ENCODE peaks, blacklist, GTF);
#        on the PATH: fastp, bowtie2, samtools (>= 1.16), bedtools, macs3;
#        R packages GenomicRanges, rtracklayer, ggplot2, patchwork.
#
# The processing (the same commands as in exercise.qmd) runs once into a work folder,
# by default ./atac-work (set ATAC_WORK to put it elsewhere; about 3 GB, 15-30 min on 8 cores).
# Re-running the script only redoes missing steps.
#
# Writes fig04-reads-per-step.png, fig08-insertions-vs-fragments-real.png,
# fig11-poisson-test.png, fig13-peaks-vs-threshold.png, fig16-rank-consistency.png,
# fig17-peaks-locus.png, fig21-promoter-fraction.png and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(ggplot2); library(patchwork)
})
D <- normalizePath("../data")
W <- Sys.getenv("ATAC_WORK", "atac-work")
dir.create(W, showWarnings = FALSE, recursive = TRUE)
W <- normalizePath(W)
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
col_rep <- c(rep1 = "#2a6f97", rep2 = "#e67e22")
col_ins <- "#c0392b"; col_frag <- "#2a6f97"; col_enc <- "#5b6b7a"; col_ok <- "#2e8b57"

# ---- 1. processing: run a bash step unless its output exists ----------------
step <- function(out, cmd) {
  if (all(file.exists(file.path(W, out)))) return(invisible())
  cat("running:", out, "\n")
  t0 <- Sys.time()
  rc <- system2("bash", c("-c", shQuote(paste0("set -euo pipefail; cd ", W, "; D=", D, "; T=", THREADS, "; ", cmd))))
  if (rc != 0) stop("step failed: ", out)
  cat(sprintf("  %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
reps <- c("rep1", "rep2")
for (r in reps) {
  step(paste0(r, ".raw.bam.bai"), sprintf(paste(
    "fastp -w $T -l 20 -i $D/gm12878_atac_%1$s_R1.fastq.gz -I $D/gm12878_atac_%1$s_R2.fastq.gz",
    "-o %1$s_R1.trim.fq.gz -O %1$s_R2.trim.fq.gz",
    "--adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT",
    "-j %1$s.fastp.json -h %1$s.fastp.html 2> %1$s.fastp.log;",
    "bowtie2 -p $T --very-sensitive -X 2000 -x $D/bowtie2/hg38_chr21_chrM",
    "-1 %1$s_R1.trim.fq.gz -2 %1$s_R2.trim.fq.gz 2> %1$s.bowtie2.log",
    "| samtools fixmate -m -u - - | samtools sort -@ 4 -o %1$s.raw.bam;",
    "samtools index %1$s.raw.bam"), r))
  step(paste0(r, ".dedup.bam.bai"), sprintf(paste(
    "samtools view -u -q 30 -f 2 -F 1804 %1$s.raw.bam chr21",
    "| samtools collate -u -O - | samtools fixmate -m -u - -",
    "| samtools view -u -f 2 - | samtools sort -u -",
    "| samtools markdup -r -s - %1$s.dedup.bam 2> %1$s.markdup.txt;",
    "samtools index %1$s.dedup.bam"), r))
  step(paste0(r, ".ins.bed"), sprintf(paste(
    "zcat $D/hg38_blacklist.bed.gz > blacklist.bed;",
    "samtools sort -n -u %1$s.dedup.bam | bedtools bamtobed -bedpe -i - 2> /dev/null",
    "| awk -v OFS='\\t' '{ s = ($2 < $5 ? $2 : $5); e = ($3 > $6 ? $3 : $6);",
    "print $1, s + 4, s + 5, \".\", \".\", \"+\"; print $1, e - 5, e - 4, \".\", \".\", \"+\" }'",
    "| bedtools intersect -v -a - -b blacklist.bed | sort -k1,1 -k2,2n > %1$s.ins.bed"), r))
  step(paste0(r, ".frag.bed"), sprintf(paste(
    "samtools sort -n -u %1$s.dedup.bam | bedtools bamtobed -bedpe -i - 2> /dev/null",
    "| awk -v OFS='\\t' '{ s = ($2 < $5 ? $2 : $5); e = ($3 > $6 ? $3 : $6); print $1, s, e }'",
    "| sort -k1,1 -k2,2n > %1$s.frag.bed"), r))
  macs <- "macs3 callpeak -g 4.0e7 --keep-dup all --outdir macs"
  cuts <- "-f BED --nomodel --shift -75 --extsize 150"
  step(sprintf("macs/%s_frag_peaks.narrowPeak", r),
       sprintf("%s -t %s.dedup.bam -f BAMPE -q 0.01 -n %s_frag 2> %s_frag.macs.log", macs, r, r, r))
  step(sprintf("macs/%s_cuts_peaks.narrowPeak", r),
       sprintf("%s -t %s.ins.bed %s -q 0.01 -n %s_cuts 2> %s_cuts.macs.log", macs, r, cuts, r, r))
  step(sprintf("macs/%s_cuts_p01_peaks.narrowPeak", r),
       sprintf("%s -t %s.ins.bed %s -p 0.01 -n %s_cuts_p01 2> %s_cuts_p01.macs.log", macs, r, cuts, r, r))
  step(sprintf("macs/%s_frag_p01_peaks.narrowPeak", r),
       sprintf("%s -t %s.dedup.bam -f BAMPE -p 0.01 -n %s_frag_p01 2> %s_frag_p01.macs.log", macs, r, r, r))
}
# threshold series for fig13 (replicate 1, both modes)
cutoffs <- c("-p 0.01", "-q 0.1", "-q 0.05", "-q 0.01", "-q 1e-3", "-q 1e-5", "-q 1e-10")
for (k in seq_along(cutoffs)) {
  step(sprintf("macs/thr%d_frag_peaks.narrowPeak", k),
       sprintf("macs3 callpeak -g 4.0e7 --keep-dup all --outdir macs -t rep1.dedup.bam -f BAMPE %s -n thr%d_frag 2> /dev/null", cutoffs[k], k))
  step(sprintf("macs/thr%d_cuts_peaks.narrowPeak", k),
       sprintf("macs3 callpeak -g 4.0e7 --keep-dup all --outdir macs -t rep1.ins.bed -f BED --nomodel --shift -75 --extsize 150 %s -n thr%d_cuts 2> /dev/null", cutoffs[k], k))
}
# pitfall: peaks called without removing chrM (and without MAPQ / dedup filters)
step("macs/rep1_nofilter_peaks.narrowPeak",
     "macs3 callpeak -t rep1.raw.bam -f BAMPE -g 4.0e7 -q 0.01 -n rep1_nofilter --outdir macs 2> rep1_nofilter.macs.log")

# ---- read counts per step (replicate 1 and 2) --------------------------------
cnt <- function(cmd) as.numeric(system2("bash", c("-c", shQuote(paste0("cd ", W, "; ", cmd))), stdout = TRUE))
funnel <- do.call(rbind, lapply(reps, function(r) {
  js <- readLines(file.path(W, paste0(r, ".fastp.json")))
  tot <- as.numeric(sub(".*: *([0-9]+).*", "\\1", grep("\"total_reads\"", js, value = TRUE)[1])) / 2
  aft <- as.numeric(sub(".*: *([0-9]+).*", "\\1", grep("\"total_reads\"", js, value = TRUE)[2])) / 2
  st <- system2("samtools", c("idxstats", file.path(W, paste0(r, ".raw.bam"))), stdout = TRUE)
  st <- read.table(text = st)
  mapped_pairs <- cnt(sprintf("samtools view -c -f 66 -F 2308 %s.raw.bam", r))     # read 1 of proper pairs
  chrM <- cnt(sprintf("samtools view -c -f 66 -F 2308 %s.raw.bam chrM", r))
  chr21_any <- cnt(sprintf("samtools view -c -f 66 -F 2308 %s.raw.bam chr21", r))
  q30 <- cnt(sprintf("samtools view -u -q 30 -f 2 -F 1804 %s.raw.bam chr21 | samtools collate -u -O - | samtools fixmate -u - - | samtools view -c -f 66 -", r))
  final <- cnt(sprintf("samtools view -c -f 64 %s.dedup.bam", r))
  ins <- cnt(sprintf("wc -l < %s.ins.bed", r))
  data.frame(rep = r,
             step = c("read pairs in", "after trimming", "properly paired", "on chr21 (not chrM)",
                      "MAPQ ≥ 30, both mates", "not duplicates"),
             n = c(tot, aft, mapped_pairs, chr21_any, q30, final), chrM = chrM, ins = ins)
}))
funnel$step <- factor(funnel$step, levels = rev(unique(funnel$step)))
print(funnel)

# ---- fig04: reads per step ----------------------------------------------------
p04 <- ggplot(funnel, aes(n / 1e6, step, fill = factor(rep, levels = rev(reps)))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f M", n / 1e6)), position = position_dodge(width = 0.8),
            hjust = -0.1, size = 4.3) +
  scale_fill_manual(values = col_rep, name = NULL, breaks = reps) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "read pairs (millions)", y = NULL) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
ggsave("fig04-reads-per-step.png", p04, width = 11, height = 4.6, dpi = 150, bg = "white")

# ---- fig11: Poisson test of k = 25 against lambda = 3 and lambda = 8 ----------
k_obs <- 25
pois <- do.call(rbind, lapply(c(3, 8), function(l) data.frame(
  lambda = sprintf("λ = %g", l), k = 0:32, p = dpois(0:32, l), tail = 0:32 >= k_obs)))
ptxt <- sapply(c(3, 8), function(l) ppois(k_obs - 1, l, lower.tail = FALSE))
lab <- data.frame(lambda = sprintf("λ = %g", c(3, 8)),
                  txt = sprintf("P(X ≥ %d) = %.1e", k_obs, ptxt))
p11 <- ggplot(pois, aes(k, p, fill = tail)) +
  geom_col(width = 0.85) +
  geom_vline(xintercept = k_obs - 0.5, colour = col_ins, linetype = 2) +
  geom_text(data = lab, aes(x = k_obs - 1.5, y = 0.1, label = txt), hjust = 1, vjust = -0.5,
            colour = col_ins, size = 5, inherit.aes = FALSE) +
  facet_wrap(~lambda, ncol = 2) +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = col_ins), guide = "none") +
  scale_y_sqrt() +
  labs(x = "count k in the window", y = "Poisson probability (sqrt scale)")
ggsave("fig11-poisson-test.png", p11, width = 11, height = 4, dpi = 150, bg = "white")
cat(sprintf("fig11: P(X>=25 | 3) = %.2g, P(X>=25 | 8) = %.2g\n", ptxt[1], ptxt[2]))

# ---- peaks ---------------------------------------------------------------------
np <- function(f) import(file.path(W, f), format = "narrowPeak")
bl <- import(file.path(D, "hg38_blacklist.bed.gz"))
enc <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
enc <- enc[seqnames(enc) == "chr21"]
pk <- list()
for (r in reps) for (m in c("frag", "cuts", "cuts_p01", "frag_p01"))
  pk[[paste(r, m, sep = "_")]] <- np(sprintf("macs/%s_%s_peaks.narrowPeak", r, m))
nofilt <- np("macs/rep1_nofilter_peaks.narrowPeak")
cat(sprintf("pitfall: unfiltered rep1 call: %d peaks, %d on chrM, %d in blacklist; top 20 by p on chrM: %d\n",
            length(nofilt), sum(seqnames(nofilt) == "chrM"), sum(overlapsAny(nofilt, bl)),
            sum(seqnames(nofilt[order(-nofilt$pValue)][1:20]) == "chrM")))
cat("peaks per call:\n"); print(sapply(pk, length))
cat("median width:\n"); print(sapply(pk, function(p) median(width(p))))
cat("in blacklist:\n"); print(sapply(pk, function(p) sum(overlapsAny(p, bl))))
cat(sprintf("ENCODE IDR peaks on chr21: %d (median width %d)\n", length(enc), median(width(enc))))
for (nm in names(pk)) cat(sprintf("%-14s %6d peaks, %5.1f%% overlap ENCODE, ENCODE recall %5.1f%%\n", nm,
  length(pk[[nm]]), 100 * mean(overlapsAny(pk[[nm]], enc)), 100 * mean(overlapsAny(enc, pk[[nm]]))))
cat(sprintf("rep1 cuts overlapping rep2 cuts: %.1f%%; rep2 in rep1: %.1f%%\n",
  100 * mean(overlapsAny(pk$rep1_cuts, pk$rep2_cuts)), 100 * mean(overlapsAny(pk$rep2_cuts, pk$rep1_cuts))))

# insertions and FRiP
ins <- lapply(setNames(reps, reps), function(r) import(file.path(W, paste0(r, ".ins.bed")), format = "bed"))
frip <- sapply(reps, function(r) c(own_cuts = mean(overlapsAny(ins[[r]], pk[[paste0(r, "_cuts")]])),
                                   own_frag = mean(overlapsAny(ins[[r]], pk[[paste0(r, "_frag")]])),
                                   encode = mean(overlapsAny(ins[[r]], enc))))
cat("FRiP:\n"); print(round(frip, 3))

# ---- helpers for browser-style panels -------------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
gtf <- gtf[seqnames(gtf) == "chr21" & gtf$gene_type == "protein_coding"]
frag <- lapply(setNames(reps, reps), function(r) import(file.path(W, paste0(r, ".frag.bed")), format = "bed"))
gene_panel <- function(region, xl, max_genes = 6) {
  tx <- subsetByOverlaps(gtf[gtf$type == "transcript"], region)
  tx <- tx[order(-width(tx))]; tx <- tx[!duplicated(tx$gene_name)]
  tx <- tx[order(start(tx))]
  if (length(tx) > max_genes) tx <- tx[seq_len(max_genes)]
  tx$row <- (seq_along(tx) - 1) %% 2 + 1
  ex <- gtf[gtf$type == "exon" & gtf$transcript_id %in% tx$transcript_id]
  exdf <- data.frame(s = pmax(start(ex), xl[1]), e = pmin(end(ex), xl[2]),
                     row = tx$row[match(ex$transcript_id, tx$transcript_id)])
  exdf <- exdf[exdf$e > exdf$s, ]
  txdf <- data.frame(s = pmax(start(tx), xl[1]), e = pmin(end(tx), xl[2]), row = tx$row,
                     g = paste0(tx$gene_name, ifelse(strand(tx) == "+", " →", " ←")))
  ggplot() +
    geom_segment(data = txdf, aes(x = s, xend = e, y = -row, yend = -row), colour = "#1b2a3a") +
    geom_rect(data = exdf, aes(xmin = s, xmax = e, ymin = -row - 0.22, ymax = -row + 0.22), fill = "#1b2a3a") +
    geom_text(data = txdf, aes(x = pmin(pmax((s + e) / 2, xl[1] + 0.05 * diff(xl)), xl[2] - 0.05 * diff(xl)),
                               y = -row - 0.5, label = g), size = 4.2, vjust = 1) +
    scale_x_continuous(limits = xl, expand = c(0, 0), labels = function(b) sprintf("%.3f Mb", b / 1e6)) +
    scale_y_continuous(breaks = NULL, limits = c(-3.1, -0.5)) +
    labs(x = "chr21", y = NULL) +
    theme_minimal(base_size = 13) + theme(panel.grid = element_blank(), plot.margin = margin(5, 30, 5, 5))
}
# per-bin CPM of insertions and of fragment coverage in a region
bin_signal <- function(gr, region, bin, total, what = c("ins", "frag")) {
  what <- match.arg(what)
  g <- subsetByOverlaps(gr, region)
  br <- seq(start(region), end(region), by = bin)
  if (what == "ins") {
    v <- tabulate(findInterval(start(g), br), nbins = length(br))
  } else {
    cv <- coverage(g)[["chr21"]][start(region):end(region)]
    v <- sapply(split(as.numeric(cv), findInterval(start(region):end(region), br)), mean)
  }
  data.frame(x = br + bin / 2, y = as.numeric(v) / total * 1e6)
}
peak_track <- function(sets, xl) {
  df <- do.call(rbind, lapply(names(sets), function(n) {
    g <- sets[[n]]; if (length(g) == 0) return(NULL)
    data.frame(set = n, s = pmax(start(g), xl[1]), e = pmin(end(g), xl[2]))
  }))
  df$set <- factor(df$set, levels = rev(names(sets)))
  ggplot(df) + geom_rect(aes(xmin = s, xmax = e, ymin = as.numeric(set) - 0.35, ymax = as.numeric(set) + 0.35,
                             fill = set)) +
    scale_y_continuous(breaks = seq_along(levels(df$set)), labels = levels(df$set), expand = c(0.1, 0.1)) +
    scale_x_continuous(limits = xl, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 12),
          legend.position = "none", plot.margin = margin(5, 30, 5, 5))
}

# ---- fig08: insertion sites vs fragment coverage at one locus (replicate 1) ------
loc8 <- GRanges("chr21", IRanges(31657501, 31662500))
n_ins <- length(ins$rep1); n_frag <- length(frag$rep1)
d8 <- rbind(cbind(track = "insertion sites\n(CPM per 10 bp)", bin_signal(ins$rep1, loc8, 10, n_ins, "ins")),
            cbind(track = "fragment coverage\n(CPM)", bin_signal(frag$rep1, loc8, 10, n_frag, "frag")))
d8$track <- factor(d8$track, levels = unique(d8$track))
xl8 <- c(start(loc8), end(loc8))
p8 <- ggplot(d8, aes(x, y, fill = track)) + geom_col(width = 10, colour = NA) +
  facet_grid(track ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(col_ins, col_frag), guide = "none") +
  scale_x_continuous(limits = xl8, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
        strip.placement = "outside", axis.text.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), plot.margin = margin(5, 30, 5, 5))
p8 <- p8 / (peak_track(list(`peaks (insertions)` = subsetByOverlaps(pk$rep1_cuts, loc8),
                           `peaks (fragments)` = subsetByOverlaps(pk$rep1_frag, loc8)), xl8) +
         scale_fill_manual(values = c(`peaks (insertions)` = col_ins, `peaks (fragments)` = col_frag))) /
  gene_panel(loc8, xl8) + plot_layout(heights = c(5, 0.9, 1.3))
ggsave("fig08-insertions-vs-fragments-real.png", p8, width = 13, height = 4.8, dpi = 150, bg = "white")

# ---- fig13: number of peaks vs threshold ------------------------------------------
thr <- do.call(rbind, lapply(seq_along(cutoffs), function(k) do.call(rbind, lapply(c("frag", "cuts"), function(m) {
  p <- np(sprintf("macs/thr%d_%s_peaks.narrowPeak", k, m))
  data.frame(cutoff = cutoffs[k], mode = m, n = length(p), width = median(width(p)),
             encode = mean(overlapsAny(p, enc)))
}))))
thr$cutoff <- factor(thr$cutoff, levels = cutoffs)
thr$mode <- factor(ifelse(thr$mode == "frag", "fragments (BAMPE)", "insertion sites"),
                   levels = c("fragments (BAMPE)", "insertion sites"))
print(thr)
p13a <- ggplot(thr, aes(cutoff, n, colour = mode, group = mode)) + geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
  geom_hline(yintercept = length(enc), linetype = 2, colour = col_enc) +
  annotate("text", x = 7.4, y = length(enc), label = sprintf("ENCODE IDR peaks: %d", length(enc)),
           hjust = 1, vjust = 1.6, colour = col_enc, size = 4.3) +
  scale_colour_manual(values = c(col_frag, col_ins), name = NULL) +
  scale_y_continuous(limits = c(0, 1900)) +
  labs(x = "MACS3 cutoff (stricter →)", y = "peaks on chr21", title = "Number of peaks") +
  theme(legend.position = "top")
p13b <- ggplot(thr, aes(cutoff, 100 * encode, colour = mode, group = mode)) + geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(col_frag, col_ins), guide = "none") +
  scale_y_continuous(limits = c(NA, 100)) +
  labs(x = "MACS3 cutoff (stricter →)", y = "% overlapping an ENCODE peak", title = "Agreement with ENCODE")
ggsave("fig13-peaks-vs-threshold.png", p13a + p13b + plot_layout(widths = c(1.2, 1)),
       width = 12.5, height = 4.6, dpi = 150, bg = "white")

# ---- fig16: replicate rank consistency (lenient calls, insertion mode) ------------
p1 <- pk$rep1_cuts_p01[order(-pk$rep1_cuts_p01$pValue)]
p2 <- pk$rep2_cuts_p01[order(-pk$rep2_cuts_p01$pValue)]
h <- findOverlaps(p1, p2)
best <- tapply(subjectHits(h), queryHits(h), min)
m <- data.frame(rank1 = seq_along(p1), s1 = p1$pValue, rank2 = NA_integer_, s2 = NA_real_)
m$rank2[as.integer(names(best))] <- best
m$s2[as.integer(names(best))] <- p2$pValue[best]
m$found <- !is.na(m$rank2)
cat_curve <- data.frame(N = seq_len(nrow(m)), frac = cumsum(m$found) / seq_len(nrow(m)))
Ns <- c(250, 500, 750, 1000, 1250, 1500, nrow(m))
cat(sprintf("lenient peaks: rep1 %d, rep2 %d; spearman of matched -log10 p %.2f\n", length(p1), length(p2),
            cor(m$s1, m$s2, use = "complete.obs", method = "spearman")))
print(data.frame(N = Ns, found = round(cat_curve$frac[Ns], 3)))
n90 <- max(which(cat_curve$frac >= 0.9 & cat_curve$N >= 100))
cat(sprintf("top-N curve stays >= 90%% up to N = %d\n", n90))
mm <- m; mm$s2[is.na(mm$s2)] <- 1.5          # unmatched: drawn on the bottom margin
mm$s2[!mm$found] <- jitter(rep(1.5, sum(!mm$found)), amount = 0.08)
p16a <- ggplot(mm, aes(s1, s2, colour = found)) + geom_point(size = 0.7, alpha = 0.35) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = col_ok, `FALSE` = "grey55"), guide = "none") +
  annotate("text", x = max(mm$s1), y = 1.5, label = "not found in rep 2", hjust = 1, vjust = -1, size = 4, colour = "grey40") +
  labs(x = "−log10 p, replicate 1", y = "−log10 p, best match in replicate 2",
       title = "Matched peak scores")
q01 <- length(pk$rep1_cuts)
p16b <- ggplot(cat_curve, aes(N, 100 * frac)) + geom_line(linewidth = 1.1, colour = col_ok) +
  geom_vline(xintercept = q01, linetype = 2, colour = col_ins) +
  annotate("text", x = q01, y = 52, label = sprintf("q < 0.01: %d peaks ", q01), hjust = 1, colour = col_ins, size = 4.3) +
  geom_vline(xintercept = length(enc), linetype = 3, colour = col_enc) +
  annotate("text", x = length(enc), y = 20, label = sprintf("ENCODE: %d ", length(enc)), hjust = 1, colour = col_enc, size = 4.3) +
  scale_y_continuous(limits = c(0, 100)) + scale_x_continuous(limits = c(0, max(nrow(m), length(enc)) * 1.02)) +
  labs(x = "top N peaks of replicate 1 (by p-value)", y = "% found in replicate 2",
       title = "Agreement at the top")
ggsave("fig16-rank-consistency.png", p16a + p16b, width = 12.5, height = 4.8, dpi = 150, bg = "white")

# ---- fig17: replicates, our peaks and ENCODE's at a locus ------------------------
loc17 <- GRanges("chr21", IRanges(44240001, 44300000))
xl17 <- c(start(loc17), end(loc17))
d17 <- rbind(cbind(track = "rep1", bin_signal(ins$rep1, loc17, 50, length(ins$rep1), "ins")),
             cbind(track = "rep2", bin_signal(ins$rep2, loc17, 50, length(ins$rep2), "ins")))
p17 <- ggplot(d17, aes(x, y, fill = track)) + geom_col(width = 50, colour = NA) +
  facet_grid(track ~ ., switch = "y") +
  scale_fill_manual(values = col_rep, guide = "none") +
  scale_x_continuous(limits = xl17, expand = c(0, 0)) +
  labs(x = NULL, y = "insertions, CPM per 50 bp") +
  theme_minimal(base_size = 13) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
        strip.placement = "outside", axis.text.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), plot.margin = margin(5, 30, 5, 5))
sets17 <- list(`rep 1 peaks` = subsetByOverlaps(pk$rep1_cuts, loc17),
               `rep 2 peaks` = subsetByOverlaps(pk$rep2_cuts, loc17),
               `ENCODE IDR` = subsetByOverlaps(enc, loc17))
p17 <- p17 / (peak_track(sets17, xl17) + scale_fill_manual(values = c(`rep 1 peaks` = col_rep[[1]],
              `rep 2 peaks` = col_rep[[2]], `ENCODE IDR` = col_enc))) /
  gene_panel(loc17, xl17) + plot_layout(heights = c(4.5, 1.3, 1.4))
ggsave("fig17-peaks-locus.png", p17, width = 13, height = 5, dpi = 150, bg = "white")

lowq <- cnt("samtools view -c -f 66 rep1.raw.bam chr21:44240000-44254000")
hiq <- cnt("samtools view -c -f 66 -q 30 rep1.raw.bam chr21:44240000-44254000")
cat(sprintf("fig17 left part chr21:44.240-44.254 Mb: %d read pairs, %d with MAPQ >= 30\n", lowq, hiq))

# ---- consensus (as in exercise R2) and fig21: promoter fraction -------------------
fixed <- function(p, sample) {
  s <- GRanges(seqnames(p), IRanges(start(p) + p$peak, width = 1))
  w <- resize(s, 500, fix = "center")
  w$score <- p$pValue / sum(p$pValue) * 1e6
  w$sample <- sample
  w
}
keep_best <- function(gr) {
  gr <- gr[order(-gr$score)]
  ov <- as.list(findOverlaps(gr, gr))
  keep <- logical(length(gr)); gone <- logical(length(gr))
  for (i in seq_along(gr)) {
    if (gone[i]) next
    keep[i] <- TRUE
    gone[ov[[i]]] <- TRUE
  }
  sort(gr[keep])
}
per_sample <- c(keep_best(fixed(pk$rep1_cuts, "rep1")), keep_best(fixed(pk$rep2_cuts, "rep2")))
cons <- keep_best(per_sample)
n_rep <- sapply(split(per_sample, per_sample$sample), function(x) overlapsAny(cons, x))
cons_rep <- cons[rowSums(n_rep) == 2]
cat(sprintf("fixed-width: per sample %s; consensus %d; reproducible %d; overlapping pairs in cons %d\n",
            paste(table(per_sample$sample), collapse = "/"), length(cons), length(cons_rep),
            sum(countOverlaps(cons, cons) > 1)))
merged <- reduce(c(granges(pk$rep1_cuts), granges(pk$rep2_cuts)))
cat(sprintf("merged: %d regions, median width %d, max %d\n", length(merged), median(width(merged)), max(width(merged))))

tx <- gtf[gtf$type == "transcript" & gtf$transcript_type == "protein_coding"]
tss <- unique(granges(resize(tx, width = 1, fix = "start")))
near <- function(g) mcols(distanceToNearest(g, tss, ignore.strand = TRUE))$distance <= 1000
dcons <- near(cons_rep)
cons_rep$decile <- cut(rank(-cons_rep$score, ties.method = "first"), 10, labels = 1:10)
dec <- data.frame(decile = 1:10, frac = tapply(dcons, cons_rep$decile, mean))
sets21 <- c(`rep 1 peaks\n(q < 0.01)` = mean(near(pk$rep1_cuts)),
            `reproducible\nconsensus` = mean(dcons), `ENCODE IDR\n(chr21)` = mean(near(enc)))
print(round(sets21, 3)); print(round(dec, 3))
cat(sprintf("TSSs with a consensus peak within 1 kb: %.1f%% of %d\n",
            100 * mean(overlapsAny(tss, cons_rep, maxgap = 1000)), length(tss)))
p21a <- ggplot(data.frame(set = factor(names(sets21), levels = names(sets21)), f = 100 * sets21),
               aes(set, f)) + geom_col(fill = c(col_ins, col_ok, col_enc), width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", f)), vjust = -0.4, size = 4.6) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "% within 1 kb of a TSS", title = "Peaks at promoters") +
  theme(axis.text.x = element_text(size = 12))
p21b <- ggplot(dec, aes(factor(decile), 100 * frac)) + geom_col(fill = col_ok, width = 0.7) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "decile of peak strength (1 = strongest)", y = "% within 1 kb of a TSS",
       title = "Consensus peaks by strength")
ggsave("fig21-promoter-fraction.png", p21a + p21b + plot_layout(widths = c(1, 1.3)),
       width = 12, height = 4.4, dpi = 150, bg = "white")
