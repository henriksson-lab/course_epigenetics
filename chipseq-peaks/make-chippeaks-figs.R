# Real-data figures for chipseq-peaks/lec.qmd.
# Run from this folder:  Rscript make-chippeaks-figs.R
#
# Needs: the CTCF and H3K36me3 read sets (data/prepare-reads.sh, data/prepare-chip.sh),
#        data/hg38_chr21_chrM.fa with its bowtie2 index, blacklist, GENCODE GTF, ENCODE CTCF peaks
#        and H3K4me3 fold change (data/download.sh);
#        on the PATH: fastp, bowtie2, samtools (>= 1.10), bedtools, macs3, idr, epic2;
#        R packages GenomicRanges, GenomicAlignments, rtracklayer, Rsamtools, ggplot2, patchwork,
#        JASPAR2024, TFBSTools, motifmatchr.
#
# Processing (the commands of exercise.qmd) runs once into ./peaks-work (set PEAKS_WORK), ~5 min on 8 cores.
# Writes fig04-lambda-real.png, fig06-threshold-motif.png, fig07-no-input.png, fig10-broad-locus.png,
# fig11-broad-genes.png, fig13-idr-real.png and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicAlignments); library(rtracklayer); library(Rsamtools)
  library(ggplot2); library(patchwork)
})
set.seed(1)
D <- normalizePath("../data")
W <- Sys.getenv("PEAKS_WORK", "peaks-work")
dir.create(W, showWarnings = FALSE, recursive = TRUE); W <- normalizePath(W)
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15)))
G21 <- 40088619

# ---- 1. processing --------------------------------------------------------------
chain <- paste(
  "chain() { local s=$1;",
  "if [ \"$3\" = 64 ]; then zcat $2 | perl -pe 'tr/@-~/!-_/ if $. % 4 == 0' | gzip > $s.fastq.gz; else cp $2 $s.fastq.gz; fi;",
  "fastp -w $T -i $s.fastq.gz -o $s.trim.fastq.gz -j $s.fastp.json -h $s.fastp.html 2> $s.fastp.log;",
  "bowtie2 -p $T -x $D/bowtie2/hg38_chr21_chrM -U $s.trim.fastq.gz 2> $s.bowtie2.log | samtools sort -@ 4 -o $s.sorted.bam -;",
  "samtools index $s.sorted.bam; samtools view -b -q 30 -F 1804 -o $s.q30.bam $s.sorted.bam;",
  "samtools markdup -s $s.q30.bam $s.markdup.bam 2> $s.markdup.txt; samtools index $s.markdup.bam;",
  "samtools view -b -F 1024 -o $s.dedup.bam $s.markdup.bam chr21;",
  "bedtools intersect -v -abam $s.dedup.bam -b blacklist.bed > $s.final.bam; samtools index $s.final.bam; };")
sh <- function(cmd) {
  rc <- system2("bash", c("-c", shQuote(paste0("set -eo pipefail; cd ", W, "; D=", D, "; T=", THREADS, "; ", chain, cmd))))
  if (rc != 0) stop("failed: ", cmd)
}
step <- function(out, cmd) if (!all(file.exists(file.path(W, out)))) { cat("running:", out, "\n"); sh(cmd) }
step("blacklist.bed", "zcat $D/hg38_blacklist.bed.gz > blacklist.bed")
smp <- data.frame(name = c("ctcf", "ctcf_rep2", "ctcf_input", "h3k36me3", "h3k36me3_rep2", "h3k36me3_input"),
                  fq = c("gm12878_ctcf_chip", "gm12878_ctcf_chip_rep2", "gm12878_ctcf_input", "gm12878_h3k36me3_chip",
                         "gm12878_h3k36me3_chip_rep2", "gm12878_h3k36me3_input"),
                  phred = c(64, 64, 64, 33, 33, 33))
for (i in seq_len(nrow(smp)))
  step(paste0(smp$name[i], ".final.bam.bai"), sprintf("chain %s $D/%s.fastq.gz %d", smp$name[i], smp$fq[i], smp$phred[i]))
step("encode_ctcf.bed", "zcat $D/gm12878_ctcf_peaks.bed.gz | awk '$1 == \"chr21\"' | sort -k1,1 -k2,2n > encode_ctcf.bed")
step("macs/ctcf_peaks.narrowPeak", "mkdir -p macs; macs3 callpeak -t ctcf.final.bam -c ctcf_input.final.bam -f BAM -g 40088619 -q 0.01 --bdg -n ctcf --outdir macs 2> macs/ctcf.log")
step("macs/ctcf_noinput_peaks.narrowPeak", "macs3 callpeak -t ctcf.final.bam -f BAM -g 40088619 -q 0.01 -n ctcf_noinput --outdir macs 2> macs/ctcf_noinput.log")
qs <- c("0.1", "0.01", "1e-5", "1e-10", "1e-20", "1e-40")
for (q in qs) step(sprintf("macs/ctcf_q%s_peaks.narrowPeak", q),
                   sprintf("macs3 callpeak -t ctcf.final.bam -c ctcf_input.final.bam -f BAM -g 40088619 -q %s -n ctcf_q%s --outdir macs 2> /dev/null", q, q))
K <- "-t h3k36me3.final.bam h3k36me3_rep2.final.bam -c h3k36me3_input.final.bam -f BAM -g 40088619 --nomodel --extsize 200"
step("macs/k36_narrow_peaks.narrowPeak", paste("macs3 callpeak", K, "-q 0.01 -n k36_narrow --outdir macs 2> macs/k36_narrow.log"))
step("macs/k36_broad_peaks.broadPeak", paste("macs3 callpeak", K, "--broad --broad-cutoff 0.1 -q 0.01 -n k36_broad --outdir macs 2> macs/k36_broad.log"))
step("k36_epic2.bed", paste("grep -w chr21 $D/hg38.chrom.sizes > chr21.sizes;",
  "epic2 -t h3k36me3.final.bam h3k36me3_rep2.final.bam -c h3k36me3_input.final.bam --chromsizes chr21.sizes",
  "--effective-genome-fraction 0.858 --fragment-size 200 --bin-size 200 --gaps-allowed 3 --false-discovery-rate-cutoff 0.05 > k36_epic2.bed 2> k36_epic2.log"))
step("ctcf_idr.txt", paste(
  "samtools merge -f ctcf_pool.final.bam ctcf.final.bam ctcf_rep2.final.bam;",
  "for s in ctcf ctcf_rep2 ctcf_pool; do macs3 callpeak -t $s.final.bam -c ctcf_input.final.bam -f BAM -g 40088619 -p 0.01 -n ${s}_p01 --outdir macs 2> macs/${s}_p01.log;",
  "sort -k8,8gr macs/${s}_p01_peaks.narrowPeak > ${s}_p01.narrowPeak; done;",
  "idr --samples ctcf_p01.narrowPeak ctcf_rep2_p01.narrowPeak --peak-list ctcf_pool_p01.narrowPeak --input-file-type narrowPeak",
  "--rank p.value --soft-idr-threshold 0.05 --use-best-multisummit-IDR --output-file ctcf_idr.txt --log-output-file ctcf_idr.log"))
step("Np.txt", paste(
  "idr_pair() { for x in $1 $2; do macs3 callpeak -t $x -c ctcf_input.final.bam -f BAM -g 40088619 -p 0.01 -n tmp --outdir macs 2> /dev/null;",
  "sort -k8,8gr macs/tmp_peaks.narrowPeak > $x.p01.narrowPeak; done;",
  "idr --samples $1.p01.narrowPeak $2.p01.narrowPeak --peak-list $3 --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.05 --use-best-multisummit-IDR --output-file $4.txt 2> /dev/null; };",
  "for s in ctcf ctcf_rep2 ctcf_pool; do samtools view -b -s 7.5 -U $s.pr2.bam -o $s.pr1.bam $s.final.bam; done;",
  "idr_pair ctcf.pr1.bam ctcf.pr2.bam ctcf_p01.narrowPeak N1; idr_pair ctcf_rep2.pr1.bam ctcf_rep2.pr2.bam ctcf_rep2_p01.narrowPeak N2;",
  "idr_pair ctcf_pool.pr1.bam ctcf_pool.pr2.bam ctcf_pool_p01.narrowPeak Np"))

f <- function(x) file.path(W, x)
np <- function(x) import(f(x), format = "narrowPeak")
rd <- function(s) granges(readGAlignments(f(paste0(s, ".final.bam"))))
reads <- lapply(setNames(smp$name, smp$name), rd)
print(sapply(reads, length))
idr_n <- function(x) sum(read.delim(f(x), header = FALSE)$V12 >= -log10(0.05))
Nt <- idr_n("ctcf_idr.txt"); N1 <- idr_n("N1.txt"); N2 <- idr_n("N2.txt"); Np <- idr_n("Np.txt")
cat("IDR < 0.05: Nt", Nt, " Np", Np, " N1", N1, " N2", N2, " rescue", round(max(Nt, Np) / min(Nt, Np), 2),
    " self-consistency", round(max(N1, N2) / min(N1, N2), 2), "\n")

# ---- annotation ---------------------------------------------------------------------
chr21 <- GRanges("chr21", IRanges(1, 46709983))
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), which = chr21)
genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]; seqlevels(genes, pruning.mode = "coarse") <- "chr21"
k4 <- import(file.path(D, "gm12878_h3k4me3_fc.bigWig"), which = chr21, as = "RleList")
genes$k4 <- mean(k4[promoters(resize(genes, 1, fix = "start"), 1000, 1001)])
genes$state <- ifelse(genes$k4 > 5, "active", ifelse(genes$k4 < 1, "silent", "other"))
tss <- resize(genes, 1, fix = "start")
enc <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak"); enc <- enc[seqnames(enc) == "chr21"]

# motif: CTCF (JASPAR MA0139.2) within 50 bp of the summit
suppressPackageStartupMessages({ library(JASPAR2024); library(TFBSTools); library(motifmatchr) })
con <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024()))
pfm <- getMatrixByID(con, "MA0139.2")
fa <- FaFile(file.path(D, "hg38_chr21_chrM.fa"))
has_motif <- function(p) {
  s <- GRanges(seqnames(p), IRanges(start(p) + p$peak - 50, width = 101))
  as.logical(motifMatches(matchMotifs(pfm, getSeq(fa, s), p.cutoff = 1e-4))[, 1])
}
rnd <- GRanges("chr21", IRanges(sample(15e6:46e6, 5000), width = 1)); rnd$peak <- 0L
bg_motif <- mean(has_motif(rnd))
cat("CTCF motif in random 101-bp windows:", round(bg_motif, 3), "\n")

# ---- fig04: the local lambda at real summits ------------------------------------------
pk <- np("macs/ctcf_peaks.narrowPeak")
summit <- GRanges("chr21", IRanges(start(pk) + pk$peak, width = 1))
trt <- import(f("macs/ctcf_treat_pileup.bdg"), format = "bedGraph")
lam <- import(f("macs/ctcf_control_lambda.bdg"), format = "bedGraph")
summit$k <- trt$score[findOverlaps(summit, trt, select = "first")]
summit$lambda <- lam$score[findOverlaps(summit, lam, select = "first")]
five_in <- resize(reads$ctcf_input, 1, fix = "start")
d <- 100; lamBG <- length(reads$ctcf_input) * d / G21
lam_w <- function(w) countOverlaps(resize(summit, w, fix = "center"), five_in, ignore.strand = TRUE) * d / w
L <- cbind(BG = lamBG, d = lam_w(d), `1 kb` = lam_w(1000), `10 kb` = lam_w(10000))
win <- colnames(L)[apply(L, 1, which.max)]
cat("lambda_BG =", round(lamBG, 3), "; lambda set by:", paste(names(table(win)), table(win)), "; cor with MACS:",
    round(cor(apply(L, 1, max), summit$lambda), 4), "\n")
summit$p_ours <- -ppois(ceiling(summit$k) - 1, summit$lambda, lower.tail = FALSE, log.p = TRUE) / log(10)
cat("p-value check, max |ours - MACS| (-log10):", signif(max(abs(summit$p_ours - pk$pValue)), 3), "\n")
cat("treatment scaling factor:", round(length(reads$ctcf_input) / length(reads$ctcf), 4), "\n")
ex <- order(-summit$lambda)[1]
cat("example summit:", start(summit)[ex], " k =", round(summit$k[ex], 1), " lambdas:", round(L[ex, ], 2),
    " -log10 p:", round(pk$pValue[ex], 1), "\n")
reg <- GRanges("chr21", IRanges(start(summit)[ex] - 3000, start(summit)[ex] + 3000))
tr <- subsetByOverlaps(trt, reg); lm <- subsetByOverlaps(lam, reg)
step_df <- function(g, what) data.frame(what = what, x = c(rbind(start(g), end(g))) / 1e3, y = rep(g$score, each = 2))
d4 <- rbind(step_df(tr, "ChIP pileup (scaled to input depth)"), step_df(lm, "local λ from the input"))
p4a <- ggplot(d4, aes(x, y, colour = what)) + geom_line(linewidth = 1.1) +
  geom_hline(yintercept = lamBG, linetype = 2, colour = "grey40") +
  annotate("text", x = max(d4$x), y = lamBG, label = "lambda[BG]", parse = TRUE, vjust = 1.4, hjust = 1, size = 4.5) +
  scale_colour_manual(values = c("#8c564b", "#2a6f97"), name = NULL) + scale_y_sqrt() +
  labs(x = "chr21 position (kb)", y = "fragments per base (sqrt scale)") + theme(legend.position = "top")
d4b <- as.data.frame(table(factor(win, colnames(L))))
p4b <- ggplot(d4b, aes(Var1, Freq)) + geom_col(fill = "#2a6f97") + geom_text(aes(label = Freq), vjust = -0.3, size = 4.5) +
  labs(x = "window that sets λ_local", y = sprintf("summits (of %d)", length(summit)), title = "Which window wins?") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
ggsave("fig04-lambda-real.png", p4a + p4b + plot_layout(widths = c(2.2, 1)), width = 12.5, height = 4.8, dpi = 150, bg = "white")

# ---- fig06: thresholds and motif by rank ------------------------------------------------
d6 <- do.call(rbind, lapply(qs, function(q) {
  p <- np(sprintf("macs/ctcf_q%s_peaks.narrowPeak", q))
  data.frame(q = q, n = length(p), enc = sum(overlapsAny(p, enc)), motif = mean(has_motif(p)))
}))
print(d6)
p6a <- ggplot(d6, aes(-log10(as.numeric(q)), n)) + geom_line(linewidth = 1.1, colour = "#8c564b") + geom_point(size = 2.5, colour = "#8c564b") +
  geom_line(aes(y = enc), linetype = 2, colour = "#5b6b7a") + geom_point(aes(y = enc), colour = "#5b6b7a", size = 2) +
  annotate("text", x = 22, y = 420, label = "all peaks", colour = "#8c564b", size = 4.5, hjust = 0) +
  annotate("text", x = 22, y = 330, label = "also in ENCODE's set", colour = "#5b6b7a", size = 4.5, hjust = 0) +
  labs(x = "q-value cutoff (-log10)", y = "CTCF peaks on chr21", title = "Fewer peaks at stricter cutoffs")
lax <- np("macs/ctcf_q0.1_peaks.narrowPeak"); lax <- lax[order(-lax$pValue)]
lax$motif <- has_motif(lax)
lax$bin <- cut(seq_along(lax), breaks = seq(0, length(lax) + 50, by = 50), labels = FALSE)
d6b <- aggregate(motif ~ bin, as.data.frame(mcols(lax)), mean)
d6b$rank <- (d6b$bin - 0.5) * 50
p6b <- ggplot(d6b, aes(rank, motif)) + geom_col(fill = "#8c564b", width = 45) +
  geom_hline(yintercept = bg_motif, linetype = 2) +
  annotate("label", x = max(d6b$rank), y = bg_motif + 0.07, label = "random DNA", hjust = 1, size = 4, label.size = 0, fill = "white") +
  geom_vline(xintercept = d6$n[d6$q == "0.01"], colour = "#c0392b") +
  annotate("text", x = d6$n[d6$q == "0.01"], y = 0.98, label = "q = 0.01", colour = "#c0392b", hjust = 1.05, size = 4.2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "peak rank (by p-value), bins of 50", y = "peaks with a CTCF motif", title = "Weaker peaks: fewer motifs")
ggsave("fig06-threshold-motif.png", p6a + p6b, width = 12.5, height = 4.8, dpi = 150, bg = "white")
cat("motif by rank bins (q < 0.1 set):", round(d6b$motif, 2), "\n")

# ---- fig07: without an input -----------------------------------------------------------------
noi <- np("macs/ctcf_noinput_peaks.narrowPeak")
extra <- noi[!overlapsAny(noi, pk)]
obs_exp <- function(p, gr) countOverlaps(GRanges("chr21", IRanges(start(p) + p$peak - 250, width = 501)), resize(gr, 1, fix = "start"),
                                          ignore.strand = TRUE) / (length(gr) * 501 / G21)
grp7 <- list(`with input` = pk, `without input` = noi, `found only without input` = extra, `ENCODE (IDR)` = enc)
d7 <- do.call(rbind, lapply(names(grp7), function(nm) { p <- grp7[[nm]]
  data.frame(set = nm, n = length(p), motif = mean(has_motif(p)), input = median(obs_exp(p, reads$ctcf_input)),
             tss = mean(overlapsAny(p, tss[tss$state == "active"], maxgap = 1000)), inENC = mean(overlapsAny(p, enc))) }))
print(d7)
d7$set <- factor(d7$set, names(grp7))
d7l <- rbind(data.frame(set = d7$set, what = "peaks with a CTCF motif", v = d7$motif),
             data.frame(set = d7$set, what = "within 1 kb of an active TSS", v = d7$tss))
p7 <- ggplot(d7l, aes(set, v, fill = set)) + geom_col(width = 0.7) + facet_wrap(~ what) +
  geom_text(aes(label = scales::percent(v, accuracy = 1)), vjust = -0.3, size = 4.3) +
  geom_text(data = d7, aes(set, -0.06, label = paste0("n = ", n)), inherit.aes = FALSE, size = 3.8, colour = "grey30") +
  scale_fill_manual(values = c("#8c564b", "#d9a58c", "#c0392b", "#5b6b7a"), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.08, 1)) +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 12), strip.text = element_text(size = 14))
ggsave("fig07-no-input.png", p7, width = 12, height = 5, dpi = 150, bg = "white")

# ---- fig10/11: broad calling of H3K36me3 -------------------------------------------------
nar <- np("macs/k36_narrow_peaks.narrowPeak")
brd <- import(f("macs/k36_broad_peaks.broadPeak"), format = "BED", extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric"))
ep <- read.delim(f("k36_epic2.bed")); epg <- GRanges(ep$X.Chromosome, IRanges(ep$Start + 1, ep$End))
calls <- list(`MACS3 narrow` = nar, `MACS3 --broad` = brd, epic2 = epg)
gb <- genes[width(genes) > 5000 & genes$state != "other"]
d11 <- do.call(rbind, lapply(names(calls), function(nm) { g <- reduce(granges(calls[[nm]]))
  data.frame(caller = nm, regions = length(calls[[nm]]), kb = sum(width(g)) / 1e3, median_w = median(width(calls[[nm]])),
             active = sum(width(intersect(g, reduce(gb[gb$state == "active"], ignore.strand = TRUE), ignore.strand = TRUE))) / sum(width(reduce(gb[gb$state == "active"], ignore.strand = TRUE))),
             silent = sum(width(intersect(g, reduce(gb[gb$state == "silent"], ignore.strand = TRUE), ignore.strand = TRUE))) / sum(width(reduce(gb[gb$state == "silent"], ignore.strand = TRUE))),
             in_genes = sum(width(intersect(g, reduce(gb[gb$state == "active"], ignore.strand = TRUE), ignore.strand = TRUE))) / sum(width(g))) }))
print(d11)
cat("epic2 domains >= 10 kb:", sum(width(epg) >= 1e4), "; longest", max(width(epg)), "\n")
win10 <- GRanges("chr21", IRanges(33220000, 33600000))
pool <- c(reads$h3k36me3, reads$h3k36me3_rep2)
b10 <- tile(win10, width = 1000)[[1]]
cpm1k <- function(gr) countOverlaps(b10, resize(gr, 200, fix = "start"), ignore.strand = TRUE) * 1e6 / length(gr)
d10 <- rbind(data.frame(track = "H3K36me3 (2 replicates)", x = (start(b10) + 500) / 1e6, v = cpm1k(pool)),
             data.frame(track = "input", x = (start(b10) + 500) / 1e6, v = cpm1k(reads$h3k36me3_input)))
p10a <- ggplot(d10, aes(x, v, fill = track)) + geom_col(width = 0.001) + facet_grid(track ~ ., switch = "y") +
  scale_fill_manual(values = c("#2e8b57", "#5b6b7a"), guide = "none") + labs(x = NULL, y = "reads per kb per M") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12, face = "bold"), strip.placement = "outside", axis.text.x = element_blank())
cdf <- do.call(rbind, lapply(names(calls), function(nm) { g <- subsetByOverlaps(granges(calls[[nm]]), win10)
  if (length(g) == 0) return(NULL); data.frame(caller = nm, s = pmax(start(g), start(win10)) / 1e6, e = pmin(end(g), end(win10)) / 1e6) }))
cdf$caller <- factor(cdf$caller, rev(names(calls)))
gw <- subsetByOverlaps(genes, win10); gw <- sort(gw[width(gw) > 2000 & !grepl("^ENSG|-", gw$gene_name)])
gdf <- data.frame(s = pmax(start(gw), start(win10)) / 1e6, e = pmin(end(gw), end(win10)) / 1e6, name = gw$gene_name, state = gw$state)
gdf$y <- rep(c(0.6, 1.4, 2.2), length.out = nrow(gdf))
p10b <- ggplot() + geom_segment(data = cdf, aes(x = s, xend = pmax(e, s + 0.0008), y = caller, yend = caller), linewidth = 6, colour = "#2e8b57") +
  labs(x = NULL, y = NULL) + theme(axis.text.y = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), panel.grid = element_blank()) +
  coord_cartesian(xlim = c(start(win10), end(win10)) / 1e6)
p10c <- ggplot(gdf) + geom_segment(aes(x = s, xend = e, y = y, yend = y, colour = state), linewidth = 3) +
  geom_text(aes(x = (s + e) / 2, y = y + 0.42, label = name), size = 3.5) +
  scale_colour_manual(values = c(active = "#2e8b57", silent = "#a8324a", other = "#9aa5b1"), guide = "none") +
  scale_y_continuous(limits = c(0.3, 2.9)) + labs(x = "chr21 position (Mb)", y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank()) + coord_cartesian(xlim = c(start(win10), end(win10)) / 1e6)
ggsave("fig10-broad-locus.png", p10a / p10b / p10c + plot_layout(heights = c(2.4, 1.1, 1)), width = 12, height = 6, dpi = 150, bg = "white")
d11l <- rbind(data.frame(caller = d11$caller, what = "active gene bodies", v = d11$active),
              data.frame(caller = d11$caller, what = "silent gene bodies", v = d11$silent))
d11l$caller <- factor(d11l$caller, names(calls))
p11 <- ggplot(d11l, aes(caller, v, fill = what)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = scales::percent(v, accuracy = 0.1)), position = position_dodge(0.8), vjust = -0.3, size = 4.2) +
  scale_fill_manual(values = c("#2e8b57", "#a8324a"), name = NULL) + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "fraction of gene-body bases covered") + theme(legend.position = "top")
ggsave("fig11-broad-genes.png", p11, width = 8.5, height = 4.8, dpi = 150, bg = "white")

# ---- fig13: IDR on the two CTCF replicates ------------------------------------------------------
idr <- read.delim(f("ctcf_idr.txt"), header = FALSE)
idr$pass <- idr$V12 >= -log10(0.05)
idr$r1 <- rank(-idr$V15); idr$r2 <- rank(-idr$V19)
cat("IDR file:", nrow(idr), "matched peaks;", sum(idr$pass), "pass IDR 0.05; Spearman of signals",
    round(cor(idr$V15, idr$V19, method = "spearman"), 2), "\n")
cat(readLines(f("ctcf_idr.log")) |> grep(pattern = "Final parameter|passing", value = TRUE), sep = "\n")
idrg <- GRanges(idr$V1, IRanges(idr$V2 + 1, idr$V3))
cat("IDR peaks in ENCODE conservative set:", sum(overlapsAny(idrg[idr$pass], enc)), "of", sum(idr$pass),
    "; ENCODE peaks recovered:", sum(overlapsAny(enc, idrg[idr$pass])), "of", length(enc), "\n")
p13a <- ggplot(idr, aes(r1, r2, colour = pass)) + geom_point(size = 1.4, alpha = 0.8) +
  scale_colour_manual(values = c(`TRUE` = "#2e8b57", `FALSE` = "#a8324a"), labels = c("IDR ≥ 0.05", "IDR < 0.05"), name = NULL) +
  labs(x = "rank in replicate 1", y = "rank in replicate 2", title = sprintf("%d matched peaks", nrow(idr))) +
  theme(legend.position = c(0.2, 0.88)) + coord_equal()
idr <- idr[order(-idr$V12), ]
idr$global <- 10^(-idr$V12)
p13b <- ggplot(idr, aes(seq_along(global), global)) + geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "#c0392b") +
  annotate("text", x = 5, y = 0.07, label = "IDR = 0.05", colour = "#c0392b", hjust = 0, size = 4.3) +
  geom_vline(xintercept = sum(idr$pass), colour = "#2e8b57") +
  annotate("text", x = sum(idr$pass), y = 0.6, label = paste(sum(idr$pass), "peaks"), colour = "#2e8b57", hjust = -0.1, size = 4.3) +
  labs(x = "peaks, most reproducible first", y = "IDR", title = "Choosing the cutoff")
ggsave("fig13-idr-real.png", p13a + p13b, width = 12, height = 5, dpi = 150, bg = "white")
