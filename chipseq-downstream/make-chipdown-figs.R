# Real-data figures for chipseq-downstream/lec.qmd.
# Run from this folder:  Rscript make-chipdown-figs.R
#
# Needs: data/h3k27ac_{gm12878,k562}_rep{1,2}.chr21.bam, data/input_{gm12878,k562}.chr21.bam and
#        data/h3k27ac_idxstats.tsv (data/prepare-chip.sh); ENCODE H3K27ac peaks of both cell lines, CTCF peaks,
#        GM12878 fold-change bigWigs, blacklist, GENCODE GTF (data/download.sh);
#        on the PATH: bedtools, deepTools (computeMatrix, plotHeatmap, bamCoverage);
#        R packages GenomicRanges, GenomicAlignments, rtracklayer, DESeq2, ChIPseeker, txdbmaker, ggplot2, patchwork.
#
# Work files go to ./down-work (set DOWN_WORK). About 3 min on 8 cores.
# Writes fig03-annotation.png, fig07-tss-heatmap.png, fig08-h3k27ac-heatmap.png, fig12-normalisation.png,
# fig13-copy-number.png, fig15-spikein-sim.png, fig16-db-locus.png and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicAlignments); library(rtracklayer); library(DESeq2)
  library(ChIPseeker); library(txdbmaker); library(ggplot2); library(patchwork)
})
D <- normalizePath("../data")
W <- Sys.getenv("DOWN_WORK", "down-work"); dir.create(W, showWarnings = FALSE, recursive = TRUE); W <- normalizePath(W)
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
theme_set(theme_minimal(base_size = 15) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15)))
col_cell <- c(GM12878 = "#2a6f97", K562 = "#a8324a")
sh <- function(cmd) { rc <- system2("bash", c("-c", shQuote(paste0("set -e; cd ", W, "; D=", D, "; T=", THREADS, "; ", cmd)))); if (rc != 0) stop(cmd) }
step <- function(out, cmd) if (!all(file.exists(file.path(W, out)))) { cat("running:", out, "\n"); sh(cmd) }
f <- function(x) file.path(W, x)

# ---- 1. command-line part (as in the exercise) ----------------------------------------
step("consensus.bed", paste(
  "zcat $D/hg38_blacklist.bed.gz > blacklist.bed;",
  "for c in gm12878 k562; do zcat $D/${c}_h3k27ac_peaks.bed.gz | awk '$1 == \"chr21\"' | sort -k1,1 -k2,2n > ${c}_peaks.bed; done;",
  "cat gm12878_peaks.bed k562_peaks.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge | bedtools intersect -v -a - -b blacklist.bed > consensus.bed;",
  "bedtools intersect -v -a consensus.bed -b k562_peaks.bed > GM12878_only.bed;",
  "bedtools intersect -v -a consensus.bed -b gm12878_peaks.bed > K562_only.bed;",
  "bedtools intersect -u -a consensus.bed -b gm12878_peaks.bed | bedtools intersect -u -a - -b k562_peaks.bed > both.bed"))
for (s in c("gm12878_rep1", "gm12878_rep2", "k562_rep1", "k562_rep2"))
  step(sprintf("h3k27ac_%s.cpm.bw", s), sprintf("bamCoverage -b $D/h3k27ac_%s.chr21.bam -o h3k27ac_%s.cpm.bw --normalizeUsing CPM --extendReads 200 --binSize 50 --blackListFileName blacklist.bed -p $T 2> /dev/null", s, s))
step("genes_chr21.bed", paste(
  "zcat $D/gencode.v46.basic.annotation.gtf.gz | awk -F'\\t' '$1 == \"chr21\" && $3 == \"gene\" && $9 ~ /gene_type \"protein_coding\"/'",
  "| awk -F'\\t' -v OFS='\\t' '{ match($9, /gene_name \"[^\"]+\"/); n = substr($9, RSTART + 11, RLENGTH - 12); print $1, $4 - 1, $5, n, 0, $7 }'",
  "| sort -k1,1 -k2,2n > genes_chr21.bed"))
step("tss_heatmap.png", paste(
  "computeMatrix reference-point --referencePoint TSS -R genes_chr21.bed",
  "-S $D/gm12878_h3k4me3_fc.bigWig $D/gm12878_h3k27ac_fc.bigWig $D/gm12878_h3k36me3_fc.bigWig $D/gm12878_h3k27me3_fc.bigWig",
  "--samplesLabel H3K4me3 H3K27ac H3K36me3 H3K27me3 -a 5000 -b 5000 --binSize 100 --missingDataAsZero --skipZeros -p $T -o tss.mat.gz 2> /dev/null;",
  "plotHeatmap -m tss.mat.gz -o tss_heatmap.png --sortUsingSamples 1 --colorMap Reds Oranges Greens Blues --zMax 30 30 6 6",
  "--heatmapHeight 12 --refPointLabel TSS --regionsLabel 'protein-coding genes'"))
step("h3k27ac_heatmap.png", paste(
  "computeMatrix reference-point --referencePoint center -R GM12878_only.bed both.bed K562_only.bed",
  "-S h3k27ac_gm12878_rep1.cpm.bw h3k27ac_gm12878_rep2.cpm.bw h3k27ac_k562_rep1.cpm.bw h3k27ac_k562_rep2.cpm.bw",
  "--samplesLabel GM12878_1 GM12878_2 K562_1 K562_2 -a 3000 -b 3000 --binSize 50 -p $T -o h3k27ac.mat.gz 2> /dev/null;",
  "plotHeatmap -m h3k27ac.mat.gz -o h3k27ac_heatmap.png --colorMap Oranges --zMax 150 --heatmapHeight 12 --refPointLabel centre",
  "--regionsLabel 'GM12878 only' 'both' 'K562 only'"))
file.copy(f("tss_heatmap.png"), "fig07-tss-heatmap.png", overwrite = TRUE)
file.copy(f("h3k27ac_heatmap.png"), "fig08-h3k27ac-heatmap.png", overwrite = TRUE)
cat("consensus:", length(readLines(f("consensus.bed"))), " GM12878 only:", length(readLines(f("GM12878_only.bed"))),
    " both:", length(readLines(f("both.bed"))), " K562 only:", length(readLines(f("K562_only.bed"))), "\n")

# ---- annotation (fig03) -------------------------------------------------------------------
chr21 <- GRanges("chr21", IRanges(1, 46709983))
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), which = chr21)
seqlevels(gtf, pruning.mode = "coarse") <- "chr21"
txdb <- suppressWarnings(makeTxDbFromGRanges(gtf[gtf$type %in% c("gene", "transcript", "exon", "CDS", "UTR")]))
np <- function(x) { p <- import(x, format = "narrowPeak"); p[seqnames(p) == "chr21"] }
sets <- list(CTCF = np(file.path(D, "gm12878_ctcf_peaks.bed.gz")),
             `H3K27ac GM12878` = np(file.path(D, "gm12878_h3k27ac_peaks.bed.gz")),
             `H3K27ac K562` = np(file.path(D, "k562_h3k27ac_peaks.bed.gz")))
ann <- lapply(sets, annotatePeak, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = FALSE)
st <- do.call(rbind, lapply(names(ann), function(n) data.frame(set = n, getAnnoStat(ann[[n]]))))
st$Feature <- as.character(st$Feature)
st$class <- ifelse(st$Feature == "Promoter", "promoter (±1 kb)", ifelse(grepl("Exon|UTR", st$Feature), "exon / UTR",
            ifelse(grepl("Intron", st$Feature), "intron", "distal intergenic")))
st2 <- aggregate(Frequency ~ set + class, st, sum)
st2$class <- factor(st2$class, c("promoter (±1 kb)", "exon / UTR", "intron", "distal intergenic"))
print(reshape(st2, idvar = "set", timevar = "class", direction = "wide"))
p3a <- ggplot(st2, aes(set, Frequency / 100, fill = class)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c("#c0392b", "#e9c46a", "#2a9d8f", "#9aa5b1"), name = NULL) +
  scale_y_continuous(labels = scales::percent) + coord_flip() + labs(x = NULL, y = "fraction of peaks") +
  theme(legend.position = "top")
dd <- do.call(rbind, lapply(names(ann), function(n) data.frame(set = n, d = as.data.frame(ann[[n]])$distanceToTSS)))
print(aggregate(d ~ set, dd, function(x) c(within1kb = mean(abs(x) <= 1000), over10kb = mean(abs(x) > 1e4), median = median(abs(x)))))
p3b <- ggplot(dd, aes(abs(d) + 1, colour = set)) + stat_ecdf(linewidth = 1.2) +
  scale_x_log10(breaks = c(1, 100, 1e3, 1e4, 1e5, 1e6), labels = c("0", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb")) +
  scale_colour_manual(values = c("#8c564b", "#2a6f97", "#a8324a"), name = NULL) +
  labs(x = "distance to the nearest TSS (any gene)", y = "cumulative fraction of peaks") + theme(legend.position = c(0.25, 0.85))
ggsave("fig03-annotation.png", p3a + p3b + plot_layout(widths = c(1.2, 1)), width = 13, height = 4.8, dpi = 150, bg = "white")

# ---- counts and DESeq2 --------------------------------------------------------------------------
bl <- import(file.path(D, "hg38_blacklist.bed.gz"))
cons <- import(f("consensus.bed"))
smp <- c(GM_1 = "h3k27ac_gm12878_rep1", GM_2 = "h3k27ac_gm12878_rep2", K562_1 = "h3k27ac_k562_rep1", K562_2 = "h3k27ac_k562_rep2",
         GM_input = "input_gm12878", K562_input = "input_k562")
centre <- function(s) { r <- granges(readGAlignments(file.path(D, paste0(s, ".chr21.bam")))); r <- r[!overlapsAny(r, bl)]
  resize(resize(r, 200, fix = "start"), 1, fix = "center") }
mid <- lapply(smp, centre)
cnt <- sapply(mid, function(m) countOverlaps(cons, m, ignore.strand = TRUE))
cat("reads:", sapply(mid, length), "\nin peaks:", colSums(cnt), "\nFRiP:", round(colSums(cnt) / sapply(mid, length), 3), "\n")
cd <- data.frame(cell = factor(c("GM12878", "GM12878", "K562", "K562")), row.names = names(smp)[1:4])
dds <- DESeq(DESeqDataSetFromMatrix(cnt[, 1:4], cd, ~ cell), quiet = TRUE)
lib <- sapply(mid[1:4], length)
bins <- tileGenome(c(chr21 = 46709983), tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
bins <- bins[start(bins) > 13e6 & !overlapsAny(bins, cons, maxgap = 5000) & !overlapsAny(bins, bl)]
bg <- sapply(mid[1:4], function(m) countOverlaps(bins, m, ignore.strand = TRUE))
sf <- cbind(`reads in peaks (RLE, DESeq2 default)` = sizeFactors(dds), `library size` = lib / exp(mean(log(lib))),
            `background bins (csaw-style)` = estimateSizeFactorsForMatrix(bg))
print(round(sf, 3))
test <- function(s) { d <- dds; sizeFactors(d) <- s; results(DESeq(d, quiet = TRUE), contrast = c("cell", "K562", "GM12878")) }
rs <- lapply(as.data.frame(sf, check.names = FALSE), test)
tab <- sapply(rs, function(r) c(up = sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
                                down = sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE), median = median(r$log2FoldChange)))
print(tab)
d12 <- do.call(rbind, lapply(names(rs), function(n) data.frame(norm = n, A = log2(rs[[n]]$baseMean), M = rs[[n]]$log2FoldChange,
                                                                sig = !is.na(rs[[n]]$padj) & rs[[n]]$padj < 0.05)))
d12$norm <- factor(d12$norm, names(rs))
lab12 <- data.frame(norm = factor(names(rs), names(rs)), lab = sprintf("K562 up %d, down %d\nmedian log2FC %.2f", tab["up", ], tab["down", ], tab["median", ]))
p12 <- ggplot(d12, aes(A, M, colour = sig)) + geom_point(size = 0.9, alpha = 0.7) + geom_hline(yintercept = 0) +
  geom_hline(data = data.frame(norm = factor(names(rs), names(rs)), m = tab["median", ]), aes(yintercept = m), linetype = 2, colour = "#c0392b") +
  geom_label(data = lab12, aes(x = 0.2, y = -9, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3.8, fill = "white", label.size = 0) +
  facet_wrap(~ norm) + scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#1b2a3a"), guide = "none") +
  coord_cartesian(ylim = c(-9, 9)) + labs(x = "log2 mean normalised count", y = "log2 fold change, K562 / GM12878") +
  theme(strip.text = element_text(size = 12.5))
ggsave("fig12-normalisation.png", p12, width = 13, height = 4.8, dpi = 150, bg = "white")

# ---- copy number from the inputs (fig13) ---------------------------------------------------------------
ix <- read.delim(file.path(D, "h3k27ac_idxstats.tsv"))
ix <- ix[ix$chrom %in% paste0("chr", c(1:22, "X")), ]
share <- with(ix, tapply(mapped, list(chrom, sample), sum)); share <- sweep(share, 2, colSums(share), "/")
cn <- data.frame(chrom = rownames(share), input = share[, "input_k562"] / share[, "input_gm12878"],
                 h3k27ac = (share[, "h3k27ac_k562_rep1"] + share[, "h3k27ac_k562_rep2"]) / (share[, "h3k27ac_gm12878_rep1"] + share[, "h3k27ac_gm12878_rep2"]))
cn$chrom <- factor(cn$chrom, paste0("chr", c(1:22, "X")))
print(cn[order(cn$chrom), ], digits = 3)
p13 <- ggplot(cn, aes(chrom, input, fill = chrom == "chr21")) + geom_col(width = 0.75) + geom_hline(yintercept = 1, linetype = 2) +
  geom_point(aes(y = h3k27ac), shape = 21, fill = "white", size = 2.5) +
  scale_fill_manual(values = c(`TRUE` = "#a8324a", `FALSE` = "#9aa5b1"), guide = "none") +
  labs(x = NULL, y = "share of reads, K562 / GM12878", title = "Bars: inputs.  Circles: H3K27ac ChIPs.") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11))
ggsave("fig13-copy-number.png", p13, width = 11, height = 4.4, dpi = 150, bg = "white")
ip <- cnt[, c("GM_input", "K562_input")]
cat("input reads in consensus peaks per million input reads: GM12878", round(sum(ip[, 1]) / length(mid$GM_input) * 1e6),
    " K562", round(sum(ip[, 2]) / length(mid$K562_input) * 1e6), "\n")

# ---- spike-in simulation (fig15) ---------------------------------------------------------------------------
set.seed(42)
n <- 2000; mu <- rlnorm(n, log(100), 1)
change <- rep(0.5, n); change[1:100] <- 1
truth <- cbind(mu, mu, mu, mu * change, mu * change, mu * change)
depth <- 1e6 / colSums(truth)
counts <- matrix(rnbinom(length(truth), mu = sweep(truth, 2, depth, "*") * 0.2, size = 20), n)
spike <- c(5000, 5000, 5000, 10000, 10000, 10000) * c(1, 1.05, 0.95, 1, 0.97, 1.03)
colnames(counts) <- paste0(c("c", "c", "c", "t", "t", "t"), 1:6)
cd2 <- data.frame(group = factor(rep(c("ctrl", "trt"), each = 3), c("ctrl", "trt")), row.names = colnames(counts))
d2 <- DESeqDataSetFromMatrix(counts, cd2, ~ group)
r_rle <- results(DESeq(d2, quiet = TRUE)); sizeFactors(d2) <- spike / exp(mean(log(spike))); r_spk <- results(DESeq(d2, quiet = TRUE))
sim <- rbind(data.frame(norm = "standard (RLE)", A = log2(r_rle$baseMean), M = r_rle$log2FoldChange, kept = seq_len(n) <= 100, sig = r_rle$padj < 0.05),
             data.frame(norm = "spike-in size factors", A = log2(r_spk$baseMean), M = r_spk$log2FoldChange, kept = seq_len(n) <= 100, sig = r_spk$padj < 0.05))
sim$norm <- factor(sim$norm, c("standard (RLE)", "spike-in size factors"))
print(aggregate(cbind(M, sig) ~ norm + kept, sim, function(x) round(mean(x), 3)))
p15 <- ggplot(sim[order(sim$kept), ], aes(A, M, colour = kept)) + geom_point(size = 0.8, alpha = 0.7) + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  facet_wrap(~ norm) + scale_colour_manual(values = c(`FALSE` = "#2a6f97", `TRUE` = "#e67e22"),
                                           labels = c("95% of peaks: truly halved", "5% of peaks: truly unchanged"), name = NULL) +
  coord_cartesian(ylim = c(-3, 2.5)) + labs(x = "log2 mean normalised count", y = "log2 fold change, treated / control") +
  theme(legend.position = "top", strip.text = element_text(size = 14))
ggsave("fig15-spikein-sim.png", p15, width = 11.5, height = 4.8, dpi = 150, bg = "white")

# ---- a differential locus (fig16) ---------------------------------------------------------------------------
r <- rs[[1]]
cons$lfc <- r$log2FoldChange; cons$sig <- !is.na(r$padj) & r$padj < 0.05
pc <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
win <- GRanges("chr21", IRanges(25300000, 26250000))
b <- tile(win, width = 2000)[[1]]
cpm <- function(m) countOverlaps(b, m, ignore.strand = TRUE) * 1e6 / length(m)
d16 <- do.call(rbind, lapply(1:4, function(i) data.frame(sample = names(smp)[i], cell = ifelse(i <= 2, "GM12878", "K562"),
                                                         x = (start(b) + 1000) / 1e6, v = cpm(mid[[i]]))))
d16$sample <- factor(d16$sample, names(smp)[1:4])
pw <- subsetByOverlaps(cons, win)
pdf16 <- data.frame(s = start(pw) / 1e6, e = end(pw) / 1e6, dir = ifelse(!pw$sig, "not significant", ifelse(pw$lfc > 0, "higher in K562", "higher in GM12878")))
gw <- subsetByOverlaps(pc, win)
p16a <- ggplot(d16, aes(x, v, fill = cell)) + geom_col(width = 0.002) + facet_grid(sample ~ ., switch = "y") +
  scale_fill_manual(values = col_cell, guide = "none") + labs(x = NULL, y = "CPM, 2-kb bins") +
  theme(strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"), strip.placement = "outside", axis.text.x = element_blank())
p16b <- ggplot() + geom_segment(data = pdf16, aes(x = s, xend = pmax(e, s + 0.002), y = 1, yend = 1, colour = dir), linewidth = 7) +
  geom_segment(data = data.frame(s = start(gw) / 1e6, e = end(gw) / 1e6), aes(x = s, xend = e, y = 0.45, yend = 0.45), linewidth = 3, colour = "#1b2a3a") +
  geom_text(data = data.frame(x = (start(gw) + end(gw)) / 2e6, n = gw$gene_name, y = rep(c(0.15, -0.15), length.out = length(gw))), aes(x, y, label = n), size = 3.8) +
  scale_colour_manual(values = c(`higher in GM12878` = "#2a6f97", `higher in K562` = "#a8324a", `not significant` = "#9aa5b1"), name = NULL) +
  scale_y_continuous(limits = c(-0.3, 1.3)) + coord_cartesian(xlim = c(start(win), end(win)) / 1e6) +
  labs(x = "chr21 position (Mb)", y = NULL) + theme(axis.text.y = element_blank(), panel.grid = element_blank(), legend.position = "top")
ggsave("fig16-db-locus.png", p16a / p16b + plot_layout(heights = c(3, 1.3)), width = 12, height = 5.6, dpi = 150, bg = "white")
cat("fig16 peaks:", nrow(pdf16), table(pdf16$dir), "genes:", gw$gene_name, "\n")
