# Real-data figures for methylation-bisulfite/lec.qmd.
# Run from this folder:  Rscript make-bisulfite-figs.R   (about 1 minute, 1 GB RAM,
# plus about 5 minutes for the Bismark run if it is not there yet)
# Needs: data/download.sh methylation-bisulfite (GM12878 WGBS reads for
#        chr21:33-35 Mb, the paired-end M-bias table, ENCODE per-CpG calls
#        for chr21, hg38 chr21); fastp and Bismark on the PATH;
#        R packages data.table, ggplot2, patchwork.
# The Bismark outputs are those of the exercise (run-bismark.sh). Set BS_RUN to
# an existing run folder to reuse it.
#
# Writes fig05, fig09, fig11, fig12, fig13, fig14, fig16 (*.png) and prints the
# numbers quoted in the slides.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(1)
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
run <- Sys.getenv("BS_RUN", file.path(tempdir(), "bismark-run"))
if (!file.exists(file.path(run, "gm_cpg.CpG_report.merged_CpG_evidence.cov")))
  system2("bash", c("run-bismark.sh", run))
f <- function(x) file.path(run, x)

# ---- fig05: base composition of bisulfite reads vs the genome ----------------------
reads <- readLines(gzfile(file.path(D, "gm12878_wgbs_reads.fastq.gz")), n = 400000)
seqs <- reads[seq(2, length(reads), 4)]
bc <- table(factor(strsplit(paste(seqs, collapse = ""), "")[[1]], levels = c("A", "C", "G", "T")))
g <- toupper(paste(readLines(gzfile(file.path(D, "hg38_chr21.fa.gz")))[-1], collapse = ""))
gc <- table(factor(strsplit(g, "")[[1]], levels = c("A", "C", "G", "T")))
comp <- rbind(data.table(source = "chr21 reference", base = names(gc), frac = as.numeric(gc / sum(gc))),
              data.table(source = "bisulfite reads (read 1)", base = names(bc), frac = as.numeric(bc / sum(bc))))
print(dcast(comp, base ~ source, value.var = "frac"))
p05 <- ggplot(comp, aes(base, frac, fill = source)) + geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * frac)), position = position_dodge(0.7), vjust = -0.3, size = 4.5) +
  scale_fill_manual(values = c(`chr21 reference` = "#5b6b7a", `bisulfite reads (read 1)` = "#2a6f97"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.55)) +
  labs(x = NULL, y = "fraction of bases", title = "Base composition: GM12878 WGBS reads vs the genome") +
  theme(legend.position = "top")
save(p05, "fig05-base-composition.png", w = 9, h = 4.6)

# ---- fig09: what happens to the reads ------------------------------------------------
fp <- jsonlite::fromJSON(f("fastp.json"))
al <- readLines(f("trimmed_bismark_bt2_SE_report.txt"))
dd <- readLines(f("trimmed_bismark_bt2.deduplication_report.txt"))
num <- function(lines, pat) as.numeric(sub(".*:\\s*([0-9]+).*", "\\1", grep(pat, lines, value = TRUE)[1]))
fate <- data.table(step = c("reads", "after trimming\n(fastp)", "aligned uniquely\n(Bismark)", "after removing\nduplicates"),
                   n = c(fp$summary$before_filtering$total_reads, fp$summary$after_filtering$total_reads,
                         num(al, "unique best hit"), num(dd, "leftover sequences")))
fate$step <- factor(fate$step, levels = fate$step)
print(fate)
p09 <- ggplot(fate, aes(step, n / 1e3)) + geom_col(fill = "#2a6f97", width = 0.65) +
  geom_text(aes(label = format(n, big.mark = " ")), vjust = -0.4, size = 5) +
  scale_y_continuous(limits = c(0, 420)) +
  labs(x = NULL, y = "thousands of reads", title = "GM12878 WGBS, chr21:33-35 Mb, single-end 125 bp")
save(p09, "fig09-read-fate.png", w = 9, h = 4.4)

# ---- fig11: M-bias -------------------------------------------------------------------
read_mbias <- function(file, label) {
  l <- readLines(file); out <- list(); cur <- NULL
  for (x in l) {
    if (grepl("context", x)) { cur <- x; next }
    if (grepl("^[0-9]+\t", x)) {
      v <- strsplit(x, "\t")[[1]]
      out[[length(out) + 1]] <- data.table(section = cur, pos = as.integer(v[1]), meth = as.numeric(v[2]),
                                          unmeth = as.numeric(v[3]))
    }
  }
  x <- rbindlist(out)
  x[, context := sub(" context.*", "", section)]
  x[, read := ifelse(grepl("R2", section), "read 2", "read 1")]
  x[, sample := label][, pct := 100 * meth / (meth + unmeth)][, cov := meth + unmeth]
  x
}
mb <- rbind(read_mbias(file.path(D, "gm12878_wgbs_pe_mbias.txt"), "ENCODE paired-end (3 M pairs)"),
            read_mbias(f("trimmed_bismark_bt2.deduplicated.M-bias.txt"), "this exercise (read 1 only)"))
mb <- mb[cov >= 2000]
print(mb[context == "CpG" & pos <= 5, .(sample, read, pos, pct = round(pct, 1), cov)])
print(mb[context == "CpG" & pos >= 30 & pos <= 100, .(plateau = round(median(pct), 1)), by = .(sample, read)])
print(mb[context == "CHH" & pos <= 3, .(sample, read, pos, pct = round(pct, 2))])
cols <- c(`read 1` = "#2a6f97", `read 2` = "#c0392b")
p11a <- ggplot(mb[context == "CpG" & grepl("ENCODE", sample)], aes(pos, pct, colour = read)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = cols, name = NULL) + coord_cartesian(ylim = c(35, 70)) +
  labs(x = "position in read (bp)", y = "% CpG methylated", title = "ENCODE GM12878, paired-end") +
  theme(legend.position = c(0.75, 0.88))
p11b <- ggplot(mb[context == "CpG" & !grepl("ENCODE", sample)], aes(pos, pct)) + geom_line(linewidth = 1, colour = cols[1]) +
  coord_cartesian(ylim = c(45, 75)) +
  labs(x = "position in read (bp)", y = "% CpG methylated", title = "This exercise, read 1")
save(p11a + p11b, "fig11-mbias.png", w = 12, h = 4.4)

# ---- fig12: methylation by context: conversion -------------------------------------
sp <- readLines(f("trimmed_bismark_bt2.deduplicated_splitting_report.txt"))
ctx <- rbindlist(lapply(c("CpG", "CHG", "CHH"), function(k)
  data.table(context = k, M = num(sp, paste0("Total methylated C's in ", k)),
             U = num(sp, paste0("Total C to T conversions in ", k)))))
ctx[, pct := 100 * M / (M + U)]
print(ctx)
cat("conversion rate from CHH:", round(100 - ctx[context == "CHH", pct], 2), "%\n")
ctx$context <- factor(ctx$context, levels = c("CpG", "CHG", "CHH"))
p12 <- ggplot(ctx, aes(context, pct)) + geom_col(fill = c("#1b2a3a", "#5b6b7a", "#5b6b7a"), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%s calls)", pct, formatC(M + U, big.mark = " ", format = "d"))), vjust = -0.15, size = 4.5) +
  scale_y_continuous(limits = c(0, 85)) +
  labs(x = "cytosine context", y = "% of calls methylated", title = "Methylation calls by context")
save(p12, "fig12-contexts.png", w = 6.5, h = 4.6)

# ---- per-CpG table: ours and ENCODE ----------------------------------------------------
ours <- fread(f("gm_cpg.CpG_report.merged_CpG_evidence.cov"), col.names = c("chr", "start", "end", "pct", "M", "U"))
ours[, N := M + U]
enc <- fread(file.path(D, "gm12878_wgbs_rep1_chr21.bed.gz"), header = FALSE, select = c(1, 2, 6, 10, 11, 13),
             col.names = c("chr", "s0", "strand", "N", "p", "call"))
enc <- enc[s0 >= 33e6 - 1 & s0 < 35e6 & call == "CG"]
enc[, start := ifelse(strand == "+", s0, s0 - 1L)][, M := round(N * p / 100)]
enc <- enc[, .(M = sum(M), N = sum(N)), by = start]
cat("our CpGs:", nrow(ours), " in 33-35 Mb:", nrow(ours[start >= 33e6 & start < 35e6]),
    " median coverage", median(ours$N), " >= 10 reads:", round(mean(ours$N >= 10), 3), "\n")

# ---- fig13: coverage and precision ------------------------------------------------------
p13a <- ggplot(ours[start >= 33e6 & start < 35e6], aes(pmin(N, 50))) + geom_histogram(binwidth = 1, fill = "#2a6f97") +
  labs(x = "reads per CpG (both strands; 50 = 50 or more)", y = "CpGs", title = "Coverage, chr21:33-35 Mb")
ci <- data.table(N = 1:60)[, half := 1.96 * sqrt(0.25 / N)]
p13b <- ggplot(ci, aes(N, half)) + geom_line(linewidth = 1.2, colour = "#c0392b") +
  geom_hline(yintercept = c(0.1, 0.2), linetype = 2, colour = "grey60") +
  annotate("point", x = 14, y = 1.96 * sqrt(0.25 / 14), size = 3) +
  annotate("text", x = 16, y = 1.96 * sqrt(0.25 / 14) + 0.04, label = "14 reads: ± 0.26", hjust = 0, size = 4.8) +
  labs(x = "reads covering the CpG", y = "95% interval half-width", title = "Precision of one CpG at β = 0.5")
save(p13a + p13b, "fig13-coverage-ci.png", w = 12, h = 4.4)

# ---- fig14: our calls vs ENCODE -------------------------------------------------------------
m <- merge(ours[, .(start, Mo = M, No = N)], enc[, .(start, Me = M, Ne = N)], by = "start")
fm <- m[No >= 10 & Ne >= 10]
cat("shared CpGs:", nrow(m), " identical counts:", round(mean(m$Mo == m$Me & m$No == m$Ne), 3),
    "  >= 10 reads in both:", nrow(fm), " Pearson r:", round(cor(fm$Mo / fm$No, fm$Me / fm$Ne), 4),
    " mean |diff|:", round(mean(abs(fm$Mo / fm$No - fm$Me / fm$Ne)), 4), "\n")
p14 <- ggplot(fm, aes(Me / Ne, Mo / No)) + geom_bin2d(bins = 40) +
  scale_fill_gradient(low = "#dfe6ec", high = "#1b2a3a", trans = "log10", name = "CpGs") +
  geom_abline(linetype = 2, colour = "#c0392b") +
  labs(x = "ENCODE (gemBS pipeline)", y = "this exercise (Bismark)",
       title = sprintf("%s CpGs with >= 10 reads in both: r = %.3f", format(nrow(fm), big.mark = " "),
                       cor(fm$Mo / fm$No, fm$Me / fm$Ne))) + coord_equal()
save(p14, "fig14-compare-encode.png", w = 7, h = 6)

# ---- fig16: SNPs at CpGs ---------------------------------------------------------------------
x <- fread(file.path(D, "gm12878_wgbs_rep1_chr21.bed.gz"), header = FALSE, select = c(2, 6, 13),
           col.names = c("s0", "strand", "call"))[s0 >= 33e6 - 1 & s0 < 35e6]
x[, start := ifelse(strand == "+", s0, s0 - 1L)]
# one class per CpG (the two strands report the same SNP: YG on + is CR on -, TG is CA)
snp <- x[call != "CG", .(calls = paste(sort(unique(call)), collapse = ",")), by = start]
snp[, grp := ifelse(grepl("TG|CA", calls), "homozygous C>T:\nno CpG at all",
              ifelse(grepl("YG|CR", calls), "heterozygous C/T:\nCpG on one allele", "other SNP\n(C/A, C/G, G/T)"))]
o <- merge(ours[N >= 8], snp[, .(start, grp)], by = "start", all.x = TRUE)
o[is.na(grp), grp := "no SNP"]
lv <- c("no SNP", "heterozygous C/T:\nCpG on one allele", "homozygous C>T:\nno CpG at all", "other SNP\n(C/A, C/G, G/T)")
s16 <- o[, .(beta = sum(M) / sum(N), n = .N), by = grp][, grp := factor(grp, levels = lv)][order(grp)]
print(s16)
p16 <- ggplot(s16, aes(grp, beta)) + geom_col(fill = c("#1b2a3a", "#c0392b", "#a8324a", "#5b6b7a"), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f\n(n = %s)", beta, formatC(n, big.mark = " ", format = "d"))), vjust = -0.15, size = 4.5) +
  scale_y_continuous(limits = c(0, 0.85)) +
  labs(x = "genotype of GM12878 at the CpG (called by ENCODE)", y = "methylation reported by Bismark",
       title = "C>T SNPs look like unmethylated CpGs (chr21:33-35 Mb, >= 8 reads)")
save(p16, "fig16-snp-data.png", w = 11, h = 4.8)
