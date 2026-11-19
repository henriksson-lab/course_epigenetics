# Real-data figures for longread-epigenomics/lec.qmd.
# Run from this folder:  Rscript make-longread-figs.R      (about 5 minutes)
# Needs: data/prepare-meth2.sh longreads (HG001 nanopore modBAM regions, trio SNPs,
#        GM12878 WGBS subset, ONT haplotype pileups of chrX/chr20, HG002 EcoGII reads)
#        and data/download.sh (hg38.fa.bgz, GENCODE v46, GM12878 ATAC peaks).
#        Command-line tools: samtools (>= 1.15), bcftools, modkit (>= 0.6), whatshap (>= 2.0);
#        set MODKIT=/path/to/modkit and WHATSHAP=/path/to/whatshap if they are not on PATH.
#        R packages: data.table, ggplot2, patchwork, ggrepel, GenomicRanges, rtracklayer.
# Intermediate files go to $LR_CACHE (default: a folder in tempdir()).
# Writes fig05, fig07, fig08, fig10, fig11, fig12, fig13 (*.png) and prints
# the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(GenomicRanges); library(rtracklayer)
})
set.seed(1)
D <- "../data"
C <- Sys.getenv("LR_CACHE", file.path(tempdir(), "lr")); dir.create(C, showWarnings = FALSE, recursive = TRUE)
MODKIT <- Sys.getenv("MODKIT", "modkit"); WHATSHAP <- Sys.getenv("WHATSHAP", "whatshap")
SAMTOOLS <- Sys.getenv("SAMTOOLS", "samtools")
run <- function(cmd) { cat("$", cmd, "\n"); stopifnot(system(cmd) == 0) }
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
colHP <- c(paternal = "#e67e22", maternal = "#2a6f97")
bedm <- c("chrom", "start", "end", "code", "score", "strand", "s2", "e2", "color",
          "nvalid", "pct", "nmod", "ncanon", "nother", "ndel", "nfail", "ndiff", "nnocall")
readbm <- function(f) { x <- fread(f, header = FALSE); setnames(x, bedm[seq_len(ncol(x))]); x }

# ---- command-line steps (cached) -------------------------------------------------
ref <- file.path(C, "hg38_lr.fa")
if (!file.exists(paste0(ref, ".fai"))) {
  run(sprintf("%s faidx %s/hg38.fa.bgz chr11 chr15 chr21 chrX > %s", SAMTOOLS, D, ref))
  run(sprintf("%s faidx %s", SAMTOOLS, ref))
}
allbed <- file.path(C, "hg001_cpg.bed")
if (!file.exists(allbed))
  run(sprintf("%s pileup %s/hg001_ont_regions.bam %s --ref %s --cpg --modified-bases 5mC 5hmC --combine-strands --include-bed %s/longread_regions.bed --sampling-frac 1 -t 8 --log-filepath %s/mk1.log",
              MODKIT, D, allbed, ref, D, C))
vcf <- file.path(C, "trio_phased.vcf")
if (!file.exists(paste0(vcf, ".gz"))) {
  awk <- paste0("awk 'BEGIN { OFS = \"\\t\"; print \"##fileformat=VCFv4.2\"; ",
    "print \"##FORMAT=<ID=GT,Number=1,Type=String,Description=\\\"Genotype paternal|maternal\\\">\"; ",
    "print \"##contig=<ID=chr11,length=135086622>\"; print \"##contig=<ID=chr15,length=101991189>\"; ",
    "print \"##contig=<ID=chr21,length=46709983>\"; print \"##contig=<ID=chrX,length=156040895>\"; ",
    "print \"#CHROM\", \"POS\", \"ID\", \"REF\", \"ALT\", \"QUAL\", \"FILTER\", \"INFO\", \"FORMAT\", \"NA12878\" } ",
    "{ gsub(/\\|/, \"/\", $5); gsub(/\\|/, \"/\", $6); gsub(/\\|/, \"/\", $7); ",
    "if ($5 != \"0/1\" && $5 != \"1/0\") next; f = $6; m = $7; pat = \"\"; ",
    "if (f == \"0/0\" || f == \"0\") pat = 0; else if (f == \"1/1\" || f == \"1\") pat = 1; ",
    "else if (m == \"0/0\") pat = 1; else if (m == \"1/1\") pat = 0; ",
    "if (pat == \"\") next; print $1, $2, \".\", $3, $4, \".\", \"PASS\", \".\", \"GT\", pat \"|\" 1 - pat }'")
  run(sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s/hg001_trio_regions.vcf.gz | %s > %s", D, awk, vcf))
  run(sprintf("bgzip -f %s && tabix -f -p vcf %s.gz", vcf, vcf))
}
tagged <- file.path(C, "hg001_tagged.bam")
if (!file.exists(paste0(tagged, ".bai"))) {
  run(sprintf("%s haplotag --ignore-read-groups --reference %s --output-haplotag-list %s/haplotags.tsv -o %s %s.gz %s/hg001_ont_regions.bam > %s/whatshap.log 2>&1",
              WHATSHAP, ref, C, tagged, vcf, D, C))
  run(sprintf("%s index %s", SAMTOOLS, tagged))
}
if (!file.exists(file.path(C, "phased", "hp1.bedmethyl")))
  run(sprintf("%s pileup %s %s/phased --phased --ref %s --cpg --modified-bases 5mC 5hmC --combine-strands --include-bed %s/longread_regions.bed --sampling-frac 1 -t 8 --log-filepath %s/mk2.log",
              MODKIT, tagged, C, ref, D, C))
h19calls <- file.path(C, "h19_calls.tsv")
if (!file.exists(h19calls))
  run(sprintf("%s extract calls %s %s --region chr11:1997000-2004000 --reference %s --mapped-only --cpg --sampling-frac 1 -t 8 --log-filepath %s/mk3.log",
              MODKIT, tagged, h19calls, ref, C))
atpcalls <- file.path(C, "atp5po_calls.tsv")
if (!file.exists(atpcalls))
  run(sprintf("%s extract calls %s/hg002_ecogii_regions.bam %s --region chr21:33913000-33919000 --reference %s --mapped-only --sampling-frac 1 -t 8 --log-filepath %s/mk4.log",
              MODKIT, D, atpcalls, ref, C))
for (x in c("ecogii", "native")) {
  f <- file.path(C, paste0("m6a_", x, ".bed"))
  if (!file.exists(f))
    run(sprintf("%s pileup %s/hg002_%s_regions.bam %s --ref %s --modified-bases 6mA --region chr21:33000000-34000000 --sampling-frac 1 -t 8 --log-filepath %s/mk_%s.log",
                MODKIT, D, x, f, ref, C, x))
}

# ---- annotation -----------------------------------------------------------------
gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), feature.type = "exon")
gtf <- gtf[seqnames(gtf) %in% c("chr11", "chr15", "chr21", "chrX")]
gene_track <- function(chr, from, to, genes) {
  e <- as.data.table(gtf[seqnames(gtf) == chr & end(gtf) >= from & start(gtf) <= to & gtf$gene_name %in% genes])
  e[, canon := grepl("Ensembl_canonical", sapply(tag, paste, collapse = ","))]
  e <- e[, if (any(canon)) .SD[canon == TRUE] else .SD[transcript_id == transcript_id[1]], by = gene_name]  # H19 has no canonical tag
  e <- e[, .(start = pmax(start, from), end = pmin(end, to), strand, gene_name)]
  tx <- e[, .(start = min(start), end = max(end), strand = strand[1]), by = gene_name]
  ggplot() + geom_segment(data = tx, aes(x = start, xend = end, y = 0, yend = 0), colour = "#1b2a3a") +
    geom_rect(data = e, aes(xmin = start, xmax = end, ymin = -0.3, ymax = 0.3), fill = "#1b2a3a") +
    geom_text(data = tx, aes(x = (start + end) / 2, y = 0.75,
                             label = paste0(gene_name, ifelse(strand == "+", "  >", "  <"))), size = 4.3) +
    scale_x_continuous(limits = c(from, to), labels = function(x) sprintf("%.3f Mb", x / 1e6)) +
    scale_y_continuous(limits = c(-0.5, 1.2)) + theme_void(base_size = 15)
}

# ---- fig05: call probabilities and the filter threshold ----------------------------
ex <- fread(h19calls)
thr <- as.numeric(sub(".*setting threshold ([0-9.]+).*", "\\1",
                      grep("setting threshold", readLines(file.path(C, "mk3.log")), value = TRUE)[1]))
cat("H19 extract: calls", nrow(ex), " threshold", thr, " filtered (fail)", round(mean(ex$fail), 3), "\n")
ex[, call := factor(c(`-` = "C (unmodified)", m = "5mC", h = "5hmC")[call_code], levels = c("5mC", "C (unmodified)", "5hmC"))]
p05 <- ggplot(ex, aes(call_prob, fill = call)) + geom_histogram(binwidth = 0.01, boundary = 0) +
  geom_vline(xintercept = thr, colour = "#c0392b", linewidth = 1, linetype = 2) +
  annotate("text", x = thr - 0.01, y = Inf, vjust = 1.5, hjust = 1, colour = "#c0392b", size = 4.5,
           label = sprintf("modkit threshold %.2f\n(lowest-confidence 10%% dropped)", thr)) +
  scale_fill_manual(values = c(`5mC` = "#1b2a3a", `C (unmodified)` = "#a9b5c0", `5hmC` = "#e67e22"), name = "called as") +
  labs(x = "probability of the called state (from ML)", y = "CpG calls",
       title = "Per-read CpG calls at the H19 locus (HG001, nanopore)")
save(p05, "fig05-call-probs.png", w = 11, h = 4.4)

# ---- fig07: nanopore vs WGBS ----------------------------------------------------
a <- readbm(allbed)
mh <- merge(a[code == "m", .(chrom, start, nvalid, m = nmod)], a[code == "h", .(chrom, start, h = nmod)], by = c("chrom", "start"))
rdw <- function(f) { w <- fread(file.path(D, f), header = FALSE); w[, pos := ifelse(V6 == "+", V2, V2 - 1)]
  w[, .(cov = sum(V10), meth = sum(V10 * V11 / 100)), by = .(chrom = V1, pos)] }
wc <- merge(rdw("gm12878_wgbs_regions.bed.gz"), rdw("gm12878_wgbs_rep2_regions.bed.gz"), by = c("chrom", "pos"), suffixes = c("1", "2"))
wc[, `:=`(wcov = cov1 + cov2, wmeth = meth1 + meth2)]
j <- merge(mh, wc, by.x = c("chrom", "start"), by.y = c("chrom", "pos"))
j[, `:=`(ont = (m + h) / nvalid, ont5mc = m / nvalid, wgbs = wmeth / wcov, rep1 = meth1 / cov1, rep2 = meth2 / cov2)]
k <- j[nvalid >= 10 & cov1 >= 10 & cov2 >= 10]
cat("nanopore: CpGs", nrow(mh), " median coverage", median(mh$nvalid), " 5mC", round(sum(mh$m) / sum(mh$nvalid), 3),
    " 5hmC", round(sum(mh$h) / sum(mh$nvalid), 4), " (pooled over calls)\n")
cat("CpGs with >= 10 calls in nanopore and in each WGBS replicate:", nrow(k), "\n")
cat("  r(nanopore 5mC+5hmC, pooled WGBS) =", round(cor(k$ont, k$wgbs), 3), "  r(5mC only) =", round(cor(k$ont5mc, k$wgbs), 3),
    "  median |diff| =", round(median(abs(k$ont - k$wgbs)), 3), "\n")
cat("  r(WGBS rep1, rep2) =", round(cor(k$rep1, k$rep2), 3), "  r(nanopore, rep1) =", round(cor(k$ont, k$rep1), 3),
    "  r(nanopore, rep2) =", round(cor(k$ont, k$rep2), 3), "\n")
hexp <- function(d, x, y, xl, yl, tt) ggplot(d, aes({{ x }}, {{ y }})) + geom_bin2d(bins = 40) +
  scale_fill_gradient(low = "#dfe6ec", high = "#1b2a3a", trans = "log10", guide = "none") +
  geom_abline(colour = "#c0392b", linetype = 2) + coord_equal() + labs(x = xl, y = yl, title = tt)
p07a <- hexp(k, wgbs, ont, "GM12878 WGBS, reps pooled", "HG001 nanopore (5mC + 5hmC)", sprintf("nanopore vs WGBS: r = %.2f", cor(k$ont, k$wgbs)))
p07c <- hexp(k, rep1, rep2, "WGBS replicate 1", "WGBS replicate 2", sprintf("WGBS rep 1 vs rep 2: r = %.2f", cor(k$rep1, k$rep2)))
tr <- j[chrom == "chr21" & start >= 33530000 & start <= 33600000]
trl <- rbind(tr[nvalid >= 5, .(start, f = ont, assay = "nanopore (5mC+5hmC)")], tr[wcov >= 5, .(start, f = wgbs, assay = "WGBS")])
p07b <- ggplot(trl, aes(start, f, colour = assay)) + geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.04, se = FALSE, linewidth = 0.9) +
  scale_colour_manual(values = c(`nanopore (5mC+5hmC)` = "#2a6f97", WGBS = "#e67e22"), name = NULL) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6)) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "chr21 (Mb)", y = "methylation", title = "Along chr21 (GART, SON, DONSON)") + theme(legend.position = "bottom")
save(p07a + p07c + p07b + plot_layout(widths = c(1, 1, 1.5)), "fig07-ont-vs-wgbs.png", w = 14, h = 5)

# ---- phased pileups -------------------------------------------------------------
hp <- merge(readbm(file.path(C, "phased", "hp1.bedmethyl"))[code == "m", .(chrom, start, n1 = nvalid, f1 = nmod / nvalid)],
            readbm(file.path(C, "phased", "hp2.bedmethyl"))[code == "m", .(chrom, start, n2 = nvalid, f2 = nmod / nvalid)],
            by = c("chrom", "start"))
ht <- fread(file.path(C, "haplotags.tsv"))
cat("haplotagged reads:", sum(ht$haplotype != "none"), "of", nrow(ht), " (H1", sum(ht$haplotype == "H1"), ", H2", sum(ht$haplotype == "H2"), ")\n")
loci <- list(H19 = list("chr11", 1994000, 2006000, c("H19")),
             KvDMR = list("chr11", 2694000, 2706000, c("KCNQ1OT1", "KCNQ1")),
             SNRPN = list("chr15", 24948000, 24966000, c("SNRPN")))
icr <- list(H19 = c("chr11", 1998000, 2003000), KvDMR = c("chr11", 2699000, 2701000), SNRPN = c("chr15", 24954000, 24957000))
mk_locus <- function(nm) {
  L <- loci[[nm]]; x <- hp[chrom == L[[1]] & start >= L[[2]] & start <= L[[3]] & n1 >= 5 & n2 >= 5]
  x[, bin := floor(start / 400) * 400 + 200]
  xb <- rbind(x[, .(f = mean(f1), n = .N, hap = "paternal"), by = bin], x[, .(f = mean(f2), n = .N, hap = "maternal"), by = bin])[n >= 3]
  e <- as.data.table(gtf[seqnames(gtf) == L[[1]] & end(gtf) >= L[[2]] & start(gtf) <= L[[3]] & gtf$gene_name %in% L[[4]]])
  e[, canon := grepl("Ensembl_canonical", sapply(tag, paste, collapse = ","))]
  e <- e[, if (any(canon)) .SD[canon == TRUE] else .SD[transcript_id == transcript_id[1]], by = gene_name]
  e[, `:=`(start = pmax(start, L[[2]]), end = pmin(end, L[[3]]))]
  tx <- e[, .(start = min(start), end = max(end), strand = strand[1]), by = gene_name]
  g <- icr[[nm]]
  ggplot(xb, aes(bin, f, colour = hap)) +
    annotate("rect", xmin = as.numeric(g[2]), xmax = as.numeric(g[3]), ymin = -Inf, ymax = Inf, fill = "#2e8b57", alpha = 0.10) +
    geom_line(linewidth = 1.1) + geom_point(size = 1.4) +
    geom_segment(data = tx, aes(x = start, xend = end, y = -0.18, yend = -0.18), colour = "#1b2a3a", inherit.aes = FALSE) +
    geom_rect(data = e, aes(xmin = start, xmax = end, ymin = -0.23, ymax = -0.13), fill = "#1b2a3a", inherit.aes = FALSE) +
    geom_text(data = tx, aes(x = (start + end) / 2, y = -0.32, label = paste0(gene_name, ifelse(strand == "+", " >", " <"))),
              size = 4, inherit.aes = FALSE) +
    scale_colour_manual(values = colHP, name = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.38, 1.02)) +
    scale_x_continuous(limits = c(L[[2]], L[[3]]), labels = function(v) sprintf("%.3f", v / 1e6)) +
    labs(x = paste(L[[1]], "(Mb)"), y = "5mC (400-bp bins)", title = paste0(nm, ifelse(nm == "SNRPN", " (PWS-IC)", " ICR"))) +
    theme(legend.position = if (nm == "H19") "top" else "none")
}
p10 <- wrap_plots(lapply(names(loci), mk_locus), nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "top")
save(p10, "fig10-imprinted-loci.png", w = 14, h = 4.8)
# ICR summaries quoted in the slides (shaded in fig10)
for (n in names(icr)) { g <- icr[[n]]
  x <- hp[chrom == g[1] & start >= as.numeric(g[2]) & start <= as.numeric(g[3]) & n1 >= 5 & n2 >= 5]
  cat(n, g[1], g[2], g[3], ": CpGs", nrow(x), " paternal", round(mean(x$f1), 2), " maternal", round(mean(x$f2), 2), "\n") }

# ---- fig08: single molecules at the H19 ICR ---------------------------------------
from <- 1997000; to <- 2004000
r <- ex[fail == FALSE & ref_position >= from & ref_position <= to & alignment_start <= 1999000 & alignment_end >= 2001500]
r <- merge(r, ht[, .(read_id = `#readname`, haplotype)], by = "read_id")[haplotype %in% c("H1", "H2")]
r[, hap := ifelse(haplotype == "H1", "paternal", "maternal")]
rd <- unique(r[, .(read_id, hap)])[order(hap, read_id)]
rd <- rbind(rd[hap == "paternal"][seq_len(min(.N, 20))], rd[hap == "maternal"][seq_len(min(.N, 20))])
rd[, y := seq_len(.N) + ifelse(hap == "maternal", 2, 0)]
r <- merge(r[read_id %in% rd$read_id], rd, by = c("read_id", "hap"))
r[, state := ifelse(call_code %in% c("m", "h"), "methylated", "unmethylated")]
cat("fig08: reads shown", nrow(rd), " fraction methylated: paternal", round(r[hap == "paternal", mean(state == "methylated")], 2),
    " maternal", round(r[hap == "maternal", mean(state == "methylated")], 2), "\n")
lab8 <- rd[, .(y = mean(y)), by = hap]
h19 <- as.data.table(gtf[seqnames(gtf) == "chr11" & gtf$transcript_name == "H19-213"])[, .(start = pmax(start, from - 300), end)]
ybot <- max(rd$y) + 2
p08 <- ggplot(r, aes(ref_position, y)) +
  annotate("rect", xmin = 1998000, xmax = 2003000, ymin = -Inf, ymax = Inf, fill = "#2e8b57", alpha = 0.08) +
  annotate("text", x = 2000500, y = 0, label = "imprinting control region", colour = "#2e8b57", size = 4.5) +
  geom_segment(data = r[, .(x0 = min(ref_position), x1 = max(ref_position)), by = .(y, hap)],
               aes(x = x0, xend = x1, y = y, yend = y, colour = hap), linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(aes(fill = state), shape = 21, size = 1.7, stroke = 0.25, colour = "#1b2a3a") +
  geom_rect(data = h19, aes(xmin = start, xmax = end, ymin = ybot - 0.4, ymax = ybot + 0.4), fill = "#1b2a3a", inherit.aes = FALSE) +
  annotate("text", x = 1997400, y = ybot + 1.3, label = "H19 <", size = 4.3) +
  annotate("text", x = from - 200, y = lab8$y, label = lab8$hap, angle = 90, colour = colHP[lab8$hap], size = 5) +
  scale_fill_manual(values = c(methylated = "#1b2a3a", unmethylated = "#ffffff"), name = NULL) +
  scale_colour_manual(values = colHP, guide = "none") +
  scale_x_continuous(limits = c(from - 300, to), labels = function(v) sprintf("%.3f", v / 1e6)) + scale_y_reverse() +
  labs(x = "chr11 (Mb)", y = NULL, title = "Single molecules across the H19 imprinting control region (HG001)") +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(), legend.position = "top")
save(p08, "fig08-h19-reads.png", w = 12, h = 5.6)

# ---- fig11: chrX vs chr20 promoter CpG islands --------------------------------------
rdh <- function(f) setnames(fread(f), c("chrom", "start", "end", "n", "m"))
h1 <- rdh(file.path(D, "hg001_hp1_chrX_chr20_5mc.bed.gz")); h2 <- rdh(file.path(D, "hg001_hp2_chrX_chr20_5mc.bed.gz"))
cgi <- fread(file.path(D, "hg38_cpgIslandExt.txt.gz"))[V2 %in% c("chrX", "chr20")]
cg <- GRanges(cgi$V2, IRanges(cgi$V3 + 1, cgi$V4))
gg <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"), feature.type = "gene")
gg <- gg[gg$gene_type == "protein_coding" & seqnames(gg) %in% c("chrX", "chr20")]
tss <- promoters(gg, upstream = 1000, downstream = 1000)
pcg <- cg[overlapsAny(cg, tss)]
o <- findOverlaps(pcg, tss); pcg$gene <- NA_character_; pcg$gene[queryHits(o)] <- tss$gene_name[subjectHits(o)]
agg <- function(h) { g <- GRanges(h$chrom, IRanges(h$start + 1, width = 1)); ov <- findOverlaps(g, pcg)
  data.table(i = subjectHits(ov), n = h$n[queryHits(ov)], m = h$m[queryHits(ov)])[n >= 5, .(f = sum(m) / sum(n), k = .N), by = i] }
xc <- merge(agg(h1), agg(h2), by = "i", suffixes = c("1", "2"))[k1 >= 5 & k2 >= 5]
xc[, `:=`(chrom = as.character(seqnames(pcg))[i], gene = pcg$gene[i])]
xc[, `:=`(hi = pmax(f1, f2), lo = pmin(f1, f2))]
sm <- xc[, .(n = .N, one_methylated = mean(hi - lo > 0.5), both_low = mean(hi < 0.2)), by = chrom]
print(sm)
cat("chr20 allelic CGIs:", paste(xc[chrom == "chr20" & hi - lo > 0.5, gene], collapse = ", "), "\n")
esc <- xc[chrom == "chrX" & hi < 0.2, gene]
cat("chrX CGIs unmethylated on both haplotypes:", length(esc), ":", paste(esc, collapse = ", "), "\n")
lab <- xc[chrom == "chr20" & hi - lo > 0.5][!duplicated(gene)]
xc[, chrom := factor(chrom, levels = c("chr20", "chrX"))]; lab[, chrom := factor(chrom, levels = c("chr20", "chrX"))]
escl <- data.table(chrom = factor("chrX", levels = c("chr20", "chrX")), x = 0.06, y = 0.03,
                   label = "both unmethylated (escape XCI):\nKDM5C, KDM6A, DDX3X, ZFX,\nSTS, PAR1 genes ...")
p11 <- ggplot(xc, aes(lo, hi)) + geom_point(aes(colour = chrom), alpha = 0.6, size = 1.8) +
  geom_abline(linetype = 2, colour = "#5b6b7a") +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3.8, min.segment.length = 0, nudge_x = 0.15) +
  geom_text(data = escl, aes(x, y, label = label), hjust = 0, vjust = 0, size = 3.6, lineheight = 0.9) +
  facet_wrap(~ chrom, labeller = labeller(chrom = c(chr20 = sprintf("chr20 (%d promoter CGIs)", sm[chrom == "chr20", n]),
                                                    chrX = sprintf("chrX (%d promoter CGIs)", sm[chrom == "chrX", n])))) +
  scale_colour_manual(values = c(chr20 = "#5b6b7a", chrX = "#a8324a"), guide = "none") + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "5mC on the less methylated haplotype", y = "5mC on the more methylated haplotype",
       title = "Promoter CpG islands, HG001 (female), haplotype-resolved")
save(p11, "fig11-chrx-cgis.png", w = 11, h = 5.4)

# ---- fig12 / fig13: m6A accessibility (HG002, EcoGII) ------------------------------
rdA <- function(f) { x <- fread(f, header = FALSE, select = c(1, 2, 10, 12)); setnames(x, c("chrom", "pos", "nvalid", "nmod")); x }
e6 <- rdA(file.path(C, "m6a_ecogii.bed")); n6 <- rdA(file.path(C, "m6a_native.bed"))
atac <- fread(file.path(D, "gm12878_atac_peaks.bed.gz"))[V1 == "chr21" & V2 >= 33e6 & V3 <= 34e6]
pk <- reduce(GRanges("chr21", IRanges(atac$V2 + 1, atac$V3)))
inp <- function(x) overlapsAny(GRanges(x$chrom, IRanges(x$pos + 1, width = 1)), pk)
fr <- function(x, i) sum(x$nmod[i]) / sum(x$nvalid[i])
ie <- inp(e6); inn <- inp(n6)
cat("m6A fraction: EcoGII in ATAC peaks", round(fr(e6, ie), 4), " outside", round(fr(e6, !ie), 4),
    " | native in peaks", round(fr(n6, inn), 4), " outside", round(fr(n6, !inn), 4), " | peaks:", length(pk), "\n")
prof <- function(x, centers, wdw = 1000, b = 25) {
  r <- rbindlist(lapply(centers, function(cc) x[pos >= cc - wdw & pos <= cc + wdw, .(d = pos - cc, nvalid, nmod)]))
  r[, bin := floor(d / b) * b + b / 2][, .(f = sum(nmod) / sum(nvalid)), by = bin][order(bin)]
}
summ <- atac$V2 + atac$V10
pr <- rbind(prof(e6, summ)[, sample := "EcoGII-treated"], prof(n6, summ)[, sample := "untreated control"])
p13a <- ggplot(pr, aes(bin, 100 * f, colour = sample)) + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c(`EcoGII-treated` = "#c0392b", `untreated control` = "#5b6b7a"), name = NULL) +
  labs(x = "distance to GM12878 ATAC summit (bp)", y = "% of A called m6A",
       title = sprintf("%d ATAC peak summits on chr21 (33-34 Mb)", nrow(atac))) + theme(legend.position = "top")
bars <- data.table(sample = rep(c("EcoGII-treated", "untreated control"), each = 2), where = rep(c("in ATAC peaks", "outside"), 2),
                   f = c(fr(e6, ie), fr(e6, !ie), fr(n6, inn), fr(n6, !inn)))
p13b <- ggplot(bars, aes(where, 100 * f, fill = sample)) + geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", 100 * f)), position = position_dodge(0.7), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c(`EcoGII-treated` = "#c0392b", `untreated control` = "#5b6b7a"), guide = "none") +
  labs(x = NULL, y = "% of A called m6A", title = "Peaks vs background")
save(p13a + p13b + plot_layout(widths = c(1.6, 1)), "fig13-m6a-peaks.png", w = 13, h = 4.8)

am <- fread(atpcalls)
from <- 33913000; to <- 33919000
am <- am[fail == FALSE & ref_position >= from & ref_position <= to & alignment_start <= from + 500 & alignment_end >= to - 500]
reads <- unique(am$read_id); am[, y := match(read_id, reads)]
cat("fig12: reads spanning ATP5PO window", length(reads), "\n")
nd <- am[canonical_base == "A" & ref_position >= 33915500 & ref_position <= 33916100]
fl <- am[canonical_base == "A" & (ref_position < 33914500 | ref_position > 33917500)]
cat("ATP5PO: m6A fraction in promoter window", round(mean(nd$call_code == "a"), 3), " in flanks", round(mean(fl$call_code == "a"), 3), "\n")
cg5 <- am[canonical_base == "C" & call_code %in% c("m", "-")]
cat("ATP5PO CpG 5mC: in CGI window", round(cg5[ref_position >= 33914800 & ref_position <= 33916500, mean(call_code == "m")], 3),
    " flanks", round(cg5[ref_position < 33914500 | ref_position > 33917500, mean(call_code == "m")], 3), "\n")
pk12 <- as.data.table(reduce(GRanges("chr21", IRanges(atac$V2 + 1, atac$V3))))[end >= from & start <= to][, .(V2 = start, V3 = end)]
atp <- as.data.table(gtf[seqnames(gtf) == "chr21" & gtf$transcript_name == "ATP5PO-201" & end(gtf) >= from])[, .(start = pmax(start, from), end)]
yb <- length(reads) + 2
p12 <- ggplot() +
  geom_rect(data = pk12, aes(xmin = V2, xmax = V3, ymin = -Inf, ymax = Inf), fill = "#2e8b57", alpha = 0.12) +
  geom_segment(data = am[, .(x0 = min(ref_position), x1 = max(ref_position)), by = y], aes(x = x0, xend = x1, y = y, yend = y),
               colour = "#dfe6ec", linewidth = 0.8) +
  geom_point(data = am[canonical_base == "C" & call_code == "m"], aes(ref_position, y), colour = "#1b2a3a", size = 0.7) +
  geom_point(data = am[canonical_base == "A" & call_code == "a"], aes(ref_position, y), colour = "#c0392b", size = 0.8) +
  annotate("text", x = mean(c(pk12$V2[1], pk12$V3[1])), y = 0, label = "GM12878 ATAC peak", colour = "#2e8b57", size = 4.5, vjust = 0.2) +
  annotate("segment", x = from, xend = max(atp$end), y = yb, yend = yb, colour = "#1b2a3a") +
  geom_rect(data = atp, aes(xmin = start, xmax = end, ymin = yb - 0.5, ymax = yb + 0.5), fill = "#1b2a3a", inherit.aes = FALSE) +
  annotate("text", x = max(atp$end) - 900, y = yb + 1.6, label = "ATP5PO <", size = 4.3) +
  scale_x_continuous(limits = c(from, to), labels = function(v) sprintf("%.3f", v / 1e6)) + scale_y_reverse() +
  labs(x = "chr21 (Mb)", y = "one read per row", title = "HG002 chromatin fibres at the ATP5PO promoter: m6A (red) and 5mC at CpGs (black)") +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
save(p12, "fig12-atp5po-fibres.png", w = 12, h = 5.8)
cat("done\n")
