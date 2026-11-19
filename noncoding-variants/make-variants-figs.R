# Real-data figures for noncoding-variants/lec.qmd ("Epigenomics and genetic variation").
# Run from this folder:  Rscript make-variants-figs.R   (about 4 minutes, 4 GB RAM)
# Needs (data/download.sh noncoding-variants; data/prepare-integ.sh for the GWAS and allele files):
#   gwas_catalog_gws.tsv.gz, gwas_autoimmune_gm12878_alleles.tsv, gm12878_atac_peaks.bed.gz,
#   k562_atac_peaks.bed.gz, hepg2_atac_peaks.bed.gz, hg38.chrom.sizes, hg38.fa.bgz (+ .fai, .gzi),
#   corces2016_atac_counts.txt.gz, hg19ToHg38.over.chain.gz, hg38_chr21_chrM.fa,
#   na12878_chr21_het_snps.vcf.gz
# and, for fig18, the filtered GM12878 ATAC BAMs rep1.dedup.bam and rep2.dedup.bam made in the
# exercise of *ATAC-seq processing and peak calling* (set BAMDIR below), plus bcftools on the PATH.
# R packages: data.table, GenomicRanges, rtracklayer, Rsamtools, Biostrings, JASPAR2024, TFBSTools,
# motifmatchr, ggplot2, patchwork, ggrepel.
#
# Writes fig02, fig05, fig07, fig08, fig09, fig14, fig15, fig18 (*.png) and prints the numbers
# quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer); library(Rsamtools); library(Biostrings)
  library(JASPAR2024); library(TFBSTools); library(motifmatchr)
  library(ggplot2); library(patchwork); library(ggrepel)
})
set.seed(1)
D <- "../data"
BAMDIR <- Sys.getenv("BAMDIR", ".")          # folder with rep1.dedup.bam, rep2.dedup.bam
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
chroms <- paste0("chr", c(1:22, "X"))

# ---- the GWAS Catalog ----------------------------------------------------------------
gw <- fread(file.path(D, "gwas_catalog_gws.tsv.gz"))
cat("associations", nrow(gw), " unique SNPs", uniqueN(gw$snp), " traits", uniqueN(gw$trait), " studies", uniqueN(gw$study), "\n")
snp <- gw[, .(chr = chr[1], pos = pos[1], context = context[1]), by = snp]
coding <- c("missense_variant", "synonymous_variant", "stop_gained", "frameshift_variant", "inframe_insertion",
            "inframe_deletion", "start_lost", "stop_lost", "stop_retained_variant", "protein_altering_variant",
            "coding_sequence_variant", "splice_donor_variant", "splice_acceptor_variant")
snp[, class := fcase(context %in% coding, "coding or splice site",
                     context %in% c("3_prime_UTR_variant", "5_prime_UTR_variant"), "UTR",
                     context == "intron_variant", "intron",
                     context == "intergenic_variant", "intergenic", default = "other non-coding")]
tab <- snp[, .N, by = class][, frac := N / sum(N)][order(-N)]
print(tab)
tab[, class := factor(class, levels = rev(class))]
p02 <- ggplot(tab, aes(class, frac)) + geom_col(fill = ifelse(tab$class == "coding or splice site", "#c0392b", "#2a6f97")) +
  geom_text(aes(label = sprintf("%.1f%%  (%s)", 100 * frac, format(N, big.mark = " "))), hjust = -0.05, size = 4.5) +
  coord_flip() + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.35))) +
  labs(x = NULL, y = "share of unique SNPs", title = sprintf("%s genome-wide significant SNPs in the GWAS Catalog",
       format(nrow(snp), big.mark = " ")))
save(p02, "fig02-context.png", w = 11, h = 4)

# ---- trait sets and independent loci ------------------------------------------------------
sets <- list(RA = "rheumatoid arthritis", SLE = "systemic lupus erythematosus", MS = "multiple sclerosis",
             IBD = c("inflammatory bowel disease", "Crohn disease", "ulcerative colitis"),
             LDL = "low density lipoprotein cholesterol measurement", schizophrenia = "schizophrenia",
             height = "body height")
sets$autoimmune <- unlist(sets[c("RA", "SLE", "MS", "IBD")])
prune <- function(x, w = 250000) {               # keep the strongest SNP, drop SNPs within w of a kept one
  x <- x[order(-mlog10p)]; keep <- logical(nrow(x))
  for (ch in unique(x$chr)) { i <- which(x$chr == ch); kept <- integer(0)
    for (j in i) if (!length(kept) || all(abs(x$pos[kept] - x$pos[j]) > w)) kept <- c(kept, j)
    keep[kept] <- TRUE }
  x[keep] }
leads <- rbindlist(lapply(names(sets), function(s) {
  x <- gw[trait %in% sets[[s]], .(chr = chr[1], pos = pos[1], mlog10p = max(mlog10p)), by = snp]
  y <- prune(x); y[, set := s]; cat(s, ": SNPs", nrow(x), " loci", nrow(y), "\n"); y }))
bg <- snp[, .(snp, chr, pos)]

# ---- annotations: 50 000 strongest peaks per cell line, 500 bp around the summit ---------
top_peaks <- function(f, n = 50000) { p <- fread(file.path(D, f))[V1 %in% chroms][order(-V7)][1:n]
  GRanges(p$V1, IRanges(p$V2 + p$V10 + 1 - 250, width = 500)) }
atac <- list(GM12878 = top_peaks("gm12878_atac_peaks.bed.gz"), K562 = top_peaks("k562_atac_peaks.bed.gz"),
             HepG2 = top_peaks("hepg2_atac_peaks.bed.gz"))
spec <- list(`GM12878 only` = atac$GM12878[!overlapsAny(atac$GM12878, c(atac$K562, atac$HepG2))],
             `K562 only`    = atac$K562[!overlapsAny(atac$K562, c(atac$GM12878, atac$HepG2))],
             `HepG2 only`   = atac$HepG2[!overlapsAny(atac$HepG2, c(atac$GM12878, atac$K562))],
             `all three`    = atac$GM12878[overlapsAny(atac$GM12878, atac$K562) & overlapsAny(atac$GM12878, atac$HepG2)])
print(sapply(spec, length))
near <- function(chr, pos, a, w = 1000) overlapsAny(GRanges(chr, IRanges(pmax(pos - w, 1), pos + w)), a)
enrich <- function(a, s, w = 1000) { l <- leads[set == s]; k <- sum(near(l$chr, l$pos, a, w)); b <- mean(near(bg$chr, bg$pos, a, w))
  data.table(set = s, n = nrow(l), k = k, frac = k / nrow(l), bg = b, enr = k / nrow(l) / b,
             p = binom.test(k, nrow(l), b, alternative = "greater")$p.value) }

# fig05: the background matters (GM12878 ATAC, +-1 kb)
cs <- fread(file.path(D, "hg38.chrom.sizes"))[V1 %in% chroms]
rp <- data.table(chr = sample(cs$V1, 1e5, TRUE, prob = cs$V2)); rp[, pos := floor(runif(.N, 1e6, cs$V2[match(chr, cs$V1)] - 1e6))]
f05 <- rbind(data.table(group = "random positions", frac = mean(near(rp$chr, rp$pos, atac$GM12878)), n = nrow(rp)),
             data.table(group = "all GWAS SNPs", frac = mean(near(bg$chr, bg$pos, atac$GM12878)), n = nrow(bg)),
             leads[set %in% c("autoimmune", "LDL", "schizophrenia", "height"),
                   .(frac = mean(near(chr, pos, atac$GM12878)), n = .N), by = .(group = paste(set, "loci"))])
print(f05)
f05[, group := factor(group, levels = c("random positions", "all GWAS SNPs", "autoimmune loci", "LDL loci", "schizophrenia loci", "height loci"))]
p05 <- ggplot(f05, aes(group, frac, fill = group)) + geom_col() +
  scale_fill_manual(values = c("#5b6b7a", "#5b6b7a", "#c0392b", "#2a6f97", "#2a6f97", "#2a6f97"), guide = "none") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * frac)), vjust = -0.3, size = 4.5) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "within 1 kb of a peak", title = "GM12878 ATAC (50 000 strongest peaks)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save(p05, "fig05-background.png", w = 10, h = 4.6)

# fig07: cell-type-specific peaks
f07 <- rbindlist(lapply(names(spec), function(a) rbindlist(lapply(c("autoimmune", "LDL", "schizophrenia", "height"),
  function(s) cbind(annot = a, enrich(spec[[a]], s))))))
print(f07[, .(annot, set, k, n, frac = round(frac, 3), bg = round(bg, 3), enr = round(enr, 2), p = signif(p, 2))])
f07[, annot := factor(annot, levels = names(spec))]; f07[, set := factor(set, levels = c("autoimmune", "LDL", "schizophrenia", "height"))]
p07 <- ggplot(f07, aes(annot, enr, fill = set)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(aes(label = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))),
            position = position_dodge(0.85), vjust = -0.2, size = 5) +
  scale_fill_manual(values = c(autoimmune = "#c0392b", LDL = "#d4a017", schizophrenia = "#6b4fa0", height = "#5b6b7a"), name = "GWAS loci") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "ATAC peaks found only in ...", y = "enrichment vs all GWAS SNPs", title = "Disease loci near cell-type-specific peaks (±1 kb)")
save(p07, "fig07-specific.png", w = 12, h = 4.8)

# fig08: 13 blood cell types (Corces et al. 2016), the 25 000 most cell-type-specific peaks each
cc <- fread(file.path(D, "corces2016_atac_counts.txt.gz"))
code2type <- c(`1` = "HSC", `2` = "MPP", `3` = "LMPP", `4` = "CMP", `5` = "GMP", `6` = "MEP", `7` = "Mono",
               `9` = "CD4", `10` = "CD8", `11` = "NK", `13` = "B", `14` = "CLP", `15` = "Ery")
healthy <- grep("^(Donor)?[0-9]+-", names(cc), value = TRUE)
re <- "^(?:Donor)?([0-9]+)-(?:[A-Za-z][A-Za-z0-9]*-)?([0-9]+)([A-C]?)(?:-[0-9]+)?$"
ct <- code2type[sapply(regmatches(healthy, regexec(re, healthy, perl = TRUE)), `[`, 3)]
m  <- as.matrix(cc[, ..healthy])
lc <- sapply(split(seq_along(ct), ct), function(i) log2(rowMeans(sweep(m[, i, drop = FALSE], 2, colSums(m[, i, drop = FALSE]), "/")) * 1e6 + 1))
z  <- lc - rowMeans(lc)
chainf <- tempfile(fileext = ".chain"); writeLines(readLines(gzfile(file.path(D, "hg19ToHg38.over.chain.gz"))), chainf)
lo <- liftOver(GRanges(cc$Chr, IRanges(cc$Start + 1, cc$End)), import.chain(chainf))
one <- which(lengths(lo) == 1); p38 <- unlist(lo[one])
f08 <- rbindlist(lapply(colnames(z), function(k) { A <- p38[order(-z[one, k])[1:25000]]
  rbindlist(lapply(c("RA", "SLE", "MS", "IBD", "LDL", "schizophrenia", "height"), function(s) cbind(celltype = k, enrich(A, s)))) }))
lv <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Ery", "Mono", "CLP", "B", "CD4", "CD8", "NK")
f08[, celltype := factor(celltype, levels = lv)]
f08[, set := factor(set, levels = rev(c("RA", "SLE", "MS", "IBD", "LDL", "schizophrenia", "height")))]
print(dcast(f08, set ~ celltype, value.var = "enr", fun.aggregate = function(x) round(x, 1)))
p08 <- ggplot(f08, aes(celltype, set, fill = log2(enr))) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f%s", enr, ifelse(p < 0.05 / nrow(f08), "*", ""))), size = 4) +
  scale_fill_gradient2(low = "#2a6f97", mid = "white", high = "#c0392b", midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                       name = "log2\nenrichment") +
  labs(x = NULL, y = NULL, title = "GWAS loci near the 25 000 most specific peaks of each blood cell type",
       subtitle = "enrichment vs all GWAS SNPs, ±1 kb; * = binomial p < 0.05 / 91")
save(p08, "fig08-corces-heatmap.png", w = 12, h = 4.8)

# fig09: window size
f09 <- rbindlist(lapply(c(0, 500, 1000, 2000, 5000, 10000, 20000), function(w)
  rbind(cbind(w = w, annot = "GM12878 only", enrich(spec[["GM12878 only"]], "autoimmune", w)),
        cbind(w = w, annot = "HepG2 only", enrich(spec[["HepG2 only"]], "LDL", w)))))
print(f09[, .(w, annot, set, k, frac = round(frac, 3), bg = round(bg, 3), enr = round(enr, 2), p = signif(p, 2))])
f09[, lab := paste(set, "loci,", annot, "peaks")]
p09a <- ggplot(f09, aes(pmax(w, 100) / 1000, enr, colour = lab)) + geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10, 20), labels = c("0", "0.5", "1", "2", "5", "10", "20")) +
  scale_colour_manual(values = c("#c0392b", "#d4a017"), name = NULL) + geom_hline(yintercept = 1, linetype = 2) +
  labs(x = "window around the lead SNP (± kb)", y = "enrichment", title = "Wider windows dilute the enrichment") +
  theme(legend.position = "bottom")
p09b <- ggplot(f09, aes(pmax(w, 100) / 1000, frac, colour = lab)) + geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
  geom_line(aes(y = bg), linetype = 2) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10, 20), labels = c("0", "0.5", "1", "2", "5", "10", "20")) +
  scale_y_continuous(labels = scales::percent) + scale_colour_manual(values = c("#c0392b", "#d4a017"), guide = "none") +
  labs(x = "window around the lead SNP (± kb)", y = "loci with a peak (dashed: background)", title = "... but catch more loci")
save(p09a + p09b, "fig09-windows.png", w = 12, h = 4.8)

# ---- motifs: alleles of autoimmune SNPs in GM12878 ATAC peaks -------------------------------
fa <- FaFile(file.path(D, "hg38.fa.bgz"))
jdb <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024()))
pfms <- getMatrixSet(jdb, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
names(pfms) <- paste(name(pfms), ID(pfms), sep = "_")
wd <- sapply(pfms, function(p) ncol(Matrix(p))); L <- max(wd)
al <- fread(file.path(D, "gwas_autoimmune_gm12878_alleles.tsv"))
risk <- gw[trait %in% sets$autoimmune | snp == "rs12740374", .(risk = paste(unique(risk_allele), collapse = ",")), by = snp]
al <- merge(al, risk, by = "snp", all.x = TRUE)
pick_alt <- function(alt, risk) { a <- strsplit(alt, ",")[[1]]; if (length(a) == 1) return(a)
  x <- intersect(strsplit(risk, ",")[[1]], a); if (length(x) == 1) x else NA_character_ }
al[, alt1 := mapply(pick_alt, alt, risk)]
cat("SNPs with alleles", nrow(al), "; with one alternative allele (or resolved by the risk allele)", sum(!is.na(al$alt1)), "\n")
al <- al[!is.na(alt1)]
allele_seqs <- function(chr, pos, ref, alt) {
  s <- getSeq(fa, GRanges(chr, IRanges(pos - L + 1, pos + L - 1)))
  stopifnot(all(as.character(subseq(s, L, L)) == ref))
  a <- s; subseq(a, L, L) <- DNAStringSet(alt); list(ref = s, alt = a) }
motif_change <- function(sq) {
  out <- lapply(split(seq_along(pfms), wd), function(i) { k <- wd[i[1]]
    r <- motifMatches(matchMotifs(pfms[i], subseq(sq$ref, L - k + 1, L + k - 1), bg = "even", p.cutoff = 1e-4))
    a <- motifMatches(matchMotifs(pfms[i], subseq(sq$alt, L - k + 1, L + k - 1), bg = "even", p.cutoff = 1e-4))
    list(lost = as.matrix(r & !a), gained = as.matrix(!r & a)) })
  list(lost = do.call(cbind, lapply(out, `[[`, "lost")), gained = do.call(cbind, lapply(out, `[[`, "gained"))) }
mc <- motif_change(allele_seqs(al$chr, al$pos, al$ref, al$alt1))
p <- fread(file.path(D, "gm12878_atac_peaks.bed.gz"))[V1 %in% paste0("chr", 1:22)][sample(.N, 2000)]
cpos <- p$V2 + p$V10 + sample(-100:100, nrow(p), TRUE)
cref <- as.character(getSeq(fa, GRanges(p$V1, IRanges(cpos, width = 1))))
ok <- cref %in% c("A", "C", "G", "T")
calt <- sapply(cref[ok], function(b) sample(setdiff(c("A", "C", "G", "T"), b), 1))
mcc <- motif_change(allele_seqs(p$V1[ok], cpos[ok], cref[ok], calt))
blood <- c("SPI1", "SPIB", "ETS1", "ELF1", "IRF4", "IRF8", "EBF1", "PAX5", "POU2F2", "RELA", "REL", "NFKB1", "NFKB2",
           "RUNX1", "RUNX3", "MEF2C", "BATF", "CTCF", "STAT1", "STAT3", "STAT4", "TBX21", "GATA3", "IKZF1", "BCL6")
isb <- sub("_MA.*", "", colnames(mc$lost)) %in% blood
cat("blood TF motifs:", sum(isb), "of", length(isb), "\n")
nch  <- function(x, cols) rowSums(x$lost[, cols, drop = FALSE] | x$gained[, cols, drop = FALSE])
f15 <- rbind(data.table(set = "autoimmune GWAS SNPs", motifs = "all 879 JASPAR", n = nch(mc, TRUE)),
             data.table(set = "random sites in the same peaks", motifs = "all 879 JASPAR", n = nch(mcc, TRUE)),
             data.table(set = "autoimmune GWAS SNPs", motifs = sprintf("%d blood TF motifs", sum(isb)), n = nch(mc, isb)),
             data.table(set = "random sites in the same peaks", motifs = sprintf("%d blood TF motifs", sum(isb)), n = nch(mcc, isb)))
print(f15[, .(sites = .N, any_change = round(mean(n > 0), 3), mean_changes = round(mean(n), 2)), by = .(set, motifs)])
sm <- f15[, .(frac = mean(n > 0), k = sum(n > 0), N = .N), by = .(set, motifs)]
for (mo in unique(sm$motifs)) { x <- sm[motifs == mo]
  cat(mo, ": Fisher p =", signif(fisher.test(matrix(c(x$k, x$N - x$k), 2))$p.value, 2), "\n") }
p15 <- ggplot(sm, aes(motifs, frac, fill = set)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), position = position_dodge(0.8), vjust = -0.3, size = 4.5) +
  scale_fill_manual(values = c("#c0392b", "#5b6b7a"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "sites that gain or lose ≥ 1 match", title = "Most variants change some motif match") +
  theme(legend.position = "top")
hits <- rbindlist(lapply(which(isb), function(j) data.table(snp = al$snp, motif = colnames(mc$lost)[j],
  change = ifelse(mc$lost[, j], "lost", ifelse(mc$gained[, j], "gained", NA)))))[!is.na(change)]
hits[, tf := sub("_MA.*", "", motif)]
print(hits[order(tf)])
p15b <- ggplot(hits[, .N, by = .(tf, change)], aes(reorder(tf, N, sum), N, fill = change)) + geom_col() + coord_flip() +
  scale_fill_manual(values = c(gained = "#2e8b57", lost = "#a8324a"), name = NULL) +
  labs(x = NULL, y = "autoimmune SNPs", title = "Blood TF motifs changed") + theme(legend.position = "top")
save(p15 + p15b + plot_layout(widths = c(1.4, 1)), "fig15-motif-scan.png", w = 12, h = 5)

# fig14: SORT1 rs12740374 and C/EBP
i <- which(al$snp == "rs12740374")
cat("rs12740374:", al$chr[i], al$pos[i], al$ref[i], ">", al$alt1[i], "; gained:", colnames(mc$gained)[mc$gained[i, ]],
    "; lost:", colnames(mc$lost)[mc$lost[i, ]], "\n")
cebpa <- pfms[["CEBPA_MA0102.5"]]
pwm <- as.matrix(toPWM(cebpa, type = "log2probratio", pseudocounts = 0.8)); w <- ncol(pwm)
s21 <- getSeq(fa, GRanges(al$chr[i], IRanges(al$pos[i] - 15, al$pos[i] + 15)))[[1]]
a21 <- s21; a21[16] <- DNAString(al$alt1[i])
win <- function(s) { st <- (16 - w + 1):16
  f <- PWMscoreStartingAt(pwm, s, st); r <- PWMscoreStartingAt(pwm, reverseComplement(s), 31 - st - w + 2)
  data.table(start = c(st, st), strand = rep(c("+", "-"), each = length(st)), score = c(f, r)) }
wr <- win(s21); wa <- win(a21)
br <- wr[which.max(score)]; ba <- wa[which.max(score)]
cat("CEBPA best window ref:", br$score, "alt:", ba$score, " (max possible", sum(apply(pwm, 2, max)), ")\n")
cat("ref  ", as.character(s21), "\nalt  ", as.character(a21), "\n")
cat("best alt window:", as.character(subseq(a21, ba$start, ba$start + w - 1)), "strand", ba$strand, "\n")
cebp <- grep("^CEBP[ABDG]_", names(pfms), value = TRUE)
d14 <- rbindlist(lapply(cebp, function(n) { pwm <<- as.matrix(toPWM(pfms[[n]], type = "log2probratio", pseudocounts = 0.8)); w <<- ncol(pwm)
  data.table(motif = sub("_", " ", n), allele = c(paste0("reference ", al$ref[i]), paste0("alternative ", al$alt1[i])),
             score = c(max(win(s21)$score), max(win(a21)$score)), best = sum(apply(pwm, 2, max))) }))
print(d14)
d14[, allele := factor(allele, levels = unique(allele))]
p14 <- ggplot(d14, aes(motif, score, fill = allele)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_point(aes(y = best), shape = 95, size = 14, colour = "#1b2a3a", show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey50") +
  scale_fill_manual(values = c("#5b6b7a", "#2e8b57"), name = NULL) +
  labs(x = NULL, y = "best score (bits)",
       title = "rs12740374 (SORT1 locus): the T allele creates a C/EBP site",
       subtitle = paste0("reference ", substr(as.character(s21), 9, 23), "   alternative ", substr(as.character(a21), 9, 23),
                         "   (black bar = best possible score)")) + theme(legend.position = "top")
save(p14, "fig14-sort1.png", w = 11, h = 4.4)

# ---- fig18: allele-specific accessibility at heterozygous SNPs, GM12878 chr21 ----------------
b1 <- file.path(BAMDIR, "rep1.dedup.bam"); b2 <- file.path(BAMDIR, "rep2.dedup.bam")
if (file.exists(b1) && file.exists(b2)) {
  pk <- fread(file.path(D, "gm12878_atac_peaks.bed.gz"))[V1 == "chr21"]
  pk <- reduce(GRanges(pk$V1, IRanges(pk$V2 + 1, pk$V3)))
  bed <- tempfile(fileext = ".bed"); fwrite(data.table(as.character(seqnames(pk)), start(pk) - 1, end(pk)), bed, sep = "\t", col.names = FALSE)
  cmd <- sprintf("bcftools mpileup -f %s -R %s -T %s -a AD -q 30 -Q 20 -d 1000 %s %s 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t[%%AD\\t]\\n'",
                 file.path(D, "hg38_chr21_chrM.fa"), bed, file.path(D, "na12878_chr21_het_snps.vcf.gz"), b1, b2)
  ad <- fread(cmd = cmd, header = FALSE, sep = "\t", select = 1:6, col.names = c("chr", "pos", "ref", "alts", "ad1", "ad2"))
  het <- fread(cmd = sprintf("bcftools query -f '%%POS\\t%%ALT\\n' %s", file.path(D, "na12878_chr21_het_snps.vcf.gz")), col.names = c("pos", "het_alt"))
  ad <- merge(ad, het, by = "pos")
  cnt <- function(adv, alts, a) mapply(function(x, al, h) { v <- as.integer(strsplit(x, ",")[[1]]); k <- match(h, strsplit(al, ",")[[1]])
    if (is.na(k)) 0L else v[k + 1] }, adv, alts, a)
  ad[, refc := as.integer(sub(",.*", "", ad1)) + as.integer(sub(",.*", "", ad2))]
  ad[, altc := cnt(ad1, alts, het_alt) + cnt(ad2, alts, het_alt)]
  ad[, n := refc + altc]; ad <- ad[n >= 10]
  ad[, p := mapply(function(k, n) binom.test(k, n, 0.5)$p.value, refc, n)]
  ad[, padj := p.adjust(p, "BH")]
  cat("het SNPs in chr21 peaks with >= 10 reads:", nrow(ad), " median depth", median(ad$n),
      " pooled reference fraction", round(sum(ad$refc) / sum(ad$n), 3), " FDR < 0.1:", sum(ad$padj < 0.1),
      " with zero reads of one allele:", sum(ad$refc == 0 | ad$altc == 0), "\n")
  print(ad[order(p)][1:10, .(pos, ref, het_alt, refc, altc, p = signif(p, 2), padj = signif(padj, 2))])
  p18a <- ggplot(ad, aes(refc / n)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "#2a6f97", colour = "white") +
    geom_vline(xintercept = 0.5, linetype = 2) + geom_vline(xintercept = sum(ad$refc) / sum(ad$n), colour = "#c0392b") +
    labs(x = "fraction of reads with the reference allele", y = "het SNPs",
         title = sprintf("%d het SNPs in chr21 ATAC peaks", nrow(ad)))
  p18b <- ggplot(ad, aes(log2((refc + 0.5) / (altc + 0.5)), -log10(p), colour = padj < 0.1)) + geom_point(size = 1.8) +
    scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#5b6b7a"), labels = c(`TRUE` = "FDR < 0.1", `FALSE` = "n.s."), name = NULL) +
    labs(x = "log2(reference / alternative reads)", y = "−log10 binomial p", title = "Allelic imbalance") + theme(legend.position = "top")
  save(p18a + p18b, "fig18-allelic.png", w = 12, h = 4.6)
} else cat("BAM files not found in", BAMDIR, ": fig18 not made\n")
