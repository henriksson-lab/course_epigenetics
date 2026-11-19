# Real-data figures for ml-interpretation-variants/lec.qmd.
# Run from this folder:  Rscript make-mlinterp-figs.R   (about 3 minutes on 4 CPU threads)
# Needs: data/prepare-ml2.sh  (ml_ctcf_dataset.rds, ml_ctcf_cnn.pt: the model of the
#        exercise "How a CNN reads DNA"; ctcf_asb_adastra.tsv: ADASTRA CTCF allele-specific
#        binding SNPs), data/hg38.fa.bgz (+ .fai, .gzi).
#        R packages torch, Rsamtools, GenomicRanges, Biostrings, JASPAR2024, TFBSTools,
#        RSQLite, motifmatchr, ggplot2, patchwork, ggseqlogo.
#
# Writes fig04-ism-map.png, fig05-ism-logos.png, fig08-grad-vs-ism.png,
# fig11-aggregate.png, fig12-core-flank.png, fig15-asb-concordance.png,
# fig16-snp-example.png and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(torch); library(Rsamtools); library(GenomicRanges); library(Biostrings)
  library(JASPAR2024); library(TFBSTools); library(RSQLite); library(motifmatchr)
  library(ggplot2); library(patchwork); library(ggseqlogo)
})
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
torch_set_num_threads(4)

# ---- the exercise's setup (R0) -------------------------------------------------
ds    <- readRDS(file.path(D, "ml_ctcf_dataset.rds"))
model <- torch_load(file.path(D, "ml_ctcf_cnn.pt"))
model$eval()
one_hot <- function(s) {
  m <- do.call(rbind, strsplit(s, ""))
  x <- array(0, c(length(s), 4, ncol(m)))
  for (k in 1:4) x[, k, ] <- m == c("A", "C", "G", "T")[k]
  x
}
logit <- function(seqs) with_no_grad(as.numeric(model(torch_tensor(one_hot(seqs), dtype = torch_float()))))
test <- ds[ds$split == "test", ]
test$cnn <- logit(test$seq)
jaspar <- dbConnect(SQLite(), db(JASPAR2024()))
ctcf_m <- getMatrixByID(jaspar, "MA0139.2")
pwm_score <- function(seqs)
  motifScores(matchMotifs(ctcf_m, DNAStringSet(seqs), out = "scores", bg = "even", p.cutoff = 1))[, 1]
cat(sprintf("test logits: peaks median %.2f, background median %.2f\n",
            median(test$cnn[test$label == 1]), median(test$cnn[test$label == 0])))

ism <- function(s) {
  b   <- strsplit(s, "")[[1]]
  mut <- expand.grid(alt = c("A", "C", "G", "T"), pos = seq_along(b), stringsAsFactors = FALSE)
  mut <- mut[mut$alt != b[mut$pos], ]
  mutants <- vapply(seq_len(nrow(mut)), function(i) {
    x <- b; x[mut$pos[i]] <- mut$alt[i]; paste(x, collapse = "") }, "")
  mut$ref   <- b[mut$pos]
  mut$delta <- logit(mutants) - logit(s)
  mut
}
imp_matrix <- function(s, r) {           # importance at the reference base, 4 x L for ggseqlogo
  imp <- -tapply(r$delta, r$pos, mean)
  m <- one_hot(s)[1, , ] * rep(imp, each = 4); rownames(m) <- c("A", "C", "G", "T"); m
}
logo_plot <- function(m, cols, title = NULL)
  ggseqlogo(m[, cols], method = "custom", seq_type = "dna") +
    scale_x_continuous(breaks = seq(1, length(cols), 10), labels = cols[seq(1, length(cols), 10)]) +
    labs(y = "importance", title = title) + theme(legend.position = "none")

# ---- fig04: mutation map of a typical test peak (exercise R1) ----------------
pos_test <- test[test$label == 1, ]
ex <- pos_test[which.min(abs(pos_test$cnn - median(pos_test$cnn))), ]
r  <- ism(ex$seq)
cat(sprintf("fig04: example peak %s:%d-%d, logit %.2f (p = %.3f); strongest mutation %s at %d: %.2f; %d of 600 mutations change z by more than 1\n",
            ex$chr, ex$start, ex$start + 199, ex$cnn, plogis(ex$cnn), r$alt[which.min(r$delta)],
            r$pos[which.min(r$delta)], min(r$delta), sum(abs(r$delta) > 1)))
win <- 85:135
imp_ex <- imp_matrix(ex$seq, r)
p4a <- ggseqlogo(imp_ex[, win], method = "custom", seq_type = "dna") +
  scale_x_continuous(limits = c(0.5, length(win) + 0.5), expand = c(0, 0), breaks = NULL) +
  labs(y = "importance", x = NULL, title = sprintf("%s:%s, a typical test peak (z = %.1f)", ex$chr,
                                                   format(ex$start, big.mark = ","), ex$cnn)) +
  theme(legend.position = "none")
rr <- r[r$pos %in% win, ]
refs <- data.frame(pos = win, alt = strsplit(ex$seq, "")[[1]][win])
p4b <- ggplot(rr, aes(pos, factor(alt, levels = c("T", "G", "C", "A")), fill = delta)) +
  geom_tile(colour = "white") +
  geom_text(data = refs, aes(pos, factor(alt, levels = c("T", "G", "C", "A")), label = alt),
            inherit.aes = FALSE, size = 3.6, colour = "grey40") +
  scale_fill_gradient2(low = "#2a6f97", mid = "white", high = "#c0392b", name = "Δ logit") +
  scale_x_continuous(limits = c(min(win) - 0.5, max(win) + 0.5), expand = c(0, 0), breaks = seq(85, 135, 10)) +
  labs(x = "position in the 200-bp window (summit at 100)", y = "mutated to") +
  theme(legend.position = "right")
p4 <- p4a / p4b + plot_layout(heights = c(1, 1.1))
ggsave("fig04-ism-map.png", p4, width = 12, height = 5.6, dpi = 150, bg = "white")

# ---- fig05: attribution logos of three more test peaks ------------------------
set.seed(5)
picks <- pos_test[pwm_score(pos_test$seq) > 12, ]
picks <- picks[sample(nrow(picks), 3), ]
p5 <- lapply(seq_len(nrow(picks)), function(i) {
  ri <- ism(picks$seq[i]); m <- imp_matrix(picks$seq[i], ri)
  centre <- as.integer(names(which.max(-tapply(ri$delta, ri$pos, mean))))
  cols <- max(1, centre - 25):min(200, centre + 25)
  logo_plot(m, cols, sprintf("%s:%s (z = %.1f)", picks$chr[i], format(picks$start[i], big.mark = ","), picks$cnn[i]))
})
ggsave("fig05-ism-logos.png", wrap_plots(p5, ncol = 1), width = 12, height = 6.6, dpi = 150, bg = "white")

# ---- fig08: gradient x input vs ISM (exercise R2, then 200 peaks) -------------
gxi_of <- function(s) {
  x <- torch_tensor(one_hot(s), dtype = torch_float(), requires_grad = TRUE)
  model(x)$backward()
  colSums(as.array(x$grad)[1, , ] * one_hot(s)[1, , ])
}
imp1 <- -tapply(r$delta, r$pos, mean); g1 <- gxi_of(ex$seq)
cat(sprintf("fig08a: correlation gradient x input vs ISM, example peak: %.3f\n", cor(g1, imp1)))

# the 200 peaks of exercise R3: ISM, aligned at the best CTCF match
set.seed(2)
peaks <- pos_test[pwm_score(pos_test$seq) > 12, ]
cat(sprintf("test peaks with best PWM score > 12: %d of %d\n", nrow(peaks), nrow(pos_test)))
peaks <- peaks[sample(nrow(peaks), 200), ]
hits  <- matchMotifs(ctcf_m, DNAStringSet(peaks$seq), out = "positions", bg = "even", p.cutoff = 1e-4)[[1]]
best  <- t(sapply(hits, function(h) { h <- h[which.max(mcols(h)$score)]
                                       c(start = start(h), minus = mcols(h)$strand == "-") }))
cors <- numeric(nrow(peaks))
prof <- t(sapply(seq_len(nrow(peaks)), function(i) {
  ri  <- ism(peaks$seq[i])
  imp <- -tapply(ri$delta, ri$pos, mean)
  cors[i] <<- cor(imp, gxi_of(peaks$seq[i]))
  idx <- best[i, "start"] + (-20):34
  if (best[i, "minus"]) idx <- rev(idx)
  imp[pmin(pmax(idx, 1), 200)] * (idx >= 1 & idx <= 200)
}))
cat(sprintf("fig08b: per-peak correlation gradient x input vs ISM over 200 peaks: median %.3f, 10th percentile %.3f\n",
            median(cors), quantile(cors, 0.1)))
cat(sprintf("        motifs on minus strand: %d of 200\n", sum(best[, "minus"])))
p8a <- ggplot(data.frame(ism = imp1, gxi = g1), aes(ism, gxi)) + geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") + geom_point(colour = "#2a6f97", size = 2.2) +
  annotate("text", x = min(imp1), y = max(g1), hjust = 0, vjust = 1, size = 5,
           label = sprintf("r = %.2f", cor(g1, imp1))) +
  labs(title = "One peak: 200 positions", x = "ISM importance (601 forward passes)", y = "gradient × input (1 backward pass)")
p8b <- ggplot(data.frame(r = cors), aes(r)) + geom_histogram(binwidth = 0.02, fill = "#6b4fa0") +
  labs(title = "200 test peaks", x = "correlation of the two attributions per peak", y = "peaks")
ggsave("fig08-grad-vs-ism.png", p8a + p8b, width = 12, height = 4.6, dpi = 150, bg = "white")

# ---- fig11 and fig12: what the model uses around the motif (exercise R3, R4) ---
mean_imp <- colMeans(prof)
p   <- (Matrix(ctcf_m) + 0.25) / colSums(Matrix(ctcf_m) + 0.25)
ic  <- 2 + colSums(p * log2(p))
core <- 21:35
cat(sprintf("fig11/12: correlation of mean importance with JASPAR IC over the 15 motif positions: %.3f\n",
            cor(mean_imp[core], ic)))
tab <- c(high_ic = mean(mean_imp[core][ic > 1]), low_ic = mean(mean_imp[core][ic <= 1]),
         flank_1_10 = mean(mean_imp[c(11:20, 36:45)]), flank_11_20 = mean(mean_imp[c(1:10, 46:55)]))
print(round(tab, 3))
cat(sprintf("        motif positions with IC > 1 bit: %s\n", paste(which(ic > 1), collapse = ", ")))
d11 <- data.frame(pos = -20:34, imp = mean_imp,
                  where = ifelse(-20:34 >= 0 & -20:34 <= 14, "CTCF motif (JASPAR MA0139.2)", "flanks"))
p11a <- ggseqlogo(as.matrix(Matrix(ctcf_m)), method = "bits", seq_type = "dna") +
  scale_x_continuous(limits = c(-20.5, 34.5) + 1, breaks = NULL) +
  labs(title = "JASPAR CTCF motif (from 913 ChIP-seq sites)", y = "bits", x = NULL) +
  theme(legend.position = "none")
p11b <- ggplot(d11, aes(pos, imp, fill = where)) + geom_col(width = 0.8) +
  scale_fill_manual(values = c("#8c564b", "grey60"), name = NULL) +
  scale_x_continuous(limits = c(-20.5, 34.5)) +
  labs(title = "What our CNN uses: mean ISM importance in 200 test peaks, aligned on the motif",
       x = "position relative to the first base of the best CTCF match (minus-strand sites flipped)",
       y = "mean importance") +
  theme(legend.position = "bottom")
ggsave("fig11-aggregate.png", p11a / p11b + plot_layout(heights = c(1, 1.6)), width = 12, height = 6, dpi = 150, bg = "white")

d12 <- data.frame(position = 1:15, ic = ic, imp = mean_imp[core],
                  consensus = rownames(p)[apply(p, 2, which.max)])
p12a <- ggplot(d12, aes(ic, imp, label = paste0(position, consensus))) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey70", linewidth = 0.8) +
  geom_point(size = 3, colour = "#8c564b") + ggrepel::geom_text_repel(size = 4.4, seed = 1) +
  annotate("text", x = 0.3, y = max(d12$imp), hjust = 0, vjust = 1, size = 5,
           label = sprintf("r = %.2f", cor(d12$ic, d12$imp))) +
  labs(title = "Motif positions", x = "information content in JASPAR (bits)", y = "mean ISM importance")
lv <- c("motif\nIC > 1 bit", "motif\nIC ≤ 1 bit", "flanks\n1–10 bp", "flanks\n11–20 bp")
d12b <- data.frame(where = factor(lv, levels = lv),
                   imp = tab)
p12b <- ggplot(d12b, aes(where, imp)) + geom_col(fill = c("#8c564b", "#c49c94", "grey60", "grey75"), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", imp)), vjust = -0.4, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Core vs flanks", x = NULL, y = "mean ISM importance")
ggsave("fig12-core-flank.png", p12a + p12b + plot_layout(widths = c(1.2, 1)), width = 12, height = 4.8, dpi = 150, bg = "white")

# ---- fig15: real SNPs with allele-specific CTCF binding (exercise R5) ----------
asb <- read.delim(file.path(D, "ctcf_asb_adastra.tsv"))
names(asb) <- c("chr", "pos", "ref", "alt", "rsid", "tfs", "cell_types",
                "fdr_ref", "fdr_alt", "es_ref", "es_alt")
cat(sprintf("ADASTRA CTCF ASB SNPs (FDR < 0.05, all chromosomes): %d\n", nrow(asb)))
asb <- asb[asb$chr %in% c("chr8", "chr9"), ]
sig_ref <- !is.na(asb$fdr_ref) & asb$fdr_ref > -log10(0.05)
sig_alt <- !is.na(asb$fdr_alt) & asb$fdr_alt > -log10(0.05)
cat(sprintf("chr8+chr9: %d SNPs; both alleles significant (dropped): %d\n", nrow(asb), sum(sig_ref & sig_alt)))
keep <- xor(sig_ref, sig_alt)
asb <- asb[keep, ]
asb$preferred <- ifelse(sig_ref[keep], "ref", "alt")
print(table(asb$preferred))
genome <- FaFile(file.path(D, "hg38.fa.bgz"))
ref_seq <- as.character(getSeq(genome, GRanges(asb$chr, IRanges(asb$pos - 99, width = 200))))
stopifnot(substr(ref_seq, 100, 100) == asb$ref)
alt_seq <- ref_seq; substr(alt_seq, 100, 100) <- asb$alt
asb$d_cnn <- logit(alt_seq) - logit(ref_seq)
asb$d_pwm <- pwm_score(alt_seq) - pwm_score(ref_seq)
observed <- ifelse(asb$preferred == "alt", 1, -1)
asb$in_peak <- overlapsAny(GRanges(asb$chr, IRanges(asb$pos, width = 1)),
                           GRanges(test$chr[test$label == 1], IRanges(test$start[test$label == 1], width = 200)))
cat(sprintf("SNPs inside a test peak window (+-100 bp of a summit): %d\n", sum(asb$in_peak)))
conc <- function(d, t) { k <- abs(d) > t; c(n = sum(k), conc = mean(sign(d[k]) == observed[k])) }
thr <- c(0, 0.25, 0.5, 1, 1.5, 2, 3)
d15 <- rbind(data.frame(model = "CNN (our model)", t = thr, t(sapply(thr, function(t) conc(asb$d_cnn, t)))),
             data.frame(model = "CTCF PWM", t = thr, t(sapply(thr, function(t) conc(asb$d_pwm, t)))))
print(d15, digits = 3)
base_rate <- mean(observed == -1)
cat(sprintf("fig15: base rate (always predict the reference allele is bound better): %.3f; cor(d_cnn, d_pwm) = %.3f\n",
            base_rate, cor(asb$d_cnn, asb$d_pwm)))
p15 <- ggplot(d15, aes(t, conc, colour = model)) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50") +
  geom_hline(yintercept = base_rate, linetype = 2, colour = "grey40") +
  annotate("text", x = 3, y = base_rate - 0.025, hjust = 1, size = 4.4, colour = "grey30",
           label = sprintf("always \"reference binds better\": %.2f", base_rate)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.6) +
  ggrepel::geom_text_repel(aes(label = n), size = 4, show.legend = FALSE, seed = 1,
                           direction = "y", nudge_y = ifelse(d15$model == "CTCF PWM", 0.04, -0.04),
                           min.segment.length = 1) +
  scale_colour_manual(values = c("#2a6f97", "#8c564b"), name = NULL) +
  scale_y_continuous(limits = c(0.45, 1.02)) +
  labs(x = "only SNPs with a predicted effect |Δ| larger than this", y = "direction agrees with ADASTRA",
       title = sprintf("%s CTCF allele-specific SNPs on chr8 + chr9 (numbers: SNPs kept)", format(nrow(asb), big.mark = " "))) +
  theme(legend.position = "bottom")
ggsave("fig15-asb-concordance.png", p15, width = 11, height = 5, dpi = 150, bg = "white")

# ---- fig16: one SNP up close (exercise R6) ------------------------------------
big <- asb[order(-abs(asb$d_cnn)), c("rsid", "chr", "pos", "ref", "alt", "preferred", "d_cnn", "d_pwm")]
print(head(big, 8), row.names = FALSE)
i <- which(asb$rsid == big$rsid[1])
rr <- ism(ref_seq[i]); ra <- ism(alt_seq[i])
z_ref <- logit(ref_seq[i]); z_alt <- logit(alt_seq[i])
cat(sprintf("fig16: %s %s:%d %s>%s, preferred %s; z_ref %.2f, z_alt %.2f, d_cnn %.2f, d_pwm %.2f\n",
            asb$rsid[i], asb$chr[i], asb$pos[i], asb$ref[i], asb$alt[i], asb$preferred[i],
            z_ref, z_alt, asb$d_cnn[i], asb$d_pwm[i]))
cols <- 75:125
p16 <- logo_plot(imp_matrix(ref_seq[i], rr), cols, sprintf("%s, reference allele %s (z = %.1f): bound more in ChIP-seq",
                                                         asb$rsid[i], asb$ref[i], z_ref)) /
  logo_plot(imp_matrix(alt_seq[i], ra), cols, sprintf("alternative allele %s (z = %.1f): Δ = %.1f",
                                                    asb$alt[i], z_alt, asb$d_cnn[i]))
p16 <- p16 & geom_vline(xintercept = 26, linetype = 2, colour = "#c0392b")
ggsave("fig16-snp-example.png", p16, width = 12, height = 5.4, dpi = 150, bg = "white")
