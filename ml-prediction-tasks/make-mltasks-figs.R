# Real-data figures for ml-prediction-tasks/lec.qmd.
# Run from this folder:  Rscript make-mltasks-figs.R
# Needs: data/download.sh ml-prediction-tasks (GM12878 ATAC and CTCF peaks,
#        K562 ATAC peaks, hg38 blacklist, hg38.fa.bgz + .fai + .gzi);
#        R packages GenomicRanges, rtracklayer, Rsamtools, Biostrings,
#        JASPAR2024, TFBSTools, RSQLite, motifmatchr, glmnet, class,
#        ggplot2, patchwork.
#
# Uses the same seeds as the exercise and writes fig11-leakage-real.png,
# fig16-negatives-real.png, fig17-motif-baseline.png, fig18-cross-celltype.png,
# and prints the numbers quoted in the slides. About 15 minutes on one core
# (k-nearest neighbours dominates). Set ML_CACHE=/some/dir to keep results.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Rsamtools); library(Biostrings)
  library(JASPAR2024); library(TFBSTools); library(RSQLite); library(motifmatchr)
  library(glmnet); library(class); library(ggplot2); library(patchwork)
})
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
cache_dir <- Sys.getenv("ML_CACHE", "")
cached <- function(name, expr) {
  if (cache_dir == "") return(expr)
  f <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr; saveRDS(x, f); x
}
save_png <- function(p, f, w, h) ggsave(f, p, width = w, height = h, dpi = 110, bg = "white")
fa   <- FaFile(file.path(D, "hg38.fa.bgz"))
auto <- paste0("chr", 1:22)
bl   <- import(file.path(D, "hg38_blacklist.bed.gz"))
auc_roc <- function(score, label) {
  o  <- order(score, decreasing = TRUE)
  s  <- score[o]; l <- label[o]; k <- c(diff(s) != 0, TRUE)
  tpr <- cumsum(l)[k] / sum(l); fpr <- cumsum(1 - l)[k] / sum(1 - l)
  sum(diff(c(0, fpr)) * (tpr + c(0, head(tpr, -1))) / 2)
}
kmers <- function(gr) oligonucleotideFrequency(getSeq(fa, gr), 4, as.prob = TRUE)
# open/closed windows: identical to R1 of the Machine learning from scratch exercise
dat <- cached("ml_atac", {
  atac <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
  pk  <- atac[seqnames(atac) %in% auto]
  pk  <- pk[order(-pk$signalValue)]
  pk  <- pk[!duplicated(paste(seqnames(pk), start(pk), end(pk)))]
  pos <- resize(GRanges(seqnames(pk), IRanges(start(pk) + pk$peak, width = 1)), 500, fix = "center")
  pos <- pos[!overlapsAny(pos, bl) & countOverlaps(pos, pos) == 1]
  set.seed(1)
  pos <- sort(pos[sample(length(pos), 20000)])
  len <- seqlengths(seqinfo(fa))[auto]
  set.seed(2)
  chr <- sample(auto, 60000, replace = TRUE, prob = len)
  neg <- GRanges(chr, IRanges(floor(runif(60000, 1, len[chr] - 500)), width = 500))
  neg <- neg[!overlapsAny(neg, atac, maxgap = 500) & !overlapsAny(neg, bl)]
  neg <- neg[letterFrequency(getSeq(fa, neg), "N")[, 1] == 0]
  neg <- sort(neg[1:20000])
  rbind(data.frame(chrom = as.character(seqnames(pos)), start = start(pos) - 1, end = end(pos), open = 1),
        data.frame(chrom = as.character(seqnames(neg)), start = start(neg) - 1, end = end(neg), open = 0))
})

# ---- A. leakage (exercise R1) ------------------------------------------------------
leak <- cached("leak", {
  set.seed(4)
  reg <- dat[c(sample(which(dat$open == 1), 3000), sample(which(dat$open == 0), 3000)), ]
  reg$region <- seq_len(nrow(reg))
  win <- do.call(rbind, lapply(c(-150, 0, 150), function(o)
    transform(reg, start = start + o, end = end + o, shift = o)))
  gr  <- GRanges(win$chrom, IRanges(win$start + 1, win$end))
  X   <- kmers(gr)
  win$gc <- letterFrequency(getSeq(fa, gr), "GC")[, 1] / width(gr)
  y   <- win$open
  set.seed(5)
  splits <- list(random     = runif(nrow(win)) < 0.2,
                 region     = win$region %in% sample(reg$region, 0.2 * nrow(reg)),
                 chromosome = win$chrom %in% c("chr8", "chr9"))
  evaluate <- function(te) {
    set.seed(1)
    nn1  <- knn(X[!te, ], X[te, ], factor(y[!te]), k = 1)
    nn25 <- knn(X[!te, ], X[te, ], factor(y[!te]), k = 25, prob = TRUE)
    p25  <- ifelse(nn25 == "1", attr(nn25, "prob"), 1 - attr(nn25, "prob"))
    rdg  <- cv.glmnet(X[!te, ], y[!te], family = "binomial", alpha = 0, nfolds = 5)
    gcm  <- glm(y ~ gc, binomial, data.frame(y = y, gc = win$gc)[!te, ])
    c(knn1_accuracy = mean(as.character(nn1) == y[te]),
      knn25_AUC = auc_roc(p25, y[te]),
      ridge_AUC = auc_roc(as.numeric(predict(rdg, X[te, ], s = "lambda.min")), y[te]),
      gc_AUC    = auc_roc(predict(gcm, data.frame(gc = win$gc[te])), y[te]))
  }
  res <- sapply(splits, evaluate)
  te <- splits$random
  nn <- knn(X[!te, ], X[te, ], factor(win$region[!te]), k = 1)
  list(res = res, same_region = mean(as.character(nn) == win$region[te]))
})
print(round(leak$res, 3)); cat("random split, nearest neighbour = shifted copy:", round(leak$same_region, 3), "\n")
lk <- as.data.frame(as.table(leak$res)); names(lk) <- c("measure", "split", "value")
lk$measure <- factor(lk$measure, levels = c("knn1_accuracy", "knn25_AUC", "ridge_AUC", "gc_AUC"),
                     labels = c("1-nearest neighbour\n(accuracy)", "25-nearest neighbours\n(ROC AUC)",
                                "ridge on 4-mers\n(ROC AUC)", "GC only\n(ROC AUC)"))
lk$split <- factor(lk$split, levels = c("random", "region", "chromosome"),
                   labels = c("random windows", "whole regions", "whole chromosomes"))
p11 <- ggplot(lk, aes(measure, value, fill = split)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", value)), position = position_dodge(0.8), vjust = -0.4, size = 4.2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_fill_manual(values = c("#c0392b", "#e67e22", "#2a6f97"), name = "test set = ") +
  coord_cartesian(ylim = c(0.4, 1.02)) + labs(x = NULL, y = "performance on the test windows") +
  theme(legend.position = "top")
save_png(p11, "fig11-leakage-real.png", 12, 5.2)

# ---- B. negatives for CTCF (exercise R2) -----------------------------------------------
negs <- cached("negs", {
  ctcf <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak")
  ctcf <- ctcf[seqnames(ctcf) %in% auto]
  pos  <- resize(GRanges(seqnames(ctcf), IRanges(start(ctcf) + ctcf$peak, width = 1)), 500, fix = "center")
  pos  <- pos[!overlapsAny(pos, bl) & countOverlaps(pos, pos) == 1]
  set.seed(1); pos <- sort(pos[sample(length(pos), 10000)])
  gc_of <- function(gr) { s <- getSeq(fa, gr); f <- letterFrequency(s, c("G", "C", "N"))
    data.frame(gc = (f[, 1] + f[, 2]) / width(s), n = f[, 3]) }
  pos$gc <- gc_of(pos)$gc
  len <- seqlengths(seqinfo(fa))[auto]
  set.seed(2)
  chr  <- sample(auto, 300000, replace = TRUE, prob = len)
  pool <- GRanges(chr, IRanges(floor(runif(300000, 1, len[chr] - 500)), width = 500))
  pool <- pool[!overlapsAny(pool, ctcf, maxgap = 500) & !overlapsAny(pool, bl)]
  g    <- gc_of(pool); pool <- pool[g$n == 0]; pool$gc <- g$gc[g$n == 0]
  rnd <- pool[1:10000]
  bin  <- function(x) cut(x, seq(0, 1, 0.02), include.lowest = TRUE)
  need <- table(bin(pos$gc))
  set.seed(3)
  idx <- unlist(lapply(names(need)[need > 0], function(b) {
    w <- which(bin(pool$gc) == b); w[sample.int(length(w), min(need[b], length(w)))] }))
  gcm <- pool[idx]
  atac <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
  atac <- atac[seqnames(atac) %in% auto]
  aw   <- resize(GRanges(seqnames(atac), IRanges(start(atac) + atac$peak, width = 1)), 500, fix = "center")
  aw   <- unique(aw[!overlapsAny(aw, ctcf, maxgap = 500) & !overlapsAny(aw, bl)])
  set.seed(4)
  acc  <- aw[sample(length(aw), 10000)]; acc$gc <- gc_of(acc)$gc
  jdb <- dbConnect(SQLite(), db(JASPAR2024()))
  ctcf_pwm <- getMatrixByID(jdb, "MA0139.2")
  feat <- function(gr) {
    s  <- getSeq(fa, gr)
    sc <- matchMotifs(ctcf_pwm, s, p.cutoff = 1, bg = "even", out = "scores")      # every window: its best score
    ht <- matchMotifs(ctcf_pwm, s, p.cutoff = 5e-5, bg = "even", out = "scores")   # strict: is there a match?
    list(chrom = as.character(seqnames(gr)), gc = letterFrequency(s, "GC")[, 1] / width(s),
         motif = as.numeric(motifScores(sc)[, 1]), hit = as.numeric(motifCounts(ht)[, 1]) > 0,
         k4 = oligonucleotideFrequency(s, 4, as.prob = TRUE))
  }
  fp <- feat(pos)
  sets <- list(random = rnd, gc_matched = gcm, accessible = acc)
  fns <- lapply(sets, feat)
  res <- sapply(names(sets), function(nm) {
    fn <- fns[[nm]]
    y  <- rep(1:0, c(length(pos), length(sets[[nm]])))
    te <- c(fp$chrom, fn$chrom) %in% c("chr8", "chr9")
    d  <- data.frame(y = y, gc = c(fp$gc, fn$gc), motif = c(fp$motif, fn$motif))
    K  <- rbind(fp$k4, fn$k4)
    m_gc  <- glm(y ~ gc, binomial, d[!te, ])
    m_mot <- glm(y ~ motif, binomial, d[!te, ])
    fold  <- match(c(fp$chrom, fn$chrom), auto) %% 5 + 1
    m_k4  <- cv.glmnet(K[!te, ], y[!te], family = "binomial", alpha = 0, foldid = fold[!te])
    c(gc_AUC    = auc_roc(predict(m_gc, d[te, ]), y[te]),
      motif_AUC = auc_roc(predict(m_mot, d[te, ]), y[te]),
      kmer_AUC  = auc_roc(as.numeric(predict(m_k4, K[te, ], s = "lambda.min")), y[te]),
      motif_in_negatives = mean(fn$hit))
  })
  list(res = res, gc = list(positives = pos$gc, random = rnd$gc, gc_matched = gcm$gc, accessible = acc$gc),
       score = list(positives = fp$motif, random = fns$random$motif, gc_matched = fns$gc_matched$motif,
                    accessible = fns$accessible$motif),
       hit = list(positives = fp$hit, random = fns$random$hit, gc_matched = fns$gc_matched$hit,
                  accessible = fns$accessible$hit),
       hit_pos = mean(fp$hit))
})
print(round(negs$res, 3)); cat("CTCF peaks with a motif match:", round(negs$hit_pos, 3), "\n")
print(sapply(negs$gc, function(x) round(c(mean = mean(x), n = length(x)), 3)))
setcol <- c("CTCF peaks" = "#8c564b", "random genome" = "#a8324a", "GC-matched" = "#2a6f97",
            "open, not CTCF" = "#e67e22")
setlab <- c(positives = "CTCF peaks", random = "random genome", gc_matched = "GC-matched", accessible = "open, not CTCF")
gd <- do.call(rbind, lapply(names(negs$gc), function(n) data.frame(set = setlab[n], gc = negs$gc[[n]])))
gd$set <- factor(gd$set, levels = setlab)
p16a <- ggplot(gd, aes(gc, colour = set, linetype = set == "CTCF peaks")) + geom_density(linewidth = 1.2, adjust = 1.3) +
  scale_colour_manual(values = setcol, name = NULL) + scale_linetype_manual(values = c(2, 1), guide = "none") +
  labs(title = "GC content of the windows", x = "GC fraction", y = "density") +
  theme(legend.position = c(0.8, 0.8))
nr <- as.data.frame(as.table(negs$res[1:3, ])); names(nr) <- c("model", "negatives", "AUC")
nr$negatives <- factor(setlab[as.character(nr$negatives)], levels = setlab[-1])
nr$model <- factor(nr$model, levels = c("gc_AUC", "kmer_AUC", "motif_AUC"),
                   labels = c("GC only", "256 4-mers (ridge)", "CTCF PWM score"))
p16b <- ggplot(nr, aes(model, AUC, fill = negatives)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", AUC)), position = position_dodge(0.8), vjust = -0.4, size = 4.2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_fill_manual(values = setcol[-1], name = "negatives") + coord_cartesian(ylim = c(0.4, 0.95)) +
  labs(title = "CTCF peaks vs. each negative set", x = NULL, y = "ROC AUC, chr8 + chr9") +
  theme(legend.position = "top")
save_png(p16a + p16b + plot_layout(widths = c(1, 1.35)), "fig16-negatives-real.png", 14, 5.4)

# ---- fig17: a single PWM as the baseline ---------------------------------------------------
f17 <- data.frame(set = factor(setlab, levels = setlab), frac = sapply(negs$hit, mean))
print(f17, digits = 3)
p17a <- ggplot(f17, aes(set, frac, fill = set)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), vjust = -0.4, size = 4.6) +
  scale_fill_manual(values = setcol, guide = "none") + scale_y_continuous(limits = c(0, 0.82)) +
  labs(title = "Windows with a match at p < 5e-5", x = NULL, y = "fraction of windows")
sd17 <- do.call(rbind, lapply(names(negs$score), function(n)
  data.frame(set = unname(setlab[n]), score = negs$score[[n]])))
sd17$set <- factor(sd17$set, levels = setlab)
cat("best CTCF PWM score (p.cutoff = 1), median by set:\n"); print(tapply(sd17$score, sd17$set, median))
p17b <- ggplot(sd17, aes(score, colour = set, linetype = set == "CTCF peaks")) +
  geom_density(linewidth = 1.2, adjust = 1.2) +
  scale_colour_manual(values = setcol, name = NULL) + scale_linetype_manual(values = c(2, 1), guide = "none") +
  labs(title = "Best match score in every window", x = "best CTCF PWM log-odds score in the window", y = "density") +
  theme(legend.position = c(0.14, 0.62))
save_png(p17a + p17b + plot_layout(widths = c(1, 1.3)), "fig17-motif-baseline.png", 13, 4.8)

# ---- C. cross-cell-type (exercise R3) ----------------------------------------------------------
cross <- cached("cross", {
  k562 <- import(file.path(D, "k562_atac_peaks.bed.gz"), format = "narrowPeak")
  gm   <- import(file.path(D, "gm12878_atac_peaks.bed.gz"), format = "narrowPeak")
  summits <- function(pk) {
    pk <- pk[seqnames(pk) %in% auto]
    pk <- pk[order(-pk$signalValue)]
    pk <- pk[!duplicated(paste(seqnames(pk), start(pk), end(pk)))]
    w  <- resize(GRanges(seqnames(pk), IRanges(start(pk) + pk$peak, width = 1)), 500, fix = "center")
    w[!overlapsAny(w, bl) & countOverlaps(w, w) == 1]
  }
  gm_w <- summits(gm); k5_w <- summits(k562)
  neg <- with(dat[dat$open == 0, ], GRanges(chrom, IRanges(start + 1, end)))
  neg <- neg[!overlapsAny(neg, k562, maxgap = 500)]
  test_chr <- c("chr8", "chr9")
  is_test  <- function(gr) as.character(seqnames(gr)) %in% test_chr
  train_model <- function(open_w) {
    set.seed(8)
    p <- open_w[!is_test(open_w)]; p <- p[sample(length(p), 15000)]
    n <- neg[!is_test(neg)]
    fold <- match(c(as.character(seqnames(p)), as.character(seqnames(n))), auto) %% 5 + 1
    cv.glmnet(rbind(kmers(p), kmers(n)), rep(1:0, c(length(p), length(n))),
              family = "binomial", alpha = 0, foldid = fold)
  }
  fit_gm <- train_model(gm_w); fit_k5 <- train_model(k5_w)
  t_gm <- gm_w[is_test(gm_w)]; t_k5 <- k5_w[is_test(k5_w)]
  tests <- list(GM12878_only = t_gm[!overlapsAny(t_gm, k562)],
                shared       = t_gm[overlapsAny(t_gm, k562)],
                K562_only    = t_k5[!overlapsAny(t_k5, gm)])
  t_neg <- neg[is_test(neg)]
  Kneg <- kmers(t_neg)
  res <- sapply(tests, function(w) {
    K <- rbind(kmers(w), Kneg); y <- rep(1:0, c(length(w), length(t_neg)))
    c(trained_on_GM12878 = auc_roc(as.numeric(predict(fit_gm, K, s = "lambda.min")), y),
      trained_on_K562    = auc_roc(as.numeric(predict(fit_k5, K, s = "lambda.min")), y))
  })
  list(res = res, n = c(sapply(tests, length), negatives = length(t_neg)))
})
print(round(cross$res, 3)); print(cross$n)
cd <- as.data.frame(as.table(cross$res)); names(cd) <- c("trained", "test", "AUC")
cd$trained <- factor(cd$trained, levels = c("trained_on_GM12878", "trained_on_K562"),
                     labels = c("trained on GM12878", "trained on K562"))
cd$test <- factor(cd$test, levels = c("GM12878_only", "shared", "K562_only"),
                  labels = sprintf("%s\n(n = %d)", c("open only in GM12878", "open in both", "open only in K562"),
                                   cross$n[1:3]))
p18 <- ggplot(cd, aes(test, AUC, fill = trained)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", AUC)), position = position_dodge(0.8), vjust = -0.4, size = 4.4) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_fill_manual(values = c("#2a6f97", "#c0392b"), name = NULL) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(x = sprintf("test windows on chr8 + chr9, each set against the same %d closed windows", cross$n[4]),
       y = "ROC AUC") + theme(legend.position = "top")
save_png(p18, "fig18-cross-celltype.png", 11, 5.2)
cat("done\n")
