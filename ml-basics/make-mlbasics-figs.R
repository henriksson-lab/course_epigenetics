# Real-data figures for ml-basics/lec.qmd.
# Run from this folder:  Rscript make-mlbasics-figs.R
# Needs: data/download.sh ml-basics (GM12878 ATAC peaks, hg38 blacklist,
#        hg38.fa.bgz + .fai + .gzi); R packages GenomicRanges, rtracklayer,
#        Rsamtools, Biostrings, JASPAR2024, TFBSTools, RSQLite, motifmatchr,
#        glmnet, ggplot2, patchwork.
#
# Builds the same open/closed data set as exercise R1 (same seeds), and writes
# fig04-features-real.png, fig05-logistic-fit.png, fig07-decision-boundary.png,
# fig08-log-loss.png, fig09-gradient-descent.png, fig11-overfitting.png,
# fig14-reg-path.png, fig16-cv-curve.png, fig18-threshold.png,
# fig20-roc-real.png, fig21-pr-real.png, fig22-imbalance.png,
# fig23-baselines.png, fig24-coefficients.png, and prints the numbers quoted in
# the slides. About 5 minutes on one core. Set ML_CACHE=/some/dir to keep the
# data set between runs.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Rsamtools); library(Biostrings)
  library(JASPAR2024); library(TFBSTools); library(RSQLite); library(motifmatchr)
  library(glmnet); library(ggplot2); library(patchwork)
})
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
colO <- "#2e8b57"; colC <- "#a8324a"
cache_dir <- Sys.getenv("ML_CACHE", "")
cached <- function(name, expr) {
  if (cache_dir == "") return(expr)
  f <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  x <- expr; saveRDS(x, f); x
}
fa <- FaFile(file.path(D, "hg38.fa.bgz"))
auto <- paste0("chr", 1:22)

# ---- the data set: identical to exercise R1 ------------------------------------
dat <- cached("ml_atac", {
  bl   <- import(file.path(D, "hg38_blacklist.bed.gz"))
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
  jdb  <- dbConnect(SQLite(), db(JASPAR2024()))
  pwms <- do.call(PFMatrixList, lapply(c("MA0139.2", "MA0080.7", "MA0154.5"),
                                       function(id) getMatrixByID(jdb, id)))
  names(pwms) <- c("ctcf", "spi1", "ebf1")
  features <- function(gr) {
    s  <- getSeq(fa, gr); L <- width(s)
    nc <- letterFrequency(s, "C")[, 1]; ng <- letterFrequency(s, "G")[, 1]
    cg <- vcountPattern("CG", s)
    mc <- as.matrix(motifCounts(matchMotifs(pwms, s, p.cutoff = 5e-5, bg = "even", out = "scores")))
    colnames(mc) <- names(pwms)
    data.frame(chrom = as.character(seqnames(gr)), start = start(gr) - 1, end = end(gr),
               gc = (nc + ng) / L, cpg_oe = cg * L / pmax(nc * ng, 1), mc)
  }
  rbind(cbind(features(pos), open = 1), cbind(features(neg), open = 0))
})
test <- dat$chrom %in% c("chr8", "chr9")
train <- dat[!test, ]; tst <- dat[test, ]
cat("windows", nrow(dat), "train", nrow(train), "test", nrow(tst), "\n")

curves <- function(score, label) {
  o  <- order(score, decreasing = TRUE); s <- score[o]; l <- label[o]
  tp <- cumsum(l); fp <- cumsum(1 - l); k <- c(diff(s) != 0, TRUE)
  data.frame(threshold = s[k], tpr = tp[k] / sum(l), fpr = fp[k] / sum(1 - l),
             precision = tp[k] / (tp[k] + fp[k]))
}
auc_roc  <- function(cv) sum(diff(c(0, cv$fpr)) * (cv$tpr + c(0, head(cv$tpr, -1))) / 2)
avg_prec <- function(cv) sum(diff(c(0, cv$tpr)) * cv$precision)
save_png <- function(p, f, w, h) ggsave(f, p, width = w, height = h, dpi = 110, bg = "white")
lab_oc <- function(x) factor(ifelse(x == 1, "open", "closed"), levels = c("open", "closed"))

# ---- fig04: the features, open vs closed ----------------------------------------
d4 <- data.frame(dat, state = lab_oc(dat$open))
p4a <- ggplot(d4, aes(gc, colour = state)) + geom_density(linewidth = 1.2) +
  scale_colour_manual(values = c(open = colO, closed = colC), name = NULL) +
  labs(title = "GC fraction", x = "GC fraction of the window", y = "density") +
  theme(legend.position = c(0.85, 0.85))
p4b <- ggplot(d4, aes(pmin(cpg_oe, 1.2), colour = state)) + geom_density(linewidth = 1.2) +
  scale_colour_manual(values = c(open = colO, closed = colC), guide = "none") +
  labs(title = "CpG observed / expected", x = "CpG o/e (capped at 1.2)", y = NULL)
mot <- do.call(rbind, lapply(c("ctcf", "spi1", "ebf1"), function(m)
  data.frame(motif = toupper(m), state = lab_oc(c(1, 0)),
             frac = c(mean(dat[[m]][dat$open == 1] > 0), mean(dat[[m]][dat$open == 0] > 0)))))
p4c <- ggplot(mot, aes(motif, frac, fill = state)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(open = colO, closed = colC), guide = "none") +
  labs(title = "Windows with a motif match", x = NULL, y = "fraction of windows")
save_png(p4a + p4b + p4c, "fig04-features-real.png", 14, 4.4)
print(aggregate(cbind(gc, cpg_oe, ctcf, spi1, ebf1) ~ open, dat, mean))
print(mot)

# ---- fig05: linear vs logistic fit on GC ------------------------------------------
m1 <- glm(open ~ gc, binomial, train); l1 <- lm(open ~ gc, train)
cat("glm(open ~ gc):", round(coef(m1), 3), " lm:", round(coef(l1), 3), "\n")
set.seed(10)
sub <- train[sample(nrow(train), 800), ]
grid <- data.frame(gc = seq(0.2, 0.85, by = 0.005))
grid$logistic <- predict(m1, grid, type = "response"); grid$linear <- predict(l1, grid)
cat("linear prediction at GC 0.85:", round(predict(l1, data.frame(gc = 0.85)), 2),
    " at 0.2:", round(predict(l1, data.frame(gc = 0.2)), 2), "\n")
# binned observed fraction
bins <- cut(train$gc, seq(0.2, 0.85, by = 0.05))
obs <- aggregate(train$open, list(bin = bins), mean)
obs$gc <- seq(0.225, 0.825, by = 0.05)[as.integer(obs$bin)]
p5 <- ggplot() +
  geom_jitter(data = sub, aes(gc, open), height = 0.04, width = 0, alpha = 0.25, size = 1.3,
              colour = ifelse(sub$open == 1, colO, colC)) +
  geom_point(data = obs, aes(gc, x), size = 3.2, shape = 21, fill = "white", colour = "#1b2a3a", stroke = 1.2) +
  geom_line(data = grid, aes(gc, linear), colour = "#5b6b7a", linewidth = 1.2, linetype = 2) +
  geom_line(data = grid, aes(gc, logistic), colour = "#2a6f97", linewidth = 1.6) +
  annotate("text", x = 0.70, y = 1.33, label = "straight line: goes above 1", hjust = 1, colour = "#5b6b7a", size = 4.6) +
  annotate("text", x = 0.75, y = 0.62, label = "logistic curve", hjust = 0, colour = "#2a6f97", size = 4.8) +
  annotate("text", x = 0.205, y = 0.62, label = "circles: observed fraction open\nin bins of GC", hjust = 0,
           colour = "#1b2a3a", size = 4.2) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.45, 1.45)) +
  labs(x = "GC fraction of the window", y = "open (1) or closed (0)")
save_png(p5, "fig05-logistic-fit.png", 11, 5)

# ---- fig07: two features and the decision boundary -------------------------------
m2 <- glm(open ~ gc + cpg_oe, binomial, train)
cat("glm(open ~ gc + cpg_oe):", round(coef(m2), 3), "\n")
set.seed(11)
s7 <- tst[sample(nrow(tst), 1500), ]
gx <- expand.grid(gc = seq(0.25, 0.8, length.out = 120), cpg_oe = seq(0, 1.1, length.out = 120))
gx$p <- predict(m2, gx, type = "response")
p7 <- ggplot(s7, aes(gc, cpg_oe)) +
  geom_point(aes(colour = lab_oc(open)), alpha = 0.45, size = 1.4) +
  geom_contour(data = gx, aes(z = p), breaks = c(0.25, 0.5, 0.75), colour = "#1b2a3a",
               linewidth = c(0.8), linetype = 1) +
  scale_colour_manual(values = c(open = colO, closed = colC), name = NULL) +
  coord_cartesian(xlim = c(0.25, 0.8), ylim = c(0, 1.1)) +
  labs(x = "GC fraction", y = "CpG o/e") + theme(legend.position = "right")
b7 <- coef(m2)
lab7 <- data.frame(cpg_oe = 1.06, l = c("P = 0.25", "P = 0.5", "P = 0.75"))
lab7$gc <- (qlogis(c(0.25, 0.5, 0.75)) - b7[1] - b7[3] * lab7$cpg_oe) / b7[2] + 0.008
p7 <- p7 + geom_text(data = lab7, aes(label = l), size = 4.3, hjust = 0)
save_png(p7, "fig07-decision-boundary.png", 10, 5.2)

# ---- fig08: log-loss ----------------------------------------------------------------
pp <- seq(0.005, 0.995, by = 0.005)
d8 <- rbind(data.frame(p = pp, loss = -log(pp), truth = "true label 1 (open):  −log p"),
            data.frame(p = pp, loss = -log(1 - pp), truth = "true label 0 (closed):  −log(1 − p)"))
pts8 <- data.frame(p = c(0.9, 0.5, 0.1), loss = -log(c(0.9, 0.5, 0.1)))
pts8$lab <- sprintf("p = %.1f: loss %.2f", pts8$p, pts8$loss)

p8 <- ggplot(d8, aes(p, loss, colour = truth)) + geom_line(linewidth = 1.5) +
  geom_point(data = pts8, aes(p, loss), inherit.aes = FALSE, size = 3, colour = colO) +
  annotate("text", x = 0.6, y = 3.6, hjust = 0, size = 4.6, colour = colO,
           label = paste(c("open window, loss = −log p:", pts8$lab), collapse = "\n")) +
  scale_colour_manual(values = c(colC, colO), name = NULL) +
  coord_cartesian(ylim = c(0, 4.2)) +
  labs(x = "predicted P(open)", y = "loss for this window") +
  theme(legend.position = "top")
save_png(p8, "fig08-log-loss.png", 9.5, 5)

# ---- fig09: gradient descent on the GC model ----------------------------------------
X <- cbind(1, train$gc); y <- train$open
sigmoid <- function(z) 1 / (1 + exp(-z))
log_loss <- function(b) { p <- sigmoid(X %*% b); -mean(y * log(p) + (1 - y) * log(1 - p)) }
b <- c(0, 0); path <- matrix(b, 1)
for (i in 1:2000) {
  b <- b - 5 * as.numeric(crossprod(X, sigmoid(X %*% b) - y) / nrow(X))
  if (i %in% c(1:20, seq(25, 2000, by = 25))) path <- rbind(path, b)
}
cat("gradient descent after 2000 steps:", round(b, 3), " loss", round(log_loss(b), 4),
    " start loss", round(log_loss(c(0, 0)), 4), "\n")
g9 <- expand.grid(b0 = seq(-6.5, 1, length.out = 90), b1 = seq(-1, 14, length.out = 90))
g9$loss <- apply(g9, 1, function(v) log_loss(v))
path <- data.frame(b0 = path[, 1], b1 = path[, 2])
p9 <- ggplot(g9, aes(b0, b1, z = loss)) +
  geom_contour_filled(breaks = c(0.585, 0.6, 0.62, 0.65, 0.69, 0.75, 0.85, 1, 1.5, 10), alpha = 0.85) +
  scale_fill_brewer(palette = "Blues", direction = -1, name = "mean log-loss") +
  geom_path(data = path, aes(b0, b1), inherit.aes = FALSE, colour = "#c0392b", linewidth = 1.1) +
  geom_point(data = path[c(1, nrow(path)), ], aes(b0, b1), inherit.aes = FALSE, colour = "#c0392b", size = 3.5) +
  annotate("label", x = 0, y = -0.6, label = "start (0, 0)", size = 4.4, fill = "white") +
  annotate("label", x = -4.6, y = 12.2, label = "end: glm's answer", size = 4.4, fill = "white") +
  labs(x = expression(b[0] ~ "(intercept)"), y = expression(b[1] ~ "(per unit GC)"))
save_png(p9, "fig09-gradient-descent.png", 10, 5.4)

# ---- fig11: overfitting with k-mer features -------------------------------------------
s4  <- getSeq(fa, GRanges(dat$chrom, IRanges(dat$start + 1, dat$end)))
kmr <- lapply(1:4, function(k) oligonucleotideFrequency(s4, k, as.prob = TRUE))
fit_k <- function(n, seed) {
  set.seed(seed); small <- sample(which(!test), n)
  do.call(rbind, lapply(1:4, function(k) {
    d <- data.frame(open = dat$open, kmr[[k]])
    m <- suppressWarnings(glm(open ~ ., binomial, d[small, ]))
    data.frame(n = n, k = k, p = ncol(kmr[[k]]),
               set = c("training", "test"),
               AUC = c(auc_roc(curves(predict(m, d[small, ]), dat$open[small])),
                       auc_roc(curves(suppressWarnings(predict(m, d[test, ])), dat$open[test]))))
  }))
}
o11 <- rbind(fit_k(500, 3), fit_k(5000, 3))
print(o11)
o11$n <- factor(paste(o11$n, "training windows"), levels = c("500 training windows", "5000 training windows"))
p11 <- ggplot(o11, aes(p, AUC, colour = set)) + geom_line(linewidth = 1.3) + geom_point(size = 3) +
  facet_wrap(~ n) + scale_x_log10(breaks = c(4, 16, 64, 256)) +
  scale_colour_manual(values = c(training = "#2a6f97", test = "#c0392b"), name = NULL) +
  labs(x = "number of features (k-mer counts, k = 1…4)", y = "ROC AUC") +
  theme(legend.position = "top", strip.text = element_text(size = 15, face = "bold"))
save_png(p11, "fig11-overfitting.png", 11, 5)

# ---- fig14 + fig16: regularisation paths and cross-validation ---------------------------
fold <- match(dat$chrom, auto) %% 5 + 1
set.seed(3); small <- sample(which(!test), 500)
X4 <- kmr[[4]]
fit_l <- glmnet(X4[small, ], dat$open[small], family = "binomial", alpha = 1)
fit_r <- glmnet(X4[small, ], dat$open[small], family = "binomial", alpha = 0)
sdx <- apply(X4[small, ], 2, sd)                     # plot coefficients per SD of each 4-mer
path_df <- function(fit, name) {
  b <- as.matrix(fit$beta) * sdx
  data.frame(kmer = rep(rownames(b), ncol(b)), loglambda = rep(log10(fit$lambda), each = nrow(b)),
             coef = as.numeric(b), model = name)
}
hl <- c("CGCG", "TTTT", "GGAA", "AATT")
pd <- rbind(path_df(fit_r, "L2 (ridge): all shrink, none reach 0"),
            path_df(fit_l, "L1 (lasso): many become exactly 0"))
pd$hl <- ifelse(pd$kmer %in% hl, pd$kmer, "other")
p14 <- ggplot(pd, aes(loglambda, coef, group = kmer)) +
  geom_line(data = pd[pd$hl == "other", ], colour = "grey75", linewidth = 0.4) +
  geom_line(data = pd[pd$hl != "other", ], aes(colour = hl), linewidth = 1.2) +
  facet_wrap(~ model, scales = "free") + scale_x_reverse() +
  scale_colour_manual(values = c(CGCG = "#c0392b", TTTT = "#2a6f97", GGAA = "#e67e22", AATT = "#6b4fa0"),
                      name = "4-mer") +
  labs(x = expression(log[10] ~ lambda ~ "(strong penalty on the left, weak on the right)"),
       y = "coefficient (per SD of the 4-mer count)") + theme(strip.text = element_text(size = 14, face = "bold"))
save_png(p14, "fig14-reg-path.png", 13, 5)

set.seed(4)
cvl <- cv.glmnet(X4[small, ], dat$open[small], family = "binomial", alpha = 1,
                 foldid = fold[small], type.measure = "auc")
test_auc <- sapply(cvl$lambda, function(l)
  auc_roc(curves(as.numeric(predict(cvl$glmnet.fit, X4[test, ], s = l)), dat$open[test])))
d16 <- rbind(data.frame(loglambda = log10(cvl$lambda), AUC = cvl$cvm, what = "cross-validation (training data only)"),
             data.frame(loglambda = log10(cvl$lambda), AUC = test_auc, what = "test chromosomes (for illustration)"))
cat("lasso n=500: lambda.min", signif(cvl$lambda.min, 3), "CV AUC", round(max(cvl$cvm), 3),
    "test AUC at lambda.min", round(test_auc[cvl$lambda == cvl$lambda.min], 3),
    "nonzero", cvl$nzero[cvl$lambda == cvl$lambda.min], "\n")
p16 <- ggplot(d16, aes(loglambda, AUC, colour = what)) + geom_line(linewidth = 1.3) +
  geom_vline(xintercept = log10(cvl$lambda.min), linetype = 2) +
  annotate("text", x = log10(cvl$lambda.min), y = 0.52, label = " λ chosen by CV", hjust = 0, size = 4.6) +
  annotate("text", x = max(log10(cvl$lambda)) - 0.12, y = 0.6, label = "strong penalty:\ntoo simple", hjust = 0, size = 4.4, colour = "#5b6b7a") +
  annotate("text", x = min(log10(cvl$lambda)), y = 0.6, label = "weak penalty:\noverfits", hjust = 1, size = 4.4, colour = "#5b6b7a") +
  scale_x_reverse() + scale_colour_manual(values = c("#2a6f97", "#c0392b"), name = NULL) +
  labs(x = expression(log[10] ~ lambda), y = "ROC AUC") + theme(legend.position = "top")
save_png(p16, "fig16-cv-curve.png", 10, 5)

# ---- the models compared on the test chromosomes ----------------------------------------
models <- list("GC" = open ~ gc, "GC + CpG" = open ~ gc + cpg_oe,
               "3 motifs" = open ~ ctcf + spi1 + ebf1, "all 5 features" = open ~ gc + cpg_oe + ctcf + spi1 + ebf1)
set.seed(7)
scores <- c(list("coin flip" = runif(nrow(tst))),
            lapply(models, function(f) predict(glm(f, binomial, train), tst)))
Xall <- cbind(X4, as.matrix(dat[, c("gc", "cpg_oe", "ctcf", "spi1", "ebf1")]))
big4 <- cached("big4", cv.glmnet(X4[!test, ], dat$open[!test], family = "binomial", alpha = 1, foldid = fold[!test]))
big  <- cached("big",  cv.glmnet(Xall[!test, ], dat$open[!test], family = "binomial", alpha = 1, foldid = fold[!test]))
scores[["256 4-mers (lasso)"]] <- as.numeric(predict(big4, X4[test, ], s = "lambda.min"))
scores[["4-mers + 5 features (lasso)"]] <- as.numeric(predict(big, Xall[test, ], s = "lambda.min"))
perf <- do.call(rbind, lapply(names(scores), function(n) {
  cv <- curves(scores[[n]], tst$open); data.frame(model = n, ROC = auc_roc(cv), PR = avg_prec(cv))
}))
print(perf, digits = 3)

# ---- fig18: the threshold ------------------------------------------------------------------
m_all <- glm(models[["all 5 features"]], binomial, train)
p_all <- predict(m_all, tst, type = "response")
print(round(summary(m_all)$coefficients, 3))
thr <- c(0.3, 0.5, 0.7)
tab <- t(sapply(thr, function(t) { call <- p_all >= t
  c(threshold = t, TPR = mean(call[tst$open == 1]), FPR = mean(call[tst$open == 0]),
    precision = mean(tst$open[call] == 1)) }))
print(round(tab, 3))
d18 <- data.frame(p = p_all, state = lab_oc(tst$open))
lab18 <- data.frame(x = thr, y = c(330, 330, 330),
                    l = sprintf("cut %.1f\nTPR %.2f\nFPR %.2f", tab[, 1], tab[, 2], tab[, 3]))
p18 <- ggplot(d18, aes(p, fill = state)) +
  geom_histogram(binwidth = 0.025, position = "identity", alpha = 0.55, boundary = 0) +
  geom_vline(xintercept = thr, linetype = 2) +
  geom_label(data = lab18, aes(x, y, label = l), inherit.aes = FALSE, size = 4, fill = "white") +
  scale_fill_manual(values = c(open = colO, closed = colC), name = "true label") +
  labs(x = "predicted P(open), all 5 features, test chromosomes", y = "windows") +
  coord_cartesian(ylim = c(0, 380))
save_png(p18, "fig18-threshold.png", 11, 5)

# ---- fig20 / fig21: ROC and PR curves -------------------------------------------------------
show <- c("coin flip", "GC", "3 motifs", "all 5 features", "4-mers + 5 features (lasso)")
cols <- c("coin flip" = "#9aa5b1", "GC" = "#e67e22", "3 motifs" = "#6b4fa0",
          "all 5 features" = "#2a6f97", "4-mers + 5 features (lasso)" = "#c0392b")
cdf <- do.call(rbind, lapply(show, function(n) {
  cv <- curves(scores[[n]], tst$open)
  data.frame(rbind(data.frame(fpr = 0, tpr = 0, precision = 1), cv[, c("fpr", "tpr", "precision")]),
             model = sprintf("%s (%.2f)", n, perf$ROC[perf$model == n]),
             model_pr = sprintf("%s (%.2f)", n, perf$PR[perf$model == n]), key = n)
}))
cm  <- setNames(cols[show], sprintf("%s (%.2f)", show, perf$ROC[match(show, perf$model)]))
cm2 <- setNames(cols[show], sprintf("%s (%.2f)", show, perf$PR[match(show, perf$model)]))
cdf$model <- factor(cdf$model, levels = names(cm)); cdf$model_pr <- factor(cdf$model_pr, levels = names(cm2))
p20 <- ggplot(cdf, aes(fpr, tpr, colour = model)) + geom_abline(linetype = 2, colour = "grey60") +
  geom_path(linewidth = 1.1) + scale_colour_manual(values = cm, name = "model (ROC AUC)") +
  coord_equal() + labs(x = "false positive rate", y = "true positive rate (recall)")
save_png(p20, "fig20-roc-real.png", 10, 5.6)
p21 <- ggplot(cdf[cdf$tpr > 0, ], aes(tpr, precision, colour = model_pr)) +
  geom_hline(yintercept = mean(tst$open), linetype = 2, colour = "grey60") +
  annotate("text", x = 0.02, y = mean(tst$open) - 0.03, label = sprintf("fraction open = %.2f", mean(tst$open)),
           hjust = 0, colour = "grey40", size = 4.3) +
  geom_path(linewidth = 1.1) + scale_colour_manual(values = cm2, name = "model (PR AUC)") +
  coord_cartesian(ylim = c(0, 1)) + labs(x = "recall", y = "precision")
save_png(p21, "fig21-pr-real.png", 10, 5.2)

# ---- fig22: class imbalance -------------------------------------------------------------------
pb <- scores[["4-mers + 5 features (lasso)"]]
set.seed(6)
keep <- which(tst$open == 0 | runif(nrow(tst)) < 0.1)
cb <- curves(pb, tst$open); ci <- curves(pb[keep], tst$open[keep])
cat("imbalanced: fraction open", round(mean(tst$open[keep]), 3), " n", length(keep),
    " ROC", round(auc_roc(ci), 3), " PR", round(avg_prec(ci), 3),
    " balanced ROC", round(auc_roc(cb), 3), " PR", round(avg_prec(cb), 3),
    " always-closed accuracy", round(mean(tst$open[keep] == 0), 3), "\n")
d22 <- rbind(data.frame(cb[, c("fpr", "tpr", "precision")], set = sprintf("46%% open (ROC %.2f, PR %.2f)", auc_roc(cb), avg_prec(cb))),
             data.frame(ci[, c("fpr", "tpr", "precision")], set = sprintf("%.0f%% open (ROC %.2f, PR %.2f)",
                                                                   100 * mean(tst$open[keep]), auc_roc(ci), avg_prec(ci))))
c22 <- setNames(c("#2a6f97", "#c0392b"), unique(d22$set))
p22a <- ggplot(d22, aes(fpr, tpr, colour = set)) + geom_abline(linetype = 2, colour = "grey60") +
  geom_path(linewidth = 1.2) + scale_colour_manual(values = c22, name = NULL) +
  coord_equal() + labs(title = "ROC: hardly changes", x = "false positive rate", y = "true positive rate") +
  theme(legend.position = "bottom", legend.direction = "vertical")
p22b <- ggplot(d22[d22$tpr > 0, ], aes(tpr, precision, colour = set)) + geom_path(linewidth = 1.2) +
  scale_colour_manual(values = c22, guide = "none") + coord_cartesian(ylim = c(0, 1)) +
  labs(title = "PR: precision collapses", x = "recall", y = "precision")
save_png(p22a + p22b, "fig22-imbalance.png", 12, 5.6)

# ---- fig23: baselines and models -----------------------------------------------------------------
perf$model <- factor(perf$model, levels = rev(perf$model))
d23 <- rbind(data.frame(model = perf$model, value = perf$ROC, metric = "ROC AUC"),
             data.frame(model = perf$model, value = perf$PR, metric = "PR AUC"))
p23 <- ggplot(d23, aes(value, model, fill = metric)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)), position = position_dodge(0.75), hjust = -0.15, size = 4.2) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  scale_fill_manual(values = c("ROC AUC" = "#2a6f97", "PR AUC" = "#e67e22"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  labs(x = "test chromosomes 8 and 9", y = NULL) + theme(legend.position = "top")
save_png(p23, "fig23-baselines.png", 11, 5.4)

# ---- fig24: what the coefficients say ------------------------------------------------------------
co <- summary(m_all)$coefficients[-1, ]
scale_by <- c(gc = 0.1, cpg_oe = 0.1, ctcf = 1, spi1 = 1, ebf1 = 1)
d24 <- data.frame(term = c("GC +0.1", "CpG o/e +0.1", "CTCF motif +1", "SPI1 (PU.1) motif +1", "EBF1 motif +1"),
                  OR = exp(co[, 1] * scale_by), lo = exp((co[, 1] - 1.96 * co[, 2]) * scale_by),
                  hi = exp((co[, 1] + 1.96 * co[, 2]) * scale_by))
print(d24, digits = 3)
d24$term <- factor(d24$term, levels = rev(d24$term))
p24 <- ggplot(d24, aes(OR, term)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 1) + geom_point(size = 3.5, colour = "#2a6f97") +
  geom_text(aes(label = sprintf("%.2f", OR)), vjust = -1, size = 4.4) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4)) +
  labs(x = "odds ratio for open (log scale), with 95% CI", y = NULL)
save_png(p24, "fig24-coefficients.png", 9.5, 4.6)
cat("done\n")
