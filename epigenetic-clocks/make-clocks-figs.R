# Real-data figures for epigenetic-clocks/lec.qmd.
# Run from this folder:  Rscript make-clocks-figs.R      (about 10 minutes, 8 GB RAM)
# Needs: data/prepare-meth2.sh clocks   (Hannum GSE40279 CpG subset and sample table,
#        Horvath 2013 coefficients, Moss 2018 reference atlas)
#        data/prepare-meth2.sh arrays   (GSE42861 IDATs and sample table: the external test set)
#        R packages data.table, glmnet, nnls, ggplot2, patchwork, minfi,
#        IlluminaHumanMethylation450kmanifest.
#
# Writes fig01, fig02, fig05, fig07, fig08, fig09, fig11, fig12 (*.png) and prints the
# numbers quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(nnls); library(ggplot2); library(patchwork)
})
set.seed(42)
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
mae <- function(p, a) median(abs(p - a))

# ---- Hannum data ----------------------------------------------------------------
b <- fread(file.path(D, "hannum_gse40279_beta_subset.tsv.gz"))
s <- fread(file.path(D, "hannum_gse40279_samples.tsv"))
X <- t(as.matrix(b[, -1])); colnames(X) <- b$cpg
stopifnot(identical(rownames(X), s$gsm), !anyNA(X))
age <- s$age
cat("Hannum:", nrow(X), "samples x", ncol(X), "CpGs; age", min(age), "-", max(age), " median", median(age), "\n")

# ---- fig01: single CpGs vs age --------------------------------------------------
r <- cor(X, age)[, 1]
neg <- names(which.min(r))
show <- c(ELOVL2 = "cg16867657", FHL2 = "cg06639320", `most negative` = neg)
cat("r: ELOVL2", round(r["cg16867657"], 2), " FHL2", round(r["cg06639320"], 2), " most negative", neg, round(r[neg], 2), "\n")
d1 <- rbindlist(lapply(seq_along(show), function(i) data.table(age = age, beta = X[, show[i]],
  cpg = sprintf("%s  %s  (r = %.2f)", show[i], names(show)[i], r[show[i]]))))
d1$cpg <- factor(d1$cpg, levels = unique(d1$cpg))
p01 <- ggplot(d1, aes(age, beta)) + geom_point(alpha = 0.35, size = 1, colour = "#2a6f97") +
  geom_smooth(method = "lm", se = FALSE, colour = "#c0392b") + facet_wrap(~ cpg, scales = "free_y") +
  labs(x = "chronological age (years)", y = "β", title = "Whole blood, 656 people (Hannum et al. 2013)")
save(p01, "fig01-age-cpgs.png", w = 13, h = 4.4)

# ---- fig02: distribution of correlations ------------------------------------------
cat("|r| > 0.5:", sum(abs(r) > 0.5), "  |r| > 0.3:", sum(abs(r) > 0.3), " of", length(r), "\n")
p02 <- ggplot(data.table(r = r), aes(r)) + geom_histogram(binwidth = 0.02, fill = "#5b6b7a") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, colour = "#c0392b") +
  labs(x = "correlation of β with age", y = "CpGs",
       title = sprintf("%s CpGs: %d with |r| > 0.5, %d with |r| > 0.3", format(length(r), big.mark = " "),
                       sum(abs(r) > 0.5), sum(abs(r) > 0.3)))
save(p02, "fig02-age-correlations.png", w = 10, h = 4)

# ---- split, cross-validated elastic net ------------------------------------------
test <- sample(nrow(X), round(0.3 * nrow(X))); train <- setdiff(seq_len(nrow(X)), test)
cat("train", length(train), " test", length(test), "\n")
cv <- cv.glmnet(X[train, ], age[train], alpha = 0.5, nfolds = 10)
cf <- coef(cv, s = "lambda.min")[, 1]; sel <- names(cf)[cf != 0][-1]
cat("lambda.min", signif(cv$lambda.min, 3), " CpGs selected", length(sel), " (lambda.1se:", sum(coef(cv, s = "lambda.1se")[-1] != 0), ")\n")
cat("  of which Horvath CpGs:", sum(sel %in% fread(file.path(D, "horvath2013_coefficients.csv"), skip = 2)$CpGmarker), "\n")
ptr <- predict(cv, X[train, ], s = "lambda.min")[, 1]; pte <- predict(cv, X[test, ], s = "lambda.min")[, 1]
cat("train: r", round(cor(ptr, age[train]), 3), " median abs error", round(mae(ptr, age[train]), 2), "\n")
cat("test:  r", round(cor(pte, age[test]), 3), " median abs error", round(mae(pte, age[test]), 2), "\n")

# ---- fig05: CV curve and coefficient path -----------------------------------------
cvd <- data.table(loglam = log(cv$lambda), mse = cv$cvm, lo = cv$cvlo, hi = cv$cvup, nz = cv$nzero)
p05a <- ggplot(cvd, aes(loglam, mse)) + geom_errorbar(aes(ymin = lo, ymax = hi), colour = "#a9b5c0", width = 0) +
  geom_point(colour = "#c0392b", size = 1.4) + geom_vline(xintercept = log(cv$lambda.min), linetype = 2) +
  annotate("text", x = log(cv$lambda.min), y = max(cvd$hi), hjust = -0.05, vjust = 1, size = 4.3,
           label = sprintf("lambda.min:\n%d CpGs", length(sel))) +
  labs(x = "log(λ)  (left: weak penalty)", y = "10-fold CV error (MSE, years²)", title = "Choosing λ")
fit <- cv$glmnet.fit
top <- names(sort(abs(cf[-1]), decreasing = TRUE))[1:12]
pth <- rbindlist(lapply(top, function(g) data.table(loglam = log(fit$lambda), coef = fit$beta[g, ], cpg = g)))
p05b <- ggplot(pth, aes(loglam, coef, group = cpg)) + geom_line(colour = "#2a6f97", alpha = 0.8) +
  geom_vline(xintercept = log(cv$lambda.min), linetype = 2) +
  labs(x = "log(λ)", y = "coefficient (years per unit β)", title = "12 largest coefficients along the path")
save(p05a + p05b, "fig05-cv-path.png", w = 13, h = 4.6)

# ---- fig07: train vs test, and overfitting with ordinary least squares --------------
d7 <- rbind(data.table(age = age[train], pred = ptr, set = sprintf("training (n = %d): median error %.1f y", length(train), mae(ptr, age[train]))),
            data.table(age = age[test], pred = pte, set = sprintf("test (n = %d): median error %.1f y", length(test), mae(pte, age[test]))))
d7[, set := factor(set, levels = unique(set))]
p07a <- ggplot(d7, aes(age, pred)) + geom_abline(linetype = 2, colour = "#5b6b7a") +
  geom_point(alpha = 0.5, size = 1.2, colour = "#2a6f97") + facet_wrap(~ set) + coord_equal() +
  labs(x = "chronological age", y = "predicted age", title = "Elastic net, Hannum data")
rt <- cor(X[train, ], age[train])[, 1]
ks <- c(5, 10, 25, 50, 100, 200, 300, 400, 440)
ols <- rbindlist(lapply(ks, function(k) {
  g <- names(sort(abs(rt), decreasing = TRUE))[1:k]
  f <- lm.fit(cbind(1, X[train, g]), age[train]); bb <- f$coefficients; bb[is.na(bb)] <- 0
  data.table(k = k, train = mae(cbind(1, X[train, g]) %*% bb, age[train]), test = mae(cbind(1, X[test, g]) %*% bb, age[test]))
}))
print(ols)
ol <- melt(ols, id.vars = "k", variable.name = "set", value.name = "err")
p07b <- ggplot(ol, aes(k, err, colour = set)) + geom_line(linewidth = 1.1) + geom_point() +
  scale_colour_manual(values = c(train = "#2a6f97", test = "#c0392b"), name = NULL) +
  labs(x = "number of CpGs in an unpenalised linear model", y = "median abs error (years)", title = "Least squares overfits") +
  theme(legend.position = "top")
save(p07a + p07b + plot_layout(widths = c(1.7, 1)), "fig07-train-test.png", w = 13.5, h = 5)

# ---- Horvath clock -----------------------------------------------------------------
hc <- fread(file.path(D, "horvath2013_coefficients.csv"), skip = 2)[, .(cpg = CpGmarker, w = CoefficientTraining)]
h0 <- hc[cpg == "(Intercept)", w]; hc <- hc[cpg != "(Intercept)"]
anti <- function(x, adult = 20) ifelse(x < 0, (1 + adult) * exp(x) - 1, (1 + adult) * x + adult)
trafo <- function(a, adult = 20) ifelse(a <= adult, log(a + 1) - log(adult + 1), (a - adult) / (adult + 1))
horvath <- function(B) anti(h0 + B[, hc$cpg, drop = FALSE] %*% hc$w)[, 1]
hz <- horvath(X)
cat("Horvath on Hannum (in his training data): r", round(cor(hz, age), 3), " median abs error", round(mae(hz, age), 2),
    " median offset", round(median(hz - age), 2), "\n")
p08a <- ggplot(data.table(a = seq(0, 100, 0.5))[, f := trafo(a)], aes(a, f)) + geom_line(linewidth = 1.2, colour = "#2a6f97") +
  geom_vline(xintercept = 20, linetype = 2, colour = "#5b6b7a") +
  annotate("text", x = 22, y = -2.5, hjust = 0, label = "log below 20 years,\nlinear above", size = 4.3) +
  labs(x = "chronological age", y = "F(age), the model's target", title = "Horvath's age transformation")
p08b <- ggplot(data.table(age, hz), aes(age, hz)) + geom_abline(linetype = 2, colour = "#5b6b7a") +
  geom_point(alpha = 0.4, size = 1, colour = "#e67e22") +
  labs(x = "chronological age", y = "Horvath DNAm age", title = "Horvath on Hannum (in-sample)",
       subtitle = sprintf("r = %.2f, median error %.1f y", cor(hz, age), mae(hz, age)))
save(p08a + p08b, "fig08-horvath.png", w = 12, h = 4.8)

# ---- GSE42861: an external test set (normalised from IDATs with noob) ----------------
suppressPackageStartupMessages(library(minfi))
g <- fread(file.path(D, "gse42861_samples.tsv"))
g$Basename <- file.path(D, "gse42861_idat", g$basename)
need <- unique(c(colnames(X), hc$cpg, fread(file.path(D, "moss2018_reference_atlas.csv"))$CpGs))
cache <- file.path(tempdir(), "gse42861_noob.rds")
if (!file.exists(cache)) {
  rg <- read.metharray.exp(targets = as.data.frame(g))
  bg <- getBeta(preprocessNoob(rg)); rm(rg); invisible(gc())
  saveRDS(bg[intersect(need, rownames(bg)), ], cache)
}
bg <- readRDS(cache)
Xg <- t(bg[colnames(X), ]); stopifnot(identical(rownames(Xg), g$basename))
if (anyNA(Xg)) { mu <- colMeans(Xg, na.rm = TRUE); ii <- which(is.na(Xg), arr.ind = TRUE); Xg[ii] <- mu[ii[, 2]] }
pg <- predict(cv, Xg, s = "lambda.min")[, 1]
hg <- horvath(t(bg[hc$cpg, ]))
cat("GSE42861 (n =", nrow(g), ", age", min(g$age), "-", max(g$age), "):\n")
cat("  our clock: r", round(cor(pg, g$age), 3), " median abs error", round(mae(pg, g$age), 2), " median offset", round(median(pg - g$age), 2), "\n")
cat("  Horvath:   r", round(cor(hg, g$age), 3), " median abs error", round(mae(hg, g$age), 2), " median offset", round(median(hg - g$age), 2), "\n")
d9 <- rbind(data.table(age = g$age, pred = pg, clock = sprintf("our clock: r = %.2f, median error %.1f y", cor(pg, g$age), mae(pg, g$age))),
            data.table(age = g$age, pred = hg, clock = sprintf("Horvath 2013: r = %.2f, median error %.1f y", cor(hg, g$age), mae(hg, g$age))))
d9[, clock := factor(clock, levels = unique(clock))]
p09 <- ggplot(d9, aes(age, pred)) + geom_abline(linetype = 2, colour = "#5b6b7a") +
  geom_point(alpha = 0.6, size = 1.4, colour = "#2e8b57") + geom_smooth(method = "lm", se = FALSE, colour = "#c0392b", linewidth = 0.8) +
  facet_wrap(~ clock) +
  labs(x = "chronological age", y = "predicted age", title = "External test: GSE42861, 151 blood samples from another lab and cohort")
save(p09, "fig09-external.png", w = 12, h = 5.2)

# ---- fig11: age acceleration = residual ---------------------------------------------
fitg <- lm(pg ~ g$age); acc <- residuals(fitg)
cat("age acceleration (our clock, GSE42861): sd", round(sd(acc), 2), " slope of DNAm age on age", round(coef(fitg)[2], 2), "\n")
tt <- t.test(acc ~ g$status); cat("  RA vs control: mean", round(tt$estimate, 2), " p", signif(tt$p.value, 2), "\n")
ttr <- t.test((pg - g$age) ~ g$status); cat("  raw difference RA vs control: mean", round(ttr$estimate, 2), " p", signif(ttr$p.value, 2), "\n")
ts <- t.test(acc[g$smoking %in% c("current", "never")] ~ g$smoking[g$smoking %in% c("current", "never")])
cat("  current vs never smokers: mean", round(ts$estimate, 2), " p", signif(ts$p.value, 2), "\n")
d11 <- data.table(age = g$age, pred = pg, acc = acc, status = g$status, fitted = fitted(fitg))
i1 <- which.max(acc); i2 <- which.min(acc)
p11a <- ggplot(d11, aes(age, pred)) + geom_abline(linetype = 2, colour = "#5b6b7a") +
  geom_smooth(method = "lm", se = FALSE, colour = "#c0392b") +
  geom_segment(data = d11[c(i1, i2)], aes(xend = age, yend = fitted), colour = "#c0392b", linewidth = 0.9) +
  geom_point(aes(colour = status), alpha = 0.7) +
  scale_colour_manual(values = c(control = "#5b6b7a", RA = "#e67e22"), name = NULL) +
  labs(x = "chronological age", y = "DNAm age (our clock)", title = "Acceleration = vertical distance to the red line") +
  theme(legend.position = "top")
p11b <- ggplot(d11, aes(status, acc, fill = status)) + geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.6) + geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c(control = "#5b6b7a", RA = "#e67e22"), guide = "none") +
  labs(x = NULL, y = "age acceleration (years)", title = sprintf("RA vs control: p = %.2f", tt$p.value))
save(p11a + p11b + plot_layout(widths = c(1.6, 1)), "fig11-age-accel.png", w = 12.5, h = 5)

# ---- fig12: deconvolution with the Moss 2018 atlas -----------------------------------
atl <- fread(file.path(D, "moss2018_reference_atlas.csv"))
cat("atlas rows", nrow(atl), "; unique CpGs", length(unique(atl$CpGs)), "\n")
atl <- atl[!duplicated(CpGs)]                     # markers chosen for two cell types appear twice
R <- as.matrix(atl[, -1]); rownames(R) <- atl$CpGs
common <- intersect(rownames(R), rownames(bg)); R <- R[common, ]
cat("atlas:", ncol(R), "cell types;", length(common), "of its CpGs on the 450k array\n")
deconv <- function(y, R) { ok <- !is.na(y); f <- nnls(R[ok, ], y[ok])$x; names(f) <- colnames(R); f / sum(f) }
W <- t(apply(bg[common, ], 2, deconv, R = R))
mw <- sort(colMeans(W), decreasing = TRUE); print(round(head(mw, 8), 3))
# in-silico cfDNA: a blood sample plus a known fraction of hepatocyte profile (+ noise)
liv <- R[, "Hepatocytes"]; base <- bg[common, 1]
spk <- rbindlist(lapply(c(0, 0.01, 0.02, 0.05, 0.1, 0.2), function(f) rbindlist(lapply(1:20, function(rep) {
  y <- (1 - f) * bg[common, sample(ncol(bg), 1)] + f * liv + rnorm(length(common), 0, 0.03)
  data.table(true = f, est = deconv(pmin(pmax(y, 0), 1), R)["Hepatocytes"])
}))))
print(spk[, .(mean_est = round(mean(est), 4), sd = round(sd(est), 4)), by = true])
wl <- melt(data.table(sample = seq_len(nrow(W)), W[, names(mw)[1:6]], other = 1 - rowSums(W[, names(mw)[1:6]])),
           id.vars = "sample", variable.name = "cell", value.name = "f")
wl[, cell := factor(gsub("_EPIC", "", cell), levels = c(gsub("_EPIC", "", names(mw)[1:6]), "other"))]
ord <- order(W[, names(mw)[1]]); wl[, sample := factor(sample, levels = ord)]
p12a <- ggplot(wl, aes(sample, f, fill = cell)) + geom_col(width = 1) +
  scale_fill_manual(values = c("#e67e22", "#2a6f97", "#2e8b57", "#6b4fa0", "#d4a017", "#8c564b", "#dfe6ec"), name = NULL) +
  labs(x = "151 whole-blood samples (GSE42861)", y = "estimated fraction", title = "Atlas deconvolution of blood") +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())
p12b <- ggplot(spk, aes(factor(100 * true), 100 * est)) + geom_boxplot(fill = "#e67e22", alpha = 0.5, outlier.size = 0.8) +
  geom_point(data = data.table(x = factor(100 * unique(spk$true)), y = 100 * unique(spk$true)), aes(x, y), colour = "#c0392b", shape = 4, size = 4, stroke = 1.5) +
  labs(x = "hepatocyte DNA added (%)", y = "estimated hepatocyte fraction (%)", title = "Simulated liver damage")
save(p12a + p12b + plot_layout(widths = c(1.5, 1)), "fig12-deconvolution.png", w = 13, h = 4.8)
cat("done\n")
