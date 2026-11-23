# Simulation figures for epigenetic-inheritance/lec.qmd.
# Run from this folder:  Rscript make-inheritance-figs.R   (about 2 minutes)
# Needs no data: every figure is a seeded simulation of a study design in which
# there is NO true inheritance (or, for fig19, a small true effect), to show how
# false "transgenerational" signals arise.
# R packages: ggplot2, patchwork.
#
# Writes fig15, fig16, fig17, fig18, fig19 (*.png) and prints the numbers
# quoted in the slides.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
rowt <- function(M, g) {                 # Welch t-test p-value for every row of M, groups g (logical)
  a <- M[, g]; b <- M[, !g]
  ma <- rowMeans(a); mb <- rowMeans(b)
  va <- apply(a, 1, var) / ncol(a); vb <- apply(b, 1, var) / ncol(b)
  t <- (ma - mb) / sqrt(va + vb)
  df <- (va + vb)^2 / (va^2 / (ncol(a) - 1) + vb^2 / (ncol(b) - 1))
  2 * pt(-abs(t), df)
}

# ---- fig16: 20 000 CpGs, no effect, three generations ----------------------------
set.seed(1)
n_cpg <- 20000; n <- 8
mu <- rbeta(n_cpg, 2, 2)
g <- rep(c(TRUE, FALSE), each = n)
p <- sapply(c("F1", "F2", "F3"), function(gen)
  rowt(matrix(rnorm(n_cpg * 2 * n, mean = mu, sd = 0.05), nrow = n_cpg), g))
hits <- colSums(p < 0.05); all3 <- sum(rowSums(p < 0.05) == 3)
fdr <- colSums(apply(p, 2, p.adjust, method = "BH") < 0.05)
cat("fig16: p<0.05 per generation", hits, " in all three:", all3,
    " expected", n_cpg * 0.05^3, " FDR<0.05:", fdr, "\n")
d17 <- data.frame(p = as.vector(p), gen = rep(colnames(p), each = n_cpg))
pA <- ggplot(d17, aes(p)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "#2a6f97", colour = "white") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  facet_wrap(~gen) + geom_hline(yintercept = n_cpg * 0.05, linetype = 2, colour = "#c0392b") +
  labs(x = "p-value, exposed vs control (8 vs 8 animals)", y = "CpGs",
       title = sprintf("No true effect: about %d CpGs with p < 0.05 in every generation", round(mean(hits))))
ov <- data.frame(k = factor(0:3), n = as.vector(table(factor(rowSums(p < 0.05), levels = 0:3))))
pB <- ggplot(ov[-1, ], aes(k, n)) + geom_col(fill = "#c0392b", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.3, size = 5) + scale_y_log10(limits = c(1, 1e4)) +
  labs(x = "generations with p < 0.05", y = "CpGs (log scale)", title = "'Persisting' hits")
save(pA + pB + plot_layout(widths = c(2.2, 1)), "fig16-multiple-testing.png", h = 4.4)

# ---- fig17: litter effects ----------------------------------------------------------
set.seed(2)
one_study <- function(n_litters = 4, pups = 6, sd_litter = 1, sd_pup = 1) {
  litter <- rep(1:(2 * n_litters), each = pups)
  group  <- rep(c("control", "exposed"), each = n_litters * pups)
  y <- rnorm(2 * n_litters, 0, sd_litter)[litter] + rnorm(length(litter), 0, sd_pup)
  m <- tapply(y, litter, mean); gl <- tapply(group, litter, `[`, 1)
  c(pups = t.test(y[group == "exposed"], y[group == "control"])$p.value,
    litters = t.test(m[gl == "exposed"], m[gl == "control"])$p.value)
}
sds <- c(0, 0.25, 0.5, 0.75, 1, 1.5)
fp <- t(sapply(sds, function(s) rowMeans(replicate(2000, one_study(sd_litter = s)) < 0.05)))
d18 <- data.frame(sd = rep(sds, 2), fpr = as.vector(fp),
                  unit = rep(c("each pup counted as independent (n = 24 vs 24)",
                               "litter means (n = 4 vs 4)"), each = length(sds)))
print(cbind(sd_litter = sds, round(fp, 3)))
p18 <- ggplot(d18, aes(sd, fpr, colour = unit)) + geom_hline(yintercept = 0.05, linetype = 2) +
  geom_line(linewidth = 1.3) + geom_point(size = 2.5) +
  scale_colour_manual(values = c("#c0392b", "#2a6f97"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = "SD between litters (SD within a litter = 1)", y = "false-positive rate",
       title = "4 exposed and 4 control litters of 6 pups, no true effect") +
  theme(legend.position = c(0.35, 0.85))
save(p18, "fig17-litter-effects.png", w = 10, h = 4.6)

# ---- fig18: genetic confounding (an meQTL that differs between regions) -------------
set.seed(3)
n <- 400
region  <- rep(c("north (famine)", "south"), each = n / 2)
snp     <- rbinom(n, 2, ifelse(region == "south", 0.25, 0.45))
meth    <- 0.50 + 0.04 * snp + rnorm(n, 0, 0.05)
exposed <- region == "north (famine)"
c1 <- summary(lm(meth ~ exposed))$coefficients; c2 <- summary(lm(meth ~ exposed + snp))$coefficients
cat("fig18: naive effect", round(c1[2, 1], 4), "p", signif(c1[2, 4], 2),
    "| adjusted for genotype", round(c2[2, 1], 4), "p", signif(c2[2, 4], 2), "\n")
d19 <- data.frame(region, snp = factor(snp), meth)
pA <- ggplot(d19, aes(region, meth)) + geom_jitter(aes(colour = snp), width = 0.2, size = 1.4, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "#1b2a3a") +
  scale_colour_manual(values = c("#9ecae1", "#2a6f97", "#1b2a3a"), name = "risk alleles") +
  labs(x = "grandchildren from", y = "methylation (beta)",
       title = sprintf("Naive: famine region +%.3f (p = %.3f)", c1[2, 1], c1[2, 4]))
pB <- ggplot(d19, aes(snp, meth, colour = region)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.6) +
  scale_colour_manual(values = c("#c0392b", "#5b6b7a"), name = NULL) +
  labs(x = "genotype (risk alleles)", y = "methylation (beta)",
       title = sprintf("Within genotype: +%.3f (p = %.2f)", c2[2, 1], c2[2, 4]))
save(pA + pB, "fig18-genetic-confounding.png", h = 4.6)

# ---- fig15: forking paths in a small historical cohort ------------------------------
set.seed(4)
one_cohort <- function(n = 300) {
  food <- matrix(rbinom(n * 4, 1, 0.3), n, 4, dimnames = list(NULL, c("PGF", "PGM", "MGF", "MGM")))
  sex  <- sample(c("M", "F"), n, replace = TRUE)
  outc <- cbind(allcause = rbinom(n, 1, 0.5), cvd = rbinom(n, 1, 0.2), diabetes = rbinom(n, 1, 0.05))
  p <- c()
  for (gp in colnames(food)) for (o in colnames(outc)) for (s in c("M", "F", "both")) {
    keep <- if (s == "both") rep(TRUE, n) else sex == s
    tab <- table(factor(food[keep, gp], 0:1), factor(outc[keep, o], 0:1))
    p[paste(gp, o, s)] <- fisher.test(tab)$p.value
  }
  p
}
res <- replicate(1000, one_cohort())
anyhit <- mean(apply(res, 2, min) < 0.05)
cat("fig15: tests per cohort", nrow(res), " cohorts with >=1 p<0.05:", anyhit,
    " mean number of hits", mean(colSums(res < 0.05)), "\n")
d20 <- data.frame(h = colSums(res < 0.05))
p20 <- ggplot(d20, aes(factor(h))) + geom_bar(fill = "#2a6f97", width = 0.7) +
  labs(x = sprintf("tests with p < 0.05 (out of %d) in one null cohort", nrow(res)),
       y = "simulated cohorts (of 1000)",
       title = sprintf("No true effect, 36 tests per cohort: %.0f%% of cohorts 'find' something", 100 * anyhit))
save(p20, "fig15-forking-paths.png", w = 11, h = 4.4)

# ---- fig19: a real but decaying effect, small samples -------------------------------
set.seed(5)
eff <- c(F1 = 1, F2 = 0.5, F3 = 0.25)
sim <- lapply(names(eff), function(gen) {
  r <- replicate(4000, { x <- rnorm(10, eff[gen]); y <- rnorm(10); c(p = t.test(x, y)$p.value, d = mean(x) - mean(y)) })
  data.frame(gen = gen, true = unname(eff[gen]), p = r["p", ], d = r["d", ])
})
sim <- do.call(rbind, sim)
pw <- tapply(sim$p < 0.05, sim$gen, mean)
infl <- tapply(sim$d[sim$p < 0.05], sim$gen[sim$p < 0.05], mean)
print(round(rbind(true = eff, power = pw[names(eff)], mean_significant_estimate = infl[names(eff)]), 3))
pA <- ggplot(data.frame(gen = names(eff), power = pw[names(eff)]), aes(gen, power)) +
  geom_col(fill = "#2a6f97", width = 0.6) + geom_text(aes(label = sprintf("%.0f%%", 100 * power)), vjust = -0.3, size = 5) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 1)) +
  labs(x = "true effect 1, 0.5, 0.25 SD", y = "power (n = 10 vs 10)", title = "Power")
pB <- ggplot(sim[sim$gen == "F3", ], aes(d, fill = p < 0.05)) + geom_histogram(bins = 50, colour = "white") +
  geom_vline(xintercept = 0.25, linetype = 2) +
  scale_fill_manual(values = c(`FALSE` = "#c9d3dd", `TRUE` = "#c0392b"), labels = c("not significant", "p < 0.05"), name = NULL) +
  labs(x = "estimated effect in F3 (SD units; true = 0.25)", y = "simulated studies",
       title = sprintf("F3: significant studies report %.2f on average", infl["F3"])) +
  theme(legend.position = c(0.82, 0.85))
save(pA + pB + plot_layout(widths = c(1, 2)), "fig19-winners-curse.png", h = 4.6)
