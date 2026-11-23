# Real-data figures for cancer-epigenetics/lec.qmd.
# Run from this folder:  Rscript make-cancer-figs.R   (about 2 minutes, 3 GB RAM)
# Needs: data/prepare-dis.sh (TCGA colon cancer and glioma 450k subsets,
#        sample tables, probe annotation)
#        R packages data.table, ggplot2, patchwork, glmnet.
#
# Writes fig02, fig05, fig07, fig09, fig10, fig14, fig22, fig24 (*.png) and
# prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(glmnet)
})
set.seed(1)
D <- "../data"
colTN <- c(normal = "#2e8b57", tumour = "#a8324a")
colCIMP <- c(`CIMP-high` = "#c0392b", `CIMP-low` = "#e9c46a", `non-CIMP` = "#2a6f97")
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")

# ---- load ---------------------------------------------------------------------
beta <- fread(file.path(D, "tcga_coad_450k_beta.tsv.gz"))
B <- as.matrix(beta[, -1]); rownames(B) <- beta$probe; rm(beta)
anno <- fread(file.path(D, "tcga_450k_probes.tsv.gz"))
probes <- anno[match(rownames(B), probe)]
samples <- fread(file.path(D, "tcga_coad_samples.tsv"))
stopifnot(identical(samples$sample, colnames(B)))
tum <- samples$type == "tumour"
cat("COAD:", nrow(B), "probes,", sum(tum), "tumours,", sum(!tum), "normals\n")

# ---- fig02: global loss, focal gain ---------------------------------------------
ctx <- c("island", "shore", "shelf", "open sea")
mb <- sapply(ctx, function(k) colMeans(B[probes$cgi == k, ], na.rm = TRUE))
print(round(apply(mb, 2, function(v) tapply(v, samples$type, median)), 3))
df <- data.frame(samples[, .(sample, patient, type)], island = mb[, "island"], opensea = mb[, "open sea"])
pA <- ggplot(df, aes(opensea, island, colour = type)) + geom_point(size = 2, alpha = 0.75) +
  scale_colour_manual(values = colTN, name = NULL) +
  labs(x = "mean beta, open-sea CpGs", y = "mean beta, CpG-island CpGs",
       title = "Each dot is one sample") +
  theme(legend.position = c(0.15, 0.88))
paired <- intersect(samples[type == "normal", patient], samples[type == "tumour", patient])
pl <- rbind(data.frame(patient = paired, type = "normal", ctx = "open sea", v = mb[paste0(paired, "-11"), "open sea"]),
            data.frame(patient = paired, type = "tumour", ctx = "open sea", v = mb[paste0(paired, "-01"), "open sea"]),
            data.frame(patient = paired, type = "normal", ctx = "CpG island", v = mb[paste0(paired, "-11"), "island"]),
            data.frame(patient = paired, type = "tumour", ctx = "CpG island", v = mb[paste0(paired, "-01"), "island"]))
pl$ctx <- factor(pl$ctx, levels = c("open sea", "CpG island"))
d_os  <- mb[paste0(paired, "-01"), "open sea"] - mb[paste0(paired, "-11"), "open sea"]
d_isl <- mb[paste0(paired, "-01"), "island"] - mb[paste0(paired, "-11"), "island"]
cat("paired patients:", length(paired), " open sea lower in tumour:", sum(d_os < 0),
    " median change", round(median(d_os), 3), "\n")
cat("  island higher in tumour:", sum(d_isl > 0), " median change", round(median(d_isl), 3), "\n")
pB <- ggplot(pl, aes(type, v, group = patient)) + geom_line(colour = "grey60", alpha = 0.6) +
  geom_point(aes(colour = type), size = 1.8) + facet_wrap(~ctx, scales = "free_y") +
  scale_colour_manual(values = colTN, guide = "none") +
  labs(x = NULL, y = "mean beta", title = sprintf("Same patient, normal vs tumour (n = %d)", length(paired)))
save(pA + pB + plot_layout(widths = c(1, 1.3)), "fig02-global-methylation.png", h = 4.8)

# ---- fig05: CpG islands of single genes ------------------------------------------
in_gene <- function(g) grepl(paste0("(^|;)", g, "(;|$)"), probes$genes) &
  probes$cgi == "island" & abs(probes$tss_dist) <= 1500
genes <- c("SFRP2", "GATA4", "SEPTIN9", "MGMT", "MLH1", "CDKN2A", "GAPDH")
gm <- sapply(genes, function(g) colMeans(B[in_gene(g), , drop = FALSE], na.rm = TRUE))
print(round(rbind(n_probes = sapply(genes, function(g) sum(in_gene(g))),
                  normal = apply(gm[!tum, ], 2, median), tumour = apply(gm[tum, ], 2, median),
                  frac_tum_gt0.3 = colMeans(gm[tum, ] > 0.3)), 3))
gl <- data.frame(type = rep(samples$type, length(genes)), gene = rep(genes, each = nrow(gm)), v = as.vector(gm))
gl$gene <- factor(gl$gene, levels = genes)
p05 <- ggplot(gl, aes(type, v, colour = type)) +
  geom_jitter(width = 0.25, size = 0.9, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#1b2a3a", width = 0.5) +
  facet_wrap(~gene, nrow = 1) + scale_colour_manual(values = colTN, guide = "none") +
  geom_hline(yintercept = 0.3, linetype = 2, colour = "grey50") +
  labs(x = NULL, y = "promoter CGI methylation\n(mean beta)",
       title = "Promoter CpG islands: gained in many tumours, but not at every gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save(p05, "fig05-gene-panel.png", h = 4.4)

# ---- fig07: MLH1 methylation, expression and MSI ---------------------------------
st <- samples[tum]; st$mlh1 <- gm[tum, "MLH1"]
print(table(MLH1_methylated = st$mlh1 > 0.3, MSI = st$msi, useNA = "ifany"))
print(tapply(log2(st$mlh1_rna + 1), st$mlh1 > 0.3, median, na.rm = TRUE))
st$msi_lab <- ifelse(is.na(st$msi), "unknown", st$msi)
p07 <- ggplot(st[!is.na(mlh1_rna)], aes(mlh1, log2(mlh1_rna + 1), colour = msi_lab)) +
  geom_vline(xintercept = 0.3, linetype = 2, colour = "grey50") +
  geom_point(size = 2.2, alpha = 0.8) +
  scale_colour_manual(values = c(`MSI-H` = "#c0392b", MSS = "#5b6b7a", unknown = "#d5dde5"), name = NULL) +
  labs(x = "MLH1 promoter CGI methylation (mean beta)", y = "MLH1 RNA, log2(RSEM + 1)",
       title = "297 colon tumours: MLH1 methylation silences the gene") +
  theme(legend.position = "right")
save(p07, "fig07-mlh1-msi.png", w = 9, h = 4.8)

# ---- fig09 / fig10: CIMP -----------------------------------------------------------
nm   <- rowMeans(B[, !tum], na.rm = TRUE)
cand <- which(probes$cgi == "island" & nm < 0.2 & rowSums(is.na(B[, tum])) == 0)
sdv  <- apply(B[cand, tum], 1, sd, na.rm = TRUE)
top  <- cand[order(sdv, decreasing = TRUE)[1:1000]]
X <- B[top, tum]
hc <- hclust(dist(t(X)), method = "ward.D2"); cl <- cutree(hc, k = 3)
cimp <- c("CIMP-high", "CIMP-low", "non-CIMP")[rank(-tapply(colMeans(X), cl, mean))][cl]
st$cimp <- factor(cimp, levels = names(colCIMP))
cat("candidate island CpGs (normal < 0.2):", length(cand), "\n"); print(table(st$cimp))
print(table(st$cimp, st$braf)); print(table(st$cimp, st$msi)); print(table(st$cimp, st$side))
print(table(st$cimp, MLH1_methylated = st$mlh1 > 0.3))
so <- order(st$cimp, -colMeans(X))
ro <- hclust(dist(X), method = "ward.D2")$order
hm <- data.frame(x = rep(seq_along(so), each = nrow(X)), y = rep(seq_len(nrow(X)), length(so)),
                 v = as.vector(X[ro, so]))
p09a <- ggplot(hm, aes(x, y, fill = v)) + geom_raster() +
  scale_fill_gradient(low = "#f4f6f8", high = "#1b2a3a", name = "beta", limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "1000 CGI CpGs\n(unmethylated in normal colon)") +
  theme(axis.text = element_blank(), panel.grid = element_blank())
ann <- data.frame(x = seq_along(so), cimp = st$cimp[so],
                  BRAF = ifelse(st$braf[so] %in% "V600E", "V600E", "no"),
                  MSI = ifelse(st$msi[so] %in% "MSI-H", "MSI-H", "no"),
                  side = ifelse(st$side[so] %in% "right", "right", "no"))
al <- rbind(data.frame(x = ann$x, row = "CIMP group", lab = as.character(ann$cimp)),
            data.frame(x = ann$x, row = "BRAF V600E", lab = ifelse(ann$BRAF == "V600E", "yes", "no")),
            data.frame(x = ann$x, row = "MSI-high", lab = ifelse(ann$MSI == "MSI-H", "yes", "no")),
            data.frame(x = ann$x, row = "right colon", lab = ifelse(ann$side == "right", "yes", "no")))
al$row <- factor(al$row, levels = rev(c("CIMP group", "BRAF V600E", "MSI-high", "right colon")))
p09b <- ggplot(al, aes(x, row, fill = lab)) + geom_tile() +
  scale_fill_manual(values = c(colCIMP, yes = "#1b2a3a", no = "#eaf1f7"), breaks = c(names(colCIMP), "yes", "no"),
                    labels = c(names(colCIMP), "yes", "no / unknown"), name = NULL) +
  scale_x_continuous(expand = c(0, 0)) + labs(x = sprintf("%d colon tumours", length(so)), y = NULL) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
save(p09a / p09b + plot_layout(heights = c(3, 1)), "fig09-cimp-heatmap.png", w = 11, h = 5.6)

fr <- st[, .(`BRAF V600E` = mean(braf %in% "V600E" & !is.na(braf)) / mean(!is.na(braf)),
             `MSI-high` = mean(msi %in% "MSI-H") / mean(!is.na(msi)),
             `right colon` = mean(side %in% "right") / mean(!is.na(side)),
             `MLH1 methylated` = mean(mlh1 > 0.3), n = .N), by = cimp]
print(fr)
frl <- melt(fr[, -"n"], id.vars = "cimp")
p10 <- ggplot(frl, aes(cimp, value, fill = cimp)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * value)), vjust = -0.3, size = 4.5) +
  facet_wrap(~variable, nrow = 1) + scale_fill_manual(values = colCIMP, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 1.05)) +
  labs(x = NULL, y = "fraction of tumours", title = "What goes with CIMP-high") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save(p10, "fig10-cimp-associations.png", h = 4.6)

# ---- fig14: G-CIMP in IDH-mutant glioma -------------------------------------------
gb <- fread(file.path(D, "tcga_lgg_450k_beta.tsv.gz"))
G <- as.matrix(gb[, -1]); rownames(G) <- gb$probe; rm(gb)
gs <- fread(file.path(D, "tcga_lgg_samples.tsv")); stopifnot(identical(gs$sample, colnames(G)))
gp <- anno[match(rownames(G), probe)]
wt_low <- gp$cgi == "island" & rowMeans(G[, gs$idh == "WT"], na.rm = TRUE) < 0.2
fh <- colMeans(G[wt_low, ] > 0.3, na.rm = TRUE)
cat("LGG: gliomas", ncol(G), table(gs$idh), " CGI CpGs unmethylated in IDH-wt:", sum(wt_low), "\n")
print(round(tapply(fh, gs$idh, median), 3))
isl <- gp$cgi == "island"
dm <- rowMeans(G[isl, gs$idh == "Mutant"], na.rm = TRUE) - rowMeans(G[isl, gs$idh == "WT"], na.rm = TRUE)
cat("island CpGs >0.2 higher in IDH-mutant:", sum(dm > 0.2), " lower:", sum(dm < -0.2), " of", sum(isl), "\n")
cc <- rowSums(is.na(G)) == 0
sel <- which(isl & cc)[order(-abs(dm[cc[isl]]))[1:400]]
og <- order(gs$idh, -colMeans(G[sel, ], na.rm = TRUE))
Y <- G[sel, og]
Y <- Y[hclust(dist(Y))$order, ]
hm2 <- data.frame(x = rep(seq_along(og), each = nrow(Y)), y = rep(seq_len(nrow(Y)), length(og)), v = as.vector(Y))
p14a <- ggplot(hm2, aes(x, y, fill = v)) + geom_raster() +
  scale_fill_gradient(low = "#f4f6f8", high = "#1b2a3a", name = "beta", limits = c(0, 1)) +
  geom_vline(xintercept = sum(gs$idh == "Mutant") + 0.5, colour = "#c0392b", linewidth = 1) +
  annotate("text", x = sum(gs$idh == "Mutant") / 2, y = nrow(Y) * 1.08, size = 5,
           label = sprintf("IDH-mutant\n(n = %d)", sum(gs$idh == "Mutant"))) +
  annotate("text", x = sum(gs$idh == "Mutant") + sum(gs$idh == "WT") / 2, y = nrow(Y) * 1.08, size = 5,
           label = sprintf("IDH-wt\n(n = %d)", sum(gs$idh == "WT"))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, nrow(Y) * 1.16)) +
  labs(x = "513 lower-grade gliomas (TCGA)", y = "400 CGI CpGs that differ most") +
  theme(axis.text = element_blank(), panel.grid = element_blank())
fd <- data.frame(idh = paste0("IDH-", ifelse(gs$idh == "Mutant", "mutant", "wt")), fh = fh)
p14b <- ggplot(fd, aes(idh, fh, colour = idh)) + geom_jitter(width = 0.25, size = 1, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#1b2a3a", width = 0.5) +
  scale_colour_manual(values = c(`IDH-mutant` = "#c0392b", `IDH-wt` = "#2a6f97"), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = NULL, y = "CGI CpGs hypermethylated\n(beta > 0.3)")
save(p14a + p14b + plot_layout(widths = c(2.2, 1)), "fig14-gcimp.png", h = 4.6)

# ---- fig22: a classifier for MSI from methylation -----------------------------------
set.seed(42)
pts <- unique(samples$patient)
test_pts <- sample(pts, round(0.3 * length(pts)))
tt <- tum & !is.na(samples$msi)
ok <- rowSums(is.na(B)) == 0
X2 <- t(B[ok, tt]); y2 <- factor(samples$msi[tt], levels = c("MSS", "MSI-H"))
is_test <- samples$patient[tt] %in% test_pts
set.seed(1)
fit <- cv.glmnet(X2[!is_test, ], y2[!is_test], family = "binomial", alpha = 0.5)
pr  <- predict(fit, X2[is_test, ], s = "lambda.1se", type = "response")[, 1]
cf  <- coef(fit, s = "lambda.1se"); nsel <- sum(cf[-1, 1] != 0)
cat("MSI classifier: train", sum(!is_test), "test", sum(is_test), "selected CpGs", nsel, "\n")
print(table(predicted = ifelse(pr > 0.5, "MSI-H", "MSS"), true = y2[is_test]))
cd <- data.frame(true = y2[is_test], p = pr,
                 mlh1 = ifelse(gm[tt, "MLH1"][is_test] > 0.3, "MLH1 methylated", "MLH1 unmethylated"))
print(table(cd$true, cd$p > 0.5, cd$mlh1))
p22 <- ggplot(cd, aes(true, p, colour = mlh1)) + geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_jitter(width = 0.2, height = 0, size = 2.6, alpha = 0.85) +
  scale_colour_manual(values = c(`MLH1 methylated` = "#c0392b", `MLH1 unmethylated` = "#5b6b7a"), name = NULL) +
  labs(x = "true MSI status (test patients)", y = "predicted probability of MSI-H",
       title = sprintf("Elastic net on %d CpGs, trained on %d tumours: %d CpGs used",
                       ncol(X2), sum(!is_test), nsel))
save(p22, "fig22-classifier.png", w = 10, h = 4.8)

# ---- fig24: screening arithmetic ------------------------------------------------------
prev <- 10^seq(-4, -1, length.out = 200)
ppv <- function(sens, spec, prev) sens * prev / (sens * prev + (1 - spec) * (1 - prev))
pp <- rbind(data.frame(prev, ppv = ppv(0.5, 0.993, prev), test = "sensitivity 50%, specificity 99.3%"),
            data.frame(prev, ppv = ppv(0.5, 0.95, prev), test = "sensitivity 50%, specificity 95%"))
cat("PPV at 1% prevalence:", round(ppv(0.5, 0.993, 0.01), 3), round(ppv(0.5, 0.95, 0.01), 3), "\n")
p24 <- ggplot(pp, aes(prev, ppv, colour = test)) + geom_line(linewidth = 1.3) +
  scale_x_log10(labels = function(x) paste0(100 * x, "%")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 1)) +
  geom_vline(xintercept = 0.01, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c("#2a6f97", "#c0392b"), name = NULL) +
  labs(x = "prevalence of cancer in the tested population (log scale)",
       y = "positive predictive value", title = "A positive test means cancer how often?") +
  theme(legend.position = c(0.3, 0.85))
save(p24, "fig24-ppv.png", w = 9, h = 4.6)
