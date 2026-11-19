# Real-data figures for methylation-arrays/lec.qmd.
# Run from this folder:  Rscript make-arrays-figs.R      (about 25 minutes, 10 GB RAM)
# Needs: data/prepare-meth2.sh arrays  (GSE42861 IDATs of 151 samples, sample table,
#        450k GRCh38 manifest and probe masks of Zhou et al. 2017)
#        R packages minfi, IlluminaHumanMethylation450kmanifest,
#        IlluminaHumanMethylation450kanno.ilmn12.hg19, FlowSorted.Blood.450k, limma,
#        data.table, ggplot2, patchwork, ggrepel.
#        If preprocessFunnorm stops with "return code from pthread_create() is 22",
#        reinstall preprocessCore with
#        BiocManager::install("preprocessCore", configure.args = "--disable-threading", force = TRUE)
#
# Writes fig04, fig06, fig07, fig09, fig12, fig13, fig14, fig15 (*.png) and prints the
# numbers quoted in the slides.

suppressPackageStartupMessages({
  library(minfi); library(limma); library(data.table); library(ggplot2); library(patchwork)
})
set.seed(1)
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save <- function(p, f, w = 12, h = 4.6) ggsave(f, p, width = w, height = h, dpi = 150, bg = "white")
colS <- c(control = "#5b6b7a", RA = "#e67e22")

s <- fread(file.path(D, "gse42861_samples.tsv"))
s$Basename <- file.path(D, "gse42861_idat", s$basename)
t0 <- Sys.time()
rg <- read.metharray.exp(targets = as.data.frame(s))
cat("read", ncol(rg), "arrays in", format(Sys.time() - t0), "\n")

# ---- fig07: detection p-values and sex check ----------------------------------------
dp <- detectionP(rg)
fail <- colMeans(dp > 0.01)
cat("failed probes per sample (p > 0.01): median", signif(median(fail), 2), " max", signif(max(fail), 2), "\n")
cat("probes failing in > 5% of samples:", sum(rowMeans(dp > 0.01) > 0.05), "of", nrow(dp), "\n")
mset <- preprocessRaw(rg)
sx <- getSex(mapToGenome(mset))
cat("reported vs predicted sex:\n"); print(table(reported = s$sex, predicted = sx$predictedSex))
q <- data.table(xMed = sx$xMed, yMed = sx$yMed, sex = s$sex, fail = fail, slide = s$slide)
p07a <- ggplot(q, aes(reorder(seq_along(fail), fail), 100 * fail)) + geom_col(fill = "#2a6f97") +
  geom_hline(yintercept = 1, linetype = 2, colour = "#c0392b") +
  annotate("text", x = 5, y = 1.15, label = "1% (a common cut-off)", hjust = 0, vjust = 0, colour = "#c0392b", size = 4.3) +
  labs(x = "151 arrays, sorted", y = "% probes with detection p > 0.01", title = "Failed probes per array") +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
p07b <- ggplot(q, aes(xMed, yMed, colour = sex)) + geom_point(size = 2.2, alpha = 0.8) +
  scale_colour_manual(values = c(F = "#a8324a", M = "#2a6f97"), name = "reported") +
  labs(x = "median log2 intensity, chrX probes", y = "median log2 intensity, chrY probes",
       title = "Sex from the array matches the records")
save(p07a + p07b, "fig07-qc.png", w = 13, h = 4.8)

# ---- normalisations on two slides (24 arrays) ----------------------------------------
two <- s$slide %in% c(7800246040, 7800246034)
rg2 <- rg[, two]
nb <- list(raw = getBeta(preprocessRaw(rg2), offset = 100), noob = getBeta(preprocessNoob(rg2)),
           SWAN = getBeta(preprocessSWAN(rg2)), funnorm = getBeta(preprocessFunnorm(rg2)))
ti <- getProbeInfo(rg2, type = "I")$Name; tii <- getProbeInfo(rg2, type = "II")$Name
dens <- rbindlist(lapply(names(nb), function(m) {
  x <- nb[[m]]; rbindlist(lapply(1:4, function(j) {
    rbind(data.table(beta = x[ti, j], type = "Type I"), data.table(beta = x[tii, j], type = "Type II"))[, `:=`(method = m, sample = j)]
  }))
}))
dens[, method := factor(method, levels = names(nb))]
pk <- dens[, .(peak = { d <- density(beta[beta > 0.6], na.rm = TRUE); d$x[which.max(d$y)] }), by = .(method, type)]
print(dcast(pk, method ~ type, value.var = "peak"))
p06 <- ggplot(dens, aes(beta, colour = type, group = interaction(type, sample))) + geom_density(linewidth = 0.6) +
  facet_wrap(~ method, nrow = 1) + scale_colour_manual(values = c(`Type I` = "#2a6f97", `Type II` = "#c0392b"), name = NULL) +
  labs(x = "β", y = "density", title = "β distributions of 4 arrays, by probe design") + theme(legend.position = "top")
save(p06, "fig06-normalisation.png", w = 13, h = 4.6)
rm(nb, dens); invisible(gc())

# ---- noob for all arrays; filtering -------------------------------------------------
noob <- preprocessNoob(rg); b <- getBeta(noob); rm(noob); invisible(gc())
man <- fread(file.path(D, "hm450_hg38_manifest.tsv.gz"), select = c("Probe_ID", "CpG_chrm", "CpG_beg"))
mask <- fread(file.path(D, "hm450_hg38_mask.tsv.gz"))
k1 <- rownames(b)[rowMeans(dp[rownames(b), ] > 0.01) <= 0.05]
k2 <- setdiff(k1, mask[M_general == TRUE, Probe_ID])
k3 <- setdiff(k2, man[CpG_chrm %in% c("chrX", "chrY") | is.na(CpG_chrm) | CpG_chrm == "", Probe_ID])
cat("probes:", nrow(b), "-> detection", length(k1), "-> M_general mask", length(k2), "-> autosomes", length(k3), "\n")

# ---- fig04: beta vs M ---------------------------------------------------------------
bk <- b[k3, ]; M <- log2(pmax(bk, 1e-3) / pmax(1 - bk, 1e-3))
sub <- sample(nrow(bk), 40000)
vm <- data.table(mean_beta = rowMeans(bk[sub, ]), sd_beta = apply(bk[sub, ], 1, sd), sd_M = apply(M[sub, ], 1, sd))
p04a <- ggplot(data.table(beta = seq(0.001, 0.999, length.out = 400))[, M := log2(beta / (1 - beta))], aes(beta, M)) +
  geom_line(linewidth = 1.2, colour = "#2a6f97") + geom_hline(yintercept = 0, linetype = 2, colour = "#5b6b7a") +
  labs(x = "β", y = "M = log2(β / (1 − β))", title = "β and M carry the same information")
p04b <- ggplot(vm, aes(mean_beta, sd_beta)) + geom_bin2d(bins = 60) + scale_fill_gradient(low = "#dfe6ec", high = "#1b2a3a", trans = "log10", guide = "none") +
  labs(x = "mean β (151 arrays)", y = "SD of β", title = "β: variance squeezed at 0 and 1")
p04c <- ggplot(vm, aes(mean_beta, sd_M)) + geom_bin2d(bins = 60) + scale_fill_gradient(low = "#dfe6ec", high = "#1b2a3a", trans = "log10", guide = "none") +
  labs(x = "mean β (151 arrays)", y = "SD of M", title = "M: more even variance")
save(p04a + p04b + p04c, "fig04-beta-m.png", w = 14, h = 4.6)

# ---- cell composition (Houseman, FlowSorted.Blood.450k) -----------------------------
cc <- estimateCellCounts(rg, compositeCellType = "Blood",
                         cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"), returnAll = TRUE)
fr <- cc$counts; ct <- cc$compTable
rm(rg, rg2); invisible(gc())
s <- cbind(s, as.data.table(fr))
print(s[, lapply(.SD, function(x) round(mean(x), 3)), by = status, .SDcols = colnames(fr)])
cat("Gran RA vs control: t-test p =", signif(t.test(Gran ~ status, s)$p.value, 2), "\n")
# fig12: reference profiles of marker CpGs. compTable has every CpG; pick, like minfi, the 50 most
# hypo- and 50 most hypermethylated CpGs per cell type (vs the mean of the others) among F-test p < 1e-8
cts <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
ref <- as.matrix(ct[, cts]); rownames(ref) <- rownames(ct); ok <- ct$p.value < 1e-8
mk <- unique(unlist(lapply(cts, function(x) { d <- (ref[, x] - rowMeans(ref[, setdiff(cts, x)]))[ok]
  c(names(sort(d))[1:50], names(sort(d, decreasing = TRUE))[1:50]) })))
cat("marker CpGs for the heat map:", length(mk), "\n")
ref <- ref[mk, ]
o <- order(apply(ref, 1, function(v) which.max(abs(v - median(v)))), -apply(ref, 1, max))
hm <- melt(data.table(cpg = seq_len(nrow(ref)), ref)[o][, cpg := factor(cpg, levels = cpg)], id.vars = "cpg", variable.name = "cell", value.name = "beta")
p12 <- ggplot(hm, aes(cell, cpg, fill = beta)) + geom_tile() +
  scale_fill_gradient(low = "#eaf1f7", high = "#1b2a3a", name = "β", limits = c(0, 1)) +
  labs(x = NULL, y = sprintf("%d marker CpGs", nrow(ref)), title = "Sorted blood cells (Reinius 2012)") +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())
save(p12, "fig12-reference.png", w = 6.5, h = 5.6)
# fig13: estimated fractions by status
fl <- melt(s[, c("status", colnames(fr)), with = FALSE], id.vars = "status", variable.name = "cell", value.name = "f")
p13 <- ggplot(fl, aes(cell, f, fill = status)) + geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = colS, name = NULL) +
  labs(x = NULL, y = "estimated fraction", title = "Houseman estimates, 151 blood samples") + theme(legend.position = "top")
save(p13, "fig13-cellfractions.png", w = 10, h = 4.8)

# ---- PCA (fig09) --------------------------------------------------------------------
v <- order(apply(M, 1, var), decreasing = TRUE)[1:20000]
pc <- prcomp(t(M[v, ]), rank. = 4)
ve <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2), 1)
pd <- data.table(PC1 = pc$x[, 1], PC2 = pc$x[, 2], status = s$status, Gran = s$Gran,
                 chips = ifelse(substr(s$slide, 1, 7) == "5730053", "5730053... chips (4)", "other chips (10)"))
cat("variance explained PC1, PC2:", ve, "\n")
cat("PC1 ~ slide R2:", round(summary(lm(pd$PC1 ~ factor(s$slide)))$r.squared, 2), " PC1 ~ chip series R2:", round(summary(lm(PC1 ~ chips, pd))$r.squared, 2),
    " PC1 vs Gran r:", round(cor(pd$PC1, pd$Gran), 2), " PC1 vs failed-probe fraction r:", round(cor(pd$PC1, fail), 2), "\n")
cat("PC2 vs Gran r:", round(cor(pd$PC2, pd$Gran), 2), " PC2 RA vs control t-test p:", signif(t.test(PC2 ~ status, pd)$p.value, 2), "\n")
p09a <- ggplot(pd, aes(PC1, PC2, colour = chips)) + geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("#a8324a", "#2a6f97"), name = NULL) +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = "PC1: which batch of BeadChips") + theme(legend.position = "top")
p09b <- ggplot(pd, aes(Gran, PC2, colour = status)) + geom_point(size = 2, alpha = 0.8) + scale_colour_manual(values = colS, name = NULL) +
  labs(x = "estimated granulocyte fraction", y = sprintf("PC2 (%.0f%%)", ve[2]), title = sprintf("PC2: cell composition (r = %.2f)", cor(pd$PC2, pd$Gran))) +
  theme(legend.position = "top")
save(p09a + p09b, "fig09-pca.png", w = 12.5, h = 4.8)

# ---- EWAS with limma (fig14, fig15) ----------------------------------------------------
s[, status := factor(status, levels = c("control", "RA"))]
s[, smoking := factor(smoking, levels = c("never", "ex", "occasional", "current"))]
s[, `:=`(slide = factor(slide), sex = factor(sex))]
keep <- !is.na(s$smoking); s <- s[keep]; M <- M[, keep]; bk <- bk[, keep]
cat("EWAS samples (smoking status known):", nrow(s), "\n")
lam <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
fitm <- function(f, coef = "statusRA") {
  X <- model.matrix(f, s); fit <- eBayes(lmFit(M, X))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  list(tt = tt, lambda = lam(tt$P.Value), n = sum(tt$adj.P.Val < 0.05), df = nrow(s) - ncol(X))
}
mods <- list(`status only` = ~ status,
             `+ age, sex, smoking` = ~ status + age + sex + smoking,
             `+ cell fractions` = ~ status + age + sex + smoking + CD8T + CD4T + NK + Bcell + Mono,
             `+ slide` = ~ status + age + sex + smoking + CD8T + CD4T + NK + Bcell + Mono + slide)
res <- lapply(mods, fitm)
for (n in names(res)) cat(sprintf("%-22s lambda %.2f  FDR<0.05: %6d  residual df %d\n", n, res[[n]]$lambda, res[[n]]$n, res[[n]]$df))
qq <- rbindlist(lapply(names(res), function(n) { p <- sort(res[[n]]$tt$P.Value); k <- length(p)
  i <- unique(c(1:2000, round(seq(1, k, length.out = 3000))))
  data.table(exp = -log10(i / (k + 1)), obs = -log10(p[i]),
             model = sprintf("%s\nλ = %.2f, %s CpGs at FDR 5%%", n, res[[n]]$lambda, format(res[[n]]$n, big.mark = " "))) }))
qq[, model := factor(model, levels = unique(model))]
p14 <- ggplot(qq, aes(exp, obs)) + geom_abline(linetype = 2, colour = "#5b6b7a") + geom_point(size = 0.7, colour = "#2a6f97") +
  facet_wrap(~ model, nrow = 1) + labs(x = "expected −log10 p", y = "observed −log10 p", title = "RA vs control, 442 415 CpGs: QQ plots as covariates are added")
save(p14, "fig14-qq.png", w = 14, h = 4.6)
tt2 <- res[["+ cell fractions"]]$tt
db <- rowMeans(bk[, s$status == "RA"]) - rowMeans(bk[, s$status == "control"])
cat("RA, cell-adjusted model: max |delta beta| among FDR hits", round(max(abs(db[tt2$adj.P.Val < 0.05])), 3),
    " median", round(median(abs(db[tt2$adj.P.Val < 0.05])), 3), "\n")

sm <- fitm(mods[["+ cell fractions"]], coef = "smokingcurrent"); st <- sm$tt
cat("smoking (current vs never): lambda", round(sm$lambda, 2), " FDR<0.05:", sm$n, "\n")
st <- st[order(st$P.Value), ]; print(head(st[, c("logFC", "P.Value", "adj.P.Val")], 8))
ahrr <- "cg05575921"
cat("AHRR", ahrr, "rank", which(rownames(st) == ahrr), " mean beta never", round(mean(bk[ahrr, s$smoking == "never"]), 3),
    " ex", round(mean(bk[ahrr, s$smoking == "ex"]), 3), " current", round(mean(bk[ahrr, s$smoking == "current"]), 3), "\n")
st$cpg <- rownames(st)
lab <- c(cg05575921 = "AHRR", cg03636183 = "F2RL3", cg21566642 = "2q37.1", cg06126421 = "6p21.33")
st$lab <- lab[st$cpg]
p15a <- ggplot(st, aes(logFC, -log10(P.Value))) + geom_point(size = 0.6, colour = "#a9b5c0") +
  geom_point(data = st[st$adj.P.Val < 0.05, ], colour = "#c0392b", size = 1.6) +
  ggrepel::geom_text_repel(data = st[!is.na(st$lab), ], aes(label = lab), size = 4.2, min.segment.length = 0, nudge_x = -0.4, seed = 1) +
  labs(x = "difference in M-value (current − never)", y = "−log10 p",
       title = sprintf("Smoking: %d CpGs at FDR 5%%, λ = %.2f", sm$n, sm$lambda))
ad <- data.table(smoking = s$smoking, beta = bk[ahrr, ])
p15b <- ggplot(ad, aes(smoking, beta)) + geom_boxplot(fill = "#e9c46a", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  labs(x = NULL, y = "β", title = "AHRR cg05575921")
save(p15a + p15b + plot_layout(widths = c(1.4, 1)), "fig15-smoking.png", w = 13, h = 4.8)
cat("done\n")
