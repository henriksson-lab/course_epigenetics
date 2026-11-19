# Real-data figures for ml-cnn-dna/lec.qmd.
# Run from this folder:  Rscript make-mlcnn-figs.R   (about 3 minutes on 4 CPU threads)
# Needs: data/prepare-ml2.sh  (builds ml_ctcf_dataset.rds, ml_ctcf_cnn.pt and
#        ml_ctcf_cnn_history.tsv with the exact code of the exercise, R1-R5)
#        R packages torch, Biostrings, JASPAR2024, TFBSTools, RSQLite, motifmatchr,
#        universalmotif, ggplot2, patchwork, ggrepel.
#
# Writes fig02-pwm-baseline.png, fig12-learning-curves.png, fig13-roc.png,
# fig14-filter-logos.png, fig15-filter-map.png, fig16-filter-use.png
# and prints the numbers quoted in the slides.

suppressPackageStartupMessages({
  library(torch); library(Biostrings); library(JASPAR2024); library(TFBSTools); library(RSQLite)
  library(motifmatchr); library(universalmotif); library(ggplot2); library(patchwork)
})
D <- "../data"
theme_set(theme_minimal(base_size = 15) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
col <- c(cnn = "#2a6f97", pwm = "#8c564b", gc = "grey60", train = "grey40", valid = "#2a6f97")
torch_set_num_threads(4)

# ---- data, model and helpers (the same code as the exercise) ------------------
ds <- readRDS(file.path(D, "ml_ctcf_dataset.rds"))
one_hot <- function(s) {
  m <- do.call(rbind, strsplit(s, ""))
  x <- array(0, c(length(s), 4, ncol(m)))
  for (k in 1:4) x[, k, ] <- m == c("A", "C", "G", "T")[k]
  x
}
auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
roc <- function(score, y) {
  o <- order(score, decreasing = TRUE); s <- score[o]; y <- y[o]
  k <- c(diff(s) != 0, TRUE)
  data.frame(fpr = c(0, cumsum(y == 0)[k] / sum(y == 0)), tpr = c(0, cumsum(y == 1)[k] / sum(y == 1)))
}
cnn <- nn_module("cnn",
  initialize = function(n_filters = 32, width = 19) {
    self$conv  <- nn_conv1d(4, n_filters, kernel_size = width)
    self$dense <- nn_linear(n_filters, 1)
  },
  forward = function(x) {
    h <- nnf_relu(self$conv(x))
    h <- torch_amax(h, dim = 3)
    self$dense(h)$squeeze(2)
  })
predict_logit <- function(model, x) { model$eval(); with_no_grad(as.numeric(model(x))) }
bce <- function(z, y) mean(y * log1p(exp(-z)) + (1 - y) * log1p(exp(z)))

test <- ds[ds$split == "test", ]
xt <- torch_tensor(one_hot(test$seq), dtype = torch_float())
model <- torch_load(file.path(D, "ml_ctcf_cnn.pt"))
test$cnn <- predict_logit(model, xt)
jdb <- dbConnect(SQLite(), db(JASPAR2024()))
ctcf_m <- getMatrixByID(jdb, "MA0139.2")
test$pwm <- motifScores(matchMotifs(ctcf_m, DNAStringSet(test$seq), out = "scores", bg = "even", p.cutoff = 1))[, 1]
A <- c(gc = auc(test$gc, test$label), pwm = auc(test$pwm, test$label), cnn = auc(test$cnn, test$label))
cat(sprintf("data: %d sequences (%d positives); test %d (%d positives)\n", nrow(ds), sum(ds$label),
            nrow(test), sum(test$label)))
print(table(ds$split, ds$label))
cat(sprintf("test AUC: GC %.3f, PWM %.3f, CNN %.3f\n", A["gc"], A["pwm"], A["cnn"]))

# ---- fig02: the PWM baseline --------------------------------------------------
d2 <- data.frame(score = test$pwm, set = factor(ifelse(test$label == 1, "CTCF peak (test)", "GC-matched background (test)"),
                                                 levels = c("GC-matched background (test)", "CTCF peak (test)")))
cat(sprintf("fig02: median best PWM score: peaks %.1f, background %.1f\n",
            median(test$pwm[test$label == 1]), median(test$pwm[test$label == 0])))
p2 <- ggplot(d2, aes(score, fill = set)) +
  geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.6, colour = NA) +
  scale_fill_manual(values = c("grey55", "#8c564b"), name = NULL) +
  annotate("text", x = min(d2$score), y = Inf, hjust = 0, vjust = 1.5, size = 5.2,
           label = sprintf("best CTCF PWM score in each 200-bp window\nROC AUC = %.2f", A["pwm"])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "best log-odds score of JASPAR MA0139.2 (both strands)", y = "sequences") +
  theme(legend.position = "bottom")
ggsave("fig02-pwm-baseline.png", p2, width = 11, height = 4.6, dpi = 150, bg = "white")

# ---- fig12: learning curves ---------------------------------------------------
h <- read.delim(file.path(D, "ml_ctcf_cnn_history.tsv"))
best <- h$epoch[which.min(h$valid_loss)]
cat(sprintf("fig12a: %d epochs, best validation loss %.4f at epoch %d (train %.4f); valid AUC %.4f\n",
            nrow(h), min(h$valid_loss), best, h$train_loss[best], h$valid_auc[best]))
d12 <- rbind(data.frame(epoch = h$epoch, loss = h$train_loss, set = "training"),
             data.frame(epoch = h$epoch, loss = h$valid_loss, set = "validation (chr7)"))
pa <- ggplot(d12, aes(epoch, loss, colour = set)) + geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_vline(xintercept = best, linetype = 2, colour = "grey30") +
  annotate("text", x = best, y = max(d12$loss), hjust = 1.05, vjust = 1, size = 4.6,
           label = "best: saved model\n(stopped 3 epochs later)") +
  scale_colour_manual(values = c("grey40", "#2a6f97"), name = NULL) +
  labs(title = sprintf("All %s training sequences", format(sum(ds$split == "train"), big.mark = " ")),
       x = "epoch", y = "loss (binary cross-entropy)") +
  theme(legend.position = "bottom")

# the same network trained on only 500 sequences (exercise R7), evaluated every epoch
train <- ds[ds$split == "train", ]; valid <- ds[ds$split == "valid", ]
xtr <- torch_tensor(one_hot(train$seq), dtype = torch_float()); ytr <- torch_tensor(train$label, dtype = torch_float())
xva <- torch_tensor(one_hot(valid$seq), dtype = torch_float()); yva <- torch_tensor(valid$label, dtype = torch_float())
loss_fn <- nn_bce_with_logits_loss()
set.seed(1); torch_manual_seed(1)
small <- sample(length(ytr), 500)
xs <- xtr[small, , ]; ys <- ytr[small]
m_small <- cnn(); opt <- optim_adam(m_small$parameters, lr = 1e-3)
curve <- data.frame()
for (epoch in 1:150) {
  m_small$train()
  for (b in split(sample(500), ceiling(1:500 / 32))) {
    opt$zero_grad(); loss_fn(m_small(xs[b, , ]), ys[b])$backward(); opt$step()
  }
  curve <- rbind(curve, data.frame(epoch,
    train = bce(predict_logit(m_small, xs), as.numeric(ys)),
    valid = bce(predict_logit(m_small, xva), as.numeric(yva))))
}
bs <- curve$epoch[which.min(curve$valid)]
cat(sprintf("fig12b: 500 sequences: min validation loss %.3f at epoch %d; at epoch 150 train %.3f, valid %.3f; test AUC %.3f\n",
            min(curve$valid), bs, curve$train[150], curve$valid[150], auc(predict_logit(m_small, xt), test$label)))
print(round(curve[curve$epoch %% 10 == 0, ], 3))
d12b <- rbind(data.frame(epoch = curve$epoch, loss = curve$train, set = "training"),
              data.frame(epoch = curve$epoch, loss = curve$valid, set = "validation (chr7)"))
pb <- ggplot(d12b, aes(epoch, loss, colour = set)) + geom_line(linewidth = 1) +
  geom_vline(xintercept = bs, linetype = 2, colour = "grey30") +
  annotate("text", x = bs + 3, y = 0.62, hjust = 0, size = 4.6, label = "validation loss\nstops improving") +
  annotate("text", x = 150, y = 0.24, hjust = 1, size = 4.6, colour = "grey30",
           label = "training loss → 0: memorised") +
  scale_colour_manual(values = c("grey40", "#2a6f97"), name = NULL) +
  labs(title = "Only 500 training sequences", x = "epoch", y = NULL) +
  theme(legend.position = "none")
p12 <- pa + pb
ggsave("fig12-learning-curves.png", p12, width = 12, height = 4.9, dpi = 150, bg = "white")

# ---- fig13: ROC on the test chromosomes ----------------------------------------
lab <- c(gc = sprintf("GC content (AUC %.2f)", A["gc"]), pwm = sprintf("CTCF PWM (AUC %.2f)", A["pwm"]),
         cnn = sprintf("CNN (AUC %.3f)", A["cnn"]))
d13 <- do.call(rbind, lapply(names(lab), function(m) cbind(roc(test[[m]], test$label), model = lab[[m]])))
d13$model <- factor(d13$model, levels = lab[c("cnn", "pwm", "gc")])
tpr5 <- sapply(c(cnn = "cnn", pwm = "pwm"), function(m) { r <- roc(test[[m]], test$label); r$tpr[which.max(r$fpr >= 0.05)] })
cat(sprintf("fig13: TPR at 5%% FPR: CNN %.3f, PWM %.3f\n", tpr5["cnn"], tpr5["pwm"]))
p13 <- ggplot(d13, aes(fpr, tpr, colour = model)) +
  geom_abline(linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0.05, linetype = 3, colour = "grey40") +
  geom_path(linewidth = 1.2) +
  annotate("text", x = 0.5, y = 0.22, hjust = 0, size = 4.6,
           label = sprintf("at 5%% false positives:\nCNN finds %.0f%% of peaks,\nPWM %.0f%%", 100 * tpr5["cnn"], 100 * tpr5["pwm"])) +
  scale_colour_manual(values = unname(col[c("cnn", "pwm", "gc")]), name = NULL) +
  coord_equal() +
  labs(x = "false positive rate (background called bound)", y = "true positive rate (peaks found)") +
  theme(legend.position = "right")
ggsave("fig13-roc.png", p13, width = 10, height = 5.2, dpi = 150, bg = "white")

# ---- filters: activation motifs and JASPAR matches -----------------------------
act <- with_no_grad(as.array(nnf_relu(model$conv(xt))))
w_out <- as.numeric(model$dense$weight)
filter_motif <- function(k) {
  a  <- act[, k, ]
  mx <- apply(a, 1, max); at <- apply(a, 1, which.max)
  keep <- mx > 0.5 * max(mx)
  sites <- DNAStringSet(substring(test$seq[keep], at[keep], at[keep] + 18))
  create_motif(consensusMatrix(sites)[c("A", "C", "G", "T"), ], type = "PCM",
               name = paste0("filter", k), nsites = sum(keep))
}
fm <- lapply(1:32, filter_motif)
pfms <- getMatrixSet(jdb, list(collection = "CORE", tax_group = "vertebrates"))
jas  <- lapply(seq_along(pfms), function(i) {
  m <- convert_motifs(pfms[[i]], "universalmotif"); m@name <- paste(ID(pfms[[i]]), name(pfms[[i]])); m })
cat(sprintf("JASPAR 2024 CORE vertebrates: %d matrices\n", length(jas)))
bm <- do.call(rbind, lapply(1:32, function(k) {
  r <- as.data.frame(compare_motifs(c(fm[k], jas), compare.to = 1, method = "PCC",
                                    min.mean.ic = 0, max.p = 1, tryRC = TRUE))
  data.frame(filter = k, w_out = w_out[k], sites = fm[[k]]@nsites,
             ic = sum(convert_type(fm[[k]], "ICM")@motif), consensus = fm[[k]]@consensus,
             match = r$target[1], PCC = r$score[1], p = r$Pval[1])
}))
bm$ctcf <- grepl("CTCF", bm$match) & bm$PCC >= 0.9
print(bm[order(-bm$w_out), c("filter", "w_out", "sites", "ic", "consensus", "match", "PCC", "p")], digits = 3)
cat(sprintf("filters matching CTCF/CTCFL (PCC >= 0.9): %d of 32; all with positive output weight: %s\n",
            sum(bm$ctcf), all(bm$w_out[bm$ctcf] > 0)))

# fig14: logos of the four filters with the largest output weight, next to JASPAR CTCF
top <- order(-w_out)[1:4]
ctcf_u <- convert_motifs(ctcf_m, "universalmotif"); ctcf_u@name <- "JASPAR MA0139.2 CTCF"
for (k in top) fm[[k]]@name <- sprintf("filter %d  (weight %+.2f; PCC %.2f with %s)", k, w_out[k],
                                       bm$PCC[k], sub("^\\S+ ", "", bm$match[k]))
p14 <- view_motifs(c(list(ctcf_u), fm[top]), use.type = "ICM", show.positions = FALSE,
                   names.pos = "top", text.size = 15) +
  theme(plot.margin = margin(4, 4, 4, 4), strip.text = element_text(size = 15))
ggsave("fig14-filter-logos.png", p14, width = 10, height = 8, dpi = 150, bg = "white")

# fig15: all 32 filters
bm$group <- ifelse(bm$ctcf, "best match CTCF / CTCFL", "other best match")
bm$lab <- ifelse(bm$ctcf | bm$w_out < 0, "", sub("^\\S+ ", "", bm$match))
neg <- bm[bm$w_out < 0, ]
p15 <- ggplot(bm, aes(ic, w_out, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(size = 3.2) +
  ggrepel::geom_text_repel(aes(label = lab), size = 4, show.legend = FALSE, seed = 1, max.overlaps = 30) +
  annotate("text", x = max(neg$ic) + 0.4, y = mean(neg$w_out), hjust = 0, size = 4.6, colour = "grey35",
           label = sprintf("%d filters: short GC-rich patterns, negative weights;\nbest matches PLAGL2, KLF, ZIC, SP, E2F ... (PCC %.2f-%.2f)",
                           nrow(neg), min(neg$PCC), max(neg$PCC))) +
  scale_colour_manual(values = c("#8c564b", "grey55"), name = NULL) +
  labs(x = "information content of the filter's motif (bits)",
       y = "output weight (dense layer)") +
  theme(legend.position = "bottom")
ggsave("fig15-filter-map.png", p15, width = 11, height = 5, dpi = 150, bg = "white")

# ---- fig16: where the strongest CTCF filter fires, and which filters matter ----
k1 <- top[1]
a1 <- act[, k1, ]; mx <- apply(a1, 1, max); at <- apply(a1, 1, which.max)
fires <- mx > 0.5 * max(mx)
centre <- at + 9 - 100.5                               # filter centre relative to the summit
cat(sprintf("fig16a: filter %d fires in %.1f%% of test peaks and %.1f%% of background; median |distance| %.0f bp\n",
            k1, 100 * mean(fires[test$label == 1]), 100 * mean(fires[test$label == 0]),
            median(abs(centre[fires & test$label == 1]))))
d16 <- data.frame(centre = centre[fires & test$label == 1])
p16a <- ggplot(d16, aes(centre)) + geom_histogram(binwidth = 5, boundary = 0, fill = "#8c564b") +
  labs(title = sprintf("Filter %d in test peaks", k1), x = "position of strongest match\nrelative to the summit (bp)", y = "peaks")
feat <- apply(act, c(1, 2), max)                       # max-pooled features, sequences x 32
b0 <- as.numeric(model$dense$bias)
z_with <- function(keep) b0 + feat[, keep, drop = FALSE] %*% w_out[keep]
abl <- c("all 32 filters" = auc(z_with(1:32), test$label),
         "only the CTCF-like filters" = auc(z_with(which(bm$ctcf)), test$label),
         "all except the CTCF-like ones" = auc(z_with(which(!bm$ctcf)), test$label),
         "JASPAR CTCF PWM" = A[["pwm"]])
print(round(abl, 3))
d16b <- data.frame(model = factor(names(abl), levels = rev(names(abl))), auc = abl)
p16b <- ggplot(d16b, aes(auc, model)) + geom_col(fill = c("#2a6f97", "#8c564b", "grey60", "#8c564b"), width = 0.6, alpha = c(1, 1, 1, 0.5)) +
  geom_text(aes(label = sprintf("%.3f", auc)), hjust = -0.15, size = 4.8) +
  coord_cartesian(xlim = c(0.5, 1.02)) +
  labs(title = "Test AUC using some of the filters", x = "ROC AUC (chr8 + chr9)", y = NULL)
p16 <- p16a + p16b + plot_layout(widths = c(1, 1.5))
ggsave("fig16-filter-use.png", p16, width = 12, height = 4.6, dpi = 150, bg = "white")
