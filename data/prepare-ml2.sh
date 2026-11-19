#!/bin/bash
# Teacher script for the machine-learning units (How a CNN reads DNA; What did
# the model learn? Variant effects). Students build the same files themselves in
# the exercises; this script makes them once so that the figure scripts, and
# students who cannot train a model, can use them.
#
#   data/prepare-ml2.sh        # ~2 min download + build, ~5-10 min training on 4 CPU threads
#
# Needs: curl, R with GenomicRanges, rtracklayer, Rsamtools, Biostrings and torch
# (install.packages("torch"); torch::install_torch()), and in data/:
# gm12878_ctcf_peaks.bed.gz, hg38_blacklist.bed.gz, hg38.fa.bgz (+ .fai, .gzi).
#
# Writes into data/:
#   ctcf_asb_adastra.tsv        CTCF allele-specific binding SNPs, FDR < 0.05, all cell
#                               types (ADASTRA release Mabel v6.1, Abramov et al. 2021)
#   ml_ctcf_dataset.rds         200-bp CTCF summit windows + GC-matched background,
#                               split by chromosome (train / valid chr7 / test chr8, chr9)
#   ml_ctcf_cnn.pt              the trained CNN (torch_save), seed 1
#   ml_ctcf_cnn_history.tsv     training and validation loss per epoch
# The R code below is the code of the exercise "How a CNN reads DNA" (R1, R2, R4,
# R5) unchanged, apart from the data path.
set -euo pipefail
cd "$(dirname "$0")"
D=$(pwd)

[ -s ctcf_asb_adastra.tsv ] || curl -fL -sS --retry 3 -o ctcf_asb_adastra.tsv \
  "https://adastra.autosome.org/api/v6/search/snps/advanced/tsv?transcription_factors=CTCF_HUMAN&fdr=0.05&es=0"

if [ -s ml_ctcf_cnn.pt ] && [ -s ml_ctcf_dataset.rds ]; then echo "have ml_ctcf_*"; exit 0; fi
TMP=$(mktemp -d -p "${TMPDIR:-/tmp}" ml2.XXXXXX)
cd "$TMP"
cat > build.R <<'REOF'
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Rsamtools); library(Biostrings); library(torch)
})
D <- Sys.getenv("D")
genome <- FaFile(file.path(D, "hg38.fa.bgz"))
chroms <- paste0("chr", c(1:22, "X"))
ctcf  <- import(file.path(D, "gm12878_ctcf_peaks.bed.gz"), format = "narrowPeak")
ctcf  <- ctcf[seqnames(ctcf) %in% chroms]
black <- import(file.path(D, "hg38_blacklist.bed.gz"))
avoid <- c(granges(ctcf) + 1000, black)          # background stays >1 kb from any CTCF peak
gc_content <- function(s) letterFrequency(s, "GC", as.prob = TRUE)[, 1]

set.seed(1)
one_chrom <- function(chr) {
  chrseq <- getSeq(genome, GRanges(chr, IRanges(1, seqlengths(genome)[chr])))[[1]]
  # positives: 200 bp centred on each summit
  p   <- ctcf[seqnames(ctcf) == chr]
  pos <- unique(GRanges(chr, IRanges(start(p) + p$peak - 99, width = 200)))
  pos <- pos[!overlapsAny(pos, black)]
  ps  <- DNAStringSet(Views(chrseq, ranges(pos)))
  ok  <- letterFrequency(ps, "N")[, 1] == 0
  pos <- pos[ok]; ps <- ps[ok]
  # candidate negatives: 20 random windows per positive, same chromosome
  cand <- GRanges(chr, IRanges(sample(length(chrseq) - 200, 20 * length(pos)), width = 200))
  cand <- cand[!overlapsAny(cand, avoid)]
  cs   <- DNAStringSet(Views(chrseq, ranges(cand)))
  ok   <- letterFrequency(cs, "N")[, 1] == 0
  cand <- cand[ok]; cs <- cs[ok]
  # GC matching: for each 2% GC bin, draw as many candidates as there are positives
  bins <- seq(0, 1, 0.02)
  pb <- findInterval(gc_content(ps), bins); cb <- findInterval(gc_content(cs), bins)
  chosen <- unlist(lapply(sort(unique(pb)), function(b) {
    have <- which(cb == b); need <- sum(pb == b)
    have[sample.int(length(have), min(need, length(have)))]
  }))
  data.frame(chr = chr, start = c(start(pos), start(cand[chosen])),
             label = rep(c(1L, 0L), c(length(pos), length(chosen))),
             seq = c(as.character(ps), as.character(cs[chosen])))
}
ds <- do.call(rbind, lapply(chroms, one_chrom))
ds$gc <- gc_content(DNAStringSet(ds$seq))
ds$split <- ifelse(ds$chr %in% c("chr8", "chr9"), "test",
                   ifelse(ds$chr == "chr7", "valid", "train"))
table(ds$split, ds$label)
tapply(ds$gc, ds$label, mean)
saveRDS(ds, "ctcf_dataset.rds")

one_hot <- function(s) {
  m <- do.call(rbind, strsplit(s, ""))              # one row per sequence, one column per base
  x <- array(0, c(length(s), 4, ncol(m)))           # sequences x 4 bases x positions
  for (k in 1:4) x[, k, ] <- m == c("A", "C", "G", "T")[k]
  x
}
one_hot("GCCACTAGG")[1, , ]

sets <- lapply(split(ds, ds$split), function(d)
  list(x = torch_tensor(one_hot(d$seq), dtype = torch_float()),
       y = torch_tensor(d$label, dtype = torch_float())))
sets$train$x$shape

auc <- function(score, y) {                         # area under the ROC curve, from ranks
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cnn <- nn_module("cnn",
  initialize = function(n_filters = 32, width = 19) {
    self$conv  <- nn_conv1d(4, n_filters, kernel_size = width)   # the filters (motif scanners)
    self$dense <- nn_linear(n_filters, 1)                          # logistic regression on top
  },
  forward = function(x) {
    h <- nnf_relu(self$conv(x))        # sequences x filters x 182 positions
    h <- torch_amax(h, dim = 3)        # global max pooling: sequences x filters
    self$dense(h)$squeeze(2)           # one logit (log-odds of "bound") per sequence
  })
torch_manual_seed(1)
model <- cnn()
model
sapply(model$parameters, function(p) p$numel())
x1 <- sets$test$x[1:5, , ]
torch_sigmoid(model(x1))

torch_set_num_threads(4)
set.seed(1); torch_manual_seed(1)
model   <- cnn()
opt     <- optim_adam(model$parameters, lr = 1e-3)
loss_fn <- nn_bce_with_logits_loss()                # sigmoid + binary cross-entropy in one step
predict_logit <- function(model, x) { model$eval(); with_no_grad(as.numeric(model(x))) }
bce <- function(z, y) mean(y * log1p(exp(-z)) + (1 - y) * log1p(exp(z)))

history <- data.frame(); best <- Inf; patience <- 0
n <- length(sets$train$y)
for (epoch in 1:20) {
  model$train()
  for (b in split(sample(n), ceiling(seq_len(n) / 128))) {   # shuffled mini-batches of 128
    opt$zero_grad()
    loss <- loss_fn(model(sets$train$x[b, , ]), sets$train$y[b])
    loss$backward()                                         # gradients by backpropagation
    opt$step()                                              # one gradient-descent step
  }
  ztr <- predict_logit(model, sets$train$x); zva <- predict_logit(model, sets$valid$x)
  ytr <- as.numeric(sets$train$y);           yva <- as.numeric(sets$valid$y)
  history <- rbind(history, data.frame(epoch,
    train_loss = bce(ztr, ytr), valid_loss = bce(zva, yva),
    train_auc = auc(ztr, ytr),  valid_auc = auc(zva, yva)))
  print(round(tail(history, 1), 4))
  if (tail(history$valid_loss, 1) < best) {                 # keep the best model so far
    best <- tail(history$valid_loss, 1); patience <- 0
    torch_save(model, "ctcf_cnn.pt")
  } else {
    patience <- patience + 1
    if (patience == 3) break                                # early stopping
  }
}
write.table(history, "ctcf_cnn_history.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
matplot(history$epoch, history[, c("train_loss", "valid_loss")], type = "b", pch = 19,
        col = c("grey40", "#2a6f97"), lty = 1, xlab = "epoch", ylab = "loss (binary cross-entropy)")
legend("topright", c("training", "validation (chr7)"), col = c("grey40", "#2a6f97"), pch = 19)
REOF
D="$D" Rscript build.R
mv ctcf_dataset.rds "$D/ml_ctcf_dataset.rds"
mv ctcf_cnn.pt "$D/ml_ctcf_cnn.pt"
mv ctcf_cnn_history.tsv "$D/ml_ctcf_cnn_history.tsv"
cd "$D" && rm -rf "$TMP"
echo "done: ml_ctcf_dataset.rds ml_ctcf_cnn.pt ml_ctcf_cnn_history.tsv ctcf_asb_adastra.tsv"
