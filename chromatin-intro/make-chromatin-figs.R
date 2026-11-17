# Real-data figures for chromatin-intro/lec.qmd: genome-browser style views of
# six histone marks in GM12878 (ENCODE fold change over input).
# Run from this folder:  Rscript make-chromatin-figs.R
# Needs: data/download.sh chromatin-intro   (about 2.6 GB of bigWigs + the GTF)
#
# Writes fig11-tracks-active-gene.png and fig12-tracks-silent.png.

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(ggplot2); library(patchwork)
})
D <- "../data"
marks <- c(H3K4me3 = "#c0392b", H3K4me1 = "#d4a017", H3K27ac = "#e67e22",
           H3K36me3 = "#2e8b57", H3K27me3 = "#2a6f97", H3K9me3 = "#6b4fa0")
bw_file <- function(m) file.path(D, sprintf("gm12878_%s_fc.bigWig", tolower(m)))

gtf <- import(file.path(D, "gencode.v46.basic.annotation.gtf.gz"))
gtf <- gtf[gtf$gene_type == "protein_coding"]

# Signal of one mark in a region, averaged in `nbin` equal bins (keeps plots light)
binned <- function(m, region, nbin = 600) {
  x <- import(bw_file(m), which = region, as = "RleList")[[as.character(seqnames(region))]]
  v <- as.numeric(x[start(region):end(region)])
  b <- cut(seq_along(v), nbin, labels = FALSE)
  data.frame(mark = m, x = start(region) + tapply(seq_along(v), b, mean) - 1,
             y = tapply(v, b, mean))
}

gene_panel <- function(region, xl, max_genes = 6, cluster_label = NULL) {
  tx <- subsetByOverlaps(gtf[gtf$type == "transcript"], region)
  tx <- tx[order(-width(tx))]; tx <- tx[!duplicated(tx$gene_name)]
  tx <- tx[order(start(tx))]
  if (length(tx) > max_genes) tx <- tx[seq_len(max_genes)]
  tx$row <- (seq_along(tx) - 1) %% 2 + 1                      # two rows, alternating
  ex <- gtf[gtf$type == "exon" & gtf$transcript_id %in% tx$transcript_id]
  exdf <- data.frame(s = start(ex), e = end(ex), row = tx$row[match(ex$transcript_id, tx$transcript_id)])
  txdf <- data.frame(s = pmax(start(tx), xl[1]), e = pmin(end(tx), xl[2]), row = tx$row,
                     g = paste0(tx$gene_name, ifelse(strand(tx) == "+", " →", " ←")))
  ggplot() +
    geom_segment(data = txdf, aes(x = s, xend = e, y = -row, yend = -row), colour = "#1b2a3a") +
    geom_rect(data = exdf, aes(xmin = s, xmax = e, ymin = -row - 0.22, ymax = -row + 0.22), fill = "#1b2a3a") +
    { if (is.null(cluster_label))
        geom_text(data = txdf, aes(x = (s + e) / 2, y = -row - 0.55, label = g), size = 4.2, vjust = 1)
      else annotate("text", x = mean(c(min(txdf$s), max(txdf$e))), y = -2.75, label = cluster_label, size = 4.2) } +
    scale_x_continuous(limits = xl, expand = c(0, 0), labels = function(b) sprintf("%.2f Mb", b / 1e6)) +
    scale_y_continuous(breaks = NULL, limits = c(-3.2, -0.5)) +
    labs(x = as.character(seqnames(region)), y = NULL) +
    theme_minimal(base_size = 13) + theme(panel.grid = element_blank(), plot.margin = margin(5, 30, 5, 5))
}

tracks <- function(region, title, ymax = NULL, cluster_label = NULL) {
  df <- do.call(rbind, lapply(names(marks), binned, region = region))
  df$mark <- factor(df$mark, levels = names(marks))
  xl <- c(start(region), end(region))
  p <- ggplot(df, aes(x, y, fill = mark)) +
    geom_area(colour = NA) +
    geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
    facet_grid(mark ~ ., switch = "y", scales = if (is.null(ymax)) "free_y" else "fixed") +
    scale_fill_manual(values = marks, guide = "none") +
    scale_x_continuous(limits = xl, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = 3) +
    labs(title = title, y = NULL, x = NULL) +
    theme_minimal(base_size = 13) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
          strip.placement = "outside", axis.text.x = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold"))
  if (!is.null(ymax)) p <- p + coord_cartesian(ylim = c(0, ymax))
  p / gene_panel(region, xl, max_genes = if (is.null(cluster_label)) 6 else 20, cluster_label = cluster_label) + plot_layout(heights = c(6, 1.3))
}

# --- fig11: an active, highly expressed B-cell gene --------------------------
active <- GRanges("chr5", IRanges(150385001, 150430000))       # CD74 and neighbours
p11 <- tracks(active, "CD74 in GM12878 (B cells): an active gene")
ggsave("fig11-tracks-active-gene.png", p11, width = 12, height = 6.2, dpi = 150, bg = "white")

# --- fig12: Polycomb-silenced HOXD cluster and an H3K9me3-marked ZNF cluster ---
hox <- GRanges("chr2",  IRanges(176060001, 176220000))
znf <- GRanges("chr19", IRanges(21700001, 22030000))
p12 <- wrap_elements(tracks(hox, "HOXD cluster: Polycomb (H3K27me3)", ymax = 10,
                            cluster_label = "HOXD13 … HOXD1 (and EVX2), all silent")) |
       wrap_elements(tracks(znf, "ZNF cluster, chr19: H3K9me3", ymax = 10))
ggsave("fig12-tracks-silent.png", p12, width = 15, height = 6.4, dpi = 150, bg = "white")

# Numbers for the slides: mean fold change per mark over each region
for (r in list(CD74 = active, HOXD = hox, ZNF = znf)) {
  print(sapply(names(marks), function(m) round(mean(binned(m, r)$y), 2)))
}
