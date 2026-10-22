# Supplementary Figure 1 (Whole versus Control: CPM Pearson correlation heatmap and mean-CPM
# scatter) and Supplementary Figures 2 and 3 (hierarchical clustering of DESeq2
# median-of-ratios normalized counts, ln(x + 1), Euclidean distance, complete linkage).
# Figure 2 is written at 3000 x 2400 px and Figure 3 at 2700 x 2400 px, both at 300 ppi;
# the PDFs use the same aspect ratios (10 x 8 in and 9 x 8 in).

source("R/common.R")
require_pins(c("DESeq2", "ggplot2", "pheatmap"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

counts <- load_counts()
c_cpm <- cpm(counts)
label <- function(libs) paste(culture_of(libs), compartment_of(libs))

# Supplementary Figure 1A: all 12 Whole and 11 Control libraries, linear CPM.
wc <- c(libraries(counts, "Whole"), libraries(counts, "Control"))
r <- cor(c_cpm[, wc])
dimnames(r) <- list(label(wc), label(wc))
heat <- function() {
  pheatmap(r, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
           number_format = "%.2f", fontsize_number = 6, number_color = "black",
           breaks = seq(0, 1, length.out = 101),
           color = colorRampPalette(c("#3B4CC0", "#DDDDDD", "#B40426"))(100),
           border_color = "white", angle_col = 90, legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
           main = "Pairwise Pearson correlation of whole-neuron transcriptomes (CPM)", silent = TRUE)
}
for (ext in c("pdf", "png")) {
  path <- file.path(FIGURES, paste0("supplementary_figure1a_correlation_heatmap.", ext))
  if (ext == "pdf") pdf(path, width = 10, height = 9.5) else png(path, width = 10, height = 9.5, units = "in", res = 200)
  grid::grid.draw(heat()$gtable)
  dev.off()
}
offdiag <- r[upper.tri(r)]
mean_r <- function(lib) mean(r[label(lib), setdiff(label(wc), label(lib))])

# Supplementary Figure 1B: log10 mean CPM per condition, genes with nonzero mean in both.
mw <- rowMeans(c_cpm[, libraries(counts, "Whole")])
mc <- rowMeans(c_cpm[, libraries(counts, "Control")])
both <- mw > 0 & mc > 0
sc <- data.frame(whole = log10(mw[both]), control = log10(mc[both]))
p <- ggplot(sc, aes(whole, control)) +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = "#E41A1C") +
  coord_fixed() +
  labs(x = expression(Log[10] ~ "mean CPM (micropatterned, Whole)"),
       y = expression(Log[10] ~ "mean CPM (standard culture, Control)")) +
  theme_bw(base_size = 12)
ggsave(file.path(FIGURES, "supplementary_figure1b_cpm_scatter.pdf"), p, width = 6, height = 6)
ggsave(file.path(FIGURES, "supplementary_figure1b_cpm_scatter.png"), p, width = 6, height = 6, dpi = 200)

# Figure 2 orders neurite libraries by culture, then soma; Figure 3 orders proximal libraries, then distal.
display <- c(Distal_Dendrites = "Distal neurites", Proximal_Dendrites = "Proximal neurites", Soma = "Soma")
colours <- c("Distal neurites" = "#1B7837", "Proximal neurites" = "#2166AC", "Soma" = "#B2182B")
# Drawing only. cluster() keeps the DESeq2 size-factor and hclust computation unchanged.
# px is the PNG size in pixels; both devices are 300 ppi, so the PDF page
# has the same aspect ratio.
draw_dendrogram <- function(hc, title, subtitle, px) {
  comp <- sub(" [0-9]+$", "", hc$labels)
  present <- colours[names(colours) %in% comp]
  labs <- hc$labels[hc$order]
  n <- length(labs)
  W <- px[1] / 300
  H <- px[2] / 300
  shrink <- function(label, cex, limit, font = 1, floor = 0.62) {
    while (cex > floor && strwidth(label, units = "inches", cex = cex, font = font) > limit)
      cex <- cex - 0.02
    cex
  }
  edge <- 0.20
  right_in <- 0.30
  title_cex <- 1.22
  sub_cex <- 0.88
  lab_cex <- 1.00
  axis_cex <- 0.90
  ylab_cex <- 0.98
  leg_cex <- 0.96
  dot_cex <- 1.20
  ylab <- "Height (Euclidean distance)"
  gap_title_sub <- 0.16
  gap_sub_plot <- 0.15
  gap_tree_lab <- 0.11
  gap_lab_dot <- 0.09
  gap_dot_leg <- 0.11
  gap_ylab_tick <- 0.16
  gap_tick_axis <- 0.10

  par(mar = c(0, 0, 0, 0), xpd = FALSE)
  plot.new()
  title_cex <- shrink(title, title_cex, W - 2 * edge, font = 2, floor = 0.85)
  title_h <- strheight(title, units = "inches", cex = title_cex, font = 2)
  tick_w <- strwidth("350", units = "inches", cex = axis_cex)
  ylab_h <- strheight(ylab, units = "inches", cex = ylab_cex)
  left_in <- edge + ylab_h + gap_ylab_tick + tick_w + gap_tick_axis
  plot_w <- W - left_in - right_in
  sub_cex <- shrink(subtitle, sub_cex, plot_w - 0.20, floor = 0.64)
  if (strwidth(subtitle, units = "inches", cex = sub_cex) > plot_w - 0.05)
    stop("subtitle does not fit on one line")
  lab_h <- max(strheight(labs, units = "inches", cex = lab_cex))
  while (lab_h > (plot_w / n) * 0.68 && lab_cex > 0.72) {
    lab_cex <- lab_cex - 0.02
    lab_h <- max(strheight(labs, units = "inches", cex = lab_cex))
  }
  if (lab_h > (plot_w / n) * 0.85) stop("leaf labels would overlap")
  sub_h <- strheight(subtitle, units = "inches", cex = sub_cex)
  lab_len <- max(strwidth(labs, units = "inches", cex = lab_cex))
  leg_h <- strheight("Ag", units = "inches", cex = leg_cex)
  dot_r <- 0.046 * dot_cex
  y_title_c <- H - edge - title_h / 2
  y_sub_c <- y_title_c - title_h / 2 - gap_title_sub - sub_h / 2
  y_plot_top <- y_sub_c - sub_h / 2 - gap_sub_plot
  below_y0 <- gap_tree_lab + lab_len + gap_lab_dot + 2 * dot_r + gap_dot_leg + leg_h + edge
  # C_dendwindow leaves a small pad below height 0. Reserve the annotations
  # as if that pad were thin, so a larger pad only adds whitespace.
  pad_frac <- 0.015
  y_plot_bot <- (below_y0 - pad_frac * y_plot_top) / (1 - pad_frac)
  plot_h <- y_plot_top - y_plot_bot
  if (plot_h < 3) stop("dendrogram plot region is too short")
  ylab_len <- strwidth(ylab, units = "inches", cex = ylab_cex)
  if (ylab_len > plot_h * 0.86) {
    ylab_cex <- ylab_cex * (plot_h * 0.86) / ylab_len
    ylab_h <- strheight(ylab, units = "inches", cex = ylab_cex)
    left_in <- edge + ylab_h + gap_ylab_tick + tick_w + gap_tick_axis
  }
  par(mai = c(y_plot_bot, left_in, H - y_plot_top, right_in))
  # The measurement plot.new() above already opened this page. Reuse it.
  par(new = TRUE)
  plot(hc, labels = FALSE, hang = -1, main = "", sub = "", xlab = "", ylab = "",
       axes = FALSE, ann = FALSE)
  axis(2, at = pretty(hc$height), cex.axis = axis_cex, las = 1, mgp = c(3, 0.6, 0), xpd = FALSE)
  usr <- par("usr")
  x_mid_in <- grconvertX(mean(usr[1:2]), from = "user", to = "inches")
  text(grconvertX(x_mid_in, "inches", "user"), grconvertY(y_title_c, "inches", "user"),
       title, cex = title_cex, font = 2, adj = c(0.5, 0.5), xpd = NA)
  text(grconvertX(x_mid_in, "inches", "user"), grconvertY(y_sub_c, "inches", "user"),
       subtitle, cex = sub_cex, adj = c(0.5, 0.5), xpd = NA)
  y0_in <- grconvertY(0, from = "user", to = "inches")
  y_bot_in <- grconvertY(usr[3], from = "user", to = "inches")
  y_top_in <- grconvertY(usr[4], from = "user", to = "inches")
  text(grconvertX(edge + ylab_h / 2, "inches", "user"),
       grconvertY((y_bot_in + y_top_in) / 2, "inches", "user"),
       ylab, srt = 90, adj = c(0.5, 0.5), cex = ylab_cex, xpd = NA)
  y_lab_top <- y0_in - gap_tree_lab
  text(seq_len(n), grconvertY(y_lab_top, "inches", "user"),
       labels = labs, srt = 90, adj = c(1, 0.5), cex = lab_cex, xpd = NA)
  y_dot <- y_lab_top - lab_len - gap_lab_dot - dot_r
  points(seq_len(n), rep(grconvertY(y_dot, "inches", "user"), n),
         pch = 19, cex = dot_cex, col = colours[comp[hc$order]], xpd = NA)
  legend(x = mean(usr[1:2]), y = grconvertY(y_dot - dot_r - gap_dot_leg, "inches", "user"),
         xjust = 0.5, yjust = 1, horiz = TRUE, legend = names(present),
         pch = 19, col = present, bty = "n", cex = leg_cex, pt.cex = dot_cex, xpd = NA)
}

cluster <- function(libs, min_count, min_samples, file, title, px) {
  dds <- DESeqDataSetFromMatrix(counts[, libs], data.frame(row.names = libs), design = ~ 1)
  dds <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ]
  dds <- estimateSizeFactors(dds)
  logc <- log(counts(dds, normalized = TRUE) + 1)
  colnames(logc) <- paste(display[compartment_of(libs)], culture_of(libs))
  hc <- hclust(dist(t(logc), method = "euclidean"), method = "complete")
  subtitle <- sprintf("ln(DESeq2-normalized count + 1), \u2265%d count%s in \u2265%d samples (%s genes); Euclidean distance, complete linkage",
                      min_count, if (min_count > 1) "s" else "", min_samples,
                      format(nrow(logc), big.mark = ",", scientific = FALSE, trim = TRUE))
  # The default pdf() device replaces U+2265 with a dot. cairo_pdf embeds a font that draws it.
  if (!nzchar(Sys.getenv("FONTCONFIG_FILE")) && file.exists("/opt/homebrew/etc/fonts/fonts.conf"))
    Sys.setenv(FONTCONFIG_FILE = "/opt/homebrew/etc/fonts/fonts.conf")
  cairo_pdf(file.path(FIGURES, paste0(file, ".pdf")), width = px[1] / 300, height = px[2] / 300, family = "sans")
  draw_dendrogram(hc, title, subtitle, px)
  dev.off()
  png(file.path(FIGURES, paste0(file, ".png")), width = px[1], height = px[2], units = "px", res = 300, bg = "white")
  draw_dendrogram(hc, title, subtitle, px)
  dev.off()
  list(genes = nrow(logc), order = hc$labels[hc$order], height = max(hc$height), hc = hc)
}

neurite <- colnames(counts)[grepl("_Dendrites$", colnames(counts))]
neurite <- neurite[order(culture_of(neurite))]
s2_libs <- c(neurite, libraries(counts, "Soma"))
s3_libs <- c(libraries(counts, "Proximal_Dendrites"), libraries(counts, "Distal_Dendrites"))
s2_title <- "Hierarchical clustering of soma, proximal and distal neurite samples"
s3_title <- "Hierarchical clustering of proximal and distal neurite samples"
f2 <- cluster(s2_libs, 3, 4, "supplementary_figure2_dendrogram", s2_title, c(3000, 2400))
f3 <- cluster(s3_libs, 1, 3, "supplementary_figure3_dendrogram", s3_title, c(2700, 2400))

leaf <- function(fig, f) data.frame(figure = fig, position = seq_along(f$order), leaf = f$order)
write_table(rbind(leaf("S2", f2), leaf("S3", f3)), "supplementary_figures_leaf_order.tsv")

write_kv("supplementary_figures", list(
  s1a_libraries = length(wc),
  s1a_min_offdiag_r = round(min(offdiag), 3),
  s1a_7_Whole_mean_r = round(mean_r("7_Whole"), 3),
  s1a_10_Whole_mean_r = round(mean_r("10_Whole"), 3),
  s1b_genes_detected_both = sum(both),
  s1b_log10_mean_cpm_r = round(cor(sc$whole, sc$control), 4),
  s2_libraries = length(f2$order),
  s2_genes = f2$genes,
  s2_max_height = round(f2$height, 1),
  s3_libraries = length(f3$order),
  s3_genes = f3$genes,
  s3_max_height = round(f3$height, 1)
))
