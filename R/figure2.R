# Figure 2, panels a-d, from the tables python/figure2.py writes to results/figure2/.
suppressMessages({library(ggplot2); library(patchwork); library(ggrepel)})

src <- file.path("results", "figure2")
RED <- "#C0392B"; BLUE <- "#5DADE2"; SLATE <- "#34495E"; TEAL <- "#0E8A6B"
th <- theme_minimal(base_size = 8) + theme(
  panel.grid.minor = element_blank(),
  plot.title = element_text(face = "bold", size = 8.2, margin = margin(b = 2)),
  plot.subtitle = element_text(size = 6.6, color = "grey30", margin = margin(b = 3)),
  axis.title = element_text(size = 7.2), axis.text = element_text(size = 6.6),
  plot.margin = margin(6, 9, 4, 6))
read_tsv <- function(name) read.table(file.path(src, name), header = TRUE, sep = "\t")
as_flag <- function(x) x %in% c("True", TRUE)

s <- read_tsv("scatter.tsv"); s$is54 <- as_flag(s$is54)
proximal_examples <- c("Gria2", "Dclk1", "Hnrnpa2b1", "Nptn")
smfish <- c("Tubb4b", "Gabarapl1", "Tnrc6a", "Ddah1")
A <- ggplot(s, aes(prox, dist)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(data = subset(s, !is54), color = "grey80", size = 0.3, alpha = 0.4) +
  geom_point(data = subset(s, is54), color = RED, size = 1.3, alpha = 0.95) +
  annotate("text", x = 2.9, y = 3.2, label = "equal abundance (y = x)", angle = 45, size = 2.1, color = "grey50") +
  annotate("text", x = 0.05, y = 3.9, hjust = 0, size = 2.4, color = "grey30",
           label = "most transcripts are lower distally\n(median \u22120.34)") +
  annotate("text", x = 3.95, y = 0.13, hjust = 1, size = 2.1, color = RED,
           label = "54 higher distally (median +1.18)") +
  geom_point(data = subset(s, gene %in% proximal_examples), color = BLUE, size = 1.4, alpha = 0.95) +
  geom_text_repel(data = subset(s, gene %in% smfish), aes(label = gene), color = RED,
                  fontface = "bold.italic", size = 1.85, box.padding = 0.7, point.padding = 0.4,
                  min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf, seed = 1) +
  geom_text_repel(data = subset(s, gene %in% proximal_examples), aes(label = gene), color = BLUE,
                  fontface = "bold.italic", size = 1.85, box.padding = 0.5, min.segment.length = 0,
                  segment.size = 0.2, max.overlaps = Inf, seed = 1) +
  coord_fixed(xlim = c(0, 4), ylim = c(0, 4)) +
  labs(title = "Distal versus proximal relative abundance",
       subtitle = "12,396 expressed genes; red = 54 distal-higher candidates; blue = proximal-higher examples",
       x = expression(log[10] ~ "(mean proximal CPM + 1)"), y = expression(log[10] ~ "(mean distal CPM + 1)")) + th

nulls <- read_tsv("nulls.tsv")
b <- data.frame(k = nulls$count, null = ifelse(nulls$null == "per_set_rate", "g", "s"))
B <- ggplot(b, aes(k, fill = null)) +
  geom_histogram(aes(y = after_stat(density) * 5000), binwidth = 1, position = "identity", alpha = 0.65, color = NA) +
  scale_fill_manual(values = c(g = "grey55", s = BLUE), guide = "none") +
  annotate("text", x = 14, y = 620, hjust = 0, size = 2.4, color = "grey35",
           label = "pair-preserving null\n(about 7 expected)") +
  annotate("text", x = 70, y = 470, hjust = 0.5, size = 2.4, color = "#2471A3",
           label = "naive symmetric null\n(inappropriate comparator;\nabout 70)") +
  geom_vline(xintercept = 54, color = RED, linewidth = 0.8) +
  annotate("point", x = 54, y = 0, color = RED, size = 1.8) +
  annotate("text", x = 52, y = Inf, vjust = 1.5, hjust = 1, label = "observed = 54", color = RED,
           fontface = "bold", size = 2.8) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 54, 70, 90), limits = c(-1, 95)) +
  labs(title = "The 54 candidates exceed the null expectation",
       subtitle = expression(italic(P) * " = 5 x 10"^-5 * " (20,000 permutations)"),
       x = "genes higher distally in \u2265 10 of 12 matched sets (by chance)", y = "null draws") + th

ab <- read_tsv("abundance.tsv")
ab$group <- factor(ab$group, levels = c("universe", "54"),
                   labels = c("all expressed\n(n = 12,396)", "54 distal-higher\ncandidates"))
meds <- tapply(ab$cpm, ab$group, median)
med <- data.frame(group = levels(ab$group), y = log10(meds + 1),
                  lab = c(sprintf("%.1f", meds[[1]]), sprintf("%.0f", meds[[2]])))
C <- ggplot(ab, aes(group, log10_cpm, fill = group)) +
  geom_violin(scale = "width", color = NA, alpha = 0.85) +
  geom_point(data = med, aes(group, y), inherit.aes = FALSE, size = 1.4) +
  geom_text(data = med, aes(group, y, label = paste("median", lab)), inherit.aes = FALSE,
            nudge_x = c(-0.36, 0.42), size = 2.4, color = "grey15") +
  scale_fill_manual(values = c("grey75", RED), guide = "none") +
  annotate("text", x = 1.5, y = 3.95, size = 2.3, color = "grey30",
           label = "detected in \u2265 11/12 distal samples (45/54);\ndepth-matching retains 20\u201342 (median 30); pseudo-distal overlap \u2264 1") +
  coord_cartesian(ylim = c(0, 4), clip = "off") +
  labs(title = "Candidate transcripts have higher normalized abundance",
       subtitle = "neurite-mean CPM (12 proximal + 12 distal)",
       x = NULL, y = expression(log[10] ~ "(neurite-mean CPM + 1)")) + th

cd <- read_tsv("categories.tsv"); cd$in_mid <- as_flag(cd$in_middleton_s4)
cd$category <- factor(cd$category, levels = c("Not assigned", "RNA-regulatory", "Cytoskeletal",
                                              "Signaling/synaptic", "Mitochondrial", "Translation/ribosome"))
set.seed(7); cd$cy <- as.numeric(cd$category) + runif(nrow(cd), -0.17, 0.17)
labc <- cd[cd$gene %in% smfish, ]; labx <- cd[cd$gene %in% c("Cox6c", "Rpl19", "Timm8b"), ]
D <- ggplot(cd, aes(log2fc, cy)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
  geom_point(aes(fill = in_mid), size = 1.7, shape = 21, color = "grey25", stroke = 0.35) +
  geom_point(data = labc, size = 2.7, shape = 21, fill = NA, color = TEAL, stroke = 1.0) +
  geom_text_repel(data = labx, aes(label = gene), size = 2, color = "grey35", fontface = "italic",
                  box.padding = 0.45, point.padding = 0.3, min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = Inf, seed = 1) +
  geom_text_repel(data = labc, aes(label = gene), size = 2.15, color = TEAL, fontface = "bold.italic",
                  box.padding = 0.6, point.padding = 0.4, min.segment.length = 0, segment.size = 0.25,
                  segment.color = TEAL, max.overlaps = Inf, seed = 2) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = SLATE), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.11))) +
  scale_y_continuous(breaks = 1:6, labels = levels(cd$category), expand = expansion(add = 0.45)) +
  labs(title = "Functional annotations of the 54 candidates",
       subtitle = expression(atop("filled = in a prior dendritic catalog (29/54; 1.4-fold, " * italic(P) * " = 0.016);",
                                  "teal rings = candidates prioritized for smFISH")),
       x = expression(log[2] ~ "[(mean distal CPM + 1) / (mean proximal CPM + 1)]"), y = NULL) +
  th + theme(axis.text.y = element_text(size = 6.8))

fig <- A + B + C + D + plot_layout(design = "AABB\nCCDD", heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 10))
dir.create("figures", showWarnings = FALSE)
ggsave(file.path("figures", "Figure2.pdf"), fig, width = 7.2, height = 6.5)
ggsave(file.path("figures", "Figure2.png"), fig, width = 7.2, height = 6.5, dpi = 300, bg = "white")
cat("figures/Figure2.pdf\n")
