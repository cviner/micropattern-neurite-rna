# Supplementary Files 1, 2 and 4: DESeq2 with apeglm shrinkage, BH-adjusted P < 0.05.

source("R/common.R")
require_pins(c("DESeq2", "apeglm", "readxl"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

counts <- load_counts()

shrunk_table <- function(dds, coef) {
  res <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  out <- data.frame(gene = rownames(res), as.data.frame(res), row.names = NULL)
  out[order(out$log2FoldChange, decreasing = TRUE), ]
}

significant <- function(tab) tab[!is.na(tab$padj) & tab$padj < 0.05, ]

# File 1: micropatterned (Whole) versus standard (Control), 7_Whole and 10_Whole excluded.
file1 <- function() {
  whole <- setdiff(libraries(counts, "Whole"), c("7_Whole", "10_Whole"))
  control <- libraries(counts, "Control")
  mat <- counts[, c(whole, control)]
  coldata <- data.frame(
    type = factor(rep(c("micropattern", "control"), c(length(whole), length(control))),
                  levels = c("control", "micropattern")),
    batch = factor(batch_of(colnames(mat)), levels = c("Batch2", "Batch1")),
    row.names = colnames(mat)
  )
  dds <- DESeqDataSetFromMatrix(mat, coldata, design = ~ batch + type)
  dds <- dds[rowSums(counts(dds) >= 5) >= 3, ]
  dds$type <- relevel(dds$type, ref = "control")
  dds$batch <- relevel(dds$batch, ref = "Batch2")
  dds <- DESeq(dds, minReplicatesForReplace = 7L, quiet = TRUE)
  shrunk_table(dds, "type_micropattern_vs_control")
}

# Files 2 and 4 use alphabetical sampleID levels (X1, X10, X11, X12, X2, ...),
# which fix the design-matrix column order.
paired_design <- function(libs, compartment, ref, prefilter_count, prefilter_n) {
  mat <- counts[, libs]
  coldata <- data.frame(
    sampleID = factor(sprintf("X%s", culture_of(libs))),
    compartment = factor(compartment),
    row.names = libs
  )
  dds <- DESeqDataSetFromMatrix(mat, coldata, design = ~ sampleID + compartment)
  dds <- dds[rowSums(counts(dds) >= prefilter_count) >= prefilter_n, ]
  dds$compartment <- relevel(dds$compartment, ref = ref)
  DESeq(dds, quiet = TRUE)
}

file2 <- function() {
  neurite <- colnames(counts)[grepl("_Dendrites$", colnames(counts))]
  neurite <- neurite[order(culture_of(neurite))]
  soma <- libraries(counts, "Soma")
  dds <- paired_design(c(neurite, soma),
                       rep(c("dendrites", "soma"), c(length(neurite), length(soma))),
                       "soma", 3, 4)
  shrunk_table(dds, "compartment_dendrites_vs_soma")
}

file4 <- function() {
  proximal <- libraries(counts, "Proximal_Dendrites")
  distal <- libraries(counts, "Distal_Dendrites")
  dds <- paired_design(c(proximal, distal),
                       rep(c("proximal", "distal"), c(length(proximal), length(distal))),
                       "proximal", 1, 3)
  shrunk_table(dds, "compartment_distal_vs_proximal")
}

f1 <- file1()
f2 <- file2()
f4 <- file4()
s1 <- significant(f1)
s2 <- significant(f2)
s4 <- significant(f4)

write_table(f1, "deseq2_micropattern_vs_control_all.tsv")
write_table(s1, "deseq2_micropattern_vs_control_significant.tsv")
write_table(f2, "deseq2_neurite_vs_soma_all.tsv")
write_table(s2, "deseq2_neurite_vs_soma_significant.tsv")
write_table(f4, "deseq2_distal_vs_proximal_all.tsv")
write_table(s4, "deseq2_distal_vs_proximal_significant.tsv")

compare_sheet <- function(file, sheet, computed) {
  pub <- read_sheet(file, sheet)
  names(pub)[1:3] <- c("gene", "lfc", "padj")
  pub$lfc <- as.numeric(pub$lfc)
  pub$padj <- as.numeric(pub$padj)
  m <- merge(pub, computed, by = "gene", suffixes = c("_pub", ""))
  data.frame(
    file = file, sheet = sheet,
    published_n = nrow(pub), computed_n = nrow(computed),
    shared = length(intersect(pub$gene, computed$gene)),
    only_published = length(setdiff(pub$gene, computed$gene)),
    only_computed = length(setdiff(computed$gene, pub$gene)),
    sets_equal = setequal(pub$gene, computed$gene),
    median_abs_lfc_diff = median(abs(m$lfc - m$log2FoldChange)),
    max_abs_lfc_diff = max(abs(m$lfc - m$log2FoldChange)),
    genes_lfc_within_1e4 = sum(abs(m$lfc - m$log2FoldChange) <= 1e-4),
    median_rel_padj_diff = median(abs(m$padj_pub - m$padj) / m$padj),
    max_rel_padj_diff = max(abs(m$padj_pub - m$padj) / m$padj)
  )
}

check <- rbind(
  compare_sheet("SuppFile1(Micropattern_Control).xlsx", "Micropatterned vs Control", s1),
  compare_sheet("SuppFile2(Dendrites_Soma).xlsx", "Neurite DE", s2[s2$log2FoldChange > 0, ]),
  compare_sheet("SuppFile2(Dendrites_Soma).xlsx", "Soma DE", s2[s2$log2FoldChange < 0, ]),
  compare_sheet("SuppFile4(Distal_Proximal).xlsx", "Distal DE", s4[s4$log2FoldChange > 0, ]),
  compare_sheet("SuppFile4(Distal_Proximal).xlsx", "Proximal DE", s4[s4$log2FoldChange < 0, ])
)
write_table(check, "deseq2_published_check.tsv")
print(check)

write_kv("deseq2", list(
  micropattern_de = nrow(s1),
  micropattern_up = sum(s1$log2FoldChange > 0),
  micropattern_down = sum(s1$log2FoldChange < 0),
  micropattern_tested = nrow(f1),
  micropattern_na_padj = sum(is.na(f1$padj)),
  neurite_soma_tested = nrow(f2),
  neurite_higher = sum(s2$log2FoldChange > 0),
  soma_higher = sum(s2$log2FoldChange < 0),
  neurite_soma_de = nrow(s2),
  distal_proximal_tested = nrow(f4),
  distal_higher = sum(s4$log2FoldChange > 0),
  proximal_higher = sum(s4$log2FoldChange < 0),
  published_sets_equal = all(check$sets_equal),
  published_max_abs_lfc_diff = max(check$max_abs_lfc_diff)
))
