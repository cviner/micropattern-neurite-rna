# Methods design check: ~ animal + compartment versus ~ batch + compartment, for neurite versus soma
# and proximal versus distal. Samples keep the count-matrix column order, animal levels are 1..12,
# compartment levels are alphabetical, and the shrunk coefficient is the last one.

source("R/common.R")
require_pins(c("DESeq2", "apeglm"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

counts <- load_counts()

sample_table <- function(pattern, neurite_label, other_label, other_pattern) {
  libs <- colnames(counts)[grepl(pattern, colnames(counts))]
  animal <- culture_of(libs)
  data.frame(
    sample = libs,
    animal = factor(animal),
    compartment = factor(ifelse(grepl(other_pattern, libs), other_label, neurite_label)),
    batch = factor(ifelse(animal <= 8, "B1", "B2")),
    stringsAsFactors = FALSE
  )
}

run <- function(samples, design, min_count, min_samples) {
  dds <- DESeqDataSetFromMatrix(counts[, samples$sample], samples, design = as.formula(design))
  dds <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ]
  dds <- DESeq(dds, quiet = TRUE)
  coef <- resultsNames(dds)[length(resultsNames(dds))]
  res <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  list(coef = coef, tested = nrow(res), genes = rownames(res)[!is.na(res$padj) & res$padj < 0.05],
       lfc = setNames(res$log2FoldChange, rownames(res)))
}

q2a <- sample_table("_(Distal_Dendrites|Proximal_Dendrites|Soma)$", "Dendrite", "Soma", "_Soma$")
q2b <- sample_table("_(Distal_Dendrites|Proximal_Dendrites)$", "Proximal", "Distal", "_Distal_Dendrites$")

fits <- list(
  neurite_soma_batch = run(q2a, "~ batch + compartment", 3, 4),
  neurite_soma_animal = run(q2a, "~ animal + compartment", 3, 4),
  proximal_distal_batch = run(q2b, "~ batch + compartment", 1, 3),
  proximal_distal_animal = run(q2b, "~ animal + compartment", 1, 3)
)

summarise_pair <- function(b, a, prefix) {
  common <- intersect(names(b$lfc), names(a$lfc))
  setNames(list(
    a$tested, length(a$genes), length(b$genes), length(intersect(a$genes, b$genes)),
    length(setdiff(a$genes, b$genes)), length(setdiff(b$genes, a$genes)),
    round(cor(b$lfc[common], a$lfc[common], use = "pairwise.complete.obs"), 4)
  ), paste0(prefix, c("_tested", "_animal_de", "_batch_de", "_overlap", "_only_animal",
                      "_only_batch", "_lfc_r")))
}

kv <- c(
  summarise_pair(fits$neurite_soma_batch, fits$neurite_soma_animal, "neurite_soma"),
  summarise_pair(fits$proximal_distal_batch, fits$proximal_distal_animal, "proximal_distal")
)

primary <- file.path(RESULTS, "deseq2_neurite_vs_soma_significant.tsv")
if (file.exists(primary)) {
  p <- read.delim(primary)$gene
  a <- fits$neurite_soma_animal$genes
  kv <- c(kv, list(
    neurite_soma_primary_de = length(p),
    neurite_soma_animal_minus_primary = length(a) - length(p),
    neurite_soma_animal_primary_shared = length(intersect(a, p)),
    neurite_soma_animal_primary_identical = setequal(a, p)
  ))
}
write_kv("design_check", kv)

genes <- do.call(rbind, lapply(names(fits), function(n) {
  data.frame(fit = n, coefficient = fits[[n]]$coef, gene = fits[[n]]$genes)
}))
write_table(genes, "design_check_de_genes.tsv")
