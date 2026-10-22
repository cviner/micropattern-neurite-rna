# Methods note on the batch-stratified micropattern DESeq2: Cook's distances of the genes whose
# BH-adjusted P is NA, micropatterned (Whole) versus control libraries.

source("R/common.R")
require_pins(c("DESeq2"))
suppressPackageStartupMessages(library(DESeq2))

counts <- load_counts()
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

res <- results(dds)
cooks <- assays(dds)[["cooks"]]
na_genes <- rownames(res)[is.na(res$padj)]
ok_genes <- rownames(res)[!is.na(res$padj)]
mp <- colnames(cooks) %in% whole
ct <- colnames(cooks) %in% control

max_sample <- colnames(cooks)[apply(cooks[na_genes, , drop = FALSE], 1, which.max)]
na_cooks_mp <- mean(cooks[na_genes, mp])
na_cooks_ct <- mean(cooks[na_genes, ct])
ok_cooks_mp <- mean(cooks[ok_genes, mp])
ok_cooks_ct <- mean(cooks[ok_genes, ct])

# NA padj arises from the Cook's outlier rule (NA pvalue) or from independent filtering.
cutoff <- qf(0.99, length(resultsNames(dds)), ncol(dds) - length(resultsNames(dds)))
cooks_flagged <- sum(is.na(res$pvalue))
filtered <- sum(!is.na(res$pvalue) & is.na(res$padj))

write_kv("cooks", list(
  tested = nrow(res),
  na_padj = length(na_genes),
  na_padj_pct = round(100 * length(na_genes) / nrow(res), 1),
  na_from_cooks_outlier = cooks_flagged,
  na_from_independent_filtering = filtered,
  cooks_cutoff_F99 = round(cutoff, 4),
  mean_cooks_na_micropattern = round(na_cooks_mp, 6),
  mean_cooks_na_control = round(na_cooks_ct, 6),
  mean_cooks_ratio = round(na_cooks_mp / na_cooks_ct, 4),
  mean_cooks_ratio_2dp = round(na_cooks_mp / na_cooks_ct, 2),
  mean_cooks_nonNA_ratio = round(ok_cooks_mp / ok_cooks_ct, 4),
  max_cooks_sample_micropattern = sum(max_sample %in% whole),
  max_cooks_sample_control = sum(max_sample %in% control),
  micropattern_libraries = sum(mp),
  control_libraries = sum(ct)
))
