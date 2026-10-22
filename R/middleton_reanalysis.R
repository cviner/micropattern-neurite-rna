# Supplementary Note 1: Middleton et al. (2019) dendrite-versus-soma counts reanalyzed with this
# study's DESeq2 settings (betaPrior = FALSE, apeglm shrinkage):
# design ~ neuronID + subcellularCompartment, genes with >= 1 count in >= 16 of 32 libraries,
# soma reference, BH-adjusted P < 0.05 and positive shrunk log2 fold change.
#
# Counts are the 3' UTR file GSE115480_combined_3utr_counts_per_gene.txt.gz.

source("R/common.R")
require_pins(c("DESeq2", "apeglm", "readxl"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

GEO <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115480/suppl/"
DOWNLOADS <- c(
  GSE115480_combined_3utr_counts_per_gene.txt.gz = "466e45601de382a1b79aa2c6e6218277e215a615ff7bf7c2b7fd10ade6130f26"
)
cache <- file.path(tempdir(), "GSE115480")
dir.create(cache, showWarnings = FALSE)
for (f in names(DOWNLOADS)) {
  path <- file.path(cache, f)
  if (!file.exists(path)) download.file(paste0(GEO, f), path, mode = "wb", quiet = TRUE)
  sha <- as.character(openssl::sha256(readBin(path, "raw", file.size(path))))
  if (sha != DOWNLOADS[[f]]) stop(f, ": sha256 ", sha, " does not match ", DOWNLOADS[[f]])
}

utr <- read.delim2(file.path(cache, "GSE115480_combined_3utr_counts_per_gene.txt.gz"),
                   header = TRUE, sep = "\t", check.names = FALSE)
utr <- utr[, nzchar(colnames(utr))]
utr <- as.matrix(utr)
storage.mode(utr) <- "integer"

s5 <- readxl::read_excel(file.path(EXTERNAL, "Middleton2019_MOESM6.xlsx"),
                         sheet = "TableS5-deDend genes", col_names = FALSE, .name_repair = "minimal")[[1]]

dendrite_higher <- function(mat) {
  stopifnot(identical(sub("^[A-Z][0-9]+", "", colnames(mat)), rep(c("D", "S"), 16)))
  coldata <- data.frame(
    subcellularCompartment = factor(rep(c("Dendrites", "Soma"), 16)),
    neuronID = factor(sub(".$", "", colnames(mat))),
    row.names = colnames(mat)
  )
  dds <- DESeqDataSetFromMatrix(mat, coldata, design = ~ neuronID + subcellularCompartment)
  dds <- dds[rowSums(counts(dds) >= 1) >= 16, ]
  dds$subcellularCompartment <- relevel(dds$subcellularCompartment, ref = "Soma")
  dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)
  res <- lfcShrink(dds, coef = "subcellularCompartment_Dendrites_vs_Soma", type = "apeglm", quiet = TRUE)
  tab <- data.frame(gene = rownames(res), as.data.frame(res), row.names = NULL)
  list(tested = nrow(dds), table = tab,
       up = tab$gene[!is.na(tab$padj) & tab$padj < 0.05 & tab$log2FoldChange > 0])
}

u <- dendrite_higher(utr)
up <- u$table[u$table$gene %in% u$up, ]
write_table(up[order(up$log2FoldChange, decreasing = TRUE), ], "middleton_reanalysis_dendrite_higher.tsv")

ours <- file.path(RESULTS, "deseq2_neurite_vs_soma_significant.tsv")
kv <- list(
  utr_libraries = ncol(utr),
  utr_genes = nrow(utr),
  utr_tested = u$tested,
  dendrite_higher = length(u$up),
  s5_genes = length(s5),
  s5_retained = length(intersect(u$up, s5)),
  s5_not_retained = length(setdiff(s5, u$up))
)
if (file.exists(ours)) {
  o <- read.delim(ours)
  neurite <- o$gene[o$log2FoldChange > 0]
  kv$overlap_with_neurite_higher <- length(intersect(u$up, neurite))
  kv$s5_overlap_with_neurite_higher <- length(intersect(s5, neurite))
}
write_kv("middleton_reanalysis", kv)
