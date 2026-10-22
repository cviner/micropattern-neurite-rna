# GO enrichment from Methods: g:Profiler ORA (gost, mmusculus, g_SCS, 0.05, GO sources),
# preranked GSEA (clusterProfiler gseGO) and the Middleton-style hypergeometric enrichGO.
# Requires the DE tables written by R/deseq2.R.
#
# g:Profiler releases and backgrounds:
#   Supplementary Files 1, 2 and 4, all annotated genes: e112_eg59_p19
#     (https://biit.cs.ut.ee/gprofiler_archive3/e112_eg59_p19), domain_scope "annotated".
#   Supplementary File 2, neuron-expressed: e112_eg59_p19, genes passing the neurite-versus-soma
#     DESeq2 prefilter (>= 3 counts in >= 4 of the 36 libraries, i.e. >= 10% of that contrast's samples).
#   Supplementary File 3: e112_eg59_p19, genes passing the proximal-versus-distal DESeq2 prefilter
#     (>= 1 count in >= 3 of 24 libraries).
#   Supplementary File 4, neuron-expressed: e114_eg62_p19
#     (https://biit.cs.ut.ee/gprofiler_archive3/e114_eg62_p19, or https://biit.cs.ut.ee/gprofiler
#     when that server reports e114_eg62_p19), >= 3 counts in >= 6 of all 59 libraries.
# gprofiler2 replaces domain_scope with "custom" whenever custom_bg is supplied.

source("R/common.R")
require_pins(c("gprofiler2", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "GO.db", "readxl"))
suppressPackageStartupMessages({
  library(gprofiler2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(BiocParallel)
})
register(SerialParam())

E112 <- "https://biit.cs.ut.ee/gprofiler_archive3/e112_eg59_p19"
E114_ARCHIVE <- "https://biit.cs.ut.ee/gprofiler_archive3/e114_eg62_p19"
LIVE <- "https://biit.cs.ut.ee/gprofiler"

need <- file.path(RESULTS, c("deseq2_micropattern_vs_control_all.tsv", "deseq2_neurite_vs_soma_all.tsv",
                             "deseq2_distal_vs_proximal_all.tsv"))
if (!all(file.exists(need))) stop("Run R/deseq2.R first")
f1 <- read.delim(need[1])
f2 <- read.delim(need[2])
f4 <- read.delim(need[3])
sig <- function(t) t[!is.na(t$padj) & t$padj < 0.05, ]
s1 <- sig(f1)
s2 <- sig(f2)
s4 <- sig(f4)

counts <- load_counts()
prox <- libraries(counts, "Proximal_Dendrites")
dist <- libraries(counts, "Distal_Dendrites")

sets <- list(
  micropattern_up = s1$gene[s1$log2FoldChange > 0],
  neurite = s2$gene[s2$log2FoldChange > 0],
  soma = s2$gene[s2$log2FoldChange < 0],
  distal = s4$gene[s4$log2FoldChange > 0],
  proximal = s4$gene[s4$log2FoldChange < 0],
  consistent_union = union(rownames(counts)[rowSums(counts[, prox] >= 1) >= 11],
                           rownames(counts)[rowSums(counts[, dist] >= 1) >= 11])
)
backgrounds <- list(
  neurite_soma_prefilter = f2$gene,
  proximal_distal_prefilter = f4$gene,
  all59_ge3_in_10pct = rownames(counts)[rowSums(counts >= 3) >= ceiling(0.1 * ncol(counts))]
)

server_version <- function(url) {
  v <- tryCatch({
    set_base_url(url)
    gost("Actb", organism = "mmusculus", sources = "GO:CC", significant = FALSE)$meta$version
  }, error = function(e) NULL)
  if (is.null(v)) NA_character_ else v
}
e114_url <- function() {
  archived <- server_version(E114_ARCHIVE)
  if (!is.na(archived) && startsWith(archived, "e114_eg62_p19")) return(E114_ARCHIVE)
  live <- server_version(LIVE)
  if (is.na(live) || !startsWith(live, "e114_eg62_p19")) {
    warning("g:Profiler e114_eg62_p19 is neither archived nor live (live: ", live,
            "); the Supplementary File 4 neuron-expressed tables will not be reproduced")
  }
  LIVE
}
E114 <- e114_url()

gp <- function(query, url, bg = NULL) {
  set_base_url(url)
  for (attempt in 1:3) {
    r <- tryCatch(
      gost(query, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE,
           exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE,
           significant = FALSE, user_threshold = 0.05, correction_method = "g_SCS",
           domain_scope = if (is.null(bg)) "annotated" else "custom", custom_bg = bg,
           sources = "GO"),
      error = function(e) NULL)
    if (!is.null(r)) break
    Sys.sleep(30)
  }
  if (is.null(r)) stop("g:Profiler query failed at ", url)
  res <- r$result[r$result$significant, c("source", "term_id", "term_name", "p_value", "term_size",
                                          "query_size", "intersection_size", "effective_domain_size")]
  list(version = r$meta$version, terms = res[order(res$source, res$p_value), ])
}

runs <- list(
  micropattern_up_annotated = list(sets$micropattern_up, E112, NULL),
  neurite_annotated = list(sets$neurite, E112, NULL),
  soma_annotated = list(sets$soma, E112, NULL),
  distal_annotated = list(sets$distal, E112, NULL),
  proximal_annotated = list(sets$proximal, E112, NULL),
  neurite_neuronbg = list(sets$neurite, E112, backgrounds$neurite_soma_prefilter),
  soma_neuronbg = list(sets$soma, E112, backgrounds$neurite_soma_prefilter),
  consistent_union_neuronbg = list(sets$consistent_union, E112, backgrounds$proximal_distal_prefilter),
  distal_neuronbg_all59_e114 = list(sets$distal, E114, backgrounds$all59_ge3_in_10pct),
  proximal_neuronbg_all59_e114 = list(sets$proximal, E114, backgrounds$all59_ge3_in_10pct)
)
gpres <- lapply(names(runs), function(n) {
  message("g:Profiler ", n)
  gp(runs[[n]][[1]], runs[[n]][[2]], runs[[n]][[3]])
})
names(gpres) <- names(runs)

terms <- do.call(rbind, lapply(names(gpres), function(n) {
  t <- gpres[[n]]$terms
  if (!nrow(t)) return(NULL)
  cbind(analysis = n, version = gpres[[n]]$version, t)
}))
write_table(terms, "go_gprofiler_terms.tsv")

published <- list(
  micropattern_up_annotated = c("SuppFile1(Micropattern_Control).xlsx", "Micropattern DE GO", "GO Term ID", "P-Value"),
  neurite_annotated = c("SuppFile2(Dendrites_Soma).xlsx", "Neurite GO (All Coding Genes)", "GO Term ID", "P-Value"),
  soma_annotated = c("SuppFile2(Dendrites_Soma).xlsx", "Soma GO (All Coding Genes)", "GO Term ID", "P-Value"),
  distal_annotated = c("SuppFile4(Distal_Proximal).xlsx", "Distal GO (All Coding Genes)", "GO Term ID", "P-Value"),
  proximal_annotated = c("SuppFile4(Distal_Proximal).xlsx", "Proximal GO (All Coding Genes)", "GO Term ID", "P-Value"),
  neurite_neuronbg = c("SuppFile2(Dendrites_Soma).xlsx", "Neurite GO (Neuron-Expressed)", "Term ID", "Adjusted P-Value"),
  soma_neuronbg = c("SuppFile2(Dendrites_Soma).xlsx", "Soma GO (Neuron-Expressed)", "Term ID", "Adjusted P-Value"),
  consistent_union_neuronbg = c("SuppFile3(ConsDistal_ConsProximal).xlsx", "ConsNeurite GO (Neuron Genes)", "GO Term ID", "P-Value"),
  distal_neuronbg_all59_e114 = c("SuppFile4(Distal_Proximal).xlsx", "Distal GO (Neuron-Expressed)", "Term ID", "Adjusted P-Value (g:SCS)"),
  proximal_neuronbg_all59_e114 = c("SuppFile4(Distal_Proximal).xlsx", "Proximal GO (Neuron-Expressed)", "Term ID", "Adjusted P-Value (g:SCS)")
)
check <- do.call(rbind, lapply(names(published), function(n) {
  spec <- published[[n]]
  pub <- read_sheet(spec[1], spec[2])
  pub <- data.frame(term_id = pub[[spec[3]]], p = suppressWarnings(as.numeric(pub[[spec[4]]])))
  pub <- pub[!is.na(pub$term_id), ]
  got <- gpres[[n]]$terms
  m <- merge(pub, got, by = "term_id")
  # Neuron-background sheets round P-values to four significant digits.
  rel <- if (nrow(m)) max(abs(m$p - m$p_value) / pmax(m$p_value, .Machine$double.xmin)) else 0
  data.frame(analysis = n, version = gpres[[n]]$version, file = spec[1], sheet = spec[2],
             published_terms = nrow(pub), computed_terms = nrow(got), shared = nrow(m),
             term_sets_equal = setequal(pub$term_id, got$term_id), max_rel_p_diff = signif(rel, 3))
}))
write_table(check, "go_gprofiler_published_check.tsv")
print(check)

nb <- gpres$neurite_neuronbg$terms
pick <- function(t, id) if (id %in% t$term_id) t$p_value[t$term_id == id] else NA

# Preranked GSEA on shrunk log2 fold changes. nPermSimple is 2,000 for the micropattern
# contrast and 10,000 for neurite versus soma and distal versus proximal. Seed is 1.
gsea <- function(tab, nperm, seed) {
  ranks <- sort(setNames(tab$log2FoldChange, tab$gene), decreasing = TRUE)
  set.seed(seed)
  g <- gseGO(geneList = ranks, ont = "ALL", OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
             nPermSimple = nperm, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
             pAdjustMethod = "BH", verbose = FALSE)
  as.data.frame(g)
}
gs <- list(
  micropattern_all_tested = gsea(f1, 2000, 1),
  neurite_soma_de_only = gsea(s2, 10000, 1),
  distal_proximal_de_only = gsea(s4, 10000, 1)
)
gsea_terms <- do.call(rbind, lapply(names(gs), function(n) {
  if (!nrow(gs[[n]])) return(NULL)
  cbind(analysis = n, gs[[n]][, c("ONTOLOGY", "ID", "Description", "setSize", "NES", "pvalue", "p.adjust")])
}))
write_table(gsea_terms, "go_gsea_terms.tsv")

# Middleton-style background: >= 1 read in at least half of the 59 libraries.
middleton_bg <- rownames(counts)[rowSums(counts >= 1) >= ceiling(ncol(counts) / 2)]
to_entrez <- function(s) unique(suppressWarnings(bitr(s, "SYMBOL", "ENTREZID", org.Mm.eg.db))$ENTREZID)
universe <- to_entrez(middleton_bg)
# DOSE 4.2.0 counts query genes outside the universe in the hypergeometric draw size, so the query
# is restricted to the background, as GOrilla's two-list mode used by Middleton et al. requires.
enrich <- function(genes, label) {
  q <- intersect(to_entrez(genes), universe)
  do.call(rbind, lapply(c("BP", "CC", "MF"), function(ont) {
    e <- as.data.frame(enrichGO(q, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
                                universe = universe, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 2000, readable = TRUE))
    if (!nrow(e)) return(NULL)
    cbind(set = label, ontology = ont,
          e[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count")])
  }))
}
eg <- rbind(enrich(sets$neurite, "neurite"), enrich(sets$soma, "soma"))
write_table(eg, "go_enrichgo_middleton_terms.tsv")
eg_n <- eg[eg$set == "neurite", ]
eg_s <- eg[eg$set == "soma", ]

write_kv("go_enrichment", list(
  gprofiler_published_version = gpres$neurite_annotated$version,
  gprofiler_e114_version = gpres$proximal_neuronbg_all59_e114$version,
  published_tables_reproduced = sum(check$term_sets_equal),
  published_tables_checked = nrow(check),
  published_max_rel_p_diff = max(check$max_rel_p_diff),
  background_neurite_soma_prefilter = length(backgrounds$neurite_soma_prefilter),
  background_proximal_distal_prefilter = length(backgrounds$proximal_distal_prefilter),
  background_all59_ge3_in_10pct = length(backgrounds$all59_ge3_in_10pct),
  neuron_bg_neurite_terms = nrow(nb),
  neuron_bg_neurite_only_GO0022626 = identical(nb$term_id, "GO:0022626"),
  neuron_bg_GO0022626_padj = pick(nb, "GO:0022626"),
  neuron_bg_soma_terms = nrow(gpres$soma_neuronbg$terms),
  consistent_union_genes = length(sets$consistent_union),
  gsea_micropattern_genes = nrow(f1),
  gsea_micropattern_significant = sum(gs$micropattern_all_tested$p.adjust < 0.05),
  gsea_micropattern_min_padj = min(gs$micropattern_all_tested$p.adjust),
  gsea_neurite_soma_genes = nrow(s2),
  gsea_neurite_soma_significant = sum(gs$neurite_soma_de_only$p.adjust < 0.05),
  gsea_neurite_soma_min_padj = min(gs$neurite_soma_de_only$p.adjust),
  gsea_distal_proximal_genes = nrow(s4),
  gsea_distal_proximal_significant = sum(gs$distal_proximal_de_only$p.adjust < 0.05),
  gsea_distal_proximal_min_padj = min(gs$distal_proximal_de_only$p.adjust),
  enrichgo_background_genes = length(middleton_bg),
  enrichgo_background_entrez = length(universe),
  enrichgo_neurite_terms = nrow(eg_n),
  enrichgo_soma_terms = nrow(eg_s),
  enrichgo_neurite_GO0022626_padj = pick(data.frame(term_id = eg_n$ID, p_value = eg_n$p.adjust), "GO:0022626"),
  enrichgo_neurite_GO0002181_padj = pick(data.frame(term_id = eg_n$ID, p_value = eg_n$p.adjust), "GO:0002181"),
  enrichgo_neurite_respiratory_chain_terms = sum(grepl("respiratory chain", eg_n$Description)),
  enrichgo_neurite_themes_recovered = all(c("GO:0022626", "GO:0002181") %in% eg_n$ID) &&
    any(grepl("respiratory chain", eg_n$Description))
))
