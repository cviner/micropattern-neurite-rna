# GOSemSim::mgoSim (Wang, BMA) between the proximal-higher and distal-higher GO term sets of
# Supplementary File 4, per ontology. The term sets are the g:Profiler sheets
# "Distal GO (All Coding Genes)" and "Proximal GO (All Coding Genes)"; R/go_enrichment.R
# regenerates both sheets exactly on g:Profiler e112_eg59_p19.

source("R/common.R")
require_pins(c("GOSemSim", "org.Mm.eg.db", "GO.db", "readxl"))
suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Mm.eg.db)
})

sheet_terms <- function(sheet) {
  x <- read_sheet("SuppFile4(Distal_Proximal).xlsx", sheet)
  data.frame(ontology = x[["GO Source"]], term_id = x[["GO Term ID"]])
}
distal <- sheet_terms("Distal GO (All Coding Genes)")
proximal <- sheet_terms("Proximal GO (All Coding Genes)")
write_table(rbind(cbind(set = "distal", distal), cbind(set = "proximal", proximal)),
            "gosemsim_term_sets.tsv")

sem <- lapply(c(BP = "BP", MF = "MF", CC = "CC"),
              function(o) suppressWarnings(godata("org.Mm.eg.db", ont = o, computeIC = FALSE)))

wang_bma <- function(d, p, sd) {
  if (exists(".SemSimCache", envir = .GlobalEnv)) rm(".SemSimCache", envir = .GlobalEnv)
  mgoSim(d, p, semData = sd, measure = "Wang", combine = "BMA")
}

kv <- list()
for (ont in c("BP", "MF", "CC")) {
  d <- distal$term_id[distal$ontology == paste0("GO:", ont)]
  p <- proximal$term_id[proximal$ontology == paste0("GO:", ont)]
  key <- tolower(ont)
  shared <- length(intersect(d, p))
  kv[[paste0(key, "_distal_terms")]] <- length(d)
  kv[[paste0(key, "_proximal_terms")]] <- length(p)
  kv[[paste0(key, "_shared_terms")]] <- shared
  kv[[paste0(key, "_similarity")]] <- wang_bma(d, p, sem[[ont]])
  kv[[paste0(key, "_similarity_pct")]] <- round(100 * kv[[paste0(key, "_similarity")]], 1)
}
info <- GO.db::GO_dbInfo()
kv$go_db_source_date <- info$value[info$name == "GOSOURCEDATE"]
write_kv("gosemsim", kv)
