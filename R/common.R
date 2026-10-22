# Shared inputs for the R analyses. Every script sources this file from the repository root.

COUNTS_TSV <- file.path("data", "Patterned_Neuron_rmdup_exon_counts.tsv")
PUBLISHED <- file.path("data", "published")
EXTERNAL <- file.path("data", "external")
RESULTS <- "results"
FIGURES <- "figures"

PINS <- c(
  DESeq2 = "1.48.2",
  apeglm = "1.30.0",
  gprofiler2 = "0.2.4",
  clusterProfiler = "4.16.0",
  GOSemSim = "2.34.0",
  org.Mm.eg.db = "3.21.0",
  AnnotationDbi = "1.70.0",
  GO.db = "3.21.0",
  readxl = "1.4.5",
  ggplot2 = "4.0.2",
  pheatmap = "1.0.13"
)

require_pins <- function(pkgs) {
  rv <- paste(R.version$major, R.version$minor, sep = ".")
  if (rv != "4.5.3") stop("R 4.5.3 is required; this is R ", rv)
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Package ", p, " is not installed")
    got <- as.character(packageVersion(p))
    if (got != PINS[[p]]) stop(p, " ", got, " is installed; ", PINS[[p]], " is required")
  }
  invisible(TRUE)
}

if (!file.exists(COUNTS_TSV)) stop("Run scripts from the repository root")

# Fall back to the system zoneinfo when R does not ship its own.
if (!nzchar(Sys.getenv("TZDIR")) && !dir.exists(file.path(R.home("share"), "zoneinfo")) &&
    dir.exists("/usr/share/zoneinfo")) {
  Sys.setenv(TZDIR = "/usr/share/zoneinfo")
}
dir.create(RESULTS, showWarnings = FALSE)
dir.create(FIGURES, showWarnings = FALSE)

load_counts <- function(ercc = FALSE) {
  df <- read.delim(COUNTS_TSV, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(df[, -1])
  rownames(m) <- df[[1]]
  if (!ercc) m <- m[!startsWith(rownames(m), "spikeIn.ERCC"), , drop = FALSE]
  storage.mode(m) <- "integer"
  m
}

culture_of <- function(libs) as.integer(sub("_.*$", "", libs))
compartment_of <- function(libs) sub("^[0-9]+_", "", libs)
batch_of <- function(libs) ifelse(culture_of(libs) <= 8L, "Batch1", "Batch2")

libraries <- function(m, comp) {
  libs <- colnames(m)[compartment_of(colnames(m)) == comp]
  libs[order(culture_of(libs))]
}

cpm <- function(m) sweep(m, 2, colSums(m), "/") * 1e6

fmt <- function(x, digits = 6) {
  if (length(x) != 1) stop("write_kv values must be scalars")
  if (is.na(x)) return("NA")
  if (is.logical(x)) return(as.character(as.integer(x)))
  if (is.numeric(x) && all(x == round(x))) return(format(x, scientific = FALSE, trim = TRUE))
  if (is.numeric(x)) return(format(signif(x, digits), scientific = abs(x) < 1e-4 & x != 0, trim = TRUE))
  as.character(x)
}

write_kv <- function(name, kv) {
  path <- file.path(RESULTS, paste0(name, ".tsv"))
  vals <- vapply(kv, fmt, character(1))
  writeLines(c("key\tvalue", paste(names(kv), vals, sep = "\t")), path)
  message("wrote ", path)
  invisible(path)
}

write_table <- function(df, file) {
  path <- file.path(RESULTS, file)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("wrote ", path)
  invisible(path)
}

read_sheet <- function(file, sheet) {
  as.data.frame(readxl::read_excel(file.path(PUBLISHED, file), sheet = sheet))
}
