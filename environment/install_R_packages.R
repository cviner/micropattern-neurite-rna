# Installs the R package versions used for the manuscript analyses.
# Run from the repository root:
#   R --vanilla -f environment/install_R_packages.R

rv <- paste(R.version$major, R.version$minor, sep = ".")
if (rv != "4.5.3") stop("R 4.5.3 is required; this is R ", rv)

lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) stop("Set R_LIBS_USER to the target library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
}
# Bioconductor 3.21 is frozen, so its binaries carry exactly these versions.
bioc <- c(
  DESeq2 = "1.48.2",
  apeglm = "1.30.0",
  AnnotationDbi = "1.70.0",
  BiocParallel = "1.42.2",
  org.Mm.eg.db = "3.21.0",
  GO.db = "3.21.0",
  GOSemSim = "2.34.0",
  clusterProfiler = "4.16.0"
)
# Paper versions. CRAN serves only the current binary, so older versions build from the archive.
cran <- c(
  gprofiler2 = "0.2.4",
  dplyr = "1.2.1",
  ggplot2 = "4.0.2",
  readxl = "1.4.5",
  pheatmap = "1.0.13",
  ggrepel = "0.9.8",
  patchwork = "1.3.2"
)

installed_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE, lib.loc = lib)) return(NA_character_)
  as.character(packageVersion(pkg, lib.loc = lib))
}

options(repos = BiocManager::repositories(version = "3.21"))
need_bioc <- names(bioc)[vapply(names(bioc), function(p) !identical(installed_version(p), bioc[[p]]), logical(1))]
if (length(need_bioc)) {
  BiocManager::install(need_bioc, version = "3.21", lib = lib, update = FALSE, ask = FALSE)
}

# When R.framework is not under /Library/Frameworks, point LIBR at the framework actually running.
framework_parent <- normalizePath(file.path(R.home(), "..", "..", "..", ".."))
if (!file.exists(file.path("/Library/Frameworks/R.framework/Resources/etc/Makeconf"))) {
  makevars <- tempfile(fileext = ".mk")
  writeLines(sprintf("LIBR = -F%s -framework R", framework_parent), makevars)
  Sys.setenv(R_MAKEVARS_USER = makevars)
}

current <- available.packages(repos = "https://cloud.r-project.org")
for (pkg in names(cran)) {
  ver <- cran[[pkg]]
  if (identical(installed_version(pkg), ver)) next
  if (pkg %in% rownames(current) && current[pkg, "Version"] == ver) {
    install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
  } else {
    install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
    src <- sprintf("https://cloud.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", pkg, pkg, ver)
    install.packages(src, lib = lib, repos = NULL, type = "source")
  }
}

want <- c(bioc, cran)
got <- vapply(names(want), installed_version, character(1))
report <- data.frame(package = names(want), required = unname(want), installed = unname(got))
report$ok <- report$required == report$installed
print(report, row.names = FALSE)
if (!all(report$ok %in% TRUE)) stop("Some package versions do not match")
