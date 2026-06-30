#!/usr/bin/env Rscript
# Install R dependencies for sc-wilms-data Phase A.
user_lib <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
if (is.na(user_lib) || !nzchar(user_lib)) {
  user_lib <- file.path(
    Sys.getenv("LOCALAPPDATA", unset = Sys.getenv("HOME")),
    "R", "win-library", paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1])
  )
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
message("[lib] ", user_lib)

pkgs_cran <- c(
  "remotes", "yaml", "dplyr", "tidyr", "ggplot2", "Matrix",
  "transport", "testthat",
  "logistf"   # Firth penalized logistic (prognostics, 16_prognostic_association.R)
)
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " ...")
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
# edgeR/limma: moderated pseudobulk DE; fgsea+msigdbr: Hallmark GSEA.
# (waddR is omitted — no current Bioconductor binary; the W1 location/size/shape
#  decomposition is computed directly in base R in 06_wasserstein_decompose.R.)
bioc <- c("SingleCellExperiment", "SummarizedExperiment", "ConsensusClusterPlus",
          "edgeR", "limma", "fgsea", "msigdbr")
for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " (Bioconductor) ...")
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}
if (!requireNamespace("ScPCAr", quietly = TRUE)) {
  message("Installing ScPCAr from GitHub (requires Rtools on Windows) ...")
  remotes::install_github("AlexsLemonade/ScPCAr", upgrade = "never")
}
if (!requireNamespace("ScPCAr", quietly = TRUE)) {
  stop(
    "ScPCAr failed to install. On Windows install Rtools 4.5:\n",
    "  winget install RProject.Rtools\n",
    "Then re-run this script. Metadata-only explore works without ScPCAr:\n",
    "  python scripts/fetch_scpca_metadata.py"
  )
}
message("[ok] R packages ready. Next: Rscript scripts/scpca_auth.R")
