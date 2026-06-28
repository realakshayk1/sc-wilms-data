#!/usr/bin/env Rscript
# Locate processed SCE from manual unzip and save to data/processed/seurat_wilms.rds
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(dirname(.script_dir), "phase1_mechanotypes", "utils.R"))

find_processed_rds <- function(root) {
  files <- list.files(root, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  merged_processed <- files[grepl("_merged.*_processed\\.rds$", files)]
  if (length(merged_processed)) return(merged_processed[1])
  merged <- files[grepl("_merged\\.rds$", files)]
  if (length(merged)) return(merged[1])
  processed <- files[grepl("_processed\\.rds$", files)]
  if (length(processed)) return(processed[1])
  NULL
}

main <- function() {
  cfg <- load_config()
  extract_root <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "scpca_downloads"))
  sce_dir <- file.path(extract_root, "single-cell-experiment")
  out_rds <- resolve_path(cfg, cfg$paths$phase_a$merged_sce_rds)
  ensure_dir(dirname(out_rds))

  if (!dir.exists(sce_dir)) {
    stop("Missing ", sce_dir, " — run scripts/ingest_manual_downloads.ps1 first")
  }

  target <- find_processed_rds(sce_dir)
  if (is.null(target)) {
    stop("No *_processed.rds under ", sce_dir)
  }

  message("[load] ", target)
  sce <- readRDS(target)

  # Attach sample-level histology (subdiagnosis) for merged objects
  sample_csv <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_samples.csv"))
  if (file.exists(sample_csv) && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    samples <- read.csv(sample_csv, stringsAsFactors = FALSE)
    meta <- as.data.frame(SummarizedExperiment::colData(sce))
    sid_col <- intersect(
      c("scpca_sample_id", "sample_id", "scpca_id"),
      colnames(meta)
    )[1]
    if (!is.na(sid_col) && "subdiagnosis" %in% colnames(samples)) {
      id_col <- intersect(c("scpca_id", "scpca_sample_id"), colnames(samples))[1]
      hist_map <- setNames(
        ifelse(grepl("anaplas", samples$subdiagnosis, ignore.case = TRUE), "anaplastic",
          ifelse(grepl("favor", samples$subdiagnosis, ignore.case = TRUE), "favorable", NA_character_)
        ),
        samples[[id_col]]
      )
      meta$histology <- unname(hist_map[as.character(meta[[sid_col]])])
      meta$subdiagnosis <- samples$subdiagnosis[match(meta[[sid_col]], samples[[id_col]])]
      SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(meta)
      message("[ok] Joined histology from ", basename(sample_csv))
    }
  }

  saveRDS(sce, out_rds)
  message("[ok] SCE -> ", out_rds)
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    message("     cells: ", ncol(sce), " genes: ", nrow(sce))
    cn <- colnames(SummarizedExperiment::colData(sce))
    message("     colData cols: ", paste(head(cn, 12), collapse = ", "), if (length(cn) > 12) " ..." else "")
  }
  invisible(out_rds)
}

if (sys.nframe() == 0) main()
