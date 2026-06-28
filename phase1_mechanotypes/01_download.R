#!/usr/bin/env Rscript
# FR-A1: Pull SCPCP000006 via ScPCAr (real data; --demo optional).
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "scpca_api.R"))

load_primary_sce <- function(file_paths) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Install SingleCellExperiment: BiocManager::install('SingleCellExperiment')")
  }
  merged <- grep("_merged.*_processed\\.rds$", file_paths, value = TRUE)
  processed <- grep("_processed\\.rds$", file_paths, value = TRUE)
  processed <- setdiff(processed, merged)
  target <- if (length(merged)) merged[1] else processed[1]
  if (!length(target)) {
    stop(
      "No *_processed.rds found.\n",
      "Files returned:\n  ", paste(head(basename(file_paths), 10), collapse = "\n  ")
    )
  }
  message("[load] ", basename(target))
  readRDS(target)
}

main <- function() {
  cfg <- load_config()
  opts <- parse_args_demo()
  set_seed_logged(cfg$features$seed, "download")

  out_rds <- resolve_path(cfg, cfg$paths$phase_a$seurat_rds)
  manifest_rds <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "download_manifest.rds"))
  ensure_dir(dirname(out_rds))
  raw_dir <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "scpca_downloads"))
  ensure_dir(raw_dir)

  if (file.exists(out_rds) && !nzchar(opts$sample_id) && !opts$spatial) {
    message("[skip] Processed object exists: ", out_rds)
    message("       Delete it or use --sample / --spatial to fetch something else.")
    return(invisible(out_rds))
  }

  if (opts$demo) {
    message("[demo] Generating synthetic snRNA-seq stand-in")
    log_provenance(cfg, "DEMO mode — synthetic snRNA-seq, not ScPCA pull")
    demo <- make_demo_seurat(cfg)
    saveRDS(demo, out_rds)
    return(invisible(out_rds))
  }

  project_id <- cfg$paths$project_id
  fmt <- if (opts$spatial) "spatial" else "sce"
  merged <- opts$merged && !opts$spatial

  log_provenance(cfg, sprintf(
    "ScPCA pull project=%s sample=%s format=%s merged=%s",
    project_id,
    if (nzchar(opts$sample_id)) opts$sample_id else "ALL",
    fmt, merged
  ))

  email <- Sys.getenv("SCPCA_NOTIFY_EMAIL", unset = "")
  if (!nzchar(email)) email <- NULL

  file_paths <- download_scpca_data(
    project_id = project_id,
    destination = raw_dir,
    sample_id = opts$sample_id,
    format = fmt,
    merged = merged,
    email = email,
    timeout_minutes = as.numeric(Sys.getenv("SCPCA_TIMEOUT_MIN", "120"))
  )

  saveRDS(list(
    project_id = project_id,
    sample_id = opts$sample_id,
    format = fmt,
    merged = merged,
    file_paths = file_paths,
    downloaded_at = Sys.time(),
    scpcar_api = if (nzchar(opts$sample_id)) "create_dataset" else "download_project"
  ), manifest_rds)

  if (opts$spatial) {
    message("[ok] Spatial download complete -> ", raw_dir)
    message("     Manifest: ", manifest_rds)
    message("     H&E: look for spatial/tissue_hires_image.png under Space Ranger folders")
    return(invisible(manifest_rds))
  }

  sce <- load_primary_sce(file_paths)
  saveRDS(sce, out_rds)
  message("[ok] SCE saved -> ", out_rds)
  message("     colnames(colData(sce)) — look for consensus_celltype_annotation")
  invisible(out_rds)
}

if (sys.nframe() == 0) main()
