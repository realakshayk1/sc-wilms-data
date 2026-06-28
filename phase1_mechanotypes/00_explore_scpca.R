#!/usr/bin/env Rscript
# Explore SCPCP000006 metadata (no download, no auth required for sample table).
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "scpca_api.R"))

main <- function() {
  cfg <- load_config()
  project_id <- cfg$paths$project_id
  ensure_scpcar()

  message("=== Project: ", project_id, " (Wilms) ===")
  proj <- ScPCAr::get_project_info(project_id)
  message("Title: ", proj$title)
  message("Modalities: ", paste(proj$modalities, collapse = ", "))
  message("Has spatial: ", proj$has_spatial_data)
  message("Has bulk RNA: ", proj$has_bulk_rna_seq)
  message("Downloadable samples: ", proj$downloadable_sample_count)

  if (!is.null(proj$diagnoses_counts)) {
    message("\n--- diagnosis counts ---")
    dc <- data.frame(count = unlist(proj$diagnoses_counts))
    dc$diagnosis <- rownames(dc)
    print(dc[order(-dc$count), c("diagnosis", "count"), drop = FALSE])
  }

  samples <- ScPCAr::get_project_samples(project_id)
  slices <- filter_wilms_samples(samples)

  out_dir <- resolve_path(cfg, cfg$paths$dirs$raw)
  ensure_dir(out_dir)
  sample_csv <- file.path(out_dir, paste0(project_id, "_samples.csv"))
  write.csv(samples, sample_csv, row.names = FALSE)
  message("\n[ok] Full sample table -> ", sample_csv)

  if (!is.null(slices$nucleus) && nrow(slices$nucleus)) {
    f <- file.path(out_dir, paste0(project_id, "_nucleus_samples.csv"))
    write.csv(slices$nucleus, f, row.names = FALSE)
    message("[ok] Nucleus (snRNA) samples (n=", nrow(slices$nucleus), ") -> ", f)
  }
  if (!is.null(slices$spot) && nrow(slices$spot)) {
    f <- file.path(out_dir, paste0(project_id, "_visium_samples.csv"))
    write.csv(slices$spot, f, row.names = FALSE)
    message("[ok] Visium spot samples (n=", nrow(slices$spot), ") -> ", f)
  }

  if ("diagnosis" %in% colnames(samples)) {
    message("\n--- diagnosis (sample-level) ---")
    print(sort(table(samples$diagnosis), decreasing = TRUE))
  }

  message("\n--- Next steps (requires get_auth) ---")
  message("  # merged snRNA exploration (recommended first download):")
  message("  Rscript phase1_mechanotypes/01_download.R")
  message("  # one nucleus sample (smaller; uses create_dataset API):")
  message("  Rscript phase1_mechanotypes/01_download.R --sample SCPCS000XXX")
  message("  # spatial / Visium + H&E:")
  message("  Rscript phase1_mechanotypes/01_download.R --spatial")
  invisible(samples)
}

if (sys.nframe() == 0) main()
