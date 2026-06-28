#!/usr/bin/env Rscript
# FR-A3/A4: Build clustering items and pairwise Wasserstein-1 per feature (1-D only).
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

wasserstein_1d <- function(x, y) {
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("Install transport: install.packages('transport')")
  }
  x <- as.numeric(x[!is.na(x)])
  y <- as.numeric(y[!is.na(y)])
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  transport::wasserstein1d(x, y)
}

build_items <- function(scores, meta, cfg) {
  min_n <- cfg$features$min_cells_per_item
  items <- list()

  for (state in cfg$features$cell_states) {
    for (hist in cfg$features$histology_groups) {
      idx <- which(meta$cell_state == state & meta$histology == hist)
      if (length(idx) >= min_n) {
        items[[length(items) + 1]] <- list(
          item_id = paste(state, hist, sep = "__"),
          cell_state = state,
          histology = hist,
          cell_idx = idx
        )
      } else if (length(idx) > 0) {
        message(sprintf(
          "[flag] Dropped %s/%s: only %d cells (< %d)",
          state, hist, length(idx), min_n
        ))
      }
    }
  }
  items
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "wasserstein_matrix")

  scores_rds <- resolve_path(cfg, cfg$paths$phase_a$scores_rds)
  items_rds <- resolve_path(cfg, cfg$paths$phase_a$items_rds)
  w_dir <- resolve_path(cfg, cfg$paths$phase_a$wasserstein_dir)
  ensure_dir(w_dir)

  if (!file.exists(scores_rds)) stop("Run 03_compute_scores.R first")
  dat <- readRDS(scores_rds)
  items <- build_items(dat$scores, dat$meta, cfg)
  saveRDS(items, items_rds)

  feature_ids <- dat$feature_ids
  for (feat in feature_ids) {
    n <- length(items)
    D <- matrix(0, n, n, dimnames = list(
      vapply(items, `[[`, "", "item_id"),
      vapply(items, `[[`, "", "item_id")
    ))
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i == j) next
        xi <- dat$scores[items[[i]]$cell_idx, feat]
        yj <- dat$scores[items[[j]]$cell_idx, feat]
        D[i, j] <- wasserstein_1d(xi, yj)
      }
    }
    out_file <- file.path(w_dir, paste0(feat, "_w1_dist.rds"))
    saveRDS(list(feature = feat, distance = D, items = items), out_file)
    message("[ok] W1 matrix: ", out_file)
  }
  invisible(w_dir)
}

if (sys.nframe() == 0) main()
