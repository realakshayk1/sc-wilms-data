#!/usr/bin/env Rscript
# FR-A6 (P1): waddR decomposition of 2-Wasserstein into location/shape/size.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
`%||%` <- function(a, b) if (!is.null(a)) a else b

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "waddR_decompose")

  if (!requireNamespace("waddR", quietly = TRUE)) {
    stop("Install waddR: BiocManager::install('waddR')")
  }

  w_dir <- resolve_path(cfg, cfg$paths$phase_a$wasserstein_dir)
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir)
  ensure_dir(out_dir)

  scores_rds <- resolve_path(cfg, cfg$paths$phase_a$scores_rds)
  items_rds <- resolve_path(cfg, cfg$paths$phase_a$items_rds)
  dat <- readRDS(scores_rds)
  items <- readRDS(items_rds)

  files <- list.files(w_dir, pattern = "_w1_dist\\.rds$", full.names = TRUE)
  rows <- list()

  for (f in files) {
    obj <- readRDS(f)
    feat <- obj$feature
    n <- length(items)
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        xi <- dat$scores[items[[i]]$cell_idx, feat]
        yj <- dat$scores[items[[j]]$cell_idx, feat]
        decomp <- tryCatch(
          waddR::wasserstein(xi, yj),
          error = function(e) NULL
        )
        if (is.null(decomp)) next
        rows[[length(rows) + 1]] <- data.frame(
          feature = feat,
          item_a = items[[i]]$item_id,
          item_b = items[[j]]$item_id,
          w2 = decomp$distance %||% NA_real_,
          location = decomp$location %||% NA_real_,
          size = decomp$size %||% NA_real_,
          shape = decomp$shape %||% NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows)) {
    df <- do.call(rbind, rows)
    out_csv <- file.path(out_dir, "waddR_decomposition.csv")
    write.csv(df, out_csv, row.names = FALSE)
    message("[ok] waddR decomposition -> ", out_csv)
  } else {
    message("[warn] No waddR decompositions computed")
  }
}

if (sys.nframe() == 0) main()
