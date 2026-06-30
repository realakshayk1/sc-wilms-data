#!/usr/bin/env Rscript
# FR-A6: location / size / shape decomposition of the 2-Wasserstein distance between
# favorable and anaplastic (and relapse vs not) score distributions, per program x
# compartment. This is the Schefzik et al. decomposition that the waddR package
# implements; waddR has no current Bioconductor binary (source-only, needs a compiler),
# so the closed-form decomposition is computed directly in base R:
#
#   d_W^2(X,Y) = (mu_X - mu_Y)^2  +  (sd_X - sd_Y)^2  +  2 sd_X sd_Y (1 - rho)
#               \___location___/    \____size_____/    \________shape________/
#
# where rho is the correlation between the quantile functions of X and Y. The
# decomposition is DESCRIPTIVE (it says *how* two distributions differ); significance
# remains the PATIENT-level permutation BH-FDR from 09_distributional_validation.R,
# which is joined in — cell-level decomposition is never used to claim significance.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

decompose_w2 <- function(x, y, ngrid = 2000) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 25 || length(y) < 25) return(NULL)
  p <- (seq_len(ngrid) - 0.5) / ngrid
  qx <- stats::quantile(x, p, names = FALSE); qy <- stats::quantile(y, p, names = FALSE)
  mux <- mean(x); muy <- mean(y); sx <- stats::sd(x); sy <- stats::sd(y)
  rho <- suppressWarnings(stats::cor(qx, qy)); if (!is.finite(rho)) rho <- 1
  location <- (mux - muy)^2
  size <- (sx - sy)^2
  shape <- 2 * sx * sy * (1 - rho)
  d2 <- location + size + shape
  if (!is.finite(d2) || d2 <= 0) return(NULL)
  data.frame(n_x = length(x), n_y = length(y),
             d_wass = sqrt(d2), location = location, size = size, shape = shape,
             perc_location = 100 * location / d2, perc_size = 100 * size / d2,
             perc_shape = 100 * shape / d2, mean_diff = mux - muy, sd_ratio = sx / sy,
             stringsAsFactors = FALSE)
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "wasserstein_decompose")
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir); ensure_dir(out_dir)
  dat <- readRDS(resolve_path(cfg, cfg$paths$phase_a$scores_rds))
  scores <- dat$scores; meta <- dat$meta
  programs <- colnames(scores)
  comps <- intersect(c("blastemal", "epithelial", "stromal"), unique(meta$cell_state))

  # group vectors per contrast
  hist <- tolower(as.character(meta$histology))
  g_hist <- ifelse(hist == "anaplastic", "pos", ifelse(hist == "favorable", "neg", NA))
  clin <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_clinical_metadata.tsv"))
  g_rel <- rep(NA_character_, nrow(meta))
  if (file.exists(clin)) {
    md <- read.delim(clin, check.names = FALSE)
    u <- unique(md[, c("scpca_sample_id", "relapse_status")])
    rmap <- setNames(as.character(u$relapse_status), u$scpca_sample_id)
    rs <- rmap[as.character(meta$sample_id)]
    g_rel <- ifelse(rs == "Yes", "pos", ifelse(rs == "No", "neg", NA))
  }
  contrasts <- list(histology = g_hist, relapse = g_rel)

  rows <- list()
  for (cn in names(contrasts)) {
    grp <- contrasts[[cn]]
    for (feat in programs) for (comp in comps) {
      sel <- meta$cell_state == comp & !is.na(grp)
      x <- scores[sel & grp == "neg", feat]   # favorable / no-relapse
      y <- scores[sel & grp == "pos", feat]   # anaplastic / relapse
      d <- tryCatch(decompose_w2(x, y), error = function(e) NULL)
      if (is.null(d)) next
      d$contrast <- cn; d$feature <- feat; d$cell_state <- comp
      d$dominant_component <- c("location", "size", "shape")[which.max(c(d$perc_location, d$perc_size, d$perc_shape))]
      rows[[length(rows) + 1]] <- d
    }
  }
  if (!length(rows)) { message("[warn] no decompositions computed"); return(invisible(NULL)) }
  res <- do.call(rbind, rows)

  # join patient-level significance from the distributional validation (if present)
  dv_files <- c(resolve_path(cfg, "results/mechanotypes/distributional_validation.csv"),
                resolve_path(cfg, "results/mechanotypes/distributional_validation_relapse.csv"))
  dv_files <- dv_files[file.exists(dv_files)]
  if (length(dv_files)) {
    dv <- do.call(rbind, lapply(dv_files, function(f) {
      d <- read.csv(f); d$contrast <- tolower(as.character(d$contrast)); d
    }))
    key <- function(d) paste(tolower(d$contrast), d$feature, d$cell_state, sep = "|")
    res$p_perm_BH <- dv$p_perm_BH[match(key(res), key(dv))]
    res$significant_BH <- !is.na(res$p_perm_BH) & res$p_perm_BH < 0.05
  }

  cols <- c("contrast", "feature", "cell_state", "n_x", "n_y", "d_wass",
            "location", "size", "shape", "perc_location", "perc_size", "perc_shape",
            "dominant_component", "mean_diff", "sd_ratio",
            intersect(c("p_perm_BH", "significant_BH"), colnames(res)))
  res <- res[, cols]
  out_csv <- file.path(out_dir, "wasserstein_decomposition.csv")
  write.csv(res, out_csv, row.names = FALSE)
  message("[ok] Wasserstein decomposition -> ", out_csv)
  message(sprintf("[ok] %d (program x compartment x axis) decompositions; dominant: %s",
                  nrow(res), paste(names(sort(table(res$dominant_component), decreasing = TRUE)),
                                   collapse = " > ")))
  if ("significant_BH" %in% colnames(res))
    message(sprintf("[ok] %d/%d patient-level BH-FDR significant",
                    sum(res$significant_BH, na.rm = TRUE), nrow(res)))
}

if (sys.nframe() == 0) main()
