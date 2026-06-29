#!/usr/bin/env Rscript
# Label-mapping sensitivity: rerun switch detection under predefined mapping variants.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "03_compute_scores.R"), local = TRUE)
source(file.path(.script_dir, "04_wasserstein_matrix.R"), local = TRUE)

map_cell_state_variant <- function(raw, variant_cfg) {
  mapped <- rep(NA_character_, length(raw))
  for (state in names(variant_cfg)) {
    patterns <- variant_cfg[[state]]
  for (pat in patterns) {
      mapped[grepl(pat, raw, ignore.case = TRUE, fixed = TRUE)] <- state
    }
  }
  mapped
}

map_histology <- function(meta) {
  raw <- as.character(meta$subdiagnosis %||% meta$histology)
  mapped <- rep(NA_character_, length(raw))
  mapped[grepl("anaplas", raw, ignore.case = TRUE)] <- "anaplastic"
  mapped[grepl("favor", raw, ignore.case = TRUE)] <- "favorable"
  mapped
}

detect_switches <- function(clusters_by_feature, cell_states) {
  rows <- list()
  for (feat in names(clusters_by_feature)) {
    clusters <- clusters_by_feature[[feat]]
    for (state in cell_states) {
      fav_id <- paste(state, "favorable", sep = "__")
      ana_id <- paste(state, "anaplastic", sep = "__")
      if (!(fav_id %in% names(clusters) && ana_id %in% names(clusters))) next
      rows[[length(rows) + 1]] <- paste(feat, state, sep = "|")
    }
  }
  rows
}

run_variant <- function(counts, meta, gene_lookup, cfg, variant_name, variant_cfg) {
  if (!"cellassign_celltype_annotation" %in% colnames(meta)) {
    stop("cellassign_celltype_annotation required for label sensitivity")
  }
  raw <- as.character(meta$cellassign_celltype_annotation)
  meta$cell_state <- map_cell_state_variant(raw, variant_cfg)
  meta$histology <- map_histology(meta)
  use <- !is.na(meta$cell_state) & !is.na(meta$histology)
  counts <- counts[, use, drop = FALSE]
  meta <- meta[use, , drop = FALSE]

  score_mat <- do.call(cbind, lapply(cfg$features$features, function(f) {
    score_feature(counts, f$genes_positive, f$genes_negative, gene_lookup)
  }))
  colnames(score_mat) <- vapply(cfg$features$features, `[[`, "", "id")

  clusters_by_feature <- list()
  for (feat in colnames(score_mat)) {
    items <- build_items(score_mat, meta, cfg)
    if (length(items) < 3) next
    n <- length(items)
    D <- matrix(0, n, n, dimnames = list(
      vapply(items, `[[`, "", "item_id"),
      vapply(items, `[[`, "", "item_id")
    ))
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i == j) next
        xi <- score_mat[items[[i]]$cell_idx, feat]
        yj <- score_mat[items[[j]]$cell_idx, feat]
        D[i, j] <- transport::wasserstein1d(xi, yj)
      }
    }
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Install cluster package")
    }
    k <- min(3L, n - 1L)
    cl <- cluster::pam(as.dist(D), k = k)$clustering
    names(cl) <- rownames(D)
    clusters_by_feature[[feat]] <- as.list(cl)
  }

  switch_keys <- list()
  for (feat in names(clusters_by_feature)) {
    clusters <- clusters_by_feature[[feat]]
    for (state in cfg$features$cell_states) {
      fav_id <- paste(state, "favorable", sep = "__")
      ana_id <- paste(state, "anaplastic", sep = "__")
      if (!(fav_id %in% names(clusters) && ana_id %in% names(clusters))) next
      if (isTRUE(clusters[[fav_id]] != clusters[[ana_id]])) {
        switch_keys[[length(switch_keys) + 1]] <- paste(feat, state, sep = "|")
      }
    }
  }
  list(
    variant = variant_name,
    n_cells = nrow(meta),
    n_switches = length(switch_keys),
    switch_keys = switch_keys
  )
}

main <- function() {
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("Install transport package")
  }
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "label_sensitivity")

  in_rds <- resolve_path(cfg, cfg$paths$phase_a$merged_sce_rds)
  out_csv <- resolve_path(cfg, cfg$paths$phase_a$label_sensitivity_csv)
  mapping_path <- file.path(cfg$root, "config", "label_mapping.yaml")
  ensure_dir(dirname(out_csv))

  if (!file.exists(in_rds)) stop("Run ingest_manual_scpca.R first")
  variants <- yaml::read_yaml(mapping_path)$variants

  raw <- readRDS(in_rds)
  if (inherits(raw, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(raw, "counts")
    if (is.null(counts)) counts <- SummarizedExperiment::assay(raw, 1)
    meta <- as.data.frame(SummarizedExperiment::colData(raw))
    gene_lookup <- NULL
    if ("gene_symbol" %in% colnames(SummarizedExperiment::rowData(raw))) {
      gene_lookup <- build_gene_lookup(
        SummarizedExperiment::rowData(raw)$gene_symbol,
        rownames(counts)
      )
    }
  } else {
    stop("Expected SingleCellExperiment from merged_sce.rds")
  }

  # QC: genes per cell
  if (requireNamespace("Matrix", quietly = TRUE) && inherits(counts, "Matrix")) {
    keep <- Matrix::colSums(counts > 0) >= 200
  } else {
    keep <- colSums(counts > 0) >= 200
  }
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  results <- lapply(names(variants), function(vname) {
    vcfg <- variants[[vname]]
    pat_cfg <- vcfg[names(vcfg) %in% cfg$features$cell_states]
    run_variant(counts, meta, gene_lookup, cfg, vname, pat_cfg)
  })

  all_keys <- unique(unlist(lapply(results, `[[`, "switch_keys")))
  rows <- list()
  for (key in all_keys) {
    parts <- strsplit(key, "\\|", fixed = FALSE)[[1]]
    feat <- parts[1]
    state <- parts[2]
    present <- vapply(results, function(r) key %in% r$switch_keys, logical(1))
    rows[[length(rows) + 1]] <- data.frame(
      feature = feat,
      cell_state = state,
      n_variants_switch = sum(present),
      n_variants_total = length(results),
      jaccard_vs_current = mean(present) / length(results),
      robust = sum(present) >= ceiling(length(results) * 0.67),
      stringsAsFactors = FALSE
    )
  }
  summary_df <- if (length(rows)) do.call(rbind, rows) else data.frame()
  variant_summary <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      variant = r$variant,
      n_cells = r$n_cells,
      n_switches = r$n_switches,
      stringsAsFactors = FALSE
    )
  }))
  out <- list(summary = summary_df, by_variant = variant_summary)
  yaml::write_yaml(out, sub("\\.csv$", ".yaml", out_csv))
  if (nrow(summary_df)) write.csv(summary_df, out_csv, row.names = FALSE)
  write.csv(variant_summary, sub("\\.csv$", "_by_variant.csv", out_csv), row.names = FALSE)
  message("[ok] Label sensitivity -> ", out_csv)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (sys.nframe() == 0) main()
