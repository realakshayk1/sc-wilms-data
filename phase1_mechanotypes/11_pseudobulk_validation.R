#!/usr/bin/env Rscript
# Pseudobulk (snRNA sample-level) and bulk RNA validation of Phase A program scores.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "03_compute_scores.R"), local = TRUE)

aggregate_sample_scores <- function(scores, meta, feature_ids, cell_states) {
  if (!"sample_id" %in% colnames(meta)) stop("sample_id required in scores metadata")
  samples <- unique(as.character(meta$sample_id))
  rows <- list()
  for (sid in samples) {
    idx <- which(meta$sample_id == sid)
    sub <- meta[idx, , drop = FALSE]
    hist <- unique(sub$histology[!is.na(sub$histology)])
    if (length(hist) != 1) next
    row <- data.frame(
      sample_id = sid,
      histology = hist[1],
      n_cells = length(idx),
      stringsAsFactors = FALSE
    )
    for (feat in feature_ids) {
      row[[feat]] <- mean(scores[idx, feat], na.rm = TRUE)
    }
    for (st in cell_states) {
      row[[paste0("frac_", st)]] <- mean(sub$cell_state == st, na.rm = TRUE)
    }
    rows[[length(rows) + 1]] <- row
  }
  do.call(rbind, rows)
}

score_bulk_matrix <- function(counts, features, gene_lookup) {
  n <- ncol(counts)
  mat <- matrix(NA_real_, nrow = n, ncol = length(features))
  colnames(mat) <- vapply(features, `[[`, "", "id")
  lib <- colSums(counts)
  lib[lib == 0] <- 1
  med <- median(lib)
  norm <- t(t(counts) / lib * med)

  for (j in seq_along(features)) {
    f <- features[[j]]
    pos <- resolve_feature_genes(f$genes_positive, gene_lookup)
    neg <- resolve_feature_genes(f$genes_negative, gene_lookup)
    pos <- intersect(pos, rownames(counts))
    neg <- intersect(neg, rownames(counts))
    pos_score <- if (length(pos)) colSums(norm[pos, , drop = FALSE]) else rep(0, n)
    neg_score <- if (length(neg)) colSums(norm[neg, , drop = FALSE]) else rep(0, n)
    mat[, j] <- log1p(pos_score) - log1p(neg_score)
  }
  mat
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "pseudobulk_validation")

  scores_rds <- resolve_path(cfg, cfg$paths$phase_a$scores_rds)
  pseudo_csv <- resolve_path(cfg, cfg$paths$phase_a$pseudobulk_csv)
  bulk_csv <- resolve_path(cfg, cfg$paths$phase_a$bulk_validation_csv)
  bulk_quant <- resolve_path(cfg, cfg$paths$phase_a$bulk_quant_tsv)
  bulk_meta <- resolve_path(cfg, cfg$paths$phase_a$bulk_metadata_tsv)
  samples_csv <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_samples.csv"))
  ensure_dir(dirname(pseudo_csv))

  dat <- readRDS(scores_rds)
  pseudo <- aggregate_sample_scores(dat$scores, dat$meta, dat$feature_ids, cfg$features$cell_states)
  write.csv(pseudo, pseudo_csv, row.names = FALSE)
  message("[ok] Sample pseudobulk -> ", pseudo_csv)

  # Wilcoxon: favorable vs anaplastic per feature at sample level
  test_rows <- list()
  for (feat in dat$feature_ids) {
    fav <- pseudo[[feat]][pseudo$histology == "favorable"]
    ana <- pseudo[[feat]][pseudo$histology == "anaplastic"]
    if (length(fav) < 3 || length(ana) < 3) next
    wt <- tryCatch(stats::wilcox.test(fav, ana)$p.value, error = function(e) NA_real_)
    test_rows[[length(test_rows) + 1]] <- data.frame(
      feature = feat,
      mean_favorable = mean(fav),
      mean_anaplastic = mean(ana),
      wilcox_p = wt,
      stringsAsFactors = FALSE
    )
  }
  if (length(test_rows)) {
    tests <- do.call(rbind, test_rows)
    write.csv(tests, sub("\\.csv$", "_histology_tests.csv", pseudo_csv), row.names = FALSE)
  }

  # Bulk RNA validation
  if (!file.exists(bulk_quant) || !file.exists(bulk_meta)) {
    message("[warn] Bulk quant/metadata missing — skip bulk validation")
    return(invisible(pseudo))
  }

  quant <- read.delim(bulk_quant, check.names = FALSE)
  gene_ids <- quant[[1]]
  lib_ids <- colnames(quant)[-1]
  counts <- as.matrix(quant[, -1, drop = FALSE])
  rownames(counts) <- gene_ids

  bmeta <- read.delim(bulk_meta, stringsAsFactors = FALSE)
  sample_hist <- NULL
  if (file.exists(samples_csv)) {
    samples <- read.csv(samples_csv, stringsAsFactors = FALSE)
    id_col <- intersect(c("scpca_id", "scpca_sample_id"), colnames(samples))[1]
    sample_hist <- setNames(
      ifelse(grepl("anaplas", samples$subdiagnosis, ignore.case = TRUE), "anaplastic",
        ifelse(grepl("favor", samples$subdiagnosis, ignore.case = TRUE), "favorable", NA_character_)
      ),
      samples[[id_col]]
    )
  }
  bmeta$histology <- unname(sample_hist[bmeta$sample_id])

  gene_lookup <- NULL
  seurat_rds <- resolve_path(cfg, cfg$paths$phase_a$seurat_rds)
  if (file.exists(seurat_rds)) {
    proc <- readRDS(seurat_rds)
    gene_lookup <- proc$gene_lookup
  }

  bulk_scores <- score_bulk_matrix(counts, cfg$features$features, gene_lookup)
  rownames(bulk_scores) <- lib_ids
  bulk_df <- data.frame(
    library_id = lib_ids,
    sample_id = bmeta$sample_id[match(lib_ids, bmeta$library_id)],
    histology = bmeta$histology[match(lib_ids, bmeta$library_id)],
    bulk_scores,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.csv(bulk_df, bulk_csv, row.names = FALSE)

  bulk_tests <- list()
  for (feat in dat$feature_ids) {
    if (!feat %in% colnames(bulk_df)) next
    fav <- bulk_df[[feat]][bulk_df$histology == "favorable"]
    ana <- bulk_df[[feat]][bulk_df$histology == "anaplastic"]
    fav <- fav[!is.na(fav)]
    ana <- ana[!is.na(ana)]
    if (length(fav) < 3 || length(ana) < 3) next
    wt <- tryCatch(stats::wilcox.test(fav, ana)$p.value, error = function(e) NA_real_)
    bulk_tests[[length(bulk_tests) + 1]] <- data.frame(
      feature = feat,
      mean_favorable = mean(fav),
      mean_anaplastic = mean(ana),
      wilcox_p = wt,
      n_favorable = length(fav),
      n_anaplastic = length(ana),
      stringsAsFactors = FALSE
    )
  }
  if (length(bulk_tests)) {
    write.csv(do.call(rbind, bulk_tests), sub("\\.csv$", "_histology_tests.csv", bulk_csv), row.names = FALSE)
  }
  message("[ok] Bulk validation -> ", bulk_csv)
  invisible(pseudo)
}

if (sys.nframe() == 0) main()
