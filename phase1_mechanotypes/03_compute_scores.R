#!/usr/bin/env Rscript
# FR-A2: Compute predefined 1-D feature scores per cell.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

score_feature <- function(counts, genes_pos, genes_neg, gene_lookup = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install Matrix: install.packages('Matrix')")
  }
  if (!is.null(gene_lookup)) {
    genes_pos <- resolve_feature_genes(genes_pos, gene_lookup)
    genes_neg <- resolve_feature_genes(genes_neg, gene_lookup)
  }
  n <- ncol(counts)
  pos <- intersect(genes_pos, rownames(counts))
  neg <- intersect(genes_neg, rownames(counts))
  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- 1
  med <- median(lib)
  pos_score <- if (length(pos)) {
    Matrix::colSums(counts[pos, , drop = FALSE]) / lib * med
  } else {
    rep(0, n)
  }
  neg_score <- if (length(neg)) {
    Matrix::colSums(counts[neg, , drop = FALSE]) / lib * med
  } else {
    rep(0, n)
  }
  as.numeric(log1p(pos_score) - log1p(neg_score))
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "compute_scores")

  in_rds <- resolve_path(cfg, cfg$paths$phase_a$seurat_rds)
  out_rds <- resolve_path(cfg, cfg$paths$phase_a$scores_rds)
  ensure_dir(dirname(out_rds))

  if (!file.exists(in_rds)) stop("Run 02_qc_normalize.R first")

  dat <- readRDS(in_rds)
  counts <- dat$counts
  meta <- dat$meta

  score_mat <- do.call(cbind, lapply(cfg$features$features, function(f) {
    score_feature(counts, f$genes_positive, f$genes_negative, dat$gene_lookup)
  }))
  colnames(score_mat) <- vapply(cfg$features$features, `[[`, "", "id")
  rownames(score_mat) <- colnames(counts)

  scores <- list(
    scores = score_mat,
    meta = meta,
    feature_ids = colnames(score_mat)
  )
  saveRDS(scores, out_rds)
  message("[ok] Feature scores -> ", out_rds)
  invisible(out_rds)
}

if (sys.nframe() == 0) main()
