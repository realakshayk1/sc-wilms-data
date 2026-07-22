#!/usr/bin/env Rscript
# WS1 / P1: per-CELL lever + extrinsic-axis scores over config/levers.yaml.
#
# Feeds coupling_core.R (Level A, cell-level partial-correlation network). Uses the SAME signed
# library-normalized score as 03_compute_scores.R (raw signed pos-neg; standardization happens
# inside coupling_core per the §2 one-convention rule). Curation-specific behaviour:
#   - emt_axis: scored ONLY in scope compartments (blastemal, epithelial); NA elsewhere
#     (de-confounds VIM = stroma vs mesenchymal-transitioned tumor cell).
#   - crowding_sensitivity: at CELL level uses the `hippo_core` subset (PIEZO1/TRPV4/AMOTL2 are
#     snRNA-dropout); the full mechanosensor panel is for BULK only.
# Output: data/processed/lever_scores_percell.rds  {scores[cells x programs], meta, provenance}.
#
# Usage:  Rscript phase1_mechanotypes/18_program_scores_percell.R
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages(library(Matrix))

# signed library-normalized score (mirrors 03_compute_scores.R::score_feature)
score_feature <- function(counts, genes_pos, genes_neg, gene_lookup) {
  pos <- resolve_feature_genes(genes_pos, gene_lookup)
  neg <- resolve_feature_genes(genes_neg, gene_lookup)
  pos <- intersect(pos, rownames(counts)); neg <- intersect(neg, rownames(counts))
  n <- ncol(counts)
  lib <- Matrix::colSums(counts); lib[lib == 0] <- 1; med <- median(lib)
  pos_s <- if (length(pos)) Matrix::colSums(counts[pos, , drop = FALSE]) / lib * med else rep(0, n)
  neg_s <- if (length(neg)) Matrix::colSums(counts[neg, , drop = FALSE]) / lib * med else rep(0, n)
  list(score = as.numeric(log1p(pos_s) - log1p(neg_s)),
       n_pos = length(pos), n_neg = length(neg))
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "program_scores_percell")
  lv <- yaml::read_yaml(file.path(cfg$root, "config", "levers.yaml"))

  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts; meta <- proc$meta; gl <- proc$gene_lookup
  cell_state <- as.character(meta$cell_state)
  message(sprintf("[data] %d genes x %d cells; %d samples",
                  nrow(counts), ncol(counts), length(unique(meta$sample_id))))

  programs <- c(lv$levers, lv$extrinsic_axes)
  ids <- vapply(programs, `[[`, "", "id")
  score_mat <- matrix(NA_real_, nrow = ncol(counts), ncol = length(programs),
                      dimnames = list(colnames(counts), ids))

  for (p in programs) {
    pos <- p$genes_positive
    # CELL-level crowding uses hippo_core (snRNA-detectable subset)
    if (identical(p$id, "crowding_sensitivity") && !is.null(p$hippo_core)) pos <- p$hippo_core
    r <- score_feature(counts, pos, p$genes_negative, gl)
    s <- r$score
    # tumor-compartment scoping (emt_axis): NA outside scope
    if (!is.null(p$scope)) s[!(cell_state %in% p$scope)] <- NA_real_
    score_mat[, p$id] <- s
    message(sprintf("[score] %-22s pos=%d neg=%d %s%s", p$id, r$n_pos, r$n_neg,
                    if (!is.null(p$scope)) sprintf("scope=%s ", paste(p$scope, collapse="/")) else "",
                    if (identical(p$id,"crowding_sensitivity")) "(hippo_core)" else ""))
  }

  keep_meta <- intersect(c("cell_id", "sample_id", "cell_state", "subdiagnosis", "histology",
                           "disease_timing", "subsets_mito_percent", "sum", "detected"),
                         colnames(meta))
  out <- list(
    scores = score_mat,
    meta = meta[, keep_meta, drop = FALSE],
    program_ids = ids,
    roles = setNames(c(rep("lever", length(lv$levers)), rep("extrinsic_axis", length(lv$extrinsic_axes))), ids),
    provenance = list(levers_yaml = "config/levers.yaml", seed = cfg$features$seed,
                      convention = "raw signed pos-neg, library-normalized; standardize in coupling_core")
  )
  out_rds <- resolve_path(cfg, "data/processed/lever_scores_percell.rds")
  ensure_dir(dirname(out_rds)); saveRDS(out, out_rds)
  message("[ok] per-cell lever scores -> ", out_rds)
  # sanity: non-NA counts per program
  nn <- colSums(!is.na(score_mat))
  print(data.frame(program = ids, n_scored = nn, row.names = NULL))
}

if (sys.nframe() == 0) main()
