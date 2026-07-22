#!/usr/bin/env Rscript
# WS2 / P2: recover the coupling network from BULK RNA-seq (same 40 SCPCP tumors), then compare to
# the single-cell result — the resolution-ladder / "what does spatial buy you" test (PLAN §4).
#
# Uses the SAME shared network function as the sc side (couplings_lib.R::tumor_partial_network), so
# any bulk-vs-sc difference is a DATA (resolution) effect, not a method difference. Bulk program
# scores use the FULL panels from levers.yaml (incl. the snRNA-dropout genes that only bulk detects;
# crowding uses its full mechanosensor set, not hippo_core). emt_axis cannot be tumor-scoped in bulk
# (no compartments) — that limitation is itself part of the answer.
#
# Outputs: results/couplings/bulk/{tumor_scores_bulk.csv, network_bulk_partial.csv,
#          concordance_bulk_vs_sc.csv, spatial_value_add.csv}
# Usage: Rscript phase1_mechanotypes/20_bulk_coupling.R   (run AFTER coupling_core.R)
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE); path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/")); dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "couplings_lib.R"))

main <- function() {
  cfg <- load_config(); set_seed_logged(cfg$features$seed, "bulk_coupling")
  lv <- yaml::read_yaml(file.path(cfg$root, "config", "levers.yaml"))
  programs <- c(lv$levers, lv$extrinsic_axes)
  lever_ids <- vapply(lv$levers, `[[`, "", "id"); axis_ids <- vapply(lv$extrinsic_axes, `[[`, "", "id")
  prog_ids <- c(lever_ids, axis_ids)

  # sc object: gene_lookup (symbol->ENSG) + counts/meta for a FAIR sc-pseudobulk reference
  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  sym2ens <- lapply(proc$gene_lookup, function(ids) unique(sub("\\..*$", "", ids)))   # for BULK (ENSG)
  sc_counts <- proc$counts
  sc_sid <- as.character(proc$meta$sample_id); sc_cs <- as.character(proc$meta$cell_state)
  # symbol(s) -> sc-counts rownames, via gene_lookup (same proven path as resolve_feature_genes / 18)
  gl <- proc$gene_lookup
  rows_for <- function(syms) intersect(unlist(gl[unlist(syms)]), rownames(sc_counts))

  # bulk quant: gene_id x library; collapse libraries -> SCPCS sample (mean), logCPM
  bq <- read.delim(resolve_path(cfg, cfg$paths$phase_a$bulk_quant_tsv), check.names = FALSE)
  bmeta <- read.delim(resolve_path(cfg, cfg$paths$phase_a$bulk_metadata_tsv), check.names = FALSE)
  lib2samp <- setNames(as.character(bmeta$sample_id), as.character(bmeta$library_id))
  M <- as.matrix(bq[, -1]); rownames(M) <- sub("\\..*$", "", bq$gene_id)
  samp <- lib2samp[colnames(M)]
  agg <- t(rowsum(t(M), group = samp, na.rm = TRUE))                     # gene x sample (sum libs)
  lib <- colSums(agg); lib[lib == 0] <- 1
  lcpm <- log2(t(t(agg) / lib * 1e6) + 1)

  score_bulk <- function(p) {                                           # signed logCPM, full panel
    pos <- unlist(p$genes_positive); neg <- unlist(p$genes_negative)
    pid <- intersect(unique(unlist(sym2ens[pos])), rownames(lcpm)); if (length(pid) < 3) return(NULL)
    s <- colMeans(lcpm[pid, , drop = FALSE])
    if (length(neg)) {
      nid <- intersect(unique(unlist(sym2ens[neg])), rownames(lcpm))
      if (length(nid) >= 1) s <- s - colMeans(lcpm[nid, , drop = FALSE])
    }
    as.numeric(scale(s))
  }
  X <- sapply(programs, function(p) { v <- score_bulk(p); if (is.null(v)) rep(NA, ncol(lcpm)) else v })
  colnames(X) <- prog_ids; rownames(X) <- colnames(lcpm)
  X <- X[, colSums(is.na(X)) == 0, drop = FALSE]                        # drop unscoreable programs
  prog_ids <- intersect(prog_ids, colnames(X))

  out <- resolve_path(cfg, "results/couplings/bulk"); ensure_dir(out)
  write.csv(data.frame(sample_id = rownames(X), round(X, 5), check.names = FALSE, row.names = NULL),
            file.path(out, "tumor_scores_bulk.csv"), row.names = FALSE)

  net <- tumor_partial_network(X, prog_ids, cfg$features$seed)
  write.csv(net$edges, file.path(out, "network_bulk_partial.csv"), row.names = FALSE)
  cat("\n===== BULK network — top |partial r| =====\n")
  print(format(head(net$edges[, c("a","b","partial_r","ci_lo","ci_hi","bh_fdr")], 10), digits = 2),
        row.names = FALSE)

  # ---- FAIR concordance: sc PSEUDOBULK (same method as bulk) vs bulk ---------------------------
  # Compute sc scores the SAME way as bulk (per-sample pseudobulk logCPM, signed) so any bulk<->sc
  # gap is a RESOLUTION effect, not an aggregation-method artifact. emt_axis is tumor-scoped on the
  # sc side (its lever definition); bulk cannot scope it — that asymmetry is part of the answer.
  suppressPackageStartupMessages(library(Matrix))
  pseudobulk_lcpm <- function(mask) {
    cols <- if (is.null(mask)) seq_along(sc_sid) else which(mask)
    ss <- sc_sid[cols]; samp <- sort(unique(ss))
    pb <- vapply(samp, function(s) as.numeric(Matrix::rowSums(sc_counts[, cols[ss == s], drop = FALSE])),
                 numeric(nrow(sc_counts)))
    rownames(pb) <- rownames(sc_counts); colnames(pb) <- samp
    l <- colSums(pb); l[l == 0] <- 1; log2(t(t(pb) / l * 1e6) + 1)
  }
  lcpm_all <- pseudobulk_lcpm(NULL); lcpm_tum <- pseudobulk_lcpm(sc_cs %in% c("blastemal", "epithelial"))
  score_sc <- function(p) {
    lc <- if (identical(p$id, "emt_axis")) lcpm_tum else lcpm_all
    pid <- intersect(rows_for(p$genes_positive), rownames(lc)); if (length(pid) < 3) return(NULL)
    s <- colMeans(lc[pid, , drop = FALSE])
    nid <- intersect(rows_for(p$genes_negative), rownames(lc))
    if (length(nid) >= 1) s <- s - colMeans(lc[nid, , drop = FALSE])
    setNames(as.numeric(scale(s)), colnames(lc))
  }
  sc_scores <- setNames(lapply(programs, function(p) {
    v <- score_sc(p); if (is.null(v)) NULL else v[colnames(lcpm_all)] }),
    vapply(programs, `[[`, "", "id"))
  common <- intersect(colnames(lcpm_all), rownames(X))
  conc <- do.call(rbind, lapply(prog_ids, function(p) {
    v <- sc_scores[[p]]; if (is.null(v)) return(NULL)
    r <- cor(v[common], X[common, p], use = "complete.obs")
    data.frame(program = p, r_bulk_sc = round(r, 3), n = length(common))
  }))
  write.csv(conc, file.path(out, "concordance_bulk_vs_sc.csv"), row.names = FALSE)
  cat("\n===== bulk<->sc per-program concordance (same tumors) =====\n"); print(conc, row.names = FALSE)

  # ---- spatial value-add: three-way tag per lever (PLAN §4.3, audit E/L) ----------------------
  drop <- read.csv(resolve_path(cfg, "results/couplings/gene_detection_sc_bulk.csv"))
  # dropout-prone program = a program whose panel has >=2 genes with sc_detect_frac < 0.02
  prog_dropout <- setNames(vapply(programs, function(p) {
    g <- c(unlist(p$genes_positive), unlist(p$genes_negative))
    sum(drop$sc_detect_frac[match(g, drop$gene)] < 0.02, na.rm = TRUE) >= 2
  }, logical(1)), vapply(programs, `[[`, "", "id"))
  tag <- sapply(prog_ids, function(p) {
    r <- conc$r_bulk_sc[conc$program == p]
    if (length(r) == 0 || is.na(r)) return("unscored")
    if (r >= 0.65) "recoverable_from_bulk"
    else if (isTRUE(prog_dropout[[p]])) "snrna_dropout(bulk_better)"
    else "needs_spatial(arrangement)"
  })
  sva <- data.frame(program = prog_ids,
                    r_bulk_sc = conc$r_bulk_sc[match(prog_ids, conc$program)],
                    tag = tag)
  write.csv(sva, file.path(out, "spatial_value_add.csv"), row.names = FALSE)
  cat("\n===== SPATIAL VALUE-ADD, per-lever (three-way, PLAN §4.3) =====\n"); print(sva, row.names = FALSE)

  # ---- network comparison: which COUPLINGS does bulk recover vs need cell resolution? ----------
  scnet <- read.csv(resolve_path(cfg, "results/couplings/network_tumorB_partial.csv"))
  key <- function(d) paste(pmin(d$a, d$b), pmax(d$a, d$b))
  scnet$k <- key(scnet); bl <- net$edges; bl$k <- key(bl)
  cmp <- merge(scnet[, c("k","a","b","partial_r","bh_fdr")],
               bl[, c("k","partial_r","bh_fdr")], by = "k", suffixes = c("_sc","_bulk"))
  cmp$sc_sig <- cmp$bh_fdr_sc < 0.10; cmp$bulk_sig <- cmp$bh_fdr_bulk < 0.10
  cmp$verdict <- ifelse(cmp$sc_sig & cmp$bulk_sig, "recovered_by_bulk",
                 ifelse(cmp$sc_sig & !cmp$bulk_sig, "NEEDS_cell_resolution",
                 ifelse(!cmp$sc_sig & cmp$bulk_sig, "bulk_only(composition?)", "neither")))
  cmp <- cmp[order(cmp$verdict, -abs(cmp$partial_r_sc)), c("a","b","partial_r_sc","partial_r_bulk","sc_sig","bulk_sig","verdict")]
  write.csv(cmp, file.path(out, "network_concordance.csv"), row.names = FALSE)
  cat("\n===== COUPLING recovery: sc (cell-scoped) vs bulk =====\n")
  print(cmp[cmp$verdict != "neither", ], row.names = FALSE)

  # ---- (§4.4) bulk -> favorable/anaplastic, LOOCV balanced accuracy --------------------------
  tmeta <- read.csv(resolve_path(cfg, "results/couplings/tumor_meta.csv"))
  y <- setNames(grepl("anapl", tolower(tmeta$subdiagnosis)), tmeta$sample_id)[rownames(X)]
  keep <- !is.na(y); Xc <- as.data.frame(X[keep, prog_ids, drop = FALSE]); yc <- as.integer(y[keep])
  loo <- vapply(seq_len(nrow(Xc)), function(i) {
    fit <- suppressWarnings(glm(yy ~ ., data = cbind(yy = yc[-i], Xc[-i, ]), family = binomial))
    as.numeric(predict(fit, newdata = Xc[i, , drop = FALSE], type = "response"))
  }, numeric(1))
  pred <- as.integer(loo >= 0.5)
  sens <- mean(pred[yc == 1] == 1); spec <- mean(pred[yc == 0] == 0)
  bal_acc <- mean(c(sens, spec))
  clf <- data.frame(n = length(yc), n_anaplastic = sum(yc), sensitivity = round(sens, 3),
                    specificity = round(spec, 3), balanced_accuracy = round(bal_acc, 3))
  write.csv(clf, file.path(out, "bulk_histology_clf.csv"), row.names = FALSE)
  cat("\n===== (§4.4) BULK -> favorable/anaplastic (LOOCV) =====\n"); print(clf, row.names = FALSE)
  cat(sprintf("\n[ok] bulk WS2 -> %s\n", out))
}

if (sys.nframe() == 0) main()
