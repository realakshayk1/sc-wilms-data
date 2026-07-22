#!/usr/bin/env Rscript
# WS2 §4.2: reference-based deconvolution of BULK RNA-seq into compartment fractions, validated
# against the paired single-cell composition (same 40 SCPCP tumors).
#
# METHOD NOTE: MuSiC/BisqueRNA have no binary build for this R 4.6.1 (2026-06-24) and need a
# source/Rtools compile; `nnls` installed cleanly, so we use NNLS reference-based deconvolution —
# the mathematical core of reference methods (CIBERSORT/Bisque): a marker signature matrix from the
# sc reference, then non-negative least squares per bulk sample. Honest limitation vs MuSiC: no
# cross-subject gene weighting; and the sc reference also defines the validation truth, so this
# tests deconvolution MECHANICS, not signature independence (plan §4.2).
#
# Output: results/couplings/bulk/deconv_fractions.csv + deconv_vs_sc.csv (per-compartment r, CCC)
# Usage: Rscript phase1_mechanotypes/25_bulk_deconv.R
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE); path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/")); dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages({library(Matrix); library(nnls)})

COMPARTMENTS <- c("blastemal", "epithelial", "stromal")
MIN_CELLS <- 100L; N_MARK <- 60L

ccc <- function(x, y) {                              # Lin's concordance correlation coefficient
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}

main <- function() {
  cfg <- load_config(); set_seed_logged(cfg$features$seed, "bulk_deconv")
  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts; cs <- as.character(proc$meta$cell_state); sid <- as.character(proc$meta$sample_id)
  types <- names(which(table(cs[!is.na(cs)]) >= MIN_CELLS))          # cell types with enough cells
  message(sprintf("[ref] %d cells, cell types used: %s", ncol(counts), paste(types, collapse = ", ")))

  # --- signature matrix: per-cell-type mean CPM, marker genes only -----------------------------
  cpm_type <- vapply(types, function(t) {
    idx <- which(cs == t & !is.na(cs))
    v <- as.numeric(Matrix::rowSums(counts[, idx, drop = FALSE])); v / sum(v) * 1e6
  }, numeric(nrow(counts)))
  rownames(cpm_type) <- rownames(counts)
  # markers: top N genes per type by log2 fold-change vs the mean of the others
  logc <- log2(cpm_type + 1)
  markers <- unique(unlist(lapply(seq_along(types), function(j) {
    fc <- logc[, j] - rowMeans(logc[, -j, drop = FALSE])
    expr <- cpm_type[, j] > 5
    names(sort(fc[expr], decreasing = TRUE))[seq_len(N_MARK)]
  })))
  markers <- intersect(markers, rownames(counts))
  C <- cpm_type[markers, , drop = FALSE]                              # signature: markers x types
  message(sprintf("[sig] %d marker genes across %d types", length(markers), length(types)))

  # --- bulk CPM over markers, NNLS per sample --------------------------------------------------
  bq <- read.delim(resolve_path(cfg, cfg$paths$phase_a$bulk_quant_tsv), check.names = FALSE)
  bmeta <- read.delim(resolve_path(cfg, cfg$paths$phase_a$bulk_metadata_tsv), check.names = FALSE)
  lib2samp <- setNames(as.character(bmeta$sample_id), as.character(bmeta$library_id))
  M <- as.matrix(bq[, -1]); rownames(M) <- sub("\\..*$", "", bq$gene_id)
  agg <- t(rowsum(t(M), group = lib2samp[colnames(M)], na.rm = TRUE))
  bcpm <- t(t(agg) / colSums(agg) * 1e6)
  gm <- intersect(sub("\\..*$", "", markers), rownames(bcpm))          # markers present in bulk (bare ENSG)
  # map sig rownames (counts ids) to bare-ENSG so rows line up with bulk
  sig_bare <- sub("\\..*$", "", rownames(C)); keep <- sig_bare %in% gm & !duplicated(sig_bare)
  Cb <- C[keep, , drop = FALSE]; rownames(Cb) <- sig_bare[keep]
  Bb <- bcpm[rownames(Cb), , drop = FALSE]

  fr <- t(apply(Bb, 2, function(b) {
    f <- nnls::nnls(Cb, b)$x; if (sum(f) == 0) return(rep(NA, ncol(Cb))); f / sum(f)
  }))
  colnames(fr) <- types
  frac <- as.data.frame(fr); frac$sample_id <- rownames(fr)
  out <- resolve_path(cfg, "results/couplings/bulk"); ensure_dir(out)
  write.csv(frac[, c("sample_id", types)], file.path(out, "deconv_fractions.csv"), row.names = FALSE)

  # --- validation vs single-cell composition (same tumors) -------------------------------------
  sc_fr <- prop.table(table(sid, cs), 1)
  cmp <- do.call(rbind, lapply(intersect(COMPARTMENTS, types), function(t) {
    common <- intersect(rownames(fr), rownames(sc_fr))
    x <- fr[common, t]; y <- as.numeric(sc_fr[common, t])
    data.frame(compartment = t, n = length(common),
               pearson_r = round(cor(x, y, use = "complete.obs"), 3),
               lin_ccc = round(ccc(x, y), 3),
               mean_bulk = round(mean(x, na.rm = TRUE), 3), mean_sc = round(mean(y), 3))
  }))
  write.csv(cmp, file.path(out, "deconv_vs_sc.csv"), row.names = FALSE)
  cat("\n===== bulk NNLS deconvolution vs single-cell composition =====\n"); print(cmp, row.names = FALSE)
  cat(sprintf("\n[ok] -> %s/deconv_fractions.csv + deconv_vs_sc.csv\n", out))
}

if (sys.nframe() == 0) main()
