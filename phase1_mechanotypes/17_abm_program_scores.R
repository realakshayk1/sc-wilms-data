#!/usr/bin/env Rscript
# A-3b: ABM PROGRAM SCORES. Augments results/mechanotypes/per_tumor_scores.csv with the
# extrinsic-axis program scores consumed by the Phase C ABM parameter mapping
# (phase2_histology_ml/17_positives_to_abm.py): EMT (epithelial / mesenchymal), contact
# inhibition (crowding), and hypoxia tolerance.
#
# Method is IDENTICAL to the proliferation/tp53 scores in 16_prognostic_association.R:
# per-sample pseudobulk logCPM, program score = mean logCPM(positive) - mean logCPM(negative),
# z-scored across samples. Gene sets are read from config/abm_programs.yaml (single source of
# truth, shared with the Python mapping). This runs AFTER 16 (which writes the base table);
# it merges new columns in and rewrites per_tumor_scores.csv, so it is non-destructive to the
# existing composition / proliferation / tp53 columns. No logistf dependency.
#
# Usage: scripts\rscript.bat phase1_mechanotypes\17_abm_program_scores.R
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages(library(Matrix))

MIN_POS_GENES <- 3L                      # need >=3 positive genes present to form a score

# per-sample pseudobulk logCPM (same as 16_prognostic_association.R::logcpm_sample)
logcpm_sample <- function(counts, samples) {
  samp <- sort(unique(samples))
  pb <- vapply(samp, function(s)
    as.numeric(Matrix::rowSums(counts[, which(samples == s), drop = FALSE])),
    numeric(nrow(counts)))
  rownames(pb) <- rownames(counts); colnames(pb) <- samp
  lib <- colSums(pb); lib[lib == 0] <- 1
  log2(t(t(pb) / lib * 1e6) + 1)
}

main <- function() {
  cfg <- load_config()
  pt_path <- resolve_path(cfg, file.path(cfg$paths$phase_a$consensus_dir, "per_tumor_scores.csv"))
  if (!file.exists(pt_path))
    stop("run 16_prognostic_association.R first (per_tumor_scores.csv missing): ", pt_path)
  programs <- yaml::read_yaml(file.path(cfg$root, "config", "abm_programs.yaml"))$programs

  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts
  id2sym <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
  if (!is.null(proc$gene_lookup))
    for (sym in names(proc$gene_lookup))
      for (id in proc$gene_lookup[[sym]]) if (id %in% names(id2sym)) id2sym[id] <- sym
  lcpm <- logcpm_sample(counts, as.character(proc$meta$sample_id))

  # signed program score: mean logCPM(pos) - mean logCPM(neg), then z across samples
  score_of <- function(pos, neg) {
    pid <- intersect(names(id2sym)[id2sym %in% pos], rownames(lcpm))
    if (length(pid) < MIN_POS_GENES) return(NULL)
    s <- colMeans(lcpm[pid, , drop = FALSE])
    nid <- if (length(neg)) intersect(names(id2sym)[id2sym %in% neg], rownames(lcpm)) else character(0)
    if (length(nid) >= 1) s <- s - colMeans(lcpm[nid, , drop = FALSE])
    setNames(as.numeric(scale(s)), names(s))
  }

  pt <- read.csv(pt_path, check.names = FALSE, stringsAsFactors = FALSE)
  added <- character(0)
  for (pg in programs) {
    col <- paste0(pg$id, "_score")
    sc <- score_of(unlist(pg$genes_positive), unlist(pg$genes_negative))
    if (is.null(sc)) {
      message(sprintf("[skip] %-22s: <%d positive genes detected", col, MIN_POS_GENES))
      next
    }
    pt[[col]] <- sc[as.character(pt$sample_id)]     # NA for samples absent from lcpm
    added <- c(added, col)
    message(sprintf("[ok]   %-22s: %d/%d samples scored (n pos genes present)",
                    col, sum(!is.na(pt[[col]])), nrow(pt)))
  }

  # --- compartment-resolved scores (cell x tumor type; addresses ABM direction 2b) -------
  # Same signed-logCPM score but on per-(sample x compartment) pseudobulk, z-scored WITHIN
  # each compartment across samples. Emits <program>_score__<compartment> columns so the ABM
  # can give e.g. blastemal vs epithelial cells compartment-specific EMT/crowding/hypoxia.
  MIN_CELLS_GRP <- 20L
  COMPARTMENTS <- c("blastemal", "epithelial", "stromal")
  dat <- readRDS(resolve_path(cfg, cfg$paths$phase_a$scores_rds))
  cs <- as.character(dat$meta$cell_state)
  sid <- as.character(proc$meta$sample_id)
  if (length(cs) == ncol(counts) && length(sid) == ncol(counts)) {
    grp <- paste(sid, cs, sep = "||")
    lcpm_g <- logcpm_sample(counts, grp)
    gsize <- table(grp)
    score_grp <- function(pos, neg) {               # raw signed score per group column
      pid <- intersect(names(id2sym)[id2sym %in% pos], rownames(lcpm_g))
      if (length(pid) < MIN_POS_GENES) return(NULL)
      s <- colMeans(lcpm_g[pid, , drop = FALSE])
      nid <- if (length(neg)) intersect(names(id2sym)[id2sym %in% neg], rownames(lcpm_g)) else character(0)
      if (length(nid) >= 1) s <- s - colMeans(lcpm_g[nid, , drop = FALSE])
      s
    }
    n_comp <- 0
    for (pg in programs) {
      sc <- score_grp(unlist(pg$genes_positive), unlist(pg$genes_negative))
      if (is.null(sc)) next
      for (comp in COMPARTMENTS) {
        cols <- grep(paste0("\\|\\|", comp, "$"), names(sc), value = TRUE)
        if (!length(cols)) next
        v <- sc[cols]
        v[as.integer(gsize[cols]) < MIN_CELLS_GRP] <- NA        # drop tiny pseudobulk groups
        z <- as.numeric(scale(v))                               # z within compartment
        names(z) <- sub(paste0("\\|\\|", comp, "$"), "", cols)  # -> sample ids
        pt[[paste0(pg$id, "_score__", comp)]] <- z[as.character(pt$sample_id)]
        n_comp <- n_comp + 1
      }
    }
    message(sprintf("[ok]   compartment-resolved: %d <program>x<compartment> score columns", n_comp))
  } else {
    message("[warn] cell_state not aligned to counts; skipped compartment-resolved scores")
  }

  write.csv(pt, pt_path, row.names = FALSE)
  message(sprintf("[ok] per_tumor_scores.csv augmented with %d tumor-level program score(s): %s",
                  length(added), paste(added, collapse = ", ")))
  message("[next] re-run phase2_histology_ml/17_positives_to_abm.py to activate the axes.")
}

main()
