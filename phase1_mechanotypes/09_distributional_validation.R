#!/usr/bin/env Rscript
# Distributional validation with PATIENT-LEVEL inference, over multiple contrasts.
#
# Rigor notes (see audit):
#  * Experimental unit is the SAMPLE (patient), not the cell. Cell-level permutation
#    is pseudoreplication and inflates significance (Squair et al. 2021). p-values are
#    built by permuting the contrast label ACROSS SAMPLES.
#  * Cliff's delta and bootstrap W1 CIs are effect-size DESCRIPTORS, not inference.
#  * BH-FDR across (feature x compartment) within each contrast.
#  * Contrasts:
#      - histology : favorable vs anaplastic (morphologic axis)
#      - relapse   : relapsed vs not (RFS axis; Yang et al. 2025 report the molecular
#                    signal tracks relapse, not simply favorable/anaplastic)
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "03_compute_scores.R"), local = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b
.str_hash <- function(s) sum(utf8ToInt(as.character(s)), na.rm = TRUE) %% 10000
isTRUE_vec <- function(x) !is.na(x) & x

wasserstein_1d <- function(x, y) {
  x <- as.numeric(x[!is.na(x)]); y <- as.numeric(y[!is.na(y)])
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  transport::wasserstein1d(x, y)
}

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  cmp <- outer(x, y, "-")
  (sum(cmp > 0) - sum(cmp < 0)) / (length(x) * length(y))
}

# Patient-level permutation: permute the (sample -> group) label, keep each cell
# tied to its true sample, recompute the cell-pooled W1. grp values in {"pos","neg"}.
sample_perm_p_w1 <- function(score, sample_id, grp, n_perm, seed) {
  obs <- wasserstein_1d(score[grp == "pos"], score[grp == "neg"])
  out0 <- list(w1 = obs, p_perm = NA_real_, n_pos_samples = NA_integer_,
               n_neg_samples = NA_integer_, null_median = NA_real_)
  if (!is.finite(obs)) return(out0)
  smap <- unique(data.frame(s = as.character(sample_id), g = as.character(grp),
                            stringsAsFactors = FALSE))
  smap <- smap[!duplicated(smap$s), , drop = FALSE]
  labs <- smap$g
  n_pos <- sum(labs == "pos"); n_neg <- sum(labs == "neg")
  out0$n_pos_samples <- n_pos; out0$n_neg_samples <- n_neg
  if (n_pos < 2L || n_neg < 2L) return(out0)
  cell_sample <- as.character(sample_id)
  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(i) {
    perm <- sample(labs)
    pos_samples <- smap$s[perm == "pos"]
    is_pos <- cell_sample %in% pos_samples
    wasserstein_1d(score[is_pos], score[!is_pos])
  }, numeric(1))
  n_ok <- sum(is.finite(null))
  out0$p_perm <- (sum(null >= obs, na.rm = TRUE) + 1) / (n_ok + 1)
  out0$null_median <- median(null, na.rm = TRUE)
  out0
}

bootstrap_w1_ci <- function(x, y, n_boot, seed) {
  set.seed(seed)
  vals <- vapply(seq_len(n_boot), function(i) {
    wasserstein_1d(sample(x, length(x), replace = TRUE), sample(y, length(y), replace = TRUE))
  }, numeric(1))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(NA_real_, NA_real_))
  stats::quantile(vals, probs = c(0.025, 0.975), names = FALSE)
}

.subsample_idx <- function(n, max_n, seed) {
  if (n <= max_n) return(seq_len(n))
  set.seed(seed); sort(sample.int(n, max_n))
}

# Run all (feature x compartment) patient-level tests for one contrast.
run_contrast <- function(scores, meta, feature_ids, cell_states, grp, cfg,
                         n_perm, n_boot, max_cells, min_cells, seed, label) {
  rows <- list()
  for (feat in feature_ids) {
    for (state in cell_states) {
      sel <- which(meta$cell_state == state & grp %in% c("pos", "neg"))
      if (!length(sel)) next
      sc <- as.numeric(scores[sel, feat]); sid <- as.character(meta$sample_id[sel]); g <- grp[sel]
      pos_i <- which(g == "pos"); neg_i <- which(g == "neg")
      if (length(pos_i) < min_cells || length(neg_i) < min_cells) next
      p_keep <- pos_i[.subsample_idx(length(pos_i), max_cells, seed + .str_hash(paste("p", label, feat, state)))]
      n_keep <- neg_i[.subsample_idx(length(neg_i), max_cells, seed + .str_hash(paste("n", label, feat, state)))]
      keep <- c(p_keep, n_keep)
      perm <- sample_perm_p_w1(sc[keep], sid[keep], g[keep], n_perm,
                               seed + .str_hash(paste(label, feat)) + .str_hash(state))
      ci <- bootstrap_w1_ci(sc[p_keep], sc[n_keep], n_boot, seed + 1 + .str_hash(paste(label, feat)))
      rows[[length(rows) + 1]] <- data.frame(
        contrast = label, feature = feat, cell_state = state,
        n_pos_cells = length(pos_i), n_neg_cells = length(neg_i),
        n_pos_samples = perm$n_pos_samples, n_neg_samples = perm$n_neg_samples,
        w1_observed = perm$w1, w1_null_median = perm$null_median, p_perm_sample = perm$p_perm,
        w1_ci_low = ci[1], w1_ci_high = ci[2],
        cliffs_delta = cliffs_delta(sc[p_keep], sc[n_keep]),
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)
  if (!is.null(df)) {
    df$p_perm_BH <- p.adjust(df$p_perm_sample, method = "BH")
    df$significant_BH <- isTRUE_vec(df$p_perm_BH < 0.05)
  }
  df
}

main <- function() {
  if (!requireNamespace("transport", quietly = TRUE)) stop("Install transport")
  cfg <- load_config()
  seed <- set_seed_logged(cfg$features$seed, "distributional_validation")
  phase_b_path <- file.path(cfg$root, "config", "phase_b.yaml")
  pb <- if (file.exists(phase_b_path)) yaml::read_yaml(phase_b_path) else list()
  val <- pb$validation %||% list()
  n_perm <- as.integer(val$permutation_n %||% 999)
  n_boot <- as.integer(val$bootstrap_n %||% 500)
  n_null <- as.integer(val$null_gene_sets %||% 30)
  max_cells <- as.integer(val$max_cells_per_test %||% 2000L)
  min_cells <- cfg$features$min_cells_per_item

  scores_rds <- resolve_path(cfg, cfg$paths$phase_a$scores_rds)
  out_csv <- resolve_path(cfg, cfg$paths$phase_a$distributional_csv)
  loso_csv <- resolve_path(cfg, cfg$paths$phase_a$loso_stability_csv)
  null_csv <- resolve_path(cfg, cfg$paths$phase_a$null_feature_csv)
  ensure_dir(dirname(out_csv))

  if (!file.exists(scores_rds)) stop("Run 03_compute_scores.R first")
  dat <- readRDS(scores_rds)
  scores <- dat$scores; meta <- dat$meta
  if (!"sample_id" %in% colnames(meta)) stop("sample_id missing - re-run ingest + QC")

  # --- Build contrast grouping vectors (values: "pos"/"neg"/NA) ---
  contrasts <- list()
  # histology: pos = anaplastic, neg = favorable
  hist <- tolower(as.character(meta$histology))
  g_hist <- ifelse(hist == "anaplastic", "pos", ifelse(hist == "favorable", "neg", NA_character_))
  contrasts[["histology"]] <- list(grp = g_hist, out = out_csv,
                                    note = "pos=anaplastic, neg=favorable")
  # relapse: join from ScPCA clinical metadata; pos = Yes, neg = No
  clin <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_clinical_metadata.tsv"))
  if (file.exists(clin)) {
    md <- read.delim(clin, check.names = FALSE)
    u <- unique(md[, c("scpca_sample_id", "relapse_status")])
    rmap <- setNames(as.character(u$relapse_status), u$scpca_sample_id)
    rs <- rmap[as.character(meta$sample_id)]
    g_rel <- ifelse(rs == "Yes", "pos", ifelse(rs == "No", "neg", NA_character_))
    contrasts[["relapse"]] <- list(grp = g_rel,
                                   out = sub("\\.csv$", "_relapse.csv", out_csv),
                                   note = "pos=relapse, neg=no relapse")
  } else {
    message("[warn] clinical metadata not found (", clin, "); skipping relapse contrast")
  }

  for (nm in names(contrasts)) {
    ct <- contrasts[[nm]]
    df <- run_contrast(scores, meta, dat$feature_ids, cfg$features$cell_states, ct$grp, cfg,
                       n_perm, n_boot, max_cells, min_cells, seed, nm)
    if (is.null(df)) { message("[warn] no testable groups for contrast ", nm); next }
    write.csv(df, ct$out, row.names = FALSE)
    message(sprintf("[ok] %s contrast (%s) -> %s", nm, ct$note, ct$out))
    message(sprintf("    %d/%d significant at BH-FDR < 0.05 | min BH p = %.3f",
                    sum(df$significant_BH, na.rm = TRUE), nrow(df), min(df$p_perm_BH, na.rm = TRUE)))
  }

  # --- Null random gene-set sensitivity (blastemal, histology contrast; control) ---
  seurat_rds <- resolve_path(cfg, cfg$paths$phase_a$seurat_rds)
  if (file.exists(seurat_rds)) {
    requireNamespace("Matrix", quietly = TRUE)
    proc <- readRDS(seurat_rds); counts <- proc$counts; lookup <- proc$gene_lookup
    state <- "blastemal"
    fav_idx <- which(meta$cell_state == state & tolower(meta$histology) == "favorable")
    ana_idx <- which(meta$cell_state == state & tolower(meta$histology) == "anaplastic")
    if (length(fav_idx) >= min_cells && length(ana_idx) >= min_cells) {
      real_w1 <- wasserstein_1d(scores[fav_idx, "blastemal_program"], scores[ana_idx, "blastemal_program"])
      set.seed(seed + 999)
      fav_sub <- sample(fav_idx, min(2000L, length(fav_idx)))
      ana_sub <- sample(ana_idx, min(2000L, length(ana_idx)))
      sub_counts <- counts[, c(fav_sub, ana_sub), drop = FALSE]; n_f <- length(fav_sub)
      score_random <- function(sd) {
        set.seed(sd)
        det <- Matrix::rowMeans(sub_counts > 0); elig <- rownames(sub_counts)[det >= 0.05]
        if (length(elig) < 12L) elig <- rownames(sub_counts)
        pos <- sample(elig, 6L); neg <- sample(setdiff(elig, pos), 6L)
        score_feature(sub_counts, pos, neg, gene_lookup = NULL)
      }
      null_w1 <- vapply(seq_len(n_null), function(i) {
        rnd <- score_random(seed + 1000 + i)
        wasserstein_1d(rnd[seq_len(n_f)], rnd[n_f + seq_along(ana_sub)])
      }, numeric(1))
      write.csv(data.frame(
        compartment = state, real_feature = "blastemal_program", real_w1 = real_w1,
        null_median_w1 = median(null_w1, na.rm = TRUE),
        null_q95_w1 = as.numeric(stats::quantile(null_w1, 0.95, na.rm = TRUE)),
        real_exceeds_null_q95 = real_w1 > stats::quantile(null_w1, 0.95, na.rm = TRUE),
        n_null_sets = n_null), null_csv, row.names = FALSE)
      message("[ok] Null gene-set sensitivity -> ", null_csv)
    }
  }
  invisible(NULL)
}

if (sys.nframe() == 0) main()
