#!/usr/bin/env Rscript
# A2/A3: Compartment & subgroup COMPOSITION analysis (the signal Phase A actually
# carries). Per-sample fractions of triphasic compartments and fetal subgroups,
# tested favorable-vs-anaplastic AND relapse-vs-no-relapse.
#
# Rigor: the SAMPLE is the unit (one fraction vector per patient). Tested with
# Wilcoxon rank-sum on raw fractions and on centred-log-ratio (CLR) values
# (compositional-data aware), BH-FDR within each (axis x representation).
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
`%||%` <- function(a, b) if (!is.null(a)) a else b

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  cmp <- outer(x, y, "-"); (sum(cmp > 0) - sum(cmp < 0)) / (length(x) * length(y))
}

# per-sample fraction matrix for a categorical label column (rows=sample, cols=levels)
sample_fractions <- function(sample_id, labels, levels_keep) {
  ok <- !is.na(labels) & labels %in% levels_keep
  samp <- unique(as.character(sample_id))
  fr <- t(vapply(samp, function(s) {
    idx <- as.character(sample_id) == s & ok
    n <- sum(idx)
    if (n == 0) return(setNames(rep(NA_real_, length(levels_keep)), levels_keep))
    tab <- table(factor(labels[idx], levels = levels_keep))
    as.numeric(tab) / n
  }, numeric(length(levels_keep))))
  colnames(fr) <- levels_keep; rownames(fr) <- samp
  fr
}

clr <- function(mat, pseudo = 1e-3) {
  m <- mat + pseudo
  lg <- log(m)
  lg - rowMeans(lg)
}

test_composition <- function(fr, grp, axis, representation, pos_lab, neg_lab) {
  fr_clr <- clr(fr)
  rows <- list()
  for (comp in colnames(fr)) {
    pos_raw <- fr[grp == "pos", comp]; neg_raw <- fr[grp == "neg", comp]
    pos_clr <- fr_clr[grp == "pos", comp]; neg_clr <- fr_clr[grp == "neg", comp]
    pos_raw <- pos_raw[!is.na(pos_raw)]; neg_raw <- neg_raw[!is.na(neg_raw)]
    pos_clr <- pos_clr[!is.na(pos_clr)]; neg_clr <- neg_clr[!is.na(neg_clr)]
    if (length(pos_raw) < 3 || length(neg_raw) < 3) next
    p_raw <- tryCatch(stats::wilcox.test(pos_raw, neg_raw)$p.value, error = function(e) NA_real_)
    p_clr <- tryCatch(stats::wilcox.test(pos_clr, neg_clr)$p.value, error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.frame(
      axis = axis, representation = representation, component = comp,
      pos_label = pos_lab, neg_label = neg_lab,
      n_pos = length(pos_raw), n_neg = length(neg_raw),
      mean_frac_pos = mean(pos_raw), mean_frac_neg = mean(neg_raw),
      cliffs_delta = cliffs_delta(pos_raw, neg_raw),
      wilcox_p_raw = p_raw, wilcox_p_clr = p_clr,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  df$p_BH_clr <- p.adjust(df$wilcox_p_clr, "BH")
  df
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "composition_analysis")
  dat <- readRDS(resolve_path(cfg, cfg$paths$phase_a$scores_rds))
  meta <- dat$meta
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir)
  ensure_dir(out_dir)

  # contrast grouping vectors
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

  comp_levels <- intersect(c("blastemal", "epithelial", "stromal"), unique(meta$cell_state))
  sub_levels <- if ("subgroup" %in% colnames(meta)) {
    intersect(c("CM_blastemal", "CM_epithelial", "PV", "UB", "fibroblast"), unique(meta$subgroup))
  } else character(0)

  representations <- list(compartment = list(col = meta$cell_state, levels = comp_levels))
  if (length(sub_levels)) representations$subgroup <- list(col = meta$subgroup, levels = sub_levels)

  axes <- list(
    histology = list(grp = g_hist, pos = "anaplastic", neg = "favorable"),
    relapse   = list(grp = g_rel,  pos = "relapse",    neg = "no_relapse")
  )

  all_rows <- list()
  for (an in names(axes)) {
    ax <- axes[[an]]
    # one grouping value per sample (sample-level) — collapse cell-level grp to sample
    samp <- unique(as.character(meta$sample_id))
    samp_grp <- vapply(samp, function(s) {
      g <- ax$grp[as.character(meta$sample_id) == s]
      g <- g[!is.na(g)]
      if (!length(g)) NA_character_ else g[1]
    }, character(1))
    for (rn in names(representations)) {
      rep_def <- representations[[rn]]
      fr <- sample_fractions(meta$sample_id, rep_def$col, rep_def$levels)
      fr <- fr[samp, , drop = FALSE]
      grp <- samp_grp[rownames(fr)]
      df <- test_composition(fr, grp, an, rn, ax$pos, ax$neg)
      if (!is.null(df)) all_rows[[length(all_rows) + 1]] <- df
    }
  }
  res <- do.call(rbind, all_rows)
  out_csv <- file.path(out_dir, "composition_analysis.csv")
  write.csv(res, out_csv, row.names = FALSE)
  message("[ok] Composition analysis -> ", out_csv)

  # report
  res2 <- res[order(res$axis, res$representation, res$wilcox_p_clr), ]
  print(format(res2[, c("axis", "representation", "component", "n_pos", "n_neg",
                        "mean_frac_pos", "mean_frac_neg", "cliffs_delta",
                        "wilcox_p_clr", "p_BH_clr")], digits = 2), row.names = FALSE)
  sig <- res[!is.na(res$p_BH_clr) & res$p_BH_clr < 0.05, ]
  message(sprintf("[ok] %d composition tests significant at BH-FDR(CLR) < 0.05", nrow(sig)))
  invisible(res)
}

if (sys.nframe() == 0) main()
