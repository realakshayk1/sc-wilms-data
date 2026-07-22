#!/usr/bin/env Rscript
# WS1 coupling core (resolution-agnostic). Level A = within-tumor, cell-level PARTIAL-correlation
# network among the curated disjoint levers/axes, combined across tumors.
#
# Design decisions baked in (from the plan audit):
#   - PARTIAL correlation (conditioning on the other programs) via the precision matrix, so an edge
#     is the direct association, not a marginal one confounded by a third lever.
#   - Covariates = TECHNICAL ONLY: log1p(depth) + mito%. NOT cell-cycle phase (audit B: proliferation
#     is itself a node; regressing cell cycle would null the prolif edges).
#   - Restricted to TUMOR compartments (blastemal + epithelial): emt_axis is defined there, and this
#     is the tumor-intrinsic coupling structure (stroma is a different population).
#   - Combine tumors by Fisher-z, weighted by df; report a stability fraction + CI. n-per-tumor is
#     large (well-powered), so this is the PRIMARY evidence (tumor-level Level B is confirmatory).
#   - Standardization happens HERE (§2 one-convention): scores come in raw signed, z-scored inside.
#
# Usage: Rscript phase1_mechanotypes/coupling_core.R
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "couplings_lib.R"))

TUMOR_COMPARTMENTS <- c("blastemal", "epithelial")
MIN_CELLS <- 100L        # min complete-case cells for a tumor to contribute an estimate
RIDGE <- 1e-3            # light diagonal shrinkage so the correlation matrix inverts cleanly

# residualize a numeric vector on technical covariates (returns residuals; raw-centered if no covars)
resid_tech <- function(y, depth, mito) {
  ok <- is.finite(y)
  df <- data.frame(y = y, d = depth, m = mito)
  fit <- tryCatch(lm(y ~ d + m, data = df, na.action = na.exclude), error = function(e) NULL)
  r <- if (is.null(fit)) y - mean(y, na.rm = TRUE) else residuals(fit)
  as.numeric(r)
}

# partial correlation matrix from a cells x p residual matrix (precision-matrix method)
partial_cor <- function(M) {
  cc <- stats::complete.cases(M); M <- M[cc, , drop = FALSE]
  if (nrow(M) < MIN_CELLS) return(NULL)
  sdv <- apply(M, 2, sd); keep <- is.finite(sdv) & sdv > 1e-8
  if (sum(keep) < 3) return(NULL)
  C <- cor(M[, keep, drop = FALSE])
  C <- (1 - RIDGE) * C + RIDGE * diag(nrow(C))     # light shrinkage -> invertible
  P <- tryCatch(solve(C), error = function(e) NULL); if (is.null(P)) return(NULL)
  d <- sqrt(diag(P)); pr <- -P / outer(d, d); diag(pr) <- 1
  full <- matrix(NA_real_, ncol(M), ncol(M), dimnames = list(colnames(M), colnames(M)))
  full[keep, keep] <- pr
  attr(full, "n") <- nrow(M); full
}

# combine per-tumor partial-corr matrices: Fisher-z weighted by df=(n-p-1)
combine_edges <- function(mats, progs) {
  ut <- t(combn(progs, 2))
  out <- lapply(seq_len(nrow(ut)), function(k) {
    a <- ut[k, 1]; b <- ut[k, 2]; zs <- c(); ws <- c()
    for (m in mats) {
      if (is.null(m) || !(a %in% rownames(m)) || is.na(m[a, b])) next
      r <- max(min(m[a, b], 0.999), -0.999)
      zs <- c(zs, atanh(r)); ws <- c(ws, max(attr(m, "n") - length(progs) - 1, 1))
    }
    if (!length(zs)) return(NULL)
    zc <- sum(ws * zs) / sum(ws); se <- 1 / sqrt(sum(ws))
    data.frame(a = a, b = b, partial_r = tanh(zc),
               ci_lo = tanh(zc - 1.96 * se), ci_hi = tanh(zc + 1.96 * se),
               n_tumors = length(zs),
               stability = mean(abs(tanh(zs)) > 0.1 & sign(zs) == sign(zc)))
  })
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

run_stratum <- function(scores, sid, depth, mito, progs, label) {
  us <- unique(sid)
  mats <- lapply(us, function(s) {
    idx <- sid == s
    M <- sapply(progs, function(p) resid_tech(scores[idx, p], depth[idx], mito[idx]))
    colnames(M) <- progs
    partial_cor(M)
  })
  mats <- mats[!vapply(mats, is.null, logical(1))]
  ed <- combine_edges(mats, progs)
  if (!is.null(ed)) { ed$stratum <- label; ed <- ed[order(-abs(ed$partial_r)), ] }
  list(edges = ed, n_tumors = length(mats))
}

# ---- Level B: TUMOR-level (between-patient) partial network + bootstrap CIs -------------------
# This is the level the ABM is seeded at (per-tumor params), and where the coupling signal lives.
# n=40 -> partial corr is estimable (40 > 8 vars) but noisy: bootstrap tumors for CIs, light ridge.
tumor_level <- function(scores, sid, progs, out_dir, seed) {
  X <- vapply(progs, function(p) tapply(scores[, p], sid, function(v) mean(v, na.rm = TRUE)),
              numeric(length(unique(sid))))
  colnames(X) <- progs
  X <- scale(X)                                   # z across tumors (standardize here, §2)
  write.csv(data.frame(sample_id = rownames(X), round(X, 5), check.names = FALSE, row.names = NULL),
            file.path(out_dir, "tumor_scores.csv"), row.names = FALSE)
  net <- tumor_partial_network(X, progs, seed, B = 2000, ridge = RIDGE)   # shared (couplings_lib.R)
  ed <- net$edges
  write.csv(ed, file.path(out_dir, "network_tumorB_partial.csv"), row.names = FALSE)
  write.csv(round(net$marginal, 3), file.path(out_dir, "network_tumorB_marginal_matrix.csv"))
  cat("\n===== Level B (TUMOR-level PARTIAL corr, bootstrap CIs, n=40) — top |r| =====\n")
  print(format(head(ed[, c("a","b","partial_r","marginal_r","ci_lo","ci_hi","bh_fdr")], 14), digits = 2),
        row.names = FALSE)
  sig <- ed[ed$bh_fdr < 0.10, ]
  cat(sprintf("\n[Level B] %d/%d edges BH-FDR<0.10\n", nrow(sig), nrow(ed)))
  invisible(ed)
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "coupling_core")
  dat <- readRDS(resolve_path(cfg, "data/processed/lever_scores_percell.rds"))
  S <- dat$scores; meta <- dat$meta
  progs <- dat$program_ids

  tumor <- as.character(meta$cell_state) %in% TUMOR_COMPARTMENTS
  message(sprintf("[data] %d tumor-compartment cells of %d; %d programs",
                  sum(tumor), nrow(S), length(progs)))
  S <- S[tumor, , drop = FALSE]; meta <- meta[tumor, , drop = FALSE]
  sid <- as.character(meta$sample_id)
  depth <- log1p(as.numeric(meta$sum))
  mito <- as.numeric(meta$subsets_mito_percent)

  out_dir <- resolve_path(cfg, "results/couplings"); ensure_dir(out_dir)
  # per-tumor metadata (first row per sample) for downstream bifurcation/transfer/emit
  tm <- meta[!duplicated(sid), intersect(c("sample_id","subdiagnosis","histology","disease_timing"), colnames(meta))]
  write.csv(tm, file.path(out_dir, "tumor_meta.csv"), row.names = FALSE)
  strata <- list(all = rep(TRUE, nrow(meta)))
  sub <- tolower(as.character(meta$subdiagnosis))
  strata$favorable  <- grepl("favor", sub)
  strata$anaplastic <- grepl("anapl", sub)

  all_edges <- list(); ntab <- list()
  for (nm in names(strata)) {
    m <- strata[[nm]]; if (sum(m) < MIN_CELLS) next
    r <- run_stratum(S[m, , drop = FALSE], sid[m], depth[m], mito[m], progs, nm)
    if (!is.null(r$edges)) all_edges[[nm]] <- r$edges
    ntab[[nm]] <- r$n_tumors
  }
  edges <- do.call(rbind, all_edges)
  write.csv(edges, file.path(out_dir, "network_cellA_partial.csv"), row.names = FALSE)

  # mean partial-r matrix (all stratum) for a heatmap
  ea <- all_edges[["all"]]
  M <- matrix(NA_real_, length(progs), length(progs), dimnames = list(progs, progs)); diag(M) <- 1
  for (i in seq_len(nrow(ea))) { M[ea$a[i], ea$b[i]] <- ea$partial_r[i]; M[ea$b[i], ea$a[i]] <- ea$partial_r[i] }
  write.csv(round(M, 3), file.path(out_dir, "network_cellA_matrix.csv"))

  cat(sprintf("\n[tumors contributing] %s\n",
              paste(sprintf("%s=%d", names(ntab), unlist(ntab)), collapse = "  ")))
  cat("\n===== Level A (cell-level PARTIAL corr), 'all' stratum — top |r| =====\n")
  show <- ea[order(-abs(ea$partial_r)), c("a","b","partial_r","ci_lo","ci_hi","n_tumors","stability")]
  print(format(head(show, 14), digits = 2), row.names = FALSE)
  cat("\n[ok] -> results/couplings/network_cellA_partial.csv (+ _matrix.csv)\n")
  cat("[note] cell-level partial-r is near-null (see values): the coupling signal is BETWEEN-tumor,\n")
  cat("       not within-tumor between-cell. Building the seeding Sigma from the tumor level below.\n")

  # ---- PRIMARY for seeding: tumor-level between-patient network -------------------------------
  tumor_level(S, sid, progs, out_dir, cfg$features$seed)

  # favorable-vs-anaplastic edge shift (confirmatory; audit: report with CIs)
  if (all(c("favorable","anaplastic") %in% names(all_edges))) {
    f <- all_edges$favorable; a <- all_edges$anaplastic
    key <- function(d) paste(d$a, d$b)
    mg <- merge(f[,c("a","b","partial_r")], a[,c("a","b","partial_r")], by = c("a","b"),
                suffixes = c("_fav","_anap"))
    mg$delta <- mg$partial_r_anap - mg$partial_r_fav
    mg <- mg[order(-abs(mg$delta)), ]
    write.csv(mg, file.path(out_dir, "network_cellA_histology_shift.csv"), row.names = FALSE)
    cat("\n===== Largest favorable->anaplastic partial-r SHIFTS (top 8) =====\n")
    print(format(head(mg, 8), digits = 2), row.names = FALSE)
  }
}

if (sys.nframe() == 0) main()
