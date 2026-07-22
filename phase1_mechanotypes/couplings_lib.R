#!/usr/bin/env Rscript
# Shared coupling machinery (the resolution-agnostic core, PLAN §2). Both coupling_core.R (sc) and
# 20_bulk_coupling.R (bulk) call tumor_partial_network() on a samples x programs z-matrix, so "the
# same analysis at a different resolution" is literally the same function.

# partial-correlation matrix from a (samples x programs) matrix via the precision matrix
.pcor_of <- function(M, ridge) {
  C <- cor(M, use = "pairwise.complete.obs"); C <- (1 - ridge) * C + ridge * diag(nrow(C))
  P <- solve(C); d <- sqrt(diag(P)); pr <- -P / outer(d, d); diag(pr) <- 1; pr
}

# Tumor/sample-level partial network with bootstrap-sample CIs + BH-FDR.
# X: numeric matrix, rows = samples, cols include `progs`. Caller z-scores columns first.
tumor_partial_network <- function(X, progs, seed, B = 2000, ridge = 1e-3) {
  X <- X[, progs, drop = FALSE]
  base <- .pcor_of(X, ridge); marg <- cor(X, use = "pairwise.complete.obs")
  set.seed(seed); n <- nrow(X); ut <- t(combn(progs, 2)); boot <- array(NA_real_, c(nrow(ut), B))
  for (b in seq_len(B)) {
    pb <- tryCatch(.pcor_of(X[sample.int(n, n, replace = TRUE), , drop = FALSE], ridge),
                   error = function(e) NULL)
    if (!is.null(pb)) boot[, b] <- pb[ut]
  }
  ed <- data.frame(a = ut[, 1], b = ut[, 2], partial_r = base[ut], marginal_r = marg[ut],
                   ci_lo = apply(boot, 1, quantile, 0.025, na.rm = TRUE),
                   ci_hi = apply(boot, 1, quantile, 0.975, na.rm = TRUE))
  ed$boot_p <- apply(boot, 1, function(z) 2 * min(mean(z <= 0, na.rm = TRUE), mean(z >= 0, na.rm = TRUE)))
  ed$bh_fdr <- p.adjust(ed$boot_p, "BH")
  list(edges = ed[order(-abs(ed$partial_r)), ], marginal = marg)
}
