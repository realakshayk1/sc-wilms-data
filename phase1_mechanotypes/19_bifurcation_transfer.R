#!/usr/bin/env Rscript
# WS1 P1 (§3.3 + §3.4 + §3.5): tumor-level bifurcation, intrinsic->extrinsic transfer, and emit
# config/joint_priors.yaml. Reads the tumor-level scores written by coupling_core.R.
#
# Dependency-free stats (diptest/mclust unavailable): bimodality = bimodality-coefficient (BC>0.555)
# AND a 1-D 2-component Gaussian-mixture EM whose BIC beats the 1-component fit with both components
# non-trivial. Transfer = standardized lm(axis ~ levers) with bootstrap-tumor CIs. n=40 -> everything
# is reported with CIs and framed as bounding, not best-fit.
#
# Usage: Rscript phase1_mechanotypes/19_bifurcation_transfer.R   (run AFTER coupling_core.R)
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE); path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/")); dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

# ---- bimodality --------------------------------------------------------------------------------
bimodality_coef <- function(x) {                 # BC>0.555 suggests bimodality (Pfister et al.)
  x <- x[is.finite(x)]; n <- length(x); m <- mean(x); s <- sd(x); if (s == 0) return(NA)
  g <- mean(((x - m) / s)^3); k <- mean(((x - m) / s)^4) - 3
  (g^2 + 1) / (k + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}
em2 <- function(x, iters = 300) {                # 1-D two-component Gaussian mixture (EM)
  x <- x[is.finite(x)]; n <- length(x)
  mu <- as.numeric(quantile(x, c(0.25, 0.75))); sg <- rep(sd(x), 2); pi <- c(.5, .5)
  for (i in seq_len(iters)) {
    d1 <- pi[1] * dnorm(x, mu[1], sg[1]); d2 <- pi[2] * dnorm(x, mu[2], sg[2])
    r <- d2 / (d1 + d2 + 1e-300)
    pi[2] <- mean(r); pi[1] <- 1 - pi[2]
    mu[1] <- sum((1 - r) * x) / sum(1 - r); mu[2] <- sum(r * x) / sum(r)
    sg[1] <- sqrt(sum((1 - r) * (x - mu[1])^2) / sum(1 - r)); sg[2] <- sqrt(sum(r * (x - mu[2])^2) / sum(r))
    sg <- pmax(sg, 1e-3)
  }
  ll2 <- sum(log(pi[1] * dnorm(x, mu[1], sg[1]) + pi[2] * dnorm(x, mu[2], sg[2]) + 1e-300))
  ll1 <- sum(dnorm(x, mean(x), sd(x), log = TRUE))
  list(bic1 = -2 * ll1 + 2 * log(n), bic2 = -2 * ll2 + 5 * log(n),
       mu = mu, sg = sg, pi = pi, split = mean(mu))
}
bifurcation <- function(S) {
  do.call(rbind, lapply(colnames(S), function(p) {
    x <- S[, p]; bc <- bimodality_coef(x); e <- em2(x)
    bimodal <- isTRUE(bc > 0.555) && (e$bic2 < e$bic1 - 2) && min(e$pi) > 0.15 &&
               abs(diff(e$mu)) > 1.0
    data.frame(program = p, bimodality_coef = round(bc, 3),
               dBIC_2vs1 = round(e$bic1 - e$bic2, 2), min_prop = round(min(e$pi), 3),
               mode_gap_sd = round(abs(diff(e$mu)), 2), split = round(e$split, 3),
               bimodal = bimodal)
  }))
}

# ---- transfer: standardized lm(axis ~ levers) + bootstrap CIs -----------------------------------
transfer <- function(S, levers, axes, seed, B = 2000) {
  set.seed(seed); n <- nrow(S)
  do.call(rbind, lapply(axes, function(ax) {
    y <- as.numeric(scale(S[, ax])); Xz <- scale(S[, levers, drop = FALSE])
    fit <- lm(y ~ Xz); b <- coef(fit)[-1]; names(b) <- levers
    boot <- sapply(seq_len(B), function(i) {
      ix <- sample.int(n, n, replace = TRUE)
      coef(lm(scale(S[ix, ax]) ~ scale(S[ix, levers, drop = FALSE])))[-1]
    })
    lo <- apply(boot, 1, quantile, .025, na.rm = TRUE); hi <- apply(boot, 1, quantile, .975, na.rm = TRUE)
    pval <- apply(boot, 1, function(z) 2 * min(mean(z <= 0), mean(z >= 0)))
    data.frame(axis = ax, lever = levers, beta = round(b, 3),
               ci_lo = round(lo, 3), ci_hi = round(hi, 3), boot_p = round(pval, 3))
  }))
}

main <- function() {
  cfg <- load_config(); set_seed_logged(cfg$features$seed, "bifurcation_transfer")
  cd <- resolve_path(cfg, "results/couplings")
  lv <- yaml::read_yaml(file.path(cfg$root, "config", "levers.yaml"))
  lever_ids <- vapply(lv$levers, `[[`, "", "id"); axis_ids <- vapply(lv$extrinsic_axes, `[[`, "", "id")

  ts <- read.csv(file.path(cd, "tumor_scores.csv"), check.names = FALSE)
  S <- as.matrix(ts[, c(lever_ids, axis_ids)]); rownames(S) <- ts$sample_id
  net <- read.csv(file.path(cd, "network_tumorB_partial.csv"))

  bif <- bifurcation(S); write.csv(bif, file.path(cd, "bifurcation.csv"), row.names = FALSE)
  cat("\n===== (3.3) BIFURCATION (tumor cohort, n=40) =====\n"); print(bif, row.names = FALSE)

  tr <- transfer(S, lever_ids, axis_ids, cfg$features$seed)
  write.csv(tr, file.path(cd, "transfer.csv"), row.names = FALSE)
  cat("\n===== (3.4) TRANSFER  extrinsic axis ~ levers (standardized, boot CI) =====\n")
  print(tr[order(tr$axis, -abs(tr$beta)), ], row.names = FALSE)

  # ---- (3.5) emit config/joint_priors.yaml ----------------------------------------------------
  fdr_edges <- net[net$bh_fdr < 0.10, ]
  qz <- function(p) as.numeric(round(quantile(S[, p], c(.05, .95), na.rm = TRUE), 3))
  levers_out <- setNames(lapply(lever_ids, function(p) {
    br <- bif[bif$program == p, ]
    list(z_p5_p95 = qz(p), bimodal = isTRUE(br$bimodal), split = if (isTRUE(br$bimodal)) br$split else NULL)
  }), lever_ids)
  couplings_out <- lapply(seq_len(nrow(fdr_edges)), function(i) {
    e <- fdr_edges[i, ]; list(a = e$a, b = e$b, partial_r = round(e$partial_r, 3),
                              ci = c(round(e$ci_lo, 3), round(e$ci_hi, 3)), bh_fdr = round(e$bh_fdr, 3))
  })
  transfer_out <- setNames(lapply(axis_ids, function(ax) {
    t <- tr[tr$axis == ax, ]
    setNames(lapply(seq_len(nrow(t)), function(i)
      list(beta = t$beta[i], ci = c(t$ci_lo[i], t$ci_hi[i]), boot_p = t$boot_p[i])), t$lever)
  }), axis_ids)
  jp <- list(
    provenance = list(resolution = "tumor-level (sc pseudobulk)", cohort = "SCPCP000006",
                      n_tumors = nrow(S), seed = cfg$features$seed,
                      generated = format(Sys.time(), "%Y-%m-%d"),
                      source_scripts = c("18_program_scores_percell.R", "coupling_core.R",
                                         "19_bifurcation_transfer.R"),
                      note = paste("Couplings are BETWEEN-tumor (cell-level near-null).",
                                   "Panels: config/levers.yaml. No fitting to outcome.")),
    levers = levers_out,
    couplings = list(fdr_threshold = 0.10, edges = couplings_out,
                     matrix_csv = "results/couplings/network_tumorB_marginal_matrix.csv"),
    transfer = list(
      measured = transfer_out,
      definitional = list(
        adhesion_motility = list(defined_by = "emt_axis", maps_to = c("adhesion_strength","migration_speed"),
                                 note = "literature prior (omics_to_params); not regressed"),
        igf_uptake = list(defined_by = "igf", maps_to = "igf_uptake_rate",
                          note = "literature prior; not regressed"))),
    sweep = list(mode = "virtual_cohort", cov = "tumor-level correlation, ridge-shrunk, FDR edges only",
                 n_draws = 256, fallback = "config/phase_c.yaml:omics_to_params (k scalars)")
  )
  out <- file.path(cfg$root, "config", "joint_priors.yaml")
  writeLines(c("# AUTO-GENERATED by phase1_mechanotypes/19_bifurcation_transfer.R — do not hand-edit.",
               "# Coupled, bounded prior for the PhysiCell virtual-cohort sweep (see PLAN §1).",
               yaml::as.yaml(jp)), out)
  cat(sprintf("\n[ok] emitted -> %s  (%d FDR couplings, %d bimodal levers)\n",
              "config/joint_priors.yaml", length(couplings_out), sum(bif$bimodal)))
}

if (sys.nframe() == 0) main()
