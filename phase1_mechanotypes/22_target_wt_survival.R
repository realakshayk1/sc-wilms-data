#!/usr/bin/env Rscript
# WS3 / P3 (§5.3): binary DNA/expression lever -> SURVIVAL (lever VALIDATION, not prognosis).
#
# A lever that stratifies outcome earns an orthogonal `survival_validated: true` annotation in
# joint_priors.yaml — it is NEVER fed back to re-weight the sweep (audit I). Time-to-event needs
# TARGET-WT (local ScPCA has only binary relapse/vital); this script is the machinery + the pull
# recipe, and self-tests on synthetic data so the Cox logic is verified even before the data lands.
#
# GATED: TARGET-WT is an external GDC pull (network + TCGAbiolinks/GenomicDataCommons). When absent,
# this runs the self-test and prints the recipe; it does NOT fabricate results.
#
# Usage: Rscript phase1_mechanotypes/22_target_wt_survival.R
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE); path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/")); dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages(library(survival))

# Cox per binary lever: univariate + adjusted (stage/histology); HR, CI, p, logrank; BH-FDR.
# surv_df: columns time (months), event (0/1), covariates, and one column per binary lever (0/1).
binary_survival <- function(surv_df, lever_cols, covars = character(0)) {
  ok_cov <- intersect(covars, colnames(surv_df))
  rows <- lapply(lever_cols, function(lv) {
    d <- surv_df[!is.na(surv_df[[lv]]) & !is.na(surv_df$time) & !is.na(surv_df$event), ]
    if (nrow(d) < 8 || length(unique(d[[lv]])) < 2 || sum(d$event) < 3) return(NULL)
    # build Surv() INSIDE the formula with data=d (never cbind a Surv object — that mangles it)
    f_uni <- as.formula(sprintf("survival::Surv(time, event) ~ `%s`", lv))
    uni <- tryCatch(survival::coxph(f_uni, data = d), error = function(e) NULL); if (is.null(uni)) return(NULL)
    su <- summary(uni)
    hr <- su$coefficients[1, "exp(coef)"]; p <- su$coefficients[1, "Pr(>|z|)"]
    ci <- su$conf.int[1, c("lower .95", "upper .95")]
    lr <- tryCatch(survival::survdiff(f_uni, data = d), error = function(e) NULL)
    lr_p <- if (is.null(lr)) NA_real_ else 1 - pchisq(lr$chisq, length(lr$n) - 1)
    hr_adj <- NA_real_; p_adj <- NA_real_
    if (length(ok_cov)) {
      f_adj <- as.formula(sprintf("survival::Surv(time, event) ~ `%s` + %s", lv,
                                  paste(sprintf("`%s`", ok_cov), collapse = " + ")))
      adj <- tryCatch(survival::coxph(f_adj, data = d), error = function(e) NULL)
      if (!is.null(adj)) { sa <- summary(adj); hr_adj <- sa$coefficients[1, "exp(coef)"]; p_adj <- sa$coefficients[1, "Pr(>|z|)"] }
    }
    data.frame(lever = lv, n = nrow(d), n_event = sum(d$event), hr = hr,
               hr_lo = ci[1], hr_hi = ci[2], p_univariate = p, logrank_p = lr_p,
               hr_adjusted = hr_adj, p_adjusted = p_adj)
  })
  res <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (!is.null(res)) { res$bh_fdr <- p.adjust(res$p_univariate, "BH"); res <- res[order(res$p_univariate), ] }
  res
}

# --- self-test: synthetic cohort where lever_true HALVES the hazard; lever_null is noise ---------
self_test <- function() {
  set_seed_logged(42, "survival_selftest")
  n <- 120
  lever_true <- rbinom(n, 1, 0.5); lever_null <- rbinom(n, 1, 0.5)
  base_haz <- 0.05 * exp(-0.7 * lever_true)                     # true lever -> lower hazard
  t_event <- rexp(n, base_haz); t_cens <- rexp(n, 0.02)
  time <- pmin(t_event, t_cens); event <- as.integer(t_event <= t_cens)
  df <- data.frame(time = time, event = event, stage = rbinom(n, 1, 0.5),
                   lever_true = lever_true, lever_null = lever_null)
  res <- binary_survival(df, c("lever_true", "lever_null"), covars = "stage")
  print(format(res[, c("lever","n","n_event","hr","hr_adjusted","p_univariate","bh_fdr")], digits = 3), row.names = FALSE)
  ht <- res$hr[res$lever == "lever_true"]; pt <- res$p_univariate[res$lever == "lever_true"]
  pn <- res$p_univariate[res$lever == "lever_null"]
  ha <- res$hr_adjusted[res$lever == "lever_true"]              # covariate-adjusted path must work too
  stopifnot(ht < 0.8, pt < 0.05, pn > 0.05, is.finite(ha), ha < 0.9)
  cat("[selftest] PASS — Cox (univariate + stage-adjusted) recovers the protective lever, ignores null.\n")
}

PULL_RECIPE <- '
TARGET-WT is an external GDC pull (not committed). To enable the real survival validation:
  1. install.packages("BiocManager"); BiocManager::install(c("TCGAbiolinks","survminer"))
  2. In R:
       library(TCGAbiolinks)
       q  <- GDCquery("TARGET-WT", "Transcriptome Profiling", workflow.type="STAR - Counts")
       GDCdownload(q); expr <- GDCprepare(q)                    # bulk RNA
       qm <- GDCquery("TARGET-WT", "Simple Nucleotide Variation",
                      data.type="Masked Somatic Mutation"); GDCdownload(qm); maf <- GDCprepare(qm)
       clin <- GDCquery_clinic("TARGET-WT","clinical")          # EFS/OS, stage, histology
     -> save to data/raw/target_wt/{expr.rds, maf.rds, clinical.tsv}
  3. Build a surv_df (time=EFS months, event, binary levers: mutation present/absent for
     CTNNB1/WT1/TP53/AMER1/DROSHA/SIX1/2/MYCN; expression hi/lo by the WS1 bifurcation split),
     then call binary_survival(surv_df, lever_cols, covars=c("stage","histology")).
  Guardrail (audit I): survival is a REPORT only; it never re-weights the joint_priors sweep.
'

main <- function() {
  cfg <- load_config()
  data_dir <- resolve_path(cfg, "data/raw/target_wt")
  have <- file.exists(file.path(data_dir, "clinical.tsv"))
  cat("=== WS3 survival machinery ===\n")
  self_test()
  if (!have) {
    cat("\n[gated] TARGET-WT not present at data/raw/target_wt/. Machinery verified above.\n")
    cat(PULL_RECIPE)
    return(invisible())
  }
  # --- real run (only when the pull has landed) -----------------------------------------------
  surv_df <- read.delim(file.path(data_dir, "surv_df.tsv"), check.names = FALSE)
  levers <- grep("^lever_", colnames(surv_df), value = TRUE)
  res <- binary_survival(surv_df, levers, covars = intersect(c("stage", "histology"), colnames(surv_df)))
  out <- resolve_path(cfg, "results/couplings/survival_levers.csv"); ensure_dir(dirname(out))
  write.csv(res, out, row.names = FALSE)
  cat(sprintf("\n[ok] survival levers -> %s (%d validated at BH-FDR<0.10)\n",
              out, sum(res$bh_fdr < 0.10, na.rm = TRUE)))
  print(res, row.names = FALSE)
}

if (sys.nframe() == 0) main()
