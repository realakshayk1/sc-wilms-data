#!/usr/bin/env Rscript
# A-3: PROGNOSTIC association (binary). Ties per-sample compartment COMPOSITION and a
# pseudobulk PROLIFERATION score to two binary outcomes — relapse_status and
# vital_status — via logistic regression (univariate + covariate-adjusted for age/sex/
# metastasis) and Fisher's exact on tertiles, with bootstrap CIs over samples.
#
# HONEST SCOPE: there is no time-to-event survival in the local ScPCA metadata (only
# binary relapse/vital flags), so this is a binary prognostic association, NOT Cox/RFS.
# Proper time-to-event validation needs TARGET-WT (flagged as Tier-3). Relapse n is
# small (~10), so estimates are wide — reported with CIs, not as definitive effects.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages(library(logistf))

PROLIF <- c("MKI67","TOP2A","PCNA","CDK1","CCNB1","CCNB2","CCNA2","AURKA","AURKB","BUB1",
            "BUB1B","CENPA","CENPE","KIF11","FOXM1","PLK1","CDC20","UBE2C","BIRC5","TYMS",
            "RRM2","CDC25C","KIF23","NUSAP1")
# p53 transcriptional targets (apoptosis/arrest effectors) — for the ABM apoptosis lever
TP53_TARGETS <- c("CDKN1A","MDM2","BAX","GADD45A","BBC3","SFN","TP53I3","RRM2B","ZMAT3",
                  "FAS","SESN1","TNFRSF10B","PMAIP1")

logcpm_sample <- function(counts, samples) {
  requireNamespace("Matrix", quietly = TRUE)
  samp <- sort(unique(samples))
  pb <- vapply(samp, function(s)
    as.numeric(Matrix::rowSums(counts[, which(samples == s), drop = FALSE])),
    numeric(nrow(counts)))
  rownames(pb) <- rownames(counts); colnames(pb) <- samp
  lib <- colSums(pb); lib[lib == 0] <- 1
  log2(t(t(pb) / lib * 1e6) + 1)
}

sample_fractions <- function(sample_id, labels, levels_keep) {
  samp <- sort(unique(as.character(sample_id)))
  fr <- t(vapply(samp, function(s) {
    idx <- as.character(sample_id) == s & !is.na(labels) & labels %in% levels_keep
    if (!sum(idx)) return(setNames(rep(NA_real_, length(levels_keep)), levels_keep))
    as.numeric(table(factor(labels[idx], levels = levels_keep))) / sum(idx)
  }, numeric(length(levels_keep))))
  colnames(fr) <- levels_keep; rownames(fr) <- samp; fr
}

# Firth penalized logistic regression (logistf). Firth's bias-reduction is the
# standard remedy for (quasi-)complete separation, which ordinary glm hits here at
# ~5-10 events / 4 params (ordinary ORs diverge to 0/Inf and Wald p->1). logistf
# returns penalized-likelihood-ratio p-values and profile-likelihood CIs that stay
# finite under separation. Effects are per-SD (predictor z-scored).
logit_assoc <- function(df, predictor, outcome, covars = character(0), seed = 42) {
  d <- df[!is.na(df[[predictor]]) & !is.na(df[[outcome]]), , drop = FALSE]
  for (cv in covars) d <- d[!is.na(d[[cv]]), , drop = FALSE]
  if (nrow(d) < 8 || length(unique(d[[outcome]])) < 2) return(NULL)
  d[[predictor]] <- as.numeric(scale(d[[predictor]]))
  rhs <- paste(c(predictor, covars), collapse = " + ")
  fit <- tryCatch(logistf::logistf(as.formula(paste(outcome, "~", rhs)), data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  i <- match(predictor, names(coef(fit)))
  if (is.na(i)) return(NULL)
  data.frame(outcome = outcome, predictor = predictor,
             model = if (length(covars)) paste0("firth_adj(", paste(covars, collapse = "+"), ")") else "firth_univariate",
             n = nrow(d), n_event = sum(d[[outcome]] == 1),
             OR_per_SD = unname(exp(coef(fit)[i])), wald_p = fit$prob[i],
             OR_lo = unname(exp(fit$ci.lower[i])), OR_hi = unname(exp(fit$ci.upper[i])),
             stringsAsFactors = FALSE)
}

fisher_tertile <- function(df, predictor, outcome) {
  d <- df[!is.na(df[[predictor]]) & !is.na(df[[outcome]]), ]
  if (nrow(d) < 9 || length(unique(d[[outcome]])) < 2) return(NULL)
  hi <- d[[predictor]] >= quantile(d[[predictor]], 2/3)
  tab <- table(factor(hi, c(FALSE, TRUE)), factor(d[[outcome]], c(0, 1)))
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL); if (is.null(ft)) return(NULL)
  data.frame(outcome = outcome, predictor = predictor, model = "fisher_top-tertile",
             n = nrow(d), n_event = sum(d[[outcome]] == 1),
             OR_per_SD = unname(ft$estimate), wald_p = ft$p.value,
             OR_lo = unname(ft$conf.int[1]), OR_hi = unname(ft$conf.int[2]),
             stringsAsFactors = FALSE)
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "prognostic_association")
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir); ensure_dir(out_dir)
  requireNamespace("Matrix", quietly = TRUE)

  dat <- readRDS(resolve_path(cfg, cfg$paths$phase_a$scores_rds))
  meta <- dat$meta
  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts
  id2sym <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
  if (!is.null(proc$gene_lookup))
    for (sym in names(proc$gene_lookup))
      for (id in proc$gene_lookup[[sym]]) if (id %in% names(id2sym)) id2sym[id] <- sym

  samples <- as.character(meta$sample_id)
  # per-sample compartment fractions
  comp_lv <- intersect(c("blastemal","epithelial","stromal"), unique(meta$cell_state))
  fr <- sample_fractions(samples, as.character(meta$cell_state), comp_lv)
  # per-sample proliferation score = mean logCPM of proliferation genes (z across samples)
  lcpm <- logcpm_sample(counts, as.character(proc$meta$sample_id))
  score_of <- function(genes) {
    ids <- intersect(names(id2sym)[id2sym %in% genes], rownames(lcpm))
    if (length(ids) < 5) return(NULL)
    colMeans(lcpm[ids, , drop = FALSE])
  }
  prolif_score <- score_of(PROLIF)
  tp53_score   <- score_of(TP53_TARGETS)

  pred <- data.frame(sample_id = rownames(fr), fr, check.names = FALSE)
  colnames(pred)[-1] <- paste0(comp_lv, "_frac")
  if (!is.null(prolif_score)) {
    ps <- data.frame(sample_id = names(prolif_score), proliferation_score = as.numeric(scale(prolif_score)))
    pred <- merge(pred, ps, by = "sample_id", all.x = TRUE)
  }
  if (!is.null(tp53_score)) {
    ts <- data.frame(sample_id = names(tp53_score), tp53_target_score = as.numeric(scale(tp53_score)))
    pred <- merge(pred, ts, by = "sample_id", all.x = TRUE)
  }

  # clinical outcomes + covariates
  clin <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_clinical_metadata.tsv"))
  md <- read.delim(clin, check.names = FALSE)
  keep <- c("scpca_sample_id","relapse_status","vital_status","age","sex","metastasis","subdiagnosis")
  md <- unique(md[, intersect(keep, colnames(md))])
  md$relapse <- ifelse(md$relapse_status == "Yes", 1L, ifelse(md$relapse_status == "No", 0L, NA))
  vs <- tolower(as.character(md$vital_status))            # ScPCA uses "Expired"/"Alive"
  md$deceased <- ifelse(grepl("decea|dead|expir", vs), 1L, ifelse(grepl("alive|living", vs), 0L, NA))
  md$age_num <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(md$age))))
  md$sex_f <- ifelse(tolower(as.character(md$sex)) %in% c("female","f"), 1L,
                     ifelse(tolower(as.character(md$sex)) %in% c("male","m"), 0L, NA))
  # metastasis is free-text (e.g. "Pulmonary metastasis"); NA/blank = none recorded.
  mraw <- tolower(trimws(as.character(md$metastasis)))
  md$mets <- ifelse(is.na(md$metastasis) | mraw %in% c("", "na", "none", "no", "absent", "m0"), 0L, 1L)

  df <- merge(pred, md, by.x = "sample_id", by.y = "scpca_sample_id", all.x = TRUE)
  message(sprintf("[data] %d samples; relapse events=%d, deceased events=%d",
                  nrow(df), sum(df$relapse == 1, na.rm = TRUE), sum(df$deceased == 1, na.rm = TRUE)))

  # per-tumor scores table — the substrate for the ABM parameter mapping (ABM-1)
  pt_cols <- intersect(c("sample_id", paste0(comp_lv, "_frac"), "proliferation_score",
                         "tp53_target_score", "subdiagnosis", "relapse", "deceased"), colnames(df))
  pt <- df[, pt_cols]
  write.csv(pt, file.path(out_dir, "per_tumor_scores.csv"), row.names = FALSE)
  message("[ok] Per-tumor scores -> ", file.path(out_dir, "per_tumor_scores.csv"))

  predictors <- c(paste0(comp_lv, "_frac"), if (!is.null(prolif_score)) "proliferation_score")
  covars <- intersect(c("age_num","sex_f","mets"), colnames(df))
  rows <- list()
  for (outcome in c("relapse","deceased")) {
    if (length(unique(na.omit(df[[outcome]]))) < 2) { message("[skip] ", outcome, ": one class only"); next }
    for (pr in predictors) {
      rows <- c(rows, list(logit_assoc(df, pr, outcome, covars = character(0), seed = cfg$features$seed)))
      rows <- c(rows, list(logit_assoc(df, pr, outcome, covars = covars, seed = cfg$features$seed)))
      rows <- c(rows, list(fisher_tertile(df, pr, outcome)))
    }
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows)) {
    res <- do.call(rbind, rows)
    res$BH_FDR <- p.adjust(res$wald_p, "BH")
    res <- res[order(res$outcome, res$wald_p), ]
    write.csv(res, file.path(out_dir, "prognostic_association.csv"), row.names = FALSE)
    message("[ok] Prognostic association -> ", file.path(out_dir, "prognostic_association.csv"))
    print(format(res[, c("outcome","predictor","model","n","n_event","OR_per_SD","wald_p","BH_FDR")],
                 digits = 3), row.names = FALSE)
    sig <- res[!is.na(res$wald_p) & res$wald_p < 0.05, ]
    message(sprintf("[ok] %d associations p<0.05 (%d at BH-FDR<0.05)",
                    nrow(sig), sum(res$BH_FDR < 0.05, na.rm = TRUE)))
  } else message("[warn] no estimable associations")
}

if (sys.nframe() == 0) main()
