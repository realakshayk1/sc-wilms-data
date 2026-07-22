#!/usr/bin/env Rscript
# WS3 real run: build the TARGET-WT OS survival table from the pulled open data and run the Cox
# lever-validation (binary_survival from 22_target_wt_survival.R). Primary endpoint = OS.
#
# Binary levers:
#   expression hi/lo (median split) for each config/levers.yaml lever, scored from TARGET bulk
#     (proliferation, tp53_target, wnt_canonical, blastemal_nephrogenic, igf, emt_axis)
#   mutation present/absent (open MAF): TP53, WT1, CTNNB1  (thin — 38 cases with WXS)
# Guardrail (audit I): survival is reported only; it never re-weights the joint_priors sweep.
#
# Usage: Rscript phase1_mechanotypes/27_target_wt_survival_run.R  (after target_wt_pull.py)
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE); path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/")); dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "22_target_wt_survival.R"))   # binary_survival() (main() is guarded)
suppressPackageStartupMessages(library(Matrix))

main <- function() {
  cfg <- load_config(); set_seed_logged(cfg$features$seed, "target_wt_survival")
  dir <- resolve_path(cfg, "data/raw/target_wt")
  clin <- read.delim(file.path(dir, "clinical.tsv"), check.names = FALSE)
  er <- read.delim(file.path(dir, "expr_counts.tsv"), check.names = FALSE)
  expr <- rowsum(as.matrix(er[, -1]), group = as.character(er[[1]]))   # collapse duplicate ENSG by sum
  mut  <- tryCatch(read.delim(file.path(dir, "mutations.tsv"), check.names = FALSE), error = function(e) NULL)
  lv   <- yaml::read_yaml(file.path(cfg$root, "config", "levers.yaml"))

  # symbol -> ENSG via the project gene_lookup (rownames of expr are bare ENSG)
  gl <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))$gene_lookup
  sym2ens <- lapply(gl, function(ids) unique(sub("\\..*$", "", ids)))
  rownames(expr) <- sub("\\..*$", "", rownames(expr))

  # OS: time (months) + event (death); Dead->days_to_death, else days_to_last_follow_up
  t_days <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_follow_up)
  os <- data.frame(case = clin$case, time = as.numeric(t_days) / 30.44,
                   event = as.integer(clin$vital_status == "Dead"))
  os <- os[is.finite(os$time) & os$time > 0 & !is.na(os$event), ]

  # expression levers: log2 CPM signed score per case, median-split to hi/lo
  lib <- colSums(expr); lib[lib == 0] <- 1; lcpm <- log2(t(t(expr) / lib * 1e6) + 1)
  score_case <- function(p) {
    pid <- intersect(unique(unlist(sym2ens[unlist(p$genes_positive)])), rownames(lcpm)); if (length(pid) < 3) return(NULL)
    s <- colMeans(lcpm[pid, , drop = FALSE])
    neg <- unlist(p$genes_negative)
    if (length(neg)) { nid <- intersect(unique(unlist(sym2ens[neg])), rownames(lcpm))
      if (length(nid) >= 1) s <- s - colMeans(lcpm[nid, , drop = FALSE]) }
    s
  }
  sdf <- os
  for (p in lv$levers) {
    sc <- score_case(p); if (is.null(sc)) next
    hi <- setNames(as.integer(sc > median(sc, na.rm = TRUE)), colnames(lcpm))
    sdf[[paste0("lever_", p$id, "_hi")]] <- hi[sdf$case]
  }
  # mutation levers (present/absent) — only defined for cases in the MAF cohort
  if (!is.null(mut) && nrow(mut)) {
    maf_cases <- unique(mut$case)
    for (g in c("TP53", "WT1", "CTNNB1")) {
      has <- sdf$case %in% mut$case[mut$gene == g]
      v <- ifelse(sdf$case %in% maf_cases, as.integer(has), NA_integer_)   # NA outside MAF cohort
      if (length(unique(na.omit(v))) == 2) sdf[[paste0("lever_", g, "_mut")]] <- v
    }
  }

  levers <- grep("^lever_", colnames(sdf), value = TRUE)
  res <- binary_survival(sdf, levers, covars = character(0))
  out <- resolve_path(cfg, "results/couplings/survival_levers.csv"); ensure_dir(dirname(out))
  write.csv(res, out, row.names = FALSE)
  write.csv(sdf, resolve_path(cfg, "results/couplings/target_wt_surv_df.csv"), row.names = FALSE)  # for KM figures
  cat(sprintf("\n=== TARGET-WT OS survival (n=%d cases, %d deaths) ===\n",
              nrow(sdf), sum(sdf$event)))
  print(format(res[, c("lever","n","n_event","hr","hr_lo","hr_hi","p_univariate","bh_fdr")], digits = 3),
        row.names = FALSE)
  cat(sprintf("\n[ok] %d/%d levers stratify OS at BH-FDR<0.10 -> %s\n",
              sum(res$bh_fdr < 0.10, na.rm = TRUE), nrow(res), out))
}

if (sys.nframe() == 0) main()
