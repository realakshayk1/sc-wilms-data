#!/usr/bin/env Rscript
# WS1 (advisor framing): use INTRINSIC state (proliferation + p53) to label tumors/cells as
# "expressive" (aggressive: high proliferation, low p53 activity) vs "sensitive" (restrained: low
# proliferation, high p53), then show the CONDITIONAL distribution of the EXTRINSIC axes
# (crowding, EMT, hypoxia) across those classes. This is the transfer (§3.4) shown as a labeled
# contrast rather than a regression coefficient — "intrinsic prolif/p53 -> lower crowding, etc."
#
# Exports (for 21_coupling_figures.py to plot):
#   results/couplings/sensitive_expressive_tumor.csv   (40 tumors: scores + class)
#   results/couplings/sensitive_expressive_cellsub.csv (subsampled cells: scores + class)
#   results/couplings/sensitive_expressive_stats.csv   (per axis x level: means, MWU p, Cliff's delta)
# Usage: Rscript phase1_mechanotypes/24_sensitive_expressive.R   (after coupling_core.R)
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE); path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/")); dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

EXTRINSIC <- c("crowding_sensitivity", "emt_axis", "hypoxia_tolerance")

cliffs_delta <- function(a, b) {                    # P(a>b) - P(a<b), in [-1,1]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]; if (!length(a) || !length(b)) return(NA_real_)
  mean(outer(a, b, ">")) - mean(outer(a, b, "<"))
}

# label rows by intrinsic (proliferation, tp53) median split
label_intrinsic <- function(prolif, p53) {
  ph <- prolif > median(prolif, na.rm = TRUE); th <- p53 > median(p53, na.rm = TRUE)
  ifelse(ph & !th, "expressive", ifelse(!ph & th, "sensitive", "intermediate"))
}

axis_stats <- function(df, level) {
  do.call(rbind, lapply(intersect(EXTRINSIC, colnames(df)), function(ax) {
    e <- df[[ax]][df$class == "expressive"]; s <- df[[ax]][df$class == "sensitive"]
    p <- tryCatch(wilcox.test(e, s)$p.value, error = function(z) NA_real_)
    data.frame(level = level, axis = ax,
               mean_expressive = round(mean(e, na.rm = TRUE), 3),
               mean_intermediate = round(mean(df[[ax]][df$class == "intermediate"], na.rm = TRUE), 3),
               mean_sensitive = round(mean(s, na.rm = TRUE), 3),
               mwu_p = signif(p, 3), cliffs_delta = round(cliffs_delta(e, s), 3),
               n_expressive = sum(df$class == "expressive"), n_sensitive = sum(df$class == "sensitive"))
  }))
}

main <- function() {
  cfg <- load_config(); set_seed_logged(cfg$features$seed, "sensitive_expressive")
  cd <- resolve_path(cfg, "results/couplings")

  # ---- TUMOR level (primary; matches per-tumor ABM seeding, where the couplings live) ----------
  ts <- read.csv(file.path(cd, "tumor_scores.csv"), check.names = FALSE)
  ts$class <- label_intrinsic(ts$proliferation, ts$tp53_target)
  ts$intrinsic_aggression <- round(ts$proliferation - ts$tp53_target, 3)   # continuous axis
  write.csv(ts[, c("sample_id", "class", "intrinsic_aggression", "proliferation", "tp53_target",
                   EXTRINSIC)], file.path(cd, "sensitive_expressive_tumor.csv"), row.names = FALSE)
  st_t <- axis_stats(ts, "tumor")

  # ---- CELL level (advisor's literal framing; expected attenuated — couplings are between-tumor) -
  st_c <- NULL
  rds <- resolve_path(cfg, "data/processed/lever_scores_percell.rds")
  if (file.exists(rds)) {
    d <- readRDS(rds); S <- as.data.frame(d$scores)
    tum <- as.character(d$meta$cell_state) %in% c("blastemal", "epithelial")
    S <- S[tum, ]
    S$class <- label_intrinsic(S$proliferation, S$tp53_target)
    keep <- intersect(c("proliferation", "tp53_target", EXTRINSIC, "class"), colnames(S))
    sub <- S[sample.int(nrow(S), min(4000L, nrow(S))), keep]      # subsample for plotting
    write.csv(sub, file.path(cd, "sensitive_expressive_cellsub.csv"), row.names = FALSE)
    st_c <- axis_stats(S, "cell")
  }
  stats <- rbind(st_t, st_c)
  write.csv(stats, file.path(cd, "sensitive_expressive_stats.csv"), row.names = FALSE)

  cat("\n=== class sizes (tumor) ===\n"); print(table(ts$class))
  cat("\n=== extrinsic axis by intrinsic class (expressive vs sensitive) ===\n")
  print(stats, row.names = FALSE)
  cat("\n[interpretation] expressive = high-proliferation/low-p53 (aggressive intrinsic state);\n",
      "  crowding_sensitivity should be LOWER in expressive tumors (less contact-inhibition brake).\n")
  cat("[ok] -> sensitive_expressive_{tumor,cellsub,stats}.csv\n")
}

if (sys.nframe() == 0) main()
