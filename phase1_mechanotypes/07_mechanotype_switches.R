#!/usr/bin/env Rscript
# FR-A7 (P1): Report cell states with mechanotype switches between histology groups.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
`%||%` <- function(a, b) if (!is.null(a)) a else b

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "mechanotype_switches")

  consensus_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir)
  out_csv <- resolve_path(cfg, cfg$paths$phase_a$switches_csv)
  ensure_dir(dirname(out_csv))

  files <- list.files(consensus_dir, pattern = "_consensus\\.rds$", full.names = TRUE)
  if (!length(files)) stop("Run 05_consensus_cluster.R first")

  rows <- list()
  for (f in files) {
    obj <- readRDS(f)
    feat <- obj$feature
    clusters <- obj$clusters
    names(clusters) <- names(obj$item_consensus)

    for (state in cfg$features$cell_states) {
      fav_id <- paste(state, "favorable", sep = "__")
      ana_id <- paste(state, "anaplastic", sep = "__")
      if (!(fav_id %in% names(clusters) && ana_id %in% names(clusters))) next
      switched <- clusters[[fav_id]] != clusters[[ana_id]]
      rows[[length(rows) + 1]] <- data.frame(
        feature = feat,
        cell_state = state,
        favorable_cluster = clusters[[fav_id]],
        anaplastic_cluster = clusters[[ana_id]],
        mechanotype_switch = switched,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) {
    message("[warn] No comparable state/histology pairs for switch analysis")
    return(invisible(NULL))
  }

  df <- do.call(rbind, rows)
  write.csv(df, out_csv, row.names = FALSE)
  n_switch <- sum(df$mechanotype_switch)
  message(sprintf("[ok] %d mechanotype switches across features -> %s", n_switch, out_csv))

  fig_dir <- resolve_path(cfg, cfg$paths$dirs$figures)
  ensure_dir(fig_dir)
  methods_file <- file.path(consensus_dir, "phase_a_methods.yaml")
  methods <- list(
    feature_scoring = "log1p(CPM_positive) - log1p(CPM_negative) per predefined gene program",
    wasserstein = "1-D Wasserstein-1 per feature on cell score distributions only",
    clustering_items = sprintf("(cell_state x histology) groups with >= %d cells", cfg$features$min_cells_per_item),
    consensus = "ConsensusClusterPlus PAM; k selected by low PAC + high Calinski-Harabasz",
    cell_state_labels = "Canonical marker-panel argmax on tumor cells (Option C; config/label_mapping.yaml marker_assignment); reference CellAssign/SingleR labels do NOT define Wilms compartments",
    primary_inference = "Patient-level: histology label permuted across samples (not cells); BH-FDR across 18 feature x compartment tests; see 09_distributional_validation.R",
    histology = "Sample subdiagnosis (Favorable vs Anaplastic)",
    seed = cfg$features$seed
  )
  yaml::write_yaml(methods, methods_file)
  message("[ok] Methods log -> ", methods_file)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_state, fill = mechanotype_switch)) +
      ggplot2::geom_bar(position = "dodge") +
      ggplot2::facet_wrap(~feature) +
      ggplot2::labs(
        title = "Mechanotype switches (favorable vs anaplastic)",
        y = "Count", x = "Cell state"
      ) +
      ggplot2::theme_bw()
    fig_file <- file.path(fig_dir, "mechanotype_switches.png")
    ggplot2::ggsave(fig_file, p, width = 10, height = 6, dpi = 150)
    message("[ok] Figure -> ", fig_file)
  }
  invisible(out_csv)
}

if (sys.nframe() == 0) main()
