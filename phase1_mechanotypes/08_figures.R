#!/usr/bin/env Rscript
# Publication-style figures for Phase A mechanotypes.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

plot_w1_heatmaps <- function(w_dir, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 for figures")
  }
  files <- list.files(w_dir, pattern = "_w1_dist\\.rds$", full.names = TRUE)
  if (!length(files)) return(invisible(NULL))

  rows <- list()
  for (f in files) {
    obj <- readRDS(f)
    D <- obj$distance
    feat <- obj$feature
    for (i in seq_len(nrow(D))) {
      for (j in seq_len(ncol(D))) {
        if (i == j) next
        rows[[length(rows) + 1]] <- data.frame(
          feature = feat,
          from = rownames(D)[i],
          to = colnames(D)[j],
          w1 = D[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df <- do.call(rbind, rows)
  df$from <- gsub("__", " / ", df$from)
  df$to <- gsub("__", " / ", df$to)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = to, y = from, fill = w1)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::facet_wrap(~feature, scales = "free", ncol = 3) +
    ggplot2::scale_fill_viridis_c(option = "magma", name = "W1") +
    ggplot2::labs(
      title = "Pairwise Wasserstein-1 distances (1-D feature scores)",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(face = "bold")
    )

  out <- file.path(fig_dir, "phase_a_w1_heatmaps.png")
  ggplot2::ggsave(out, p, width = 12, height = 8, dpi = 180)
  message("[fig] ", out)
}

plot_switch_heatmap <- function(switches_csv, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!file.exists(switches_csv)) return(invisible(NULL))

  df <- read.csv(switches_csv, stringsAsFactors = FALSE)
  df$switch_label <- ifelse(df$mechanotype_switch, "Switch", "Same mechanotype")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_state, y = feature, fill = switch_label)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::scale_fill_manual(
      values = c("Switch" = "#E64B35", "Same mechanotype" = "#4DBBD5"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Mechanotype switches: favorable vs anaplastic histology",
      subtitle = "Per feature × Wilms compartment (blastemal / epithelial / stromal)",
      x = "Cell compartment", y = "Feature program"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  out <- file.path(fig_dir, "phase_a_mechanotype_switch_heatmap.png")
  ggplot2::ggsave(out, p, width = 8, height = 5.5, dpi = 180)
  message("[fig] ", out)
}

plot_score_distributions <- function(scores_rds, fig_dir, features = c("wt1_activity", "blastemal_program", "proliferation")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!file.exists(scores_rds)) return(invisible(NULL))

  dat <- readRDS(scores_rds)
  meta <- dat$meta
  scores <- dat$scores

  rows <- list()
  for (feat in intersect(features, colnames(scores))) {
    rows[[length(rows) + 1]] <- data.frame(
      feature = feat,
      score = scores[, feat],
      cell_state = meta$cell_state,
      histology = meta$histology,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  df <- df[!is.na(df$cell_state) & !is.na(df$histology), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_state, y = score, fill = histology)) +
    ggplot2::geom_violin(scale = "width", alpha = 0.7, position = ggplot2::position_dodge(0.8)) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.3, position = ggplot2::position_dodge(0.8)) +
    ggplot2::facet_wrap(~feature, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = c(favorable = "#00A087", anaplastic = "#E64B35")) +
    ggplot2::labs(
      title = "1-D feature score distributions by compartment and histology",
      x = NULL, y = "Score (log CPM contrast)", fill = "Histology"
    ) +
    ggplot2::theme_bw(base_size = 10)

  out <- file.path(fig_dir, "phase_a_score_distributions.png")
  ggplot2::ggsave(out, p, width = 8, height = 9, dpi = 180)
  message("[fig] ", out)
}

plot_consensus_metrics <- function(consensus_dir, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!requireNamespace("tidyr", quietly = TRUE)) return(invisible(NULL))

  files <- list.files(consensus_dir, pattern = "_consensus\\.rds$", full.names = TRUE)
  if (!length(files)) return(invisible(NULL))

  rows <- list()
  for (f in files) {
    obj <- readRDS(f)
    for (k in names(obj$pac_by_k)) {
      rows[[length(rows) + 1]] <- data.frame(
        feature = obj$feature,
        k = as.integer(k),
        PAC = obj$pac_by_k[k],
        CHI = obj$chi_by_k[k],
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)
  long <- tidyr::pivot_longer(df, c(PAC, CHI), names_to = "metric", values_to = "value")

  p <- ggplot2::ggplot(long, ggplot2::aes(x = k, y = value, color = metric, group = metric)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~feature, scales = "free_y", ncol = 3) +
    ggplot2::scale_color_manual(values = c(PAC = "#3C5488", CHI = "#F39B7F")) +
    ggplot2::labs(
      title = "Consensus clustering metrics by k",
      subtitle = "Low PAC + high Calinski–Harabasz guides k selection",
      x = "Number of clusters (k)", y = NULL, color = NULL
    ) +
    ggplot2::theme_bw(base_size = 10)

  out <- file.path(fig_dir, "phase_a_consensus_metrics.png")
  ggplot2::ggsave(out, p, width = 11, height = 7, dpi = 180)
  message("[fig] ", out)
}

plot_distributional_validation <- function(dist_csv, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  # Combine both contrasts (histology + relapse) into one panel; show the null.
  files <- c(dist_csv, sub("\\.csv$", "_relapse.csv", dist_csv))
  dfs <- lapply(files[file.exists(files)], read.csv, stringsAsFactors = FALSE)
  if (!length(dfs)) return(invisible(NULL))
  df <- do.call(rbind, dfs)
  if (!"significant_BH" %in% colnames(df)) df$significant_BH <- df$p_perm_BH < 0.05
  df$significant_BH[is.na(df$significant_BH)] <- FALSE

  p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_state, y = w1_observed, fill = significant_BH)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_point(ggplot2::aes(y = w1_null_median), shape = 95, size = 4, color = "black") +
    ggplot2::facet_grid(contrast ~ feature, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                               name = "BH-FDR < 0.05") +
    ggplot2::labs(
      title = "Within-compartment program W1: observed (bar) vs permutation null (–)",
      subtitle = "Patient-level permutation; 0/18 significant on both axes (method-robust negative)",
      x = "Compartment", y = "Wasserstein-1"
    ) +
    ggplot2::theme_bw(base_size = 9)

  out <- file.path(fig_dir, "phase_a_distributional_validation.png")
  ggplot2::ggsave(out, p, width = 13, height = 6, dpi = 170)
  message("[fig] ", out)
}

plot_composition <- function(comp_csv, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!file.exists(comp_csv)) return(invisible(NULL))
  df <- read.csv(comp_csv, stringsAsFactors = FALSE)
  df <- df[df$representation == "compartment", , drop = FALSE]
  if (!nrow(df)) return(invisible(NULL))
  df$sig <- !is.na(df$p_BH_clr) & df$p_BH_clr < 0.05
  long <- rbind(
    data.frame(axis = df$axis, component = df$component, group = df$pos_label,
               frac = df$mean_frac_pos, sig = df$sig),
    data.frame(axis = df$axis, component = df$component, group = df$neg_label,
               frac = df$mean_frac_neg, sig = df$sig)
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(x = component, y = frac, fill = group)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~axis) +
    ggplot2::labs(
      title = "Compartment composition by clinical axis (the Phase A signal)",
      subtitle = "Epithelial up in anaplastic / stromal up in favorable (BH-FDR<0.05, CLR Wilcoxon)",
      x = "Compartment", y = "Mean per-sample fraction", fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10)
  out <- file.path(fig_dir, "phase_a_composition.png")
  ggplot2::ggsave(out, p, width = 10, height = 5, dpi = 180)
  message("[fig] ", out)
}

plot_loso_stability <- function(loso_csv, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!file.exists(loso_csv)) return(invisible(NULL))

  df <- read.csv(loso_csv, stringsAsFactors = FALSE)
  if (!all(c("switch_frequency", "stable_switch") %in% colnames(df))) return(invisible(NULL))  # deprecated schema
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_state, y = switch_frequency, fill = stable_switch)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~feature, ncol = 3) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#00A087", "FALSE" = "#B09C85"), name = "Stable (>=80%)") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      title = "LOSO mechanotype switch stability",
      x = "Compartment", y = "Switch frequency"
    ) +
    ggplot2::theme_bw(base_size = 10)

  out <- file.path(fig_dir, "phase_a_loso_stability.png")
  ggplot2::ggsave(out, p, width = 11, height = 7, dpi = 180)
  message("[fig] ", out)
}

plot_waddr_decomposition <- function(waddr_csv, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!file.exists(waddr_csv)) return(invisible(NULL))

  df <- read.csv(waddr_csv, stringsAsFactors = FALSE)
  rows <- list()
  for (state in c("blastemal", "epithelial", "stromal")) {
    fav <- paste(state, "favorable", sep = "__")
    ana <- paste(state, "anaplastic", sep = "__")
    sub <- df[(df$item_a == fav & df$item_b == ana) | (df$item_a == ana & df$item_b == fav), , drop = FALSE]
    if (!nrow(sub)) next
    for (i in seq_len(nrow(sub))) {
      rows[[length(rows) + 1]] <- data.frame(
        feature = sub$feature[i],
        cell_state = state,
        location = sub$location[i],
        size = sub$size[i],
        shape = sub$shape[i],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(invisible(NULL))
  long <- do.call(rbind, rows)
  if (!requireNamespace("tidyr", quietly = TRUE)) return(invisible(NULL))
  long <- tidyr::pivot_longer(long, c(location, size, shape), names_to = "component", values_to = "value")

  p <- ggplot2::ggplot(long, ggplot2::aes(x = cell_state, y = value, fill = component)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::facet_wrap(~feature, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = "waddR decomposition: favorable vs anaplastic (per compartment)",
      x = "Compartment", y = "2-Wasserstein component", fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10)

  out <- file.path(fig_dir, "phase_a_waddr_decomposition.png")
  ggplot2::ggsave(out, p, width = 11, height = 7, dpi = 180)
  message("[fig] ", out)
}

main <- function() {
  cfg <- load_config()
  fig_dir <- resolve_path(cfg, cfg$paths$dirs$figures)
  ensure_dir(fig_dir)

  plot_w1_heatmaps(
    resolve_path(cfg, cfg$paths$phase_a$wasserstein_dir),
    fig_dir
  )
  plot_switch_heatmap(
    resolve_path(cfg, cfg$paths$phase_a$switches_csv),
    fig_dir
  )
  plot_score_distributions(
    resolve_path(cfg, cfg$paths$phase_a$scores_rds),
    fig_dir
  )
  plot_consensus_metrics(
    resolve_path(cfg, cfg$paths$phase_a$consensus_dir),
    fig_dir
  )
  if (!is.null(cfg$paths$phase_a$distributional_csv)) {
    plot_distributional_validation(
      resolve_path(cfg, cfg$paths$phase_a$distributional_csv),
      fig_dir
    )
  }
  if (!is.null(cfg$paths$phase_a$loso_stability_csv)) {
    plot_loso_stability(
      resolve_path(cfg, cfg$paths$phase_a$loso_stability_csv),
      fig_dir
    )
  }
  plot_composition(
    file.path(resolve_path(cfg, cfg$paths$phase_a$consensus_dir), "composition_analysis.csv"),
    fig_dir
  )
  waddr_csv <- file.path(resolve_path(cfg, cfg$paths$phase_a$consensus_dir), "waddR_decomposition.csv")
  plot_waddr_decomposition(waddr_csv, fig_dir)
  message("[ok] Phase A figures -> ", fig_dir)
}

if (sys.nframe() == 0) main()
