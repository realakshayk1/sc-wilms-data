#!/usr/bin/env Rscript
# FR-A5: Consensus clustering per feature; log k, PAC, CHI; flag boundary items.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

select_k <- function(pac_by_k, chi_by_k, max_k) {
  ks <- as.integer(names(pac_by_k))
  if (length(ks) == 1L) return(ks[1L])
  pac <- pac_by_k[as.character(ks)]
  chi <- chi_by_k[as.character(ks)]
  pac[!is.finite(pac)] <- Inf
  chi[!is.finite(chi)] <- -Inf
  pac_sc <- if (length(unique(pac)) > 1L) as.numeric(scale(-pac)) else rep(0, length(pac))
  chi_sc <- if (length(unique(chi)) > 1L) as.numeric(scale(chi)) else rep(0, length(chi))
  score <- pac_sc + chi_sc
  ks[which.max(score)]
}

compute_pac <- function(consensus_mat, k) {
  # Proportion of ambiguous cluster assignments (mid-range consensus)
  cm <- consensus_mat
  if (is.null(cm)) return(NA_real_)
  ambiguous <- cm > 0.1 & cm < 0.9
  mean(ambiguous)
}

compute_chi <- function(dist_mat, clusters) {
  # Calinski-Harabasz via between/within cluster dispersion on distance embedding
  if (length(unique(clusters)) < 2) return(NA_real_)
  n <- length(clusters)
  k <- length(unique(clusters))
  # Use classical MDS embedding of distance matrix for CHI proxy
  if (!requireNamespace("stats", quietly = TRUE)) return(NA_real_)
  embed <- tryCatch(
    cmdscale(as.dist(dist_mat), k = min(5, n - 1)),
    error = function(e) NULL
  )
  if (is.null(embed)) return(NA_real_)
  overall <- colMeans(embed)
  between <- 0
  within <- 0
  for (cl in unique(clusters)) {
    idx <- which(clusters == cl)
    center <- colMeans(embed[idx, , drop = FALSE])
    between <- between + length(idx) * sum((center - overall)^2)
    within <- within + sum(apply(embed[idx, , drop = FALSE], 1, function(r) sum((r - center)^2)))
  }
  if (within == 0) return(NA_real_)
  (between / (k - 1)) / (within / (n - k))
}

run_consensus <- function(dist_mat, cfg) {
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("Install ConsensusClusterPlus via BiocManager")
  }
  reps <- cfg$features$consensus$reps
  p_item <- cfg$features$consensus$p_item
  max_k <- min(
    cfg$features$consensus$max_k,
    nrow(dist_mat) - 1L,
    max(2L, floor(nrow(dist_mat) * p_item) - 1L)
  )
  if (max_k < 2L) {
    stop("Need at least 3 clustering items for consensus; got ", nrow(dist_mat))
  }
  message("[consensus] n_items=", nrow(dist_mat), " maxK=", max_k)

  cc <- ConsensusClusterPlus::ConsensusClusterPlus(
    as.dist(dist_mat),
    maxK = max_k,
    reps = reps,
    pItem = p_item,
    pFeature = 1,
    clusterAlg = "pam",
    distance = "euclidean",
    seed = cfg$features$seed,
    plot = NULL
  )

  ks <- 2:max_k
  pac_by_k <- chi_by_k <- numeric(length(ks))
  names(pac_by_k) <- names(chi_by_k) <- as.character(ks)

  cluster_assignments <- list()
  for (k in ks) {
    res <- cc[[k]]
    clusters <- res$consensusClass
    pac_by_k[as.character(k)] <- compute_pac(res$consensusMatrix, k)
    chi_by_k[as.character(k)] <- compute_chi(dist_mat, clusters)
    cluster_assignments[[as.character(k)]] <- clusters
  }

  best_k <- select_k(pac_by_k, chi_by_k, max_k)
  list(
    cc = cc,
    pac_by_k = pac_by_k,
    chi_by_k = chi_by_k,
    best_k = best_k,
    clusters = cluster_assignments[[as.character(best_k)]]
  )
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "consensus_cluster")

  w_dir <- resolve_path(cfg, cfg$paths$phase_a$wasserstein_dir)
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir)
  ensure_dir(out_dir)

  files <- list.files(w_dir, pattern = "_w1_dist\\.rds$", full.names = TRUE)
  if (!length(files)) stop("Run 04_wasserstein_matrix.R first")

  summary_rows <- list()
  boundary_thr <- cfg$features$consensus$boundary_threshold

  for (f in files) {
    obj <- readRDS(f)
    feat <- obj$feature
    D <- obj$distance
  res <- run_consensus(D, cfg)

    # Item consensus from ConsensusClusterPlus best-k matrix
    cm <- res$cc[[res$best_k]]$consensusMatrix
    item_ids <- rownames(D)
    item_consensus <- setNames(
      vapply(seq_along(item_ids), function(i) {
        same <- which(res$clusters == res$clusters[i])
        mean(cm[i, same])
      }, numeric(1)),
      item_ids
    )
    boundary <- names(item_consensus)[item_consensus < boundary_thr]

    out <- list(
      feature = feat,
      best_k = res$best_k,
      pac_by_k = res$pac_by_k,
      chi_by_k = res$chi_by_k,
      clusters = res$clusters,
      item_consensus = item_consensus,
      boundary_items = boundary
    )
    out_file <- file.path(out_dir, paste0(feat, "_consensus.rds"))
    saveRDS(out, out_file)

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      feature = feat,
      k = res$best_k,
      PAC = res$pac_by_k[as.character(res$best_k)],
      CHI = res$chi_by_k[as.character(res$best_k)],
      n_boundary = length(boundary),
      stringsAsFactors = FALSE
    )
    message(sprintf(
      "[ok] %s: k=%d PAC=%.3f CHI=%.3f boundary=%d -> %s",
      feat, res$best_k,
      res$pac_by_k[as.character(res$best_k)],
      res$chi_by_k[as.character(res$best_k)],
      length(boundary), out_file
    ))
  }

  summary <- do.call(rbind, summary_rows)
  summary_file <- file.path(out_dir, "consensus_summary.csv")
  write.csv(summary, summary_file, row.names = FALSE)
  message("[ok] Summary -> ", summary_file)
  invisible(out_dir)
}

if (sys.nframe() == 0) main()
