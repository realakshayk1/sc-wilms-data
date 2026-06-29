#!/usr/bin/env Rscript
# FR-A1 continued: QC, normalize, retain cell-state annotations.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
source(file.path(.script_dir, "03_compute_scores.R"), local = TRUE)  # score_feature()

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Option C: assign each cell to the compartment whose canonical marker panel scores
# highest, after gating out confidently non-neoplastic microenvironment cells.
assign_compartments_markers <- function(counts, gene_lookup, meta, mcfg) {
  panels <- mcfg$assign_panels
  states <- names(panels)
  # pos-only marker score per panel (same normalization as the test programs)
  score_mat <- vapply(states, function(s) {
    score_feature(counts, panels[[s]], character(0), gene_lookup)
  }, numeric(ncol(counts)))
  colnames(score_mat) <- states

  top <- max.col(score_mat, ties.method = "first")
  sorted <- t(apply(score_mat, 1, function(r) sort(r, decreasing = TRUE)))
  margin <- sorted[, 1] - sorted[, 2]
  compartment <- states[top]
  compartment[margin <= (mcfg$tie_min_margin %||% 0)] <- NA_character_

  # Gate out non-tumor microenvironment cells
  drop <- rep(FALSE, ncol(counts))
  terms <- tolower(unlist(mcfg$nontumor_consensus_terms %||% character(0)))
  if (length(terms) && "consensus_celltype_annotation" %in% colnames(meta)) {
    cons <- tolower(as.character(meta$consensus_celltype_annotation))
    for (t in terms) drop <- drop | grepl(t, cons, fixed = TRUE)
  }
  if (isTRUE(mcfg$exclude_infercnv_reference) && "is_infercnv_reference" %in% colnames(meta)) {
    drop <- drop | (meta$is_infercnv_reference %in% TRUE)
  }
  compartment[drop] <- NA_character_
  list(values = compartment, n_dropped_nontumor = sum(drop), score_mat = score_mat)
}

# Gene-wise z-scored module scores (AddModuleScore analogue, no Seurat dependency):
# normalize -> log1p -> z-score each gene across cells -> mean within signature.
# This stops high-baseline genes (WT1, DST) from dominating the argmax.
scaled_signature_scores <- function(counts, gene_lookup, sig_list) {
  requireNamespace("Matrix", quietly = TRUE)
  resolved <- lapply(sig_list, function(g) {
    rg <- if (!is.null(gene_lookup)) resolve_feature_genes(g, gene_lookup) else as.character(unlist(g))
    intersect(rg, rownames(counts))
  })
  allg <- unique(unlist(resolved))
  lib <- as.numeric(Matrix::colSums(counts)); lib[lib == 0] <- 1; med <- median(lib)
  sub <- as.matrix(counts[allg, , drop = FALSE])          # genes x cells, dense (~60 genes)
  norm <- log1p(sweep(sub, 2, lib, "/") * med)            # CP-median, log1p (base R)
  mu <- rowMeans(norm); sds <- apply(norm, 1, stats::sd); sds[sds == 0 | is.na(sds)] <- 1
  z <- (norm - mu) / sds
  out <- vapply(names(sig_list), function(s) {
    gs <- resolved[[s]]
    if (!length(gs)) return(rep(0, ncol(counts)))
    colMeans(z[gs, , drop = FALSE])
  }, numeric(ncol(counts)))
  colnames(out) <- names(sig_list)
  out
}

# A1: assign tumor cells to fetal-kidney subgroups (CM/UB/PV/fibroblast/neural),
# split CM into blastemal/epithelial, collapse to the triphasic axis, gate non-tumor.
assign_compartments_fetal <- function(counts, gene_lookup, meta, scfg) {
  sigs <- scfg$fetal_signatures
  sig_names <- names(sigs)
  all_panels <- c(sigs, list(.cm_blast = scfg$cm_split$blastemal, .cm_epi = scfg$cm_split$epithelial))
  sm <- scaled_signature_scores(counts, gene_lookup, all_panels)
  score_mat <- sm[, sig_names, drop = FALSE]

  subgroup <- sig_names[max.col(score_mat, ties.method = "first")]
  # confidence margin: NA out cells whose top signature barely beats the runner-up
  top2 <- apply(score_mat, 1, function(r) { s <- sort(r, decreasing = TRUE); s[1] - s[2] })
  subgroup[top2 < (scfg$assign_min_margin %||% 0)] <- NA_character_
  cm <- !is.na(subgroup) & subgroup == "CM"
  subgroup[cm] <- ifelse(sm[cm, ".cm_blast"] >= sm[cm, ".cm_epi"], "CM_blastemal", "CM_epithelial")

  map <- scfg$subgroup_to_compartment
  compartment <- vapply(subgroup, function(s) {
    if (is.na(s)) return(NA_character_)
    v <- map[[s]]
    if (is.null(v) || (length(v) == 1 && is.na(v)) || identical(as.character(v), "NA")) {
      NA_character_
    } else as.character(v)
  }, character(1), USE.NAMES = FALSE)

  drop <- rep(FALSE, ncol(counts))
  terms <- tolower(unlist(scfg$nontumor_consensus_terms %||% character(0)))
  if (length(terms) && "consensus_celltype_annotation" %in% colnames(meta)) {
    cons <- tolower(as.character(meta$consensus_celltype_annotation))
    for (t in terms) drop <- drop | grepl(t, cons, fixed = TRUE)
  }
  if (isTRUE(scfg$exclude_infercnv_reference) && "is_infercnv_reference" %in% colnames(meta)) {
    drop <- drop | (meta$is_infercnv_reference %in% TRUE)
  }
  compartment[drop] <- NA_character_
  subgroup[drop] <- NA_character_
  list(values = compartment, subgroup = subgroup, n_dropped_nontumor = sum(drop))
}

normalize_counts <- function(counts) {
  lib_size <- colSums(counts)
  lib_size[lib_size == 0] <- 1
  log1p(t(t(counts) / lib_size * median(lib_size)))
}

map_cell_state <- function(meta) {
  # Wilms SCPCP000006: cellassign has compartment-relevant labels on this project.
  if ("cellassign_celltype_annotation" %in% colnames(meta)) {
    raw <- as.character(meta$cellassign_celltype_annotation)
    mapped <- rep(NA_character_, length(raw))
    mapped[grepl("Kidney progenitor|Hemangioblast|Trophoblast", raw, ignore.case = TRUE)] <- "blastemal"
    mapped[grepl("Podocyte|Juxtaglomerular", raw, ignore.case = TRUE)] <- "epithelial"
    mapped[grepl("Macrophage|Pericyte|Endothelial", raw, ignore.case = TRUE)] <- "stromal"
    return(list(column = "cellassign_celltype_annotation", values = mapped))
  }
  candidates <- c(
    "cell_state", "consensus_celltype_annotation",
    "celltype_annotation", "singler_celltype_annotation"
  )
  src <- intersect(candidates, colnames(meta))[1]
  if (is.na(src)) return(NULL)

  raw <- as.character(meta[[src]])
  mapped <- rep(NA_character_, length(raw))
  mapped[grepl("blastem", raw, ignore.case = TRUE)] <- "blastemal"
  mapped[grepl("epith", raw, ignore.case = TRUE)] <- "epithelial"
  mapped[grepl("strom|mesench", raw, ignore.case = TRUE)] <- "stromal"
  list(column = src, values = mapped)
}

map_histology <- function(meta) {
  # Wilms SCPCP000006: subdiagnosis is Favorable vs Anaplastic (sample-level).
  candidates <- c("histology", "subdiagnosis", "diagnosis", "molecular_characteristics")
  src <- intersect(candidates, colnames(meta))[1]
  if (is.na(src)) return(NULL)

  raw <- as.character(meta[[src]])
  mapped <- rep(NA_character_, length(raw))
  mapped[grepl("anaplas", raw, ignore.case = TRUE)] <- "anaplastic"
  mapped[grepl("favor", raw, ignore.case = TRUE)] <- "favorable"
  list(column = src, values = mapped)
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "qc_normalize")

  in_rds <- resolve_path(cfg, cfg$paths$phase_a$merged_sce_rds)
  out_rds <- resolve_path(cfg, cfg$paths$phase_a$seurat_rds)
  if (!file.exists(in_rds)) {
    stop("Run ingest_manual_scpca.R first. Missing: ", in_rds)
  }

  raw <- readRDS(in_rds)

  if (!is.null(raw$is_demo) && raw$is_demo) {
    counts <- raw$counts
    meta <- raw$meta
  } else if (
    requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    inherits(raw, "SingleCellExperiment")
  ) {
    sce <- raw
    counts <- SummarizedExperiment::assay(sce, "counts")
    if (is.null(counts)) counts <- SummarizedExperiment::assay(sce, 1)
    meta <- as.data.frame(SummarizedExperiment::colData(sce))
    rownames(meta) <- colnames(sce)
    raw_sce <- sce
  } else if (requireNamespace("Seurat", quietly = TRUE) && inherits(raw, "Seurat")) {
    counts <- Seurat::GetAssayData(raw, layer = "counts")
    meta <- raw@meta.data
  } else {
    stop("Unsupported input — expected demo list, SingleCellExperiment, or Seurat")
  }

  if (inherits(counts, "Matrix") || inherits(counts, "dgCMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Install Matrix: install.packages('Matrix')")
    }
    n_genes <- Matrix::colSums(counts > 0)
  } else {
    n_genes <- colSums(counts > 0)
  }
  keep <- n_genes >= 200
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  # Build gene_lookup (symbol -> id) early; marker assignment needs it.
  gene_lookup <- NULL
  if (
    exists("raw_sce", inherits = FALSE) &&
    "gene_symbol" %in% colnames(SummarizedExperiment::rowData(raw_sce))
  ) {
    gene_lookup <- build_gene_lookup(
      SummarizedExperiment::rowData(raw_sce)$gene_symbol,
      rownames(counts)
    )
    message("[ok] gene_symbol lookup for feature scoring")
  }

  # Compartment assignment: A1 fetal-kidney signatures > Option C markers > keyword map.
  sig_path <- file.path(cfg$root, "config", "cell_signatures.yaml")
  scfg <- if (file.exists(sig_path)) yaml::read_yaml(sig_path) else NULL
  mapping_path <- file.path(cfg$root, "config", "label_mapping.yaml")
  mcfg <- if (file.exists(mapping_path)) yaml::read_yaml(mapping_path)$marker_assignment else NULL
  is_demo_obj <- !is.null(raw$is_demo) && raw$is_demo

  if (!is_demo_obj && !is.null(gene_lookup) && !is.null(scfg) &&
      identical(scfg$method, "fetal_signature")) {
    asg <- assign_compartments_fetal(counts, gene_lookup, meta, scfg)
    meta$cell_state <- asg$values
    meta$subgroup <- asg$subgroup
    message(sprintf("[map] cell_state from FETAL-KIDNEY SIGNATURES (A1); dropped %d non-tumor cells",
                    asg$n_dropped_nontumor))
    print(table(subgroup = meta$subgroup, useNA = "ifany"))
    print(table(compartment = meta$cell_state, useNA = "ifany"))
  } else if (isTRUE(mcfg$enabled) && !is_demo_obj && !is.null(gene_lookup)) {
    asg <- assign_compartments_markers(counts, gene_lookup, meta, mcfg)
    meta$cell_state <- asg$values
    message(sprintf("[map] cell_state from CANONICAL MARKER PANELS (Option C); dropped %d non-tumor cells",
                    asg$n_dropped_nontumor))
    print(table(meta$cell_state, useNA = "ifany"))
  } else {
    state_map <- map_cell_state(meta)
    if (is.null(state_map)) {
      stop(
        "No cell-state column found in colData. Inspect colnames(colData(sce)) and ",
        "add a mapping in 02_qc_normalize.R map_cell_state()."
      )
    }
    meta$cell_state <- state_map$values
    message("[map] cell_state from ", state_map$column)
  }

  hist_map <- map_histology(meta)
  if (is.null(hist_map) && "histology" %in% colnames(meta) && any(!is.na(meta$histology))) {
    message("[map] histology already on object")
  } else if (is.null(hist_map)) {
    warning("No histology column mapped — add sample-level join in ingest_manual_scpca.R.")
    meta$histology <- NA_character_
  } else {
    meta$histology <- hist_map$values
    message("[map] histology from ", hist_map$column)
  }

  n_missing_state <- sum(is.na(meta$cell_state))
  if (n_missing_state > 0) {
    message("[flag] ", n_missing_state, " cells lack blastemal/epithelial/stromal mapping")
  }

  # Keep only cells usable for mechanotyping (state + histology); stay sparse — no dense matrix.
  use <- !is.na(meta$cell_state) & !is.na(meta$histology)
  counts <- counts[, use, drop = FALSE]
  meta <- meta[use, , drop = FALSE]
  message(sprintf("[filter] %d cells with cell_state + histology", ncol(counts)))

  processed <- list(
    counts = counts,
    normalized = NULL,
    meta = meta,
    gene_lookup = gene_lookup,
    is_demo = !is.null(raw$is_demo) && raw$is_demo
  )

  saveRDS(processed, out_rds)
  message(sprintf("[ok] QC/normalize: %d cells retained -> %s", ncol(counts), out_rds))
  invisible(out_rds)
}

if (sys.nframe() == 0) main()
