#!/usr/bin/env Rscript
# FR-A1 continued: QC, normalize, retain cell-state annotations.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))

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

  state_map <- map_cell_state(meta)
  if (is.null(state_map)) {
    stop(
      "No cell-state column found in colData. Inspect colnames(colData(sce)) and ",
      "add a mapping in 02_qc_normalize.R map_cell_state()."
    )
  }
  meta$cell_state <- state_map$values
  message("[map] cell_state from ", state_map$column)

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
    gene_lookup = NULL,
    is_demo = !is.null(raw$is_demo) && raw$is_demo
  )

  if (
    exists("raw_sce", inherits = FALSE) &&
    "gene_symbol" %in% colnames(SummarizedExperiment::rowData(raw_sce))
  ) {
    processed$gene_lookup <- build_gene_lookup(
      SummarizedExperiment::rowData(raw_sce)$gene_symbol,
      rownames(counts)
    )
    message("[ok] gene_symbol lookup for feature scoring")
  }

  saveRDS(processed, out_rds)
  message(sprintf("[ok] QC/normalize: %d cells retained -> %s", ncol(counts), out_rds))
  invisible(out_rds)
}

if (sys.nframe() == 0) main()
