# Shared utilities for Phase A (mechanotypes)
suppressPackageStartupMessages({
  library(yaml)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }
  normalizePath(".", winslash = "/")
}

REPO_ROOT <- function() {
  # Walk up from script location to find repo root (contains config/paths.yaml)
  here <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  candidates <- c(
    here,
    normalizePath(file.path(here, ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(here, "../.."), winslash = "/", mustWork = FALSE)
  )
  for (c in candidates) {
    if (file.exists(file.path(c, "config", "paths.yaml"))) return(c)
  }
  stop("Could not locate repo root (config/paths.yaml). Run from repo root.")
}

load_config <- function() {
  root <- REPO_ROOT()
  paths <- yaml::read_yaml(file.path(root, "config", "paths.yaml"))
  features <- yaml::read_yaml(file.path(root, "config", "features.yaml"))
  list(root = root, paths = paths, features = features)
}

resolve_path <- function(cfg, rel) {
  file.path(cfg$root, rel)
}

set_seed_logged <- function(seed, label = "global") {
  set.seed(seed)
  message(sprintf("[seed] %s = %d", label, seed))
  invisible(seed)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_provenance <- function(cfg, event) {
  log_file <- resolve_path(cfg, cfg$paths$provenance$access_log)
  ensure_dir(dirname(log_file))
  line <- sprintf("%s | %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), event)
  cat(line, file = log_file, append = file.exists(log_file))
  message("[provenance] ", trimws(line))
}

parse_args_demo <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  sample_id <- sub("^--sample=", "", args[grep("^--sample=", args)])
  if (length(sample_id) == 0) sample_id <- ""
  list(
    demo = "--demo" %in% args || Sys.getenv("WILMS_DEMO", "0") == "1",
    sample_id = sample_id,
    merged = !("--per-sample" %in% args),
    spatial = "--spatial" %in% args
  )
}

get_scpca_token <- function() {
  # ScPCAr stores token in SCPCA_AUTH_TOKEN; we also accept SCPCA_TOKEN alias.
  token <- Sys.getenv("SCPCA_AUTH_TOKEN", unset = "")
  if (!nzchar(token)) token <- Sys.getenv("SCPCA_TOKEN", unset = "")
  token
}

build_gene_lookup <- function(symbols, gene_ids) {
  ok <- !is.na(symbols) & nzchar(symbols)
  split(gene_ids[ok], symbols[ok])
}

resolve_feature_genes <- function(gene_names, lookup) {
  if (!length(gene_names)) return(character(0))
  gene_names <- as.character(unlist(gene_names))
  hits <- lookup[gene_names]
  hits <- hits[!vapply(hits, is.null, logical(1))]
  unique(unlist(hits, use.names = FALSE))
}

make_demo_seurat <- function(cfg) {
  set_seed_logged(cfg$features$seed, "demo_seurat")
  n_cells <- 800
  n_genes <- 200
  states <- rep(cfg$features$cell_states, length.out = n_cells)
  histology <- sample(cfg$features$histology_groups, n_cells, replace = TRUE)
  sample_id <- paste0("DEMO_", sample(1:8, n_cells, replace = TRUE))

  gene_names <- unique(unlist(lapply(cfg$features$features, function(f) {
    c(f$genes_positive, f$genes_negative)
  })))
  gene_names <- unique(c(gene_names, paste0("GENE_", seq_len(max(0, n_genes - length(gene_names))))))
  gene_names <- gene_names[seq_len(n_genes)]

  counts <- matrix(
    rpois(n_cells * n_genes, lambda = 3),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(gene_names, paste0("cell_", seq_len(n_cells)))
  )

  meta <- data.frame(
    cell_id = colnames(counts),
    cell_state = states,
    histology = histology,
    sample_id = sample_id,
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$cell_id

  list(counts = counts, meta = meta, is_demo = TRUE)
}
