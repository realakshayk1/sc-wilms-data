#!/usr/bin/env Rscript
# Phase A omics positive: PSEUDOBULK differential expression + gene-set enrichment.
# Aggregate snRNA counts per sample (and per sample x compartment), logCPM, then a
# vectorized Welch t-test favorable-vs-anaplastic / relapse-vs-not + BH-FDR, and a
# rank-based enrichment of anaplasia-relevant pathways (proliferation/TP53/MYC).
# (DESeq2/edgeR/limma unavailable here; pseudobulk + Welch across the 41 patients is a
# valid, conventional sample-level DE. Enrichment recovers signal individual-gene FDR misses.)
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
`%||%` <- function(a, b) if (!is.null(a)) a else b

GENE_SETS <- list(
  proliferation = c("MKI67","TOP2A","PCNA","CDK1","CCNB1","CCNB2","CCNA2","AURKA","AURKB",
                    "BUB1","BUB1B","CENPA","CENPE","KIF11","FOXM1","PLK1","CDC20","UBE2C",
                    "BIRC5","TYMS","RRM2","CDC25C","KIF23","NUSAP1"),
  tp53_targets  = c("CDKN1A","MDM2","BAX","GADD45A","BBC3","SFN","TP53I3","RRM2B","ZMAT3",
                    "FAS","SESN1","TNFRSF10B","PMAIP1"),
  myc_mycn      = c("MYCN","MYC","MAX","ODC1","NPM1","NCL","NOP56","NOP58","FBL","MRTO4"),
  nephron_prog  = c("SIX2","CITED1","EYA1","MEOX1","SALL1","GDNF","DAPL1","CRABP2"))

pseudobulk <- function(counts, groups) {
  requireNamespace("Matrix", quietly = TRUE)
  lev <- sort(unique(groups[!is.na(groups)]))
  mat <- vapply(lev, function(g) as.numeric(Matrix::rowSums(counts[, which(groups == g), drop = FALSE])),
                numeric(nrow(counts)))
  colnames(mat) <- lev; rownames(mat) <- rownames(counts); mat
}

logcpm <- function(mat) {
  lib <- colSums(mat); lib[lib == 0] <- 1
  log2(t(t(mat) / lib * 1e6) + 1)
}

.rowVars <- function(M) {
  n <- ncol(M); mu <- rowMeans(M)
  (rowSums(M * M) - n * mu * mu) / (n - 1)
}

de_welch <- function(lcpm, grp, pos, neg, id2sym) {
  keep <- rowMeans(lcpm > 0) >= 0.2
  lcpm <- lcpm[keep, , drop = FALSE]
  pi <- which(grp == pos); ni <- which(grp == neg)
  if (length(pi) < 4 || length(ni) < 4) return(NULL)
  P <- lcpm[, pi, drop = FALSE]; N <- lcpm[, ni, drop = FALSE]
  m1 <- rowMeans(P); m2 <- rowMeans(N); v1 <- .rowVars(P); v2 <- .rowVars(N)
  n1 <- length(pi); n2 <- length(ni)
  se <- sqrt(v1 / n1 + v2 / n2); se[se == 0] <- NA
  t <- (m1 - m2) / se
  df <- (v1 / n1 + v2 / n2)^2 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
  p <- 2 * stats::pt(-abs(t), df)
  res <- data.frame(gene_id = rownames(lcpm), symbol = id2sym[rownames(lcpm)],
                    log2FC = m1 - m2, t = t, p = p, stringsAsFactors = FALSE)
  res$FDR <- p.adjust(res$p, "BH")
  res[order(res$p), ]
}

enrich <- function(res, scope, contrast) {
  # rank by signed significance; Mann-Whitney of set genes vs background
  res <- res[!is.na(res$symbol) & !is.na(res$t), ]
  score <- res$t
  rows <- list()
  for (gs in names(GENE_SETS)) {
    idx <- which(toupper(res$symbol) %in% GENE_SETS[[gs]])
    if (length(idx) < 4) next
    inset <- score[idx]; bg <- score[-idx]
    p <- tryCatch(stats::wilcox.test(inset, bg)$p.value, error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.frame(
      scope = scope, contrast = contrast, gene_set = gs, n_in_set = length(idx),
      median_t_set = median(inset), median_t_bg = median(bg),
      direction = ifelse(median(inset) > median(bg), "up_in_pos", "down_in_pos"),
      enrich_p = p, stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "pseudobulk_de")
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir); ensure_dir(out_dir)
  requireNamespace("Matrix", quietly = TRUE)  # register S4 dimnames/[/rowSums for dgCMatrix
  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts; meta <- proc$meta
  id2sym <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
  if (!is.null(proc$gene_lookup)) {
    for (sym in names(proc$gene_lookup)) for (id in proc$gene_lookup[[sym]]) if (id %in% names(id2sym)) id2sym[id] <- sym
  }
  hist <- tolower(as.character(meta$histology))
  clin <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_clinical_metadata.tsv"))
  relapse <- rep(NA_character_, nrow(meta))
  if (file.exists(clin)) {
    md <- read.delim(clin, check.names = FALSE); u <- unique(md[, c("scpca_sample_id", "relapse_status")])
    rmap <- setNames(as.character(u$relapse_status), u$scpca_sample_id); rs <- rmap[as.character(meta$sample_id)]
    relapse <- ifelse(rs == "Yes", "relapse", ifelse(rs == "No", "no_relapse", NA))
  }
  samples <- as.character(meta$sample_id)
  samp_hist <- tapply(hist, samples, function(x) x[!is.na(x)][1])
  samp_rel  <- tapply(relapse, samples, function(x) x[!is.na(x)][1])

  de_sum <- list(); en_sum <- list()
  run_de <- function(pb, samp_groups, contrast, pos, neg, scope) {
    lcpm <- logcpm(pb); grp <- samp_groups[colnames(lcpm)]
    res <- de_welch(lcpm, grp, pos, neg, id2sym); if (is.null(res)) return(invisible(NULL))
    write.csv(res, file.path(out_dir, sprintf("de_%s_%s.csv", scope, contrast)), row.names = FALSE)
    nsig <- sum(res$FDR < 0.05, na.rm = TRUE)
    en <- enrich(res, scope, contrast)
    de_sum[[length(de_sum) + 1]] <<- data.frame(scope = scope, contrast = contrast, n_genes = nrow(res),
      n_raw_p05 = sum(res$p < 0.05, na.rm = TRUE), n_fdr05 = nsig,
      top = paste(head(res$symbol[res$FDR < 0.05 & !is.na(res$symbol)], 10), collapse = ","),
      stringsAsFactors = FALSE)
    if (!is.null(en)) {
      en_sum[[length(en_sum) + 1]] <<- en
      for (i in seq_len(nrow(en))) message(sprintf("[enrich] %-10s %-22s %-14s %s p=%.2e",
        scope, contrast, en$gene_set[i], en$direction[i], en$enrich_p[i]))
    }
    message(sprintf("[de] %-10s %-22s FDR<0.05=%d", scope, contrast, nsig))
  }
  pb_all <- pseudobulk(counts, samples)
  run_de(pb_all, samp_hist, "favorable_vs_anaplastic", "anaplastic", "favorable", "overall")
  run_de(pb_all, samp_rel,  "relapse_vs_norelapse",    "relapse",    "no_relapse", "overall")
  for (comp in c("blastemal", "epithelial", "stromal")) {
    sel <- which(meta$cell_state == comp); if (length(sel) < 100) next
    pb <- pseudobulk(counts[, sel, drop = FALSE], samples[sel])
    run_de(pb, samp_hist, "favorable_vs_anaplastic", "anaplastic", "favorable", comp)
    run_de(pb, samp_rel,  "relapse_vs_norelapse",    "relapse",    "no_relapse", comp)
  }
  if (length(de_sum)) write.csv(do.call(rbind, de_sum), file.path(out_dir, "de_summary.csv"), row.names = FALSE)
  if (length(en_sum)) {
    es <- do.call(rbind, en_sum); es$enrich_FDR <- p.adjust(es$enrich_p, "BH")
    write.csv(es, file.path(out_dir, "de_enrichment.csv"), row.names = FALSE)
    message("[ok] Enrichment -> ", file.path(out_dir, "de_enrichment.csv"))
  }
}

if (sys.nframe() == 0) main()
