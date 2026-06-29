#!/usr/bin/env Rscript
# A-1: Full MSigDB Hallmark GSEA (fgsea, preranked) — replaces the 4 hand-curated
# gene sets in 13_pseudobulk_de.R with all 50 Hallmark pathways. Genes are ranked
# by the limma-voom MODERATED t-statistic of the pseudobulk contrast (the standard,
# variance-stable preranked GSEA input), then fgsea computes enrichment + BH-FDR
# across pathways. Contrasts: histology (anaplastic vs favorable) and relapse,
# overall + per compartment.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages({ library(edgeR); library(limma); library(fgsea); library(msigdbr) })

pseudobulk <- function(counts, groups) {
  requireNamespace("Matrix", quietly = TRUE)
  lev <- sort(unique(groups[!is.na(groups)]))
  mat <- vapply(lev, function(g)
    as.numeric(Matrix::rowSums(counts[, which(groups == g), drop = FALSE])),
    numeric(nrow(counts)))
  colnames(mat) <- lev; rownames(mat) <- rownames(counts); mat
}

# limma-voom moderated t for one 2-level pseudobulk contrast -> named vector by symbol
ranking_stat <- function(pb, grp, pos, neg, id2sym) {
  keep_s <- !is.na(grp) & grp %in% c(pos, neg)
  pb <- pb[, keep_s, drop = FALSE]; grp <- factor(grp[keep_s], levels = c(neg, pos))
  if (sum(grp == pos) < 3 || sum(grp == neg) < 3) return(NULL)
  dge <- DGEList(counts = round(pb))
  dge <- dge[filterByExpr(dge, group = grp), , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  design <- model.matrix(~grp)
  v <- voom(dge, design)
  lf <- eBayes(lmFit(v, design))
  tt <- topTable(lf, coef = 2, number = Inf, sort.by = "none")
  sym <- id2sym[rownames(tt)]
  ok <- !is.na(sym) & is.finite(tt$t)
  stat <- tt$t[ok]; sym <- sym[ok]
  # collapse duplicate symbols to the most extreme t (keeps the strongest signal)
  o <- order(abs(stat), decreasing = TRUE)
  stat <- stat[o]; sym <- sym[o]
  stat <- stat[!duplicated(sym)]; names(stat) <- sym[!duplicated(sym)]
  stat[is.finite(stat)]
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "hallmark_gsea")
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir); ensure_dir(out_dir)
  requireNamespace("Matrix", quietly = TRUE)
  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts; meta <- proc$meta
  id2sym <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
  if (!is.null(proc$gene_lookup))
    for (sym in names(proc$gene_lookup))
      for (id in proc$gene_lookup[[sym]]) if (id %in% names(id2sym)) id2sym[id] <- sym

  # Hallmark pathways (50), symbol-based
  hm <- msigdbr(species = "Homo sapiens", collection = "H")
  pathways <- split(hm$gene_symbol, hm$gs_name)
  message(sprintf("[msigdbr] %d Hallmark pathways", length(pathways)))

  hist <- tolower(as.character(meta$histology))
  clin <- resolve_path(cfg, file.path(cfg$paths$dirs$raw, "SCPCP000006_clinical_metadata.tsv"))
  relapse <- rep(NA_character_, nrow(meta))
  if (file.exists(clin)) {
    md <- read.delim(clin, check.names = FALSE)
    u <- unique(md[, c("scpca_sample_id", "relapse_status")])
    rmap <- setNames(as.character(u$relapse_status), u$scpca_sample_id)
    rs <- rmap[as.character(meta$sample_id)]
    relapse <- ifelse(rs == "Yes", "relapse", ifelse(rs == "No", "no_relapse", NA))
  }
  samples <- as.character(meta$sample_id)
  samp_hist <- tapply(hist, samples, function(x) x[!is.na(x)][1])
  samp_rel  <- tapply(relapse, samples, function(x) x[!is.na(x)][1])

  all_res <- list()
  run_gsea <- function(pb, samp_groups, contrast, pos, neg, scope) {
    grp <- samp_groups[colnames(pb)]
    stat <- tryCatch(ranking_stat(pb, grp, pos, neg, id2sym), error = function(e) NULL)
    if (is.null(stat) || length(stat) < 50) return(invisible(NULL))
    set.seed(cfg$features$seed)
    fg <- fgsea(pathways = pathways, stats = stat, minSize = 10, maxSize = 500, eps = 0)
    fg <- as.data.frame(fg)
    fg$scope <- scope; fg$contrast <- contrast
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(head(x, 12), collapse = ","), character(1))
    all_res[[length(all_res) + 1]] <<- fg
    nsig <- sum(fg$padj < 0.05, na.rm = TRUE)
    message(sprintf("[gsea] %-10s %-22s pathways FDR<0.05 = %d", scope, contrast, nsig))
    if (nsig > 0) {
      top <- fg[order(fg$padj), ][seq_len(min(5, nsig)), c("pathway", "NES", "padj")]
      for (i in seq_len(nrow(top)))
        message(sprintf("        %-32s NES=%+.2f padj=%.2e", top$pathway[i], top$NES[i], top$padj[i]))
    }
  }
  pb_all <- pseudobulk(counts, samples)
  run_gsea(pb_all, samp_hist, "favorable_vs_anaplastic", "anaplastic", "favorable", "overall")
  run_gsea(pb_all, samp_rel,  "relapse_vs_norelapse",    "relapse",    "no_relapse", "overall")
  for (comp in c("blastemal", "epithelial", "stromal")) {
    sel <- which(meta$cell_state == comp); if (length(sel) < 100) next
    pb <- pseudobulk(counts[, sel, drop = FALSE], samples[sel])
    run_gsea(pb, samp_hist, "favorable_vs_anaplastic", "anaplastic", "favorable", comp)
    run_gsea(pb, samp_rel,  "relapse_vs_norelapse",    "relapse",    "no_relapse", comp)
  }
  if (length(all_res)) {
    res <- do.call(rbind, all_res)
    cols <- c("scope", "contrast", "pathway", "NES", "ES", "pval", "padj", "size", "leadingEdge")
    res <- res[, cols]
    write.csv(res, file.path(out_dir, "hallmark_gsea.csv"), row.names = FALSE)
    message("[ok] Hallmark GSEA -> ", file.path(out_dir, "hallmark_gsea.csv"))
    message(sprintf("[ok] total pathway-contrasts FDR<0.05 = %d", sum(res$padj < 0.05, na.rm = TRUE)))
  }
}

if (sys.nframe() == 0) main()
