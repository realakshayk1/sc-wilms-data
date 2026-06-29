#!/usr/bin/env Rscript
# A-4: MODERATED pseudobulk differential expression (edgeR-QLF + limma-voom).
# Replaces the Welch+BH stopgap in 13_pseudobulk_de.R now that edgeR/limma are
# installed. Both methods borrow variance across genes (empirical-Bayes moderation),
# which recovers single-gene FDR hits that an unmoderated Welch test misses at
# n~20/group. Pseudobulk (sum counts per sample) keeps the SAMPLE as the unit of
# inference (no pseudoreplication). Contrasts: histology (anaplastic vs favorable)
# and relapse (relapse vs no_relapse), overall + per compartment.
.script_dir <- local({
  cmd <- commandArgs(trailingOnly = FALSE)
  path <- sub("--file=", "", cmd[grep("--file=", cmd)])
  if (length(path) == 0) return(normalizePath(".", winslash = "/"))
  dirname(normalizePath(path, winslash = "/"))
})
source(file.path(.script_dir, "utils.R"))
suppressPackageStartupMessages({ library(edgeR); library(limma) })

pseudobulk <- function(counts, groups) {
  requireNamespace("Matrix", quietly = TRUE)
  lev <- sort(unique(groups[!is.na(groups)]))
  mat <- vapply(lev, function(g)
    as.numeric(Matrix::rowSums(counts[, which(groups == g), drop = FALSE])),
    numeric(nrow(counts)))
  colnames(mat) <- lev; rownames(mat) <- rownames(counts); mat
}

# edgeR quasi-likelihood F-test + limma-voom on one pseudobulk matrix / 2-level group
moderated_de <- function(pb, grp, pos, neg, id2sym, scope, contrast) {
  keep_s <- !is.na(grp) & grp %in% c(pos, neg)
  pb <- pb[, keep_s, drop = FALSE]; grp <- factor(grp[keep_s], levels = c(neg, pos))
  if (sum(grp == pos) < 3 || sum(grp == neg) < 3) return(NULL)
  dge <- DGEList(counts = round(pb))
  keep_g <- filterByExpr(dge, group = grp)
  dge <- dge[keep_g, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)                      # TMM
  design <- model.matrix(~grp)
  ids <- rownames(dge); syms <- id2sym[ids]

  # --- edgeR QLF ---
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)                 # the pos-vs-neg term
  tt <- topTags(qlf, n = Inf)$table
  edger <- data.frame(method = "edgeR_QLF", scope = scope, contrast = contrast,
                      gene_id = rownames(tt), symbol = id2sym[rownames(tt)],
                      log2FC = tt$logFC, PValue = tt$PValue, FDR = tt$FDR,
                      stringsAsFactors = FALSE)

  # --- limma-voom ---
  v <- voom(dge, design)
  lf <- eBayes(lmFit(v, design))
  lt <- topTable(lf, coef = 2, number = Inf, sort.by = "P")
  voom <- data.frame(method = "limma_voom", scope = scope, contrast = contrast,
                     gene_id = rownames(lt), symbol = id2sym[rownames(lt)],
                     log2FC = lt$logFC, PValue = lt$P.Value, FDR = lt$adj.P.Val,
                     stringsAsFactors = FALSE)
  list(table = rbind(edger, voom),
       summary = data.frame(
         scope = scope, contrast = contrast, n_pos = sum(grp == pos), n_neg = sum(grp == neg),
         n_genes_tested = nrow(tt),
         edgeR_fdr05 = sum(edger$FDR < 0.05, na.rm = TRUE),
         voom_fdr05 = sum(voom$FDR < 0.05, na.rm = TRUE),
         top_edgeR = paste(head(edger$symbol[edger$FDR < 0.05 & !is.na(edger$symbol)], 8), collapse = ","),
         stringsAsFactors = FALSE))
}

main <- function() {
  cfg <- load_config()
  set_seed_logged(cfg$features$seed, "moderated_de")
  out_dir <- resolve_path(cfg, cfg$paths$phase_a$consensus_dir); ensure_dir(out_dir)
  requireNamespace("Matrix", quietly = TRUE)
  proc <- readRDS(resolve_path(cfg, cfg$paths$phase_a$seurat_rds))
  counts <- proc$counts; meta <- proc$meta
  id2sym <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
  if (!is.null(proc$gene_lookup))
    for (sym in names(proc$gene_lookup))
      for (id in proc$gene_lookup[[sym]]) if (id %in% names(id2sym)) id2sym[id] <- sym

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

  tables <- list(); summaries <- list()
  run <- function(pb, samp_groups, contrast, pos, neg, scope) {
    grp <- samp_groups[colnames(pb)]
    r <- tryCatch(moderated_de(pb, grp, pos, neg, id2sym, scope, contrast),
                  error = function(e) { message("[warn] ", scope, " ", contrast, ": ", conditionMessage(e)); NULL })
    if (is.null(r)) return(invisible(NULL))
    tables[[length(tables) + 1]] <<- r$table
    summaries[[length(summaries) + 1]] <<- r$summary
    message(sprintf("[mod-de] %-10s %-22s edgeR FDR05=%d  voom FDR05=%d",
                    scope, contrast, r$summary$edgeR_fdr05, r$summary$voom_fdr05))
  }
  pb_all <- pseudobulk(counts, samples)
  run(pb_all, samp_hist, "favorable_vs_anaplastic", "anaplastic", "favorable", "overall")
  run(pb_all, samp_rel,  "relapse_vs_norelapse",    "relapse",    "no_relapse", "overall")
  for (comp in c("blastemal", "epithelial", "stromal")) {
    sel <- which(meta$cell_state == comp); if (length(sel) < 100) next
    pb <- pseudobulk(counts[, sel, drop = FALSE], samples[sel])
    run(pb, samp_hist, "favorable_vs_anaplastic", "anaplastic", "favorable", comp)
    run(pb, samp_rel,  "relapse_vs_norelapse",    "relapse",    "no_relapse", comp)
  }
  if (length(tables)) {
    allt <- do.call(rbind, tables)
    write.csv(allt, file.path(out_dir, "moderated_de.csv"), row.names = FALSE)
    write.csv(do.call(rbind, summaries), file.path(out_dir, "moderated_de_summary.csv"), row.names = FALSE)
    message("[ok] Moderated DE -> ", file.path(out_dir, "moderated_de.csv"))
    print(do.call(rbind, summaries)[, c("scope","contrast","n_pos","n_neg","edgeR_fdr05","voom_fdr05")],
          row.names = FALSE)
  }
}

if (sys.nframe() == 0) main()
