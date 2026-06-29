# Learnings

Accumulated lessons for this repo. Newest first.

### [better-approach] Negatives were the wrong instrument — match analysis to biology — 2026-06-29
- What happened: within-compartment distributional tests (Phase A) and cross-modal composition
  regression (Phase B) were both null. Reframing to the analyses the biology actually supports
  turned both into significant positives:
    * Phase A omics: compartment COMPOSITION (epithelial↑ anaplastic, FDR<0.05) + pseudobulk
      PATHWAY enrichment (proliferation↑/TP53↓ in relapse, FDR≤1e-9). Pathway enrichment >>
      single-gene FDR at n=20/group without DESeq2/edgeR.
    * Phase B spatial: classify HISTOLOGY (anaplasia = nuclear atypia, H&E's native signal) from
      Phikon embeddings → tumor-level AUC 0.72, perm p=0.006. Watershed-morphology AUC=0.39 (the
      weak tool); the FM embedding bypasses segmentation.
- Rule: before concluding "no signal," ask whether the *instrument* matches where the biology
  lives. Composition/pathway for omics; histology classification (not composition regression) for
  H&E. Use a pretrained pathology FM (Phikon, ungated, PyTorch — torch is installed, no TF needed).
- See also: 12_composition_analysis.R, 13_pseudobulk_de.R, 14_phase_b_positives.py.

### [gotcha] Demo manifest can clobber the real 260k-tile manifest — 2026-06-29
- What happened: tiles_manifest.json was overwritten by a 12-entry demo manifest; Phase B FM run
  loaded 12 tiles / 1 tumor and merged 0.
- Rule: if Phase B loads suspiciously few tiles, run `00_recover_manifest.py --force` to rebuild
  from spot_signatures + the surviving PNGs.

### [gotcha] Reference cell-type annotations cannot define WT compartments — 2026-06-29
- What happened: `cellassign_celltype_annotation` labels WT tumor cells as
  "Hemangioblasts"/"Trophoblast"/"Macrophage"; `consensus_*` calls ~55k tumor cells
  "Unknown"; SingleR calls 26k "neuron". The old keyword map turned these into
  blastemal/epithelial/stromal — biologically invalid (the whole Phase A premise).
- Rule: define WT compartments from **fetal-kidney signatures** (CM/UB/PV/fibroblast,
  `config/cell_signatures.yaml`) on tumor cells, not reference labels. Gate non-tumor via
  consensus immune/endothelial terms + `is_infercnv_reference`.
- See also: phase1_mechanotypes/02_qc_normalize.R; Yang et al. 2025 Front Immunol.

### [gotcha] Per-cell signature argmax is dominated by high-baseline / long genes — 2026-06-29
- What happened: naive marker-mean argmax over fetal signatures over-called CM (WT1) and
  "neural" (DST/NRXN1 — long genes accumulate intronic reads in snRNA), 52k false neural.
- Rule: z-score each gene across cells before averaging (AddModuleScore analogue), drop
  artifact-prone signatures from the argmax, and require a top-vs-runner-up margin (else NA).
- See also: scaled_signature_scores() in 02_qc_normalize.R.

### [correction] Inference unit is the PATIENT, not the cell (pseudoreplication) — 2026-06-29
- What happened: Phase A permuted/tested at the cell level over ~40k cells → every p=0.001
  (anti-conservative). Cells within a sample are not independent (Squair et al. 2021).
- Rule: permute the histology/relapse label ACROSS SAMPLES (n≈40); report effect sizes
  (Cliff's δ, bootstrap CI) + BH-FDR. Cell-level perm is never valid for condition contrasts.
- See also: sample_perm_p_w1() in 09_distributional_validation.R; tests/test_inference.R.

### [knowledge-gap] WT signal is COMPOSITIONAL, not within-compartment distributional — 2026-06-29
- What happened: with correct labels + patient-level stats, within-compartment program
  distributions show 0/18 FDR-significant on both histology AND relapse axes. But compartment
  *fractions* differ (epithelial↑ anaplastic, stromal↑ favorable, FDR<0.05). Matches Yang 2025.
- Rule: Phase A's result is compartment composition; the distributional-mechanotype framing
  is a (clean, method-robust) negative for WT. Phase B (read composition from H&E) is the
  high-value contribution.
- See also: 12_composition_analysis.R; results/mechanotypes/composition_analysis.csv.

### [gotcha] R needs Matrix namespace loaded for dgCMatrix `[`/`t`/`as.matrix` — 2026-06-29
- What happened: subsetting/transposing a sparse `dgCMatrix` errored ("object of type S4
  is not subsettable" / "not a matrix") because Matrix's S4 methods weren't registered.
- Rule: call `requireNamespace("Matrix")` before sparse ops, or densify the small submatrix
  first. `Matrix::colSums` works namespace-qualified but `[`/`t` need the methods loaded.

### [gotcha] WILMS_DEMO=1 forces synthetic Phase B — 2026-06-28
- What happened: Shell had `WILMS_DEMO=1` from demo tests; real Visium extraction skipped.
- Rule for next time: Use `scripts/run_phase_b.bat` or `WILMS_DEMO=0` before Phase B.
- See also: phase2_histology_ml/utils.py

### [gotcha] Scanpy X is cells × genes — 2026-06-28
- What happened: Gene scoring used row indexing; IndexError on Visium program scores.
- Rule for next time: Use `X[:, gene_idx]` and `X.sum(axis=1)` for AnnData.
- See also: phase2_histology_ml/spatial_utils.py

### [knowledge-gap] ScPCAr API June 2026 — 2026-06-28
- What happened: `download_sample()` is deprecated; computed-files endpoint removed.
  Single-sample pulls use `create_dataset()` → `download_dataset(await_processing=TRUE)`.
  Spatial format is `"spatial"` (not `spaceranger`). `auth_token` is no longer the 2nd
  positional arg — use `SCPCA_AUTH_TOKEN` from `get_auth()`.
- Rule for next time: Full project → `download_project()`; one sample → dataset API;
  Wilms SCPCP000006 snRNA uses `modality=SINGLE_CELL` filter inside `download_project`.
- See also: phase1_mechanotypes/scpca_api.R; https://alexslemonade.github.io/ScPCAr/

### [gotcha] Wasserstein must run on 1-D scores only — 2026-06-28