# Learnings

Accumulated lessons for this repo. Newest first.

### [gotcha] Phase C: StarDist counts undercount cells/spot — 2026-07-06
- What happened: nuclei-per-spot from the StarDist hires-tile run come out ~1 (median),
  far below the true ~10-15 cells in a 55 µm Visium spot. Using it as
  `N_cells_per_location` gives implausibly sparse ABMs.
- Rule for next time: StarDist here is for morphology, not absolute counts (resolution
  ceiling). Default `config/phase_c.yaml density.source: fixed` (literature prior);
  reserve `stardist` for relative variation once WSI-resolution segmentation exists.
- See also: phase3_abm/01_spot_density.py, 02_place_agents.py

### [knowledge] Phase C: Visium px→µm affine — 2026-07-06
- Full-res pixel spot centres → microns via `spot_diameter_fullres` from
  scalefactors_json (µm/px = 55 / spot_diameter_fullres); isotropic, recentre to origin.
- See also: phase3_abm/abm_utils.py:um_per_pixel, load_spot_coords_um

### [better-approach] StarDist fixes the watershed morphology failure; ensemble adds nothing — 2026-06-29
- What happened: watershed nuclear-morphology -> histology AUC was 0.39 (worse than chance).
  StarDist '2D_versatile_he' at prob_thresh=0.4 (default 0.69 detects ~0 nuclei on 96px hires
  tiles) gives morphology AUC 0.687 (perm p=0.021) — segmentation, not the hypothesis, was the
  bottleneck. BUT ensemble (morphology + phikon-v2 embedding) = 0.719 vs embedding-only 0.714,
  paired DeLong p=0.57: no orthogonal gain — both read the same nuclear atypia, same ~0.73 ceiling.
- Rule: for H&E nuclear morphology on these tiles use StarDist with prob_thresh~0.4 + median
  imputation (few-nuclei tumors -> NaN feats). Don't expect morphology+embedding ensembling to
  beat embedding alone here. Median nuclei/tumor stays low (~14) — Visium-hires resolution limit.
- See also: 16_stardist_morphology.py; results/classifier/stardist_morphology.json.

### [gotcha] StarDist from_pretrained needs a Windows symlink (admin) — use a junction — 2026-06-29
- What happened: `StarDist2D.from_pretrained("2D_versatile_he")` downloaded fine but crashed at
  `OSError [WinError 1314] A required privilege is not held` — csbdeep publishes the extracted
  model via `os.symlink`, which needs admin/Developer Mode on Windows.
- Rule: create a directory JUNCTION (no admin needed) once, then from_pretrained skips the
  symlink: `New-Item -ItemType Junction -Path <...>\StarDist2D\2D_versatile_he\2D_versatile_he
  -Target <...>\2D_versatile_he_extracted`. Persists across runs. (TF backend itself is fine.)

### [knowledge-gap] Phase B histology AUC is at a ~0.73 RESOLUTION ceiling — 2026-06-29
- What happened: scaling Phikon spots 60->200, upgrading v1(ViT-B,768) -> v2(ViT-L,1024), and
  swapping mean-pool -> gated attention-MIL moved tumor-level histology AUC only 0.724 -> 0.748.
  Paired DeLong MIL vs mean-pool: delta +0.014, p=0.83 (indistinguishable). Both still beat
  chance (perm p~0.003-0.006), so the signal is real but capped.
- Rule: ~0.73-0.75 looks like a genuine ceiling from Visium-HIRES tiles (not true WSI). Don't
  burn more effort on encoder/pooling tweaks for Phase B histology; the lever to lift it is
  higher-resolution input (gated FMs on real WSI, or XMAG on 5x), which needs external data.
  Report scale/MIL as a clean null-improvement, not a win.
- See also: 15_phase_b_mil.py; results/classifier/phase_b_mil_phikon-v2.json.

### [gotcha] Background jobs die on Claude session teardown; checkpoint long jobs — 2026-06-29
- What happened: the phikon-v2 embedding (8002 ViT-L tiles, ~40min CPU) was killed ~5x by
  process teardowns; the original all-or-nothing cache lost everything each time.
- Rule: any multi-10-min background compute in this repo must checkpoint incrementally.
  Pattern that worked: embed in batches, write APPEND-ONLY chunk parquets to a *_partial/ dir
  every ~8 batches (atomic tmp+rename), resume by unioning chunk spot_ids and skipping done.
  O(chunk) I/O, not O(all-so-far). On completion concat -> final cache, delete chunks.

### [gotcha] ScPCA vital_status is "Expired"/"Alive", and small-n logistic separates — 2026-06-29
- What happened: a regex `grepl("decea|dead")` reported 0 OS events because ScPCA encodes
  death as **"Expired"** (8 of 43 samples). Also, covariate-adjusted ordinary logistic for
  relapse (~10 events / 4 params) hit complete separation — ORs collapsed to exactly 1.0
  with Wald p=1, or exploded to CI [2e-321, 5e27]. Neither is real signal.
- Rule: map vital_status with `grepl("decea|dead|expir")`. For any small-n / rare-event
  logistic in this repo use **Firth penalized regression** (`logistf`, installed) — it gives
  finite ORs + profile-likelihood CIs under separation. Report univariate AND adjusted; at
  EPV~2.5 the adjusted CI is wide by design (say so, don't dress it up).
- Net Phase A prognostic (16_prognostic_association.R): proliferation_score -> relapse is a
  NOMINAL positive (Firth univariate OR~4/SD p=0.013; Fisher OR 7.6 p=0.017) that triangulates
  the GSEA/DE relapse proliferation signal, but does not survive BH-FDR or covariate adjustment.
- See also: 15_hallmark_gsea.R (E2F/G2M/MYC up in relapse, padj~1e-29); 14_moderated_de.R.

### [gotcha] R is installed but NOT on PATH — `which Rscript` lies — 2026-06-29
- What happened: an agent ran `which Rscript` / `ls "/c/Program Files/R"`, got nothing, and
  wrongly concluded "R is not installed → Phase A / edgeR blocked." R 4.6.1 is at
  `C:\Program Files (x86)\R\R-4.6.1\bin\x64\Rscript.exe` (note the **(x86)** dir) and produced
  every committed Phase A result this session.
- Rule: never decide R is absent from PATH alone. Call it by full path or via
  `scripts/rscript.bat` (which already encodes the fallback). From bash:
  `RSCRIPT="/c/Program Files (x86)/R/R-4.6.1/bin/x64/Rscript.exe"`. BiocManager is present;
  limma/edgeR install via `BiocManager::install(..., update=FALSE, ask=FALSE)`.

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