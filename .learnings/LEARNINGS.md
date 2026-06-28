# Learnings

Accumulated lessons for this repo. Newest first.

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