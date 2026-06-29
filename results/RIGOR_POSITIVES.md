# Rigor upside-lever results (branch: `rigor-audit-positives`)

Execution of the Tier-1/Tier-2 upside levers + the ABM payoff. Every claim below is
held-out and effect-size/CI reported; nulls are reported as nulls. Data: ScPCA
`SCPCP000006` (41 Visium tumors / ~40 snRNA samples). Inference unit = patient.

## Phase A (omics) — 3/3 positives, mutually consistent

| Lever | Method | Result | Acceptance |
|---|---|---|---|
| **A-4** moderated DE | edgeR-QLF + limma-voom on pseudobulk | **130 genes FDR<0.05** (histology), 39 (relapse); top: NOTCH2, PODXL, PTPRO, DACT3 | ✅ real single-gene FDR hits (Welch stopgap found ~0) |
| **A-1** Hallmark GSEA | fgsea preranked (limma-voom *t*), 50 Hallmark sets | **166 sig pathway-contrasts**; relapse axis E2F_TARGETS **padj 9e-29**, G2M, MYC — replicated in epithelial (4e-28) and stromal (6e-30) | ✅ ≥5 pathways FDR<0.05 |
| **A-3** prognostic | Firth logistic + Fisher, bootstrap/profile CI | proliferation→relapse **OR 3.97/SD, p=0.013** (Fisher OR 7.6, p=0.017) | ⚠️ nominal only |

**Convergent conclusion:** proliferation / cell-cycle (E2F/G2M/MYC) is the relapse-associated
axis in Wilms — triangulated independently at gene level (DE), pathway level (GSEA), and
patient level (prognostics).

**Honest caveats (A-3):** the proliferation→relapse association does **not** survive BH-FDR
(0.20) or covariate adjustment (at ~10 events / 4 params, even Firth gives an uninformative
CI). Overall survival is **unmodelable** locally (only 5 deaths; `vital_status` is `Expired`).
Composition fractions are *not* prognostic for relapse — composition is the histology-axis
signal (epithelial↑ anaplastic), consistent with the prior composition positive.

## Phase B (spatial) — 1 positive, 2 honest nulls (resolution ceiling)

| Lever | Method | Result |
|---|---|---|
| baseline | Phikon-v1 mean-pool, 60 spots | AUC 0.724 |
| **B-1** scale+encoder+MIL | Phikon-v2, 200 spots, attention-MIL | mean-pool **0.733** [0.55,0.86] p=0.006; MIL **0.748** [0.57,0.87] p=0.003; **paired DeLong MIL vs mean-pool p=0.83** |
| **B-3** StarDist morphology | StarDist `2D_versatile_he`, prob_thresh 0.4 | **0.687** [0.50,0.83] **p=0.021** (vs watershed **0.393**) |
| **B-2** ensemble | morphology + embedding | 0.719 [0.53,0.85] p=0.009; **paired vs embedding-only (0.714) p=0.57** |

**Conclusions:**
- The tumor-level histology (anaplasia) signal is **real and stable at ~0.73–0.75** (all models
  beat chance, perm p≈0.003–0.02).
- **B-1 underdelivered:** 3.3× more spots + a ViT-L encoder + attention-MIL together moved AUC
  only +0.024, and MIL is statistically indistinguishable from flat mean-pooling. This is a
  **resolution ceiling** — Visium-*hires* tiles, not true WSI, cap what any encoder extracts.
- **B-3 fixed the segmentation failure** (watershed 0.39 → StarDist 0.687, significant) — the
  hypothesis was right; watershed was the wrong tool. But morphology and the FM embedding read
  the **same** nuclear atypia, so the ensemble adds nothing over the embedding alone (p=0.57).
- Hard ceiling: median ~14 StarDist nuclei/tumor on hires tiles — lifting Phase B needs
  higher-resolution input (gated WSI FMs / XMAG), i.e. external data (Tier-3).

## Validation hardening (V-1)
DeLong variance CIs + 10k-permutation p on all AUCs (`phase_b_stats.py`); Firth
profile-likelihood CIs for small-n logistic; BH-FDR across all DE/GSEA/composition tests;
paired DeLong for model-vs-model comparisons. Seeds logged.

## ABM payoff (ABM-1)
`17_positives_to_abm.py` maps the positives → per-tumor PhysiCell inputs
(`results/abm/positives_to_physicell.yaml`, `per_tumor_abm_params.csv`):
- composition → **initial cell-type fractions**;
- proliferation score → **proliferation_rate** multiplier (bounded `1+0.6·z`);
- p53-target score → **apoptosis_rate** multiplier;
- H&E anaplasia probability → **high-grade regime** (extra proliferation, less adhesion).

Directional check (encodes the biology, not fit to it): mean proliferation multiplier
**1.40 in relapse vs 0.98 in non-relapse**.

## Infrastructure / fixes
- Installs (all verified): edgeR 4.10.1 + limma 3.68.4, TensorFlow 2.21 + StarDist, fgsea +
  msigdbr, logistf 1.26.1.
- Bugs fixed: R-not-on-PATH (`Program Files (x86)`); `vital_status` = `"Expired"`; small-n
  logistic separation → Firth; StarDist Windows symlink → directory junction.
- Resilience: background jobs die on session teardown — long compute now checkpoints to
  append-only chunk parquets and resumes (the phikon-v2 embedding survived ~5 teardowns).

See `.learnings/LEARNINGS.md` for the durable rules. Scripts:
`phase1_mechanotypes/14_moderated_de.R`, `15_hallmark_gsea.R`, `16_prognostic_association.R`;
`phase2_histology_ml/15_phase_b_mil.py`, `16_stardist_morphology.py`, `17_positives_to_abm.py`,
`phase_b_stats.py`.
