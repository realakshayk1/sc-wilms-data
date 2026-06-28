# PRD: Wilms Tumor Mechanotyping & Histology-Informed Spatial ABM

| | |
|---|---|
| **Status** | Implemented (pilot) — Phase A complete; Phase B pilot on 480 Visium spots |
| **Owner** | Akshay (PhD, Radhakrishnan Lab, UPenn) |
| **Repo** | `sc-wilms-data` |
| **Last updated** | 2026-06-28 |

### Change history
- v2.1 (2026-06-28): Phase A + B pilot implemented; figure scripts 08/07; comprehensive README.
- v2 (2026-06-28): Restructured to lead with problem/outcomes; added testable
  acceptance criteria (Given/When/Then), non-functional requirements, and P0/P1/P2
  scope tiers. Made self-contained (no external plan dependency).
- v1 (2026-06-28): Initial draft.

---

## 1. Overview (problem & why now)

The Radhakrishnan lab has a published framework that couples cell-intrinsic genomic
heterogeneity with extrinsic microenvironmental constraints in a multiscale model
(cellular signaling + spatial agent-based model). Two reusable methods from that work
are underexploited for pediatric cancer:

1. **Distributional "mechanotypes":** comparing whole *distributions* of a feature
   across (cell-group × condition) pairs via Wasserstein-1 distance + consensus
   clustering, instead of comparing medians.
2. **A spatial agent-based model (PhysiCell)** whose initial conditions currently rely
   on uniform or deconvolution-based assumptions rather than directly observed tissue
   morphology.

**Why Wilms tumor, why now:** Wilms tumor (nephroblastoma) is the most common pediatric
kidney cancer, has morphologically distinct compartments (blastemal / epithelial /
stromal), a clinically dominant favorable-vs-anaplastic histology axis, and — critically
— one public dataset (`SCPCP000006`, ScPCA Portal) that carries snRNA-seq, Visium
spatial transcriptomics with paired H&E images, and bulk RNA-seq for the same project.
No published work has applied the lab's distributional-mechanotype method to Wilms
single-nuclei data, nor used H&E-derived morphology to parameterize a Wilms ABM. The
data and methods both already exist; the contribution is connecting them.

## 2. Goals & success metrics

| ID | Goal | Success metric |
|---|---|---|
| G1 | Reproducible molecular mechanotypes for Wilms cell states across histology | ≥1 interpretable mechanotype-switching pattern between favorable & anaplastic, with statistical support |
| G2 | Nucleus-level cell-state classifier from H&E → per-spot fractions | Per-class balanced accuracy reported on held-out spots; fractions correlate with transcriptomic deconvolution |
| G3 | Spatially-resolved fractions parameterize the ABM | ≥1 PhysiCell run initialized from H&E-derived fractions, reproducible from config |
| G4 | Clean, reproducible GitHub repo | Fresh clone reproduces all figures from pinned environment + documented commands |

## 3. Users & use cases
- **Primary:** the author and lab collaborators (advisor R. Radhakrishnan; Stephanie).
- **Secondary:** future lab members reusing the skeleton for other pediatric tumors.
- **Core use case:** "Given Wilms multi-modal data, identify distributional cell-state
  classes and feed observed spatial composition into the lab's ABM."

## 4. Data

| Modality | Source | Use |
|---|---|---|
| snRNA-seq, 40 samples (23 favorable, 22 anaplastic) | ScPCA `SCPCP000006` | Mechanotypes; cell-state signatures |
| Visium spatial (~100 slides) + paired H&E | ScPCA `SCPCP000006` | Image ML + spatial validation |
| Bulk RNA-seq | ScPCA `SCPCP000006` | Pseudobulk sanity checks |
| Cell-type annotations (blastemal/epithelial/stromal) | OpenScPCA module (already done) | Labels/signatures |
| Optional external bulk: GSE31403, GSE10320, TARGET-WT | GEO / GDC | Validation only |

Access is via the **`ScPCAr`** R package API (token-based). Raw data is never committed;
the Portal access date is recorded because the Portal versions its data.

## 5. Functional requirements (what it must do)

Priority: **P0** = must-have for a defensible result · **P1** = important · **P2** = nice-to-have.

### Phase A — Mechanotypes
- **P0 · FR-A1.** Pull `SCPCP000006` snRNA-seq via `ScPCAr`; QC and normalize; retain
  existing cell-state annotations.
- **P0 · FR-A2.** Reduce each cell to a small, **predefined** set of interpretable 1-D
  feature scores (e.g. blastemal program, epithelial/differentiation program,
  proliferation, WT1, Wnt/β-catenin target score).
- **P0 · FR-A3.** Form items as (cell state × sample/histology) groups; include a group
  for a feature only if it has ≥25 cells.
- **P0 · FR-A4.** Compute pairwise Wasserstein-1 distance **per feature** on the 1-D
  score distributions.
- **P0 · FR-A5.** Run consensus clustering on each feature's distance matrix; select
  cluster count via low PAC + high Calinski–Harabasz; flag boundary items
  (item-consensus < 0.8).
- **P1 · FR-A6.** Decompose the 2-Wasserstein distance into location/shape/size terms
  (waddR) to interpret *why* groups differ.
- **P1 · FR-A7.** Identify and report cell states that switch mechanotype between
  favorable and anaplastic histology.

### Phase B — Histology ML → ABM
- **P0 · FR-B1.** Extract H&E tiles aligned to Visium spot coordinates; stain-normalize.
- **P0 · FR-B2.** Segment nuclei per tile (StarDist default; Cellpose fallback).
- **P0 · FR-B3.** Compute per-nucleus morphology features (area, eccentricity, solidity,
  texture, hematoxylin intensity, neighbor density).
- **P0 · FR-B4.** Derive labels via weak/spatial supervision (spot dominant state from
  Phase A signatures → propagated to nuclei). No manual pathology annotation.
- **P0 · FR-B5.** Train an interpretable classifier (random forest / gradient boosting
  first; CNN only if tabular underperforms); handle class imbalance.
- **P0 · FR-B6.** Aggregate nucleus predictions to per-spot cell-state fractions;
  validate against transcriptomic deconvolution.
- **P1 · FR-B7.** Map fractions to PhysiCell parameters (blastemal→high proliferation/low
  adhesion; stromal→ECM stiffness; epithelial→differentiated) and run ≥1 simulation.

## 6. Non-functional requirements
- **Reproducibility:** single pinned environment file; every stochastic step takes an
  explicit, logged seed; figures regenerable from saved intermediates.
- **Provenance:** record ScPCA access date and dataset version with every data pull.
- **Privacy/compliance:** no raw or patient-identifiable data committed; data is
  research-use-only per ScPCA terms.
- **Portability:** scripts run on a standard Linux workstation / lab cluster; no
  hard-coded absolute paths.
- **Maintainability:** numbered, standalone, idempotent scripts; R for
  access/clustering, Python for imaging/ML/ABM.

## 7. Acceptance criteria (testable)

- **AC1 (G1):** *Given* the snRNA-seq for `SCPCP000006`, *when* the Phase A pipeline is
  run end-to-end, *then* it outputs per-feature mechanotype assignments and figures
  without manual intervention, and logs k, PAC, and CHI for each feature.
- **AC2 (G1):** *Given* mechanotype assignments, *when* favorable and anaplastic groups
  are compared, *then* at least one cell state shows a documented mechanotype switch,
  supported by a waddR location/shape/size decomposition.
- **AC3 (G2):** *Given* a held-out set of H&E tiles, *when* segmentation runs, *then*
  nuclei masks pass a visual-overlay check recorded in the repo.
- **AC4 (G2):** *Given* the trained classifier, *when* evaluated on held-out spots,
  *then* per-class balanced accuracy is reported (not raw accuracy), and per-spot
  fractions show a reported correlation with deconvolution.
- **AC5 (G3):** *Given* H&E-derived fractions for one tumor, *when* the mapping script
  runs, *then* a PhysiCell simulation initializes and completes, reproducibly from config.
- **AC6 (G4):** *Given* a fresh clone and the pinned environment, *when* the documented
  commands are run, *then* all headline figures regenerate.

## 8. Assumptions, dependencies, risks

| Type | Item | Mitigation |
|---|---|---|
| Assumption | ≥25 cells per (state × group) for key features | If not, pool at histology level; use external bulk for validation only |
| Dependency | `ScPCAr` API availability + auth token | Cache downloaded objects locally with recorded version |
| Risk | Wasserstein misused on high-dim data | Hard rule: W1 only on predefined 1-D scores (benchmarks show multivariate W1 underperforms on scRNA-seq) |
| Risk | Feature choice drives clustering | Predefine feature list; report sensitivity to it |
| Risk | H&E segmentation quality varies (esp. kidney tissue) | Start with cleanest tiles; stain-normalize; treat poor tiles as out-of-scope |
| Risk | Spot = mixture of ~10–20 cells | Classify at nucleus level; use spot expression only for labels/validation |
| Risk | Imaging ML is a new skill area | Phase A stands alone as a result; ABM coupling (FR-B7) may move to future work |

## 9. Out of scope (this version)
- Full ABM parameter sweeps / cohort-level tumor-evolution simulation.
- Clinical or survival modeling and prognostic claims.
- Foundation-model segmentation (CellViT++, Cellpose-SAM) — noted upgrade path only.
- Generalization to other ScPCA tumor types.

## 10. Appendix — references informing methodology
- Lab pan-cancer mechanobiology paper (Wasserstein + ConsensusClusterPlus + maotai; n≥25; PAC/CHI).
- Lab prostate CRPC multiscale manuscript (intrinsic–extrinsic coupling; PhysiCell).
- Benchmark: multivariate Wasserstein underperforms on high-dim scRNA-seq → use 1-D scores.
- waddR: 2-Wasserstein decomposition (location/shape/size) for scRNA-seq.
- StarDist for H&E nuclei (CoNIC); segmentation tool comparisons; kidney-pathology caveats.
- ScPCA Portal / `SCPCP000006`; `ScPCAr` API (cite ScPCA preprint DOI: 10.1101/2024.04.19.590243).