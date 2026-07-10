# Phase C — omics-initialized spatial agent-based model (PhysiCell)

Turns the Phase A molecular scores and Phase B histology signals into per-tumor,
spatially-resolved PhysiCell inputs. Each Visium tumor becomes an agent-based model whose
initial tissue is placed from the real spot coordinates, whose cell identities come from
spot deconvolution, and whose behavior follows a cell-behavior grammar (Johnson et al.,
*Cell* 2025) anchored to that tumor's Phase A/B-derived rates.

Authoring runs on CPU (this module). The simulations themselves run on a cluster with a
grammar-enabled PhysiCell (≥1.14.1).

## Pipeline

| Stage | Script | Output (per tumor, `results/abm/<sample_id>/`) |
|---|---|---|
| 1 | `01_spot_density.py` | `results/abm/{spot,tumor}_density.csv` — StarDist nuclei/spot (diagnostic; see note) |
| 2 | `02_place_agents.py` | `cells.csv` — agents at Visium coords (µm), identities from deconvolution; **patch mode** emits one dir per (tumor, FOV patch) |
| 3 | `03_emit_rules.py` | `rules.csv` (+ `rules_annotated.csv`) — grammar; oxygen-necrosis + pressure half-maxes are **per-tumor** (hypoxia / crowding programs) |
| 4 | `04_build_model.py` | `PhysiCell_settings.xml` + `provenance.json` — domain + oxygen/IGF2/ECM substrates + cell defs (per-tumor motility + secretion) |

Downstream: `05_uq.py` sensitivity sweep, `06_run_cohort.sh` SLURM array (cluster),
`07_validate.py` emergent (patient-level growth/invasion) + spatial QoIs.

### Omics-determined parameters (v1)

Beyond the initial fractions + proliferation/apoptosis multipliers, `17_positives_to_abm.py`
derives, per tumor, from snRNA program scores, with **transparent bounded transforms (no
fitting)** and a **neutral fallback** when a score is absent. The scores are produced by
`phase1_mechanotypes/17_abm_program_scores.R` (run after `16_prognostic_association.R`; same
z-scored pseudobulk-logCPM method), which augments `per_tumor_scores.csv` from the gene sets
in `config/abm_programs.yaml`:

| Program | PhysiCell parameter | direction |
|---|---|---|
| EMT (epithelial vs mesenchymal) | `adhesion_strength` + `migration_speed` (reciprocal) | mesenchymal ↑ → adhesion ↓, motility ↑ |
| Contact-inhibition (crowding) | `pressure→cycle-entry` **half-max** | higher → brake at lower pressure |
| Hypoxia tolerance | `oxygen→necrosis` **half-max** | higher → necrosis at lower O₂ |
| IGF (IGF2/IGF1R; 11p15 LOI) **[v1.1]** | per-cell **IGF2 uptake_rate** | higher → more IGF2 consumed |

Half-maxes are the omics-determined knob because they dominate ABM QoIs (Johnson et al.,
*Cell* 2025, Fig 2E), unlike base-rate perturbations. Bounds/slopes live in
`config/phase_c.yaml → omics_to_params`. New geometry QoIs in `07_validate.py`: **clustering
index** (homotypic self-segregation) and **invasiveness** (radial projections, per the *Cell*
2025 STAR Methods), computed on the initial `cells.csv` now as the t0 baseline.

**v1.1 substrates** (`config/phase_c.yaml → substrates`; `04_build_model.py` writes the
BioFVM fields + per-cell secretion): **IGF2** — a diffusible growth factor (the Wilms analog of
the CRPC androgen-uptake axis; tumor cells consume it at the IGF-scaled rate, coupling clones
through local depletion) with an `IGF2→cycle-entry` rule; **ECM** — a near-static matrix
secreted by stromal cells, with an `ECM→migration-speed` (decrease) rule encoding the Johnson
et al. fibroblast-barrier effect. Oxygen now has a per-cell uptake so hypoxic gradients form.

## Spatial validation QoIs (`07_validate.py`)

Compartment architecture that the simulated tissue must reproduce, computed on the real
Visium now and on simulated endpoints later:

- **Neighbourhood enrichment** — permutation z-score of compartment adjacency on a spatial
  kNN graph (one number per compartment pair).
- **Co-occurrence curve** — co-location enrichment **vs distance** (a curve, not a scalar):
  `P(neighbour is B | centre is A, at range r) / P(B)`. >1 co-locate, <1 segregate.
- **Ripley's L** (squidpy only) — multi-scale clustering vs dispersion.

Both the enrichment and the co-occurrence curve are implemented in-house (scipy kNN +
`cKDTree.count_neighbors`), so the default install needs no extra dependency and they are
always computed. Set `spatial_qoi.backend: squidpy` (or `auto`, which uses squidpy when
importable) to *also* emit squidpy's battle-tested nhood z-scores and **Ripley's L** as a
cross-check (`pip install squidpy`). This is opt-in because Ripley is ~O(n²) over
full-resolution spots; the in-house path stays the fast default.

Observed baseline (real Visium): blastemal↔epithelial **co-locate** at short range
(enrichment z ≈ +2.5; near-range co-occurrence ≈ 1.3), stromal **segregates** from both
(z ≈ −8; co-occurrence < 1) — the nephrogenic transition vs the stromal compartment.

## Compute scale — representative FOV patches

Simulating the whole ~6 mm slide is ~70k agents / ~80k voxels **per tumor** — hours per run
and ~1 TB of output over the cohort. The CRPC-lab methodology instead runs small
data-initialised tissues, so `02_place_agents.py` defaults to **patch mode**
(`config/phase_c.yaml → patch`): a few compositionally-diverse ~700 µm windows per tumor
(~1.5k agents, ~2k voxels each), one model dir per `(tumor, patch)`. Patches of a tumor share
its parameters and are aggregated to the patient in `07_validate.py` (the unit of inference).

| | whole-slide | patch (default) |
|---|---|---|
| agents / run | ~70k | **~1.5k** (48× less) |
| voxels / run | ~80k | **~2.1k** (38× less) |
| cohort runs | 41 × 20 = 820 | 111 × 10 = 1,110 |
| est. core-hours (SU) | ~13k–79k | **~0.7k–1.8k** |
| output | ~1 TB | **~7 GB** |
| replicates | 20 | 10 (matches the CRPC paper) |
| UQ | all tumors (~2k runs) | representative (`05_uq.py`, ~150 runs) |

Set `patch.enabled: false` for the full whole-slide run. Absolute per-run wallclock still
needs one timed calibration run on the cluster, but the agent/voxel reductions are structural.

## Run

```bash
# authoring (CPU)
python phase3_abm/01_spot_density.py
python phase3_abm/02_place_agents.py            # all tumors (or --sample SCPCS000168)
python phase3_abm/03_emit_rules.py
python phase3_abm/04_build_model.py
pytest tests/test_phase_c.py -q

# cluster (Bridges-2), after building the grammar-enabled PhysiCell binary:
PHYSICELL_BIN=/path/to/PhysiCell/project bash phase3_abm/calibrate.sh   # go/no-go + timing
PHYSICELL_BIN=/path/to/PhysiCell/project REPLICATES=10 \
  sbatch --array=0-$(( $(wc -l < results/abm/model_manifest.txt)*10 - 1 ))%200 \
         phase3_abm/06_run_cohort.sh                                    # cohort (account/partition baked in)
```

`calibrate.sh` times a median + max patch, sanity-checks growth (flags explosion/collapse),
and extrapolates cohort core-hours — run it **before** the array and only proceed if both
verdicts are OK.

Configuration is in [`config/phase_c.yaml`](../config/phase_c.yaml) (spot geometry, density
prior, deconvolution backend, domain, simulation). Base cell rates come from
`config/physicell.yaml`; per-tumor scaling comes from
`results/abm/positives_to_physicell.yaml` (emitted by `phase2_histology_ml/17_positives_to_abm.py`).

## Scope & honest limitations

- **Intrinsic side is the grammar paradigm, not an MHS signaling model.** The CRPC lab's
  "intrinsic" is a mechanistic AR/PI3K/p53 ODE-Boolean network whose species are set by DEGs.
  Phase C instead uses **reduced-form** rate multipliers from program scores + Johnson-grammar
  environmental rules. Cell↔environment coupling *is* present (oxygen/IGF2/pressure/ECM
  dynamically modulate each cell), but there is no intracellular network computing fate. A
  Wilms signaling model would be a separate (v3) effort.
- **Compartment-resolved DE (cell × tumor type).** EMT adhesion/motility is now set
  per-compartment from `<program>_score__<compartment>` (17_abm_program_scores.R), so
  blastemal / epithelial / stromal cells differ within a tumor; crowding/hypoxia half-maxes
  remain tumor-level.
- **Uncalibrated priors → relative, not absolute.** Base rates, half-maxes, k-slopes,
  diffusion/uptake constants are literature/biophysical priors (no fitting). Read outputs as
  **relative contrasts** (favorable vs anaplastic), not absolute predictions; the `05_uq.py`
  sweep quantifies sensitivity to every uncertain knob.
- **mRNA ≠ pathway activity.** YAP/TAZ and PI3K are set post-translationally, so the crowding
  half-max rests partly on weak proxies — trust the cell-cycle effectors (CCND1/MKI67/CDKN1x)
  in the `contact_inhibition` program.
- **Histology → positions is partly deferred.** Agents are placed from Visium spot coordinates
  + deconvolution (omics positions); **necrotic-territory seeding is on** — spots with
  anomalously high mitochondrial fraction (cohort z(pct_mito) > 1.5, the necrosis readout
  validated in the regressive pilot) are seeded as inert `necrotic` tissue, so the model starts
  with the tumor's real necrosis geography. Sub-spot placement from StarDist nucleus centroids
  remains WSI-gated (the documented Visium-hires nuclei ceiling).

## Design notes

- **Coordinates.** Full-res pixel spot centres → microns via the Visium
  `spot_diameter_fullres` scalefactor (isotropic µm/px), recentred to the origin + margin.
- **Density.** A 55 µm capture spot in dense small-round-blue-cell tumor holds ~10–15 cells;
  the default `N_cells_per_location` is that literature prior. StarDist on hires spot tiles
  undercounts absolute nuclei (the documented resolution ceiling) so it is reported as a
  diagnostic and is only used for relative variation (`density.source: stardist`) once
  WSI-resolution segmentation exists.
- **Deconvolution backend** is a switch (`program_softmax` default; `nnls` / `cell2location`
  for the v2 fetal-kidney-subtype granularity).
- **Rules** share one environmental response shape (oxygen→cycle-entry, oxygen→necrosis,
  pressure→contact-inhibition); relapse/anaplastic tumors differ through their base rates
  (higher proliferation, lower adhesion), not through bespoke rules.
- **Provenance & reproducibility.** Seeded placement; each model dir records seed, sources,
  backend, and the target PhysiCell version.

The generated `PhysiCell_settings.xml` is scaffolding: validate it against the target
PhysiCell sample project on the cluster before the first run (PRD AC5).
