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
| 2 | `02_place_agents.py` | `cells.csv` — agents at Visium coords (µm), identities from deconvolution |
| 3 | `03_emit_rules.py` | `rules.csv` (+ `rules_annotated.csv`) — grammar, saturation anchored to base rates |
| 4 | `04_build_model.py` | `PhysiCell_settings.xml` + `provenance.json` — domain + oxygen + cell defs |

Downstream (cluster / follow-up): `05_uq.py` sensitivity, `06_run_cohort.sh` SLURM array,
`07_validate.py` emergent (patient-level growth/invasion) + spatial (squidpy) QoIs.

## Run

```bash
python phase3_abm/01_spot_density.py
python phase3_abm/02_place_agents.py            # all tumors (or --sample SCPCS000168)
python phase3_abm/03_emit_rules.py
python phase3_abm/04_build_model.py
pytest tests/test_phase_c.py -q
```

Configuration is in [`config/phase_c.yaml`](../config/phase_c.yaml) (spot geometry, density
prior, deconvolution backend, domain, simulation). Base cell rates come from
`config/physicell.yaml`; per-tumor scaling comes from
`results/abm/positives_to_physicell.yaml` (emitted by `phase2_histology_ml/17_positives_to_abm.py`).

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
