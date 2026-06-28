# AGENTS.md

Wilms tumor computational biology repo: cluster distributions of single-nuclei feature
scores into "mechanotypes" (Phase A), and classify cell state from Visium H&E images to
parameterize a PhysiCell agent-based model (Phase B). Data: ScPCA Portal `SCPCP000006`.

## Setup & commands

```bash
# environment (conda)
conda env create -f environment.yml && conda activate sc-wilms-data

Rscript scripts/install_r_packages.R
Rscript scripts/scpca_auth.R              # sets SCPCA_AUTH_TOKEN via get_auth()

# Phase A (R): explore metadata (no token) -> download -> mechanotypes
Rscript phase1_mechanotypes/00_explore_scpca.R
Rscript phase1_mechanotypes/01_download.R        # uses SCPCA_AUTH_TOKEN
Rscript phase1_mechanotypes/04_wasserstein_matrix.R
Rscript phase1_mechanotypes/05_consensus_cluster.R
Rscript phase1_mechanotypes/08_figures.R

# Phase B (Python): H&E -> cell state -> ABM inputs
python phase2_histology_ml/02_segment_nuclei.py
python phase2_histology_ml/04_train_classifier.py
python phase2_histology_ml/07_figures.py

# All figures
scripts/run_figures.bat   # Windows

# tests
pytest -q
Rscript -e 'testthat::test_dir("tests")'
```

R handles ScPCA access + clustering (`ScPCAr`, `ConsensusClusterPlus`, `maotai`, `waddR`).
Python handles imaging/ML/ABM (`scanpy`, `stardist`, `cellpose`, `scikit-image`, `scikit-learn`).

## Hard rules (these are the non-obvious ones — follow exactly)

- **Wasserstein distance runs ONLY on 1-D feature scores, never on the gene × cell
  matrix.** Multivariate Wasserstein underperforms on high-dimensional scRNA-seq. If a
  change would pass high-dim vectors into a Wasserstein call, stop and flag it.
- **Never commit raw or patient data.** It is pulled via the `ScPCAr` API into `data/raw/`
  (git-ignored). Log the ScPCA access date with each pull — the Portal versions its data.
- **The cell-state classifier operates at the nucleus level.** A Visium spot is ~55 µm and
  holds ~10–20 cells, so spot expression is a *mixture* — use it only for weak-supervision
  labels and for validation against deconvolution, never as a per-cell label.
- **Cell-state labels come from weak/spatial supervision, not manual annotation.** If a
  task seems to need hand-labeled pathology, flag it instead of adding an annotation step.
- **The Phase A feature list is fixed in `config/` before clustering.** Do not add or drop
  features to improve clusters; report sensitivity to the feature set instead.
- **Inclusion rule:** a (cell state × group) item enters clustering only with ≥25 cells for
  that feature.
- **Seed every stochastic step** (clustering subsampling, train/test split, model init) and
  log the seed.

## Conventions

- Numbered scripts (`NN_name.{R,py}`) are standalone and idempotent; outputs go to
  `data/processed/` or `results/`. Re-running must not corrupt prior outputs.
- No hard-coded absolute paths — use `config/` or CLI args.
- Notebooks in `notebooks/` are for exploration only; anything reproducible becomes a
  numbered script.
- Update `environment.yml` when adding a dependency.

## Validation before declaring a task done

- Phase A: pipeline runs API → figures; k / PAC / CHI logged.
- Phase B: segmentation overlays visually checked; classifier reports **per-class balanced
  accuracy** (raw accuracy is misleading under class imbalance); spot fractions checked
  against deconvolution.
- ABM: at least one PhysiCell run reproducible from config.

## Flag to a human (don't decide silently)

- Any Wasserstein call on multi-dimensional input.
- Subgroups dropping below the ≥25-cell rule after stratification.
- Switching segmentation tool (StarDist → Cellpose) — record why in the commit.
- Anything that would commit raw/identifiable data.

## Self-improvement protocol

This file is meant to get better over time. Lessons live in `.learnings/LEARNINGS.md`
(create it from the template there if missing) so this file stays short.

1. **Read first.** At the start of a non-trivial task, read `.learnings/LEARNINGS.md` and
   apply anything relevant.
2. **Log when any of these happen** — append a dated entry to `.learnings/LEARNINGS.md`:
   - a command, tool, or pipeline step failed in a non-obvious way;
   - the user corrected you ("no", "actually", "that's wrong", redo);
   - a needed capability/file/config didn't exist;
   - you found a better approach for a recurring task;
   - you hit a data/methodology gotcha worth remembering.
   Use one entry per lesson: `category | date | what happened | the rule to follow next time`.
   Categories: `correction` | `gotcha` | `better-approach` | `knowledge-gap`.
3. **Don't bloat this file.** Default to logging in `.learnings/`. Only propose editing
   AGENTS.md directly for a rule that is universal and will apply every session.
4. **Promote periodically.** When a lesson has recurred or is clearly durable, propose
   moving it from `.learnings/` into the "Hard rules" section above — as a concise line —
   and remove it from the log. Surface the proposed diff to the human; don't silently
   rewrite the rules.
5. **Keep it lean.** AGENTS.md should stay well under ~150 lines; long or situational
   guidance belongs in `.learnings/` or `PRD.md`, not here.

Methodology mirrors the lab's published mechanobiology and prostate CRPC frameworks; keep
that lineage visible. See `PRD.md` for full requirements and acceptance criteria.
