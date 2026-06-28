#!/usr/bin/env python3
"""FR-B7 (P1): Map H&E-derived fractions to PhysiCell parameters; run simulation stub."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging


def fractions_to_initial_cells(frac_row: pd.Series, abm_cfg: dict) -> list[dict]:
    """Convert spot fractions into initial cell placements for PhysiCell XML/config."""
    domain = abm_cfg["domain"]
    nx = int(domain["x_max_um"] / domain["dx_um"])
    ny = int(domain["y_max_um"] / domain["dy_um"])
    n_cells = min(nx * ny, 500)

    cells = []
    types = ["blastemal", "epithelial", "stromal"]
    weights = [frac_row.get(f"frac_{t}", 0.0) for t in types]
    if sum(weights) == 0:
        weights = [1 / 3, 1 / 3, 1 / 3]
    weights = [w / sum(weights) for w in weights]

    rng = __import__("numpy").random.default_rng(abm_cfg["seed"])
    for i in range(n_cells):
        ct = rng.choice(types, p=weights)
        params = abm_cfg["cell_types"][ct]
        mapping = abm_cfg["fraction_to_parameter"].get(ct, {})
        cells.append(
            {
                "id": i,
                "type": ct,
                "x": float(rng.uniform(0, domain["x_max_um"])),
                "y": float(rng.uniform(0, domain["y_max_um"])),
                "proliferation_rate": params["proliferation_rate"]
                * mapping.get("proliferation_scale", 1.0)
                * float(frac_row.get(f"frac_{ct}", 0.33)),
                "adhesion_strength": params["adhesion_strength"]
                * mapping.get("adhesion_scale", 1.0),
            }
        )
    return cells


def write_physicell_config(run_dir: Path, abm_cfg: dict, cells: list[dict]) -> None:
    """Write reproducible run config — PhysiCell binary invoked separately if installed."""
    ensure_dir(run_dir)
    with open(run_dir / "physicell_settings.yaml", "w") as f:
        yaml.dump(abm_cfg, f, default_flow_style=False)
    with open(run_dir / "initial_cells.json", "w") as f:
        json.dump({"cells": cells, "n_cells": len(cells)}, f, indent=2)

    # Stub simulation output for reproducibility without PhysiCell binary
    sim = abm_cfg["simulation"]
    timeline = list(range(0, sim["max_time_min"] + 1, sim["save_interval_min"]))
    summary = {
        "status": "stub_complete",
        "message": "Replace stub with PhysiCell binary when available on cluster",
        "timepoints_min": timeline,
        "final_cell_count": len(cells),
        "seed": abm_cfg["seed"],
    }
    with open(run_dir / "simulation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


def main() -> None:
    setup_logging()
    cfg = load_config()
    root = cfg["root"]
    set_seed_logged(cfg["features"]["seed"], "map_to_physicell")

    with open(root / "config" / "physicell.yaml") as f:
        abm_cfg = yaml.safe_load(f)

    frac_csv = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_fractions_csv"])
    if not frac_csv.exists():
        raise SystemExit("Run 05_spot_fractions.py first")

    frac_df = pd.read_csv(frac_csv)
    tumor_row = frac_df.iloc[0]

    out_base = resolve_path(cfg, cfg["paths"]["abm"]["output_dir"])
    run_dir = ensure_dir(out_base / abm_cfg["simulation"]["output_subdir"])

    if (run_dir / "simulation_summary.json").exists() and "--force" not in __import__("sys").argv:
        print(f"[skip] ABM run exists: {run_dir}")
        return

    cells = fractions_to_initial_cells(tumor_row, abm_cfg)
    write_physicell_config(run_dir, abm_cfg, cells)

    # Copy source config for provenance
    shutil.copy(root / "config" / "physicell.yaml", run_dir / "physicell_source.yaml")
    print(f"[ok] PhysiCell init ({len(cells)} cells) -> {run_dir}")


if __name__ == "__main__":
    main()
