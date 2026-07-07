#!/usr/bin/env python3
"""Stage 3: emit PhysiCell cell-behavior grammar rules per tumor.

Encodes the microenvironment responses in the Johnson et al. 2025 plain-language grammar
("In [cell type], [signal] increases/decreases [behavior]") as PhysiCell rules.csv rows.
Per-tumor saturation values are anchored to that tumor's base rates from
results/abm/positives_to_physicell.yaml (Phase A/B derived), so relapse/anaplastic tumors
inherit their higher proliferation / lower adhesion through the base phenotype while the
environmental response shape stays shared.

PhysiCell rules.csv is headerless with 8 columns:
  cell_type, signal, response, behavior, saturation_value, half_max, hill_power, apply_to_dead
We also write rules_annotated.csv (with header + plain-language text) for humans.

Writes results/abm/<sample_id>/rules.csv (+ rules_annotated.csv).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import COMPARTMENTS, ensure_dir, load_config, resolve_path, setup_logging  # noqa: E402

HEADER = ["cell_type", "signal", "response", "behavior",
          "saturation_value", "half_max", "hill_power", "apply_to_dead"]


def rules_for(cell_type: str, rates: dict) -> list[dict]:
    """Grammar rows for one (tumor, cell_type). saturation anchored to base rates."""
    prolif = float(rates["proliferation_rate"])
    return [
        # oxygen increases cycle entry (proliferation up to ~1.5x base at high pO2)
        dict(cell_type=cell_type, signal="oxygen", response="increases",
             behavior="cycle entry", saturation_value=round(1.5 * prolif, 6),
             half_max=10.0, hill_power=4, apply_to_dead=0,
             text=f"In {cell_type}, oxygen increases cycle entry"),
        # oxygen decreases necrosis (hypoxic core dies)
        dict(cell_type=cell_type, signal="oxygen", response="decreases",
             behavior="necrosis", saturation_value=0.0,
             half_max=5.0, hill_power=8, apply_to_dead=0,
             text=f"In {cell_type}, oxygen decreases necrosis"),
        # pressure decreases cycle entry (contact inhibition)
        dict(cell_type=cell_type, signal="pressure", response="decreases",
             behavior="cycle entry", saturation_value=0.0,
             half_max=1.0, hill_power=4, apply_to_dead=0,
             text=f"In {cell_type}, pressure decreases cycle entry"),
    ]


def main() -> None:
    setup_logging()
    cfg = load_config()
    yml = resolve_path(cfg, "results/abm/positives_to_physicell.yaml")
    if not yml.exists():
        raise FileNotFoundError(f"run 17_positives_to_abm.py first: {yml}")
    abm = yaml.safe_load(yml.read_text())
    out_dir = ensure_dir(resolve_path(cfg, "results/abm"))

    n = 0
    for sid, tumor in abm.get("tumors", {}).items():
        rows = []
        for ct in COMPARTMENTS:
            rates = tumor["cell_types"].get(ct)
            if rates is None:
                continue
            rows.extend(rules_for(ct, rates))
        if not rows:
            continue
        df = pd.DataFrame(rows)
        d = ensure_dir(out_dir / sid)
        df[HEADER].to_csv(d / "rules.csv", index=False, header=False)
        df.to_csv(d / "rules_annotated.csv", index=False)
        n += 1
    print(f"[ok] grammar rules written for {n} tumors -> results/abm/<sample_id>/rules.csv")


if __name__ == "__main__":
    main()
