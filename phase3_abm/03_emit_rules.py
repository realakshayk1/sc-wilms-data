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


# fallback half-maxes when a tumor predates the per-tumor `half_max` block (17_positives_to_abm)
DEFAULT_HM = {"oxygen_cycle": 10.0, "oxygen_necrosis": 5.0, "pressure_cycle": 1.0}


def rules_for(cell_type: str, rates: dict, half_max: dict | None = None,
              substrates: set | None = None) -> list[dict]:
    """Grammar rows for one (tumor, cell_type). Saturation anchored to base rates; the
    oxygen-necrosis and pressure half-maxes are PER-TUMOR (hypoxia / contact-inhibition
    programs, via 17_positives_to_abm.py). Half-max is where ABM sensitivity concentrates
    (Johnson et al. Cell 2025, Fig 2E), so it is the omics-determined knob. IGF2 / ECM rules
    are added only when those v1.1 substrates are configured; they layer additively onto the
    same behaviours via PhysiCell's multi-rule bilinear law."""
    prolif = float(rates["proliferation_rate"])
    hm = {**DEFAULT_HM, **(half_max or {})}
    substrates = substrates or set()
    rows = [
        # oxygen increases cycle entry (proliferation up to ~1.5x base at high pO2)
        dict(cell_type=cell_type, signal="oxygen", response="increases",
             behavior="cycle entry", saturation_value=round(1.5 * prolif, 6),
             half_max=round(float(hm["oxygen_cycle"]), 4), hill_power=4, apply_to_dead=0,
             text=f"In {cell_type}, oxygen increases cycle entry"),
        # oxygen decreases necrosis (hypoxic core dies) — half-max from hypoxia-tolerance program
        dict(cell_type=cell_type, signal="oxygen", response="decreases",
             behavior="necrosis", saturation_value=0.0,
             half_max=round(float(hm["oxygen_necrosis"]), 4), hill_power=8, apply_to_dead=0,
             text=f"In {cell_type}, oxygen decreases necrosis"),
        # pressure decreases cycle entry (contact inhibition) — half-max from crowding program
        dict(cell_type=cell_type, signal="pressure", response="decreases",
             behavior="cycle entry", saturation_value=0.0,
             half_max=round(float(hm["pressure_cycle"]), 4), hill_power=4, apply_to_dead=0,
             text=f"In {cell_type}, pressure decreases cycle entry"),
    ]
    if "IGF2" in substrates:
        # IGF2 increases cycle entry (Wilms growth-factor uptake axis; 11p15 LOI)
        rows.append(dict(cell_type=cell_type, signal="IGF2", response="increases",
                         behavior="cycle entry", saturation_value=round(1.5 * prolif, 6),
                         half_max=0.5, hill_power=4, apply_to_dead=0,
                         text=f"In {cell_type}, IGF2 increases cycle entry"))
    if "ECM" in substrates:
        # dense ECM decreases migration speed (Johnson et al. 2025 fibroblast-barrier effect)
        rows.append(dict(cell_type=cell_type, signal="ECM", response="decreases",
                         behavior="migration speed", saturation_value=0.0,
                         half_max=0.5, hill_power=4, apply_to_dead=0,
                         text=f"In {cell_type}, ECM decreases migration speed"))
    return rows


def main() -> None:
    setup_logging()
    cfg = load_config()
    yml = resolve_path(cfg, "results/abm/positives_to_physicell.yaml")
    if not yml.exists():
        raise FileNotFoundError(f"run 17_positives_to_abm.py first: {yml}")
    abm = yaml.safe_load(yml.read_text())
    out_dir = ensure_dir(resolve_path(cfg, "results/abm"))
    substrates = set(cfg["phase_c"].get("substrates", {}))   # v1.1 IGF2 / ECM if configured

    n = 0
    for sid, tumor in abm.get("tumors", {}).items():
        rows = []
        half_max = tumor.get("half_max")            # per-tumor (17); None -> DEFAULT_HM
        for ct in COMPARTMENTS:
            rates = tumor["cell_types"].get(ct)
            if rates is None:
                continue
            rows.extend(rules_for(ct, rates, half_max, substrates))
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
