#!/usr/bin/env python3
"""Stage 5: uncertainty-quantification / sensitivity sweep manifest (per tumor).

Emits a one-at-a-time (OAT) sensitivity design: each uncertain parameter is perturbed by
+/- each pct level while the others hold at base, plus one baseline run. The manifest lists
every run the cluster should execute; the run wrapper (06_run_cohort) reads it and applies
the override to that tumor's rules.csv / PhysiCell_settings.xml before launching.

The design generation is pure and unit-tested here; only the execution is cluster-side.

Writes results/abm/<sample_id>/uq/sweep_manifest.csv.

Usage: python 05_uq.py [--sample SCPCS000168 ...]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import ensure_dir, load_config, resolve_path, setup_logging  # noqa: E402


def generate_sweep(sample_id: str, base_values: dict[str, float],
                   pct_levels: list[float]) -> pd.DataFrame:
    """OAT design: baseline + (param x +/-pct) rows. Deterministic, order-stable."""
    rows = [{"run_id": f"{sample_id}__base", "sample_id": sample_id, "param": "(baseline)",
             "pct": 0.0, "factor": 1.0, "base_value": None, "perturbed_value": None}]
    for param, base in base_values.items():
        for pct in pct_levels:
            for sign in (-1, 1):
                factor = 1.0 + sign * pct / 100.0
                rows.append({
                    "run_id": f"{sample_id}__{param}__{'m' if sign < 0 else 'p'}{pct:g}",
                    "sample_id": sample_id, "param": param, "pct": sign * pct,
                    "factor": round(factor, 4), "base_value": base,
                    "perturbed_value": round(base * factor, 6)})
    return pd.DataFrame(rows)


def representative_samples(tumors: dict, n: int) -> list[str]:
    """A few tumors spanning the high-grade axis (deterministic), so a sensitivity sweep
    covers both regimes without running the whole cohort (the ~49-runs-each trap)."""
    hi = sorted(s for s, t in tumors.items() if t.get("high_grade_regime"))
    lo = sorted(s for s, t in tumors.items() if not t.get("high_grade_regime"))
    picks, i = [], 0
    while len(picks) < n and (hi or lo):
        pool = (lo, hi)[i % 2] or (hi, lo)[i % 2]
        if pool:
            picks.append(pool.pop(0))
        i += 1
    return picks


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=None, help="explicit sample_id(s)")
    ap.add_argument("--all", action="store_true", help="sweep every tumor (expensive)")
    args = ap.parse_args()

    setup_logging()
    cfg = load_config()
    uq = cfg["phase_c"]["uq"]
    pct_levels = list(uq["pct_levels"])
    base_params = dict(uq["params"])
    out_dir = resolve_path(cfg, "results/abm")

    abm = yaml.safe_load(
        resolve_path(cfg, "results/abm/positives_to_physicell.yaml").read_text())
    tumors = abm.get("tumors", {})
    if args.sample:
        samples = args.sample
    elif args.all:
        samples = list(tumors)
    else:
        samples = representative_samples(tumors, int(uq.get("n_representative", 3)))
        print(f"[info] representative UQ on {samples} (use --all to sweep every tumor)")

    total = 0
    for sid in samples:
        tumor = abm["tumors"].get(sid)
        if tumor is None:
            continue
        vals = dict(base_params)
        # per-tumor adhesion base if available (blastemal as representative)
        adh = tumor["cell_types"].get("blastemal", {}).get("adhesion_strength")
        if adh is not None:
            vals["adhesion_strength"] = float(adh)
        df = generate_sweep(sid, vals, pct_levels)
        d = ensure_dir(out_dir / sid / "uq")
        df.to_csv(d / "sweep_manifest.csv", index=False)
        total += len(df)
    n_per = 1 + len(base_params) * len(pct_levels) * 2
    print(f"[ok] UQ manifests for {len(samples)} tumors "
          f"({n_per} runs/tumor incl. replicates-to-be, {total} rows total)")
    print(f"[info] pct levels {pct_levels} x {len(base_params)} params (OAT + baseline)")


if __name__ == "__main__":
    main()
