#!/usr/bin/env python3
"""Quantities of interest (QoIs) from PhysiCell runs — used by Stage 6/7.

The scientific QoI computation (`compute_qoi`) is a pure function over a time census + a
final-frame position table, so it is unit-tested on CPU without PhysiCell. The I/O that
turns a PhysiCell output directory into those tables (`load_run`) uses pcdl (PhysiCell Data
Loader) when available; it is exercised on the cluster where real MCDS output exists.

Modes:
  --run-dir DIR --out FILE   one run  -> qoi.csv
  --aggregate SAMPLE_DIR     mean/median across replicates/rep_*/qoi.csv -> output/qoi.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def compute_qoi(census: pd.DataFrame, final_positions: pd.DataFrame) -> dict[str, float]:
    """QoIs from a total-count time census and the final live-cell positions.

    census: columns time, count (total live cells), optional per-compartment count cols.
    final_positions: columns x, y (live cells at the final frame).
    """
    census = census.sort_values("time")
    t = census["time"].to_numpy(float)
    c = census["count"].to_numpy(float)
    qoi: dict[str, float] = {}
    qoi["final_total_cells"] = float(c[-1])
    qoi["fold_growth"] = float(c[-1] / c[0]) if c[0] > 0 else float("nan")
    # growth AUC normalised by duration (mean population over the run)
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy>=2 renamed trapz
    qoi["growth_auc"] = float(_trapz(c, t) / (t[-1] - t[0])) if t[-1] > t[0] else float(c[0])

    if len(final_positions):
        xy = final_positions[["x", "y"]].to_numpy(float)
        cen = xy.mean(0)
        r = np.linalg.norm(xy - cen, axis=1)
        qoi["radial_extent_um"] = float(np.percentile(r, 95))       # invasive front
        qoi["radius_median_um"] = float(np.median(r))
        # invasion: how heavy-tailed the radial spread is (front vs bulk)
        qoi["invasion_index"] = float(np.percentile(r, 95) / (np.median(r) + 1e-9))
    for col in [c for c in census.columns if c.endswith("_count") and c != "count"]:
        comp = col[:-6]
        qoi[f"final_{comp}_frac"] = float(census[col].iloc[-1] / c[-1]) if c[-1] > 0 else 0.0
    return qoi


def load_run(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """PhysiCell output dir -> (census, final_positions) via pcdl. Cluster-side."""
    try:
        import pcdl  # type: ignore
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "pcdl (PhysiCell Data Loader) required to parse PhysiCell output; "
            "pip install pcdl. compute_qoi() itself needs no PhysiCell.") from e
    mcds_ts = pcdl.TimeSeries(str(run_dir / "output"))
    rows, final_xy = [], pd.DataFrame(columns=["x", "y"])
    for mcds in mcds_ts.get_mcds_list():
        cell = mcds.get_cell_df()
        live = cell[cell.get("dead", False) == False] if "dead" in cell else cell  # noqa: E712
        row = {"time": mcds.get_time(), "count": len(live)}
        if "cell_type" in live:
            for ct, n in live["cell_type"].value_counts().items():
                row[f"{ct}_count"] = int(n)
        rows.append(row)
        final_xy = live[["position_x", "position_y"]].rename(
            columns={"position_x": "x", "position_y": "y"})
    return pd.DataFrame(rows).fillna(0), final_xy


def aggregate(sample_dir: Path) -> pd.DataFrame:
    reps = sorted(sample_dir.glob("replicates/rep_*/qoi.csv"))
    if not reps:
        raise FileNotFoundError(f"no replicate qoi.csv under {sample_dir}/replicates")
    df = pd.concat([pd.read_csv(r) for r in reps], ignore_index=True)
    agg = df.mean(numeric_only=True).add_suffix("_mean")
    med = df.median(numeric_only=True).add_suffix("_median")
    out = pd.concat([agg, med]).to_frame().T
    out.insert(0, "n_replicates", len(reps))
    (sample_dir / "output").mkdir(parents=True, exist_ok=True)
    out.to_csv(sample_dir / "output" / "qoi.csv", index=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--aggregate", type=Path)
    args = ap.parse_args()

    if args.aggregate:
        out = aggregate(args.aggregate)
        print(f"[ok] aggregated QoIs -> {args.aggregate/'output'/'qoi.csv'}")
        print(out.to_string(index=False))
    elif args.run_dir and args.out:
        census, final_xy = load_run(args.run_dir)
        pd.DataFrame([compute_qoi(census, final_xy)]).to_csv(args.out, index=False)
        print(f"[ok] QoI -> {args.out}")
    else:
        ap.error("use --run-dir/--out or --aggregate")


if __name__ == "__main__":
    main()
