#!/usr/bin/env python3
"""Stage 1: N_cells_per_location prior from StarDist nuclei counts.

Counts segmented nuclei per Visium spot (StarDist, prob_thresh 0.4 morphology run) to
give each location an image-grounded initial cell number. Segmentation covers only a
sampled subset of spots, so we also emit a per-tumor median that backfills unsegmented
spots at placement time. Bounded by config/phase_c.yaml density.{min_n,max_n}.

Writes results/abm/spot_density.csv (per spot) + results/abm/tumor_density.csv (per tumor).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import ensure_dir, load_config, resolve_path, setup_logging  # noqa: E402


def main() -> None:
    setup_logging()
    cfg = load_config()
    dens = cfg["phase_c"]["density"]
    lo, hi = int(dens["min_n"]), int(dens["max_n"])

    sd_path = resolve_path(cfg, "data/processed/nucleus_features_stardist_80_pt40.parquet")
    if not sd_path.exists():
        raise FileNotFoundError(f"StarDist features missing: {sd_path}")
    sd = pd.read_parquet(sd_path, columns=["spot_id", "sample_id", "nucleus_id"])

    per_spot = (
        sd.groupby(["sample_id", "spot_id"], observed=True)["nucleus_id"]
        .count().rename("n_cells").reset_index()
    )
    per_spot["n_cells"] = per_spot["n_cells"].clip(lo, hi)

    per_tumor = (
        per_spot.groupby("sample_id", observed=True)["n_cells"]
        .agg(median_n="median", mean_n="mean", n_spots_segmented="count").reset_index()
    )
    per_tumor["median_n"] = per_tumor["median_n"].round().clip(lo, hi).astype(int)

    out_dir = ensure_dir(resolve_path(cfg, "results/abm"))
    per_spot.to_csv(out_dir / "spot_density.csv", index=False)
    per_tumor.to_csv(out_dir / "tumor_density.csv", index=False)
    print(f"[ok] per-spot density  -> {out_dir/'spot_density.csv'}  ({len(per_spot)} spots)")
    print(f"[ok] per-tumor density -> {out_dir/'tumor_density.csv'}  ({len(per_tumor)} tumors)")
    print(f"[info] median nuclei/spot across tumors: {int(np.median(per_tumor['median_n']))}")
    print(per_tumor.to_string(index=False))


if __name__ == "__main__":
    main()
