#!/usr/bin/env python3
"""Stage 2: place agents from Visium coordinates -> PhysiCell cells.csv (per tumor).

For each in-tissue spot: convert its full-res pixel centre to microns, look up the
N_cells_per_location prior (StarDist per-spot, else per-tumor median, else fixed), seed
that many agents whose compartment identities are drawn from the spot's deconvolved
fractions, and jitter them within the capture-spot radius. Coordinates are recentred so
the tumor sits at the origin with a configurable margin. Deterministic given the seed.

Writes results/abm/<sample_id>/cells.csv (x,y,z,cell_type) + placement_summary.csv.

Usage: python 02_place_agents.py [--sample SCPCS000168 ...]   (default: all tumors)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import (  # noqa: E402
    COMPARTMENTS, discover_library_dirs, ensure_dir, load_config,
    load_spot_coords_um, resolve_path, rng_for, setup_logging,
)


def _density_lookup(cfg) -> tuple[dict, dict]:
    per_spot, per_tumor = {}, {}
    sp = resolve_path(cfg, "results/abm/spot_density.csv")
    tp = resolve_path(cfg, "results/abm/tumor_density.csv")
    if sp.exists():
        d = pd.read_csv(sp)
        per_spot = dict(zip(d["spot_id"], d["n_cells"]))
    if tp.exists():
        d = pd.read_csv(tp)
        per_tumor = dict(zip(d["sample_id"], d["median_n"]))
    return per_spot, per_tumor


def place_tumor(cfg, sample_id, lib, sig, rng) -> pd.DataFrame | None:
    dens = cfg["phase_c"]["density"]
    spot = cfg["phase_c"]["spot"]
    frac_cols = cfg["phase_c"]["deconvolution"]["frac_cols"]
    per_spot_n, per_tumor_n = _density_lookup(cfg)
    fixed_n = int(dens["fixed_n"])
    use_stardist = dens.get("source") == "stardist"
    default_n = int(per_tumor_n.get(sample_id, fixed_n)) if use_stardist else fixed_n
    radius = float(spot["place_radius_um"])

    coords = load_spot_coords_um(lib["library_dir"], float(spot["diameter_um"]))
    s = sig[(sig["sample_id"] == sample_id) & (sig["in_tissue"] == 1)].copy()
    s = s.join(coords[["x_um", "y_um"]], on="barcode")
    s = s.dropna(subset=["x_um", "y_um"])
    if s.empty:
        return None

    fr = s[frac_cols].to_numpy(float)
    fr = np.clip(np.nan_to_num(fr), 0, None)
    rowsum = fr.sum(axis=1, keepdims=True)
    fr = np.where(rowsum > 0, fr / rowsum, 1.0 / len(COMPARTMENTS))

    rec = []
    for i, (_, row) in enumerate(s.iterrows()):
        n = int(per_spot_n.get(row["spot_id"], default_n)) if use_stardist else default_n
        n = max(int(dens["min_n"]), min(int(dens["max_n"]), n))
        types = rng.choice(len(COMPARTMENTS), size=n, p=fr[i])
        ang = rng.uniform(0, 2 * np.pi, n)
        rad = radius * np.sqrt(rng.uniform(0, 1, n))   # uniform over the disc
        for k in range(n):
            rec.append((row["x_um"] + rad[k] * np.cos(ang[k]),
                        row["y_um"] + rad[k] * np.sin(ang[k]),
                        COMPARTMENTS[types[k]]))
    cells = pd.DataFrame(rec, columns=["x", "y", "cell_type"])

    # recentre to origin + margin, add z
    margin = float(cfg["phase_c"]["domain"]["margin_um"])
    cells["x"] = cells["x"] - cells["x"].min() + margin
    cells["y"] = cells["y"] - cells["y"].min() + margin
    cells["z"] = float(cfg["phase_c"]["domain"]["z_um"])
    return cells[["x", "y", "z", "cell_type"]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=None, help="sample_id(s); default all")
    args = ap.parse_args()

    setup_logging()
    cfg = load_config()
    seed = int(cfg["phase_c"]["seed"])
    sig = pd.read_parquet(
        resolve_path(cfg, "data/processed/spot_signatures.parquet"),
        columns=["spot_id", "sample_id", "barcode", "in_tissue",
                 *cfg["phase_c"]["deconvolution"]["frac_cols"]],
    )
    libs = discover_library_dirs(resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"]))
    samples = args.sample or sorted(set(sig["sample_id"]) & set(libs))

    summ, out_dir = [], ensure_dir(resolve_path(cfg, "results/abm"))
    for sid in samples:
        if sid not in libs:
            print(f"[skip] {sid}: no Visium library on disk")
            continue
        rng = rng_for(seed, f"place:{sid}")
        cells = place_tumor(cfg, sid, libs[sid], sig, rng)
        if cells is None or cells.empty:
            print(f"[skip] {sid}: no placeable spots")
            continue
        d = ensure_dir(out_dir / sid)
        cells.to_csv(d / "cells.csv", index=False)
        comp = cells["cell_type"].value_counts().to_dict()
        summ.append({"sample_id": sid, "n_cells": len(cells),
                     "x_span_um": round(cells["x"].max(), 1),
                     "y_span_um": round(cells["y"].max(), 1),
                     **{f"n_{c}": comp.get(c, 0) for c in COMPARTMENTS}})
        print(f"[ok] {sid}: {len(cells)} agents -> {d/'cells.csv'}")

    if summ:
        pd.DataFrame(summ).to_csv(out_dir / "placement_summary.csv", index=False)
        print(f"[ok] summary -> {out_dir/'placement_summary.csv'}")


if __name__ == "__main__":
    main()
