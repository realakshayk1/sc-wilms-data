#!/usr/bin/env python3
"""Stage 2: place agents from Visium coordinates -> PhysiCell cells.csv (per tumor).

For each in-tissue spot: convert its full-res pixel centre to microns, look up the
N_cells_per_location prior (StarDist per-spot, else per-tumor median, else fixed), seed
that many agents whose compartment identities are drawn from the spot's deconvolved
fractions, and jitter them within the capture-spot radius. Coordinates are recentred so
the tumor sits at the origin with a configurable margin. Deterministic given the seed.

Two modes (config/phase_c.yaml -> patch):
  * whole-slide (patch.enabled=false): one cells.csv per tumor over ALL in-tissue spots
    (~tens of thousands of agents; expensive to simulate).
  * PATCH (default): sample a few representative field-of-view windows per tumor and emit
    one model dir per (tumor, patch) -> results/abm/<sample_id>__p<k>/cells.csv. This
    matches the CRPC-lab scale (small data-initialised tissues, ~1e2-1e3 cells) at a large
    compute saving with no cost to the tumor-level results; patches are chosen for
    compositional diversity so within-tumor heterogeneity is still represented.

Writes results/abm/<run_id>/cells.csv (+ placement_summary.csv / patch_manifest.csv).

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


def _tumor_spots(cfg, sample_id, lib, sig) -> pd.DataFrame | None:
    """In-tissue spots for a tumor with micron coords + normalised fraction columns."""
    spot = cfg["phase_c"]["spot"]
    frac_cols = cfg["phase_c"]["deconvolution"]["frac_cols"]
    coords = load_spot_coords_um(lib["library_dir"], float(spot["diameter_um"]))
    s = sig[(sig["sample_id"] == sample_id) & (sig["in_tissue"] == 1)].copy()
    s = s.join(coords[["x_um", "y_um"]], on="barcode").dropna(subset=["x_um", "y_um"])
    if s.empty:
        return None
    fr = np.clip(np.nan_to_num(s[frac_cols].to_numpy(float)), 0, None)
    rowsum = fr.sum(axis=1, keepdims=True)
    fr = np.where(rowsum > 0, fr / rowsum, 1.0 / len(COMPARTMENTS))
    for j, c in enumerate(COMPARTMENTS):
        s[f"_p_{c}"] = fr[:, j]
    return s


def _place_from_spots(cfg, s, rng) -> pd.DataFrame:
    """Seed jittered agents for a set of spots (raw microns, not recentred)."""
    dens, spot = cfg["phase_c"]["density"], cfg["phase_c"]["spot"]
    per_spot_n, per_tumor_n = _density_lookup(cfg)
    fixed_n = int(dens["fixed_n"])
    use_stardist = dens.get("source") == "stardist"
    radius = float(spot["place_radius_um"])
    default_n = fixed_n
    p_cols = [f"_p_{c}" for c in COMPARTMENTS]
    rec = []
    for _, row in s.iterrows():
        n = int(per_spot_n.get(row["spot_id"], default_n)) if use_stardist else default_n
        n = max(int(dens["min_n"]), min(int(dens["max_n"]), n))
        p = row[p_cols].to_numpy(float)
        p = p / p.sum() if p.sum() > 0 else np.full(len(p), 1.0 / len(p))  # guard fp drift
        types = rng.choice(len(COMPARTMENTS), size=n, p=p)
        ang = rng.uniform(0, 2 * np.pi, n)
        rad = radius * np.sqrt(rng.uniform(0, 1, n))       # uniform over the disc
        for k in range(n):
            rec.append((row["x_um"] + rad[k] * np.cos(ang[k]),
                        row["y_um"] + rad[k] * np.sin(ang[k]),
                        COMPARTMENTS[types[k]]))
    return pd.DataFrame(rec, columns=["x", "y", "cell_type"])


def _recentre(cells: pd.DataFrame, margin: float, z: float) -> pd.DataFrame:
    cells = cells.copy()
    cells["x"] = cells["x"] - cells["x"].min() + margin
    cells["y"] = cells["y"] - cells["y"].min() + margin
    cells["z"] = z
    return cells[["x", "y", "z", "cell_type"]]


def select_patches(s: pd.DataFrame, size_um: float, n_patches: int,
                   min_spots: int) -> list[np.ndarray]:
    """Choose up to n_patches non-overlapping FOV tiles, maximising compositional
    diversity. Deterministic: densest tile seeds the set, then greedy farthest-point in
    mean-fraction space (ties broken by tile index). Returns positional-index arrays into s."""
    xy = s[["x_um", "y_um"]].to_numpy(float)
    p = s[[f"_p_{c}" for c in COMPARTMENTS]].to_numpy(float)
    gx = np.floor((xy[:, 0] - xy[:, 0].min()) / size_um).astype(int)
    gy = np.floor((xy[:, 1] - xy[:, 1].min()) / size_um).astype(int)
    tiles: dict[tuple[int, int], list[int]] = {}
    for i, t in enumerate(zip(gx.tolist(), gy.tolist())):
        tiles.setdefault(t, []).append(i)
    tiles = {t: idx for t, idx in tiles.items() if len(idx) >= min_spots}
    if not tiles:
        return []
    keys = sorted(tiles, key=lambda t: (-len(tiles[t]), t))     # densest first, deterministic
    mean_frac = {t: p[tiles[t]].mean(axis=0) for t in keys}
    if len(keys) <= n_patches:
        chosen = keys
    else:
        chosen = [keys[0]]
        while len(chosen) < n_patches:
            best, best_d = None, -1.0
            for t in keys:
                if t in chosen:
                    continue
                d = min(float(np.linalg.norm(mean_frac[t] - mean_frac[c])) for c in chosen)
                if d > best_d:
                    best, best_d = t, d
            chosen.append(best)
    return [np.array(tiles[t], dtype=int) for t in chosen]


def place_tumor(cfg, sample_id, lib, sig, rng) -> pd.DataFrame | None:
    """Whole-slide placement: all in-tissue spots -> one recentred cells frame."""
    s = _tumor_spots(cfg, sample_id, lib, sig)
    if s is None:
        return None
    cells = _place_from_spots(cfg, s, rng)
    if cells.empty:
        return None
    return _recentre(cells, float(cfg["phase_c"]["domain"]["margin_um"]),
                     float(cfg["phase_c"]["domain"]["z_um"]))


def place_tumor_patches(cfg, sample_id, lib, sig, rng) -> list[tuple]:
    """Sample representative FOV patches -> [(patch_id, (cx,cy), cells), ...]."""
    pc = cfg["phase_c"]["patch"]
    s = _tumor_spots(cfg, sample_id, lib, sig)
    if s is None:
        return []
    masks = select_patches(s, float(pc["size_um"]), int(pc["n_patches"]), int(pc["min_spots"]))
    margin = float(cfg["phase_c"]["domain"]["margin_um"])
    zval = float(cfg["phase_c"]["domain"]["z_um"])
    out = []
    for k, idx in enumerate(masks):
        sub = s.iloc[idx]
        cells = _place_from_spots(cfg, sub, rng)
        if cells.empty:
            continue
        out.append((k, (float(sub["x_um"].mean()), float(sub["y_um"].mean())),
                    _recentre(cells, margin, zval)))
    return out


def _summary_row(run_id, sample_id, cells, **extra):
    comp = cells["cell_type"].value_counts().to_dict()
    return {"run_id": run_id, "sample_id": sample_id, "n_cells": len(cells),
            "x_span_um": round(cells["x"].max(), 1), "y_span_um": round(cells["y"].max(), 1),
            **{f"n_{c}": comp.get(c, 0) for c in COMPARTMENTS}, **extra}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=None, help="sample_id(s); default all")
    args = ap.parse_args()

    setup_logging()
    cfg = load_config()
    seed = int(cfg["phase_c"]["seed"])
    patch_cfg = cfg["phase_c"].get("patch", {})
    patch_on = bool(patch_cfg.get("enabled", False))
    sig = pd.read_parquet(
        resolve_path(cfg, "data/processed/spot_signatures.parquet"),
        columns=["spot_id", "sample_id", "barcode", "in_tissue",
                 *cfg["phase_c"]["deconvolution"]["frac_cols"]],
    )
    libs = discover_library_dirs(resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"]))
    samples = args.sample or sorted(set(sig["sample_id"]) & set(libs))

    summ, patch_rows, out_dir = [], [], ensure_dir(resolve_path(cfg, "results/abm"))
    for sid in samples:
        if sid not in libs:
            print(f"[skip] {sid}: no Visium library on disk")
            continue
        rng = rng_for(seed, f"place:{sid}")
        if patch_on:
            patches = place_tumor_patches(cfg, sid, libs[sid], sig, rng)
            if not patches:
                print(f"[skip] {sid}: no patch met min_spots")
                continue
            for k, (cx, cy), cells in patches:
                run_id = f"{sid}__p{k}"
                (ensure_dir(out_dir / run_id) / "cells.csv").write_text(cells.to_csv(index=False))
                patch_rows.append(_summary_row(run_id, sid, cells, patch_id=k,
                                                center_x_um=round(cx, 1), center_y_um=round(cy, 1)))
                print(f"[ok] {run_id}: {len(cells)} agents")
        else:
            cells = place_tumor(cfg, sid, libs[sid], sig, rng)
            if cells is None or cells.empty:
                print(f"[skip] {sid}: no placeable spots")
                continue
            (ensure_dir(out_dir / sid) / "cells.csv").write_text(cells.to_csv(index=False))
            summ.append(_summary_row(sid, sid, cells))
            print(f"[ok] {sid}: {len(cells)} agents")

    if patch_on and patch_rows:
        df = pd.DataFrame(patch_rows)
        df.to_csv(out_dir / "patch_manifest.csv", index=False)
        df.to_csv(out_dir / "placement_summary.csv", index=False)
        print(f"[ok] {len(df)} patches over {df.sample_id.nunique()} tumors "
              f"(median {int(df.n_cells.median())} agents/patch) -> patch_manifest.csv")
    elif summ:
        pd.DataFrame(summ).to_csv(out_dir / "placement_summary.csv", index=False)
        print(f"[ok] summary -> {out_dir/'placement_summary.csv'}")


if __name__ == "__main__":
    main()
