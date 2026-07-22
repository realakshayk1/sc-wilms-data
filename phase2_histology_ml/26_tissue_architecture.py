#!/usr/bin/env python3
"""WS4 / P4 (§6.1): per-tumor tissue-ARCHITECTURE descriptors from the Visium deconvolution map.

Turns the per-spot compartment map (data/processed/spot_signatures.parquet: deconv_* fractions +
dominant_state) + spot coordinates into macro-shape statistics that describe how the tumor is
organised — nodular vs diffuse, nodule size/geometry, compartment interface. These parameterize the
ABM initial *shape* (§6.2) and feed the architecture->invasion question (§6.3, ABM-gated).

Honest resolution (audit F): this is the compartment map at ~55um Visium-spot resolution — the same
resolution as its labels. Phikon tile embeddings would add coverage beyond the spot grid, not sub-spot
detail; that refinement is deferred (needs the per-tile embedding table, a heavier extraction).

Pure metric functions (morans_i, detect_nodules, architecture_metrics) are unit-tested.
Output: results/spatial/tissue_architecture.csv (per library + per-sample aggregate).
Usage: python phase2_histology_ml/26_tissue_architecture.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "phase3_abm"))
from abm_utils import (ensure_dir, load_config, load_spot_coords_um,  # noqa: E402
                       resolve_path)

COMPARTMENTS = ["blastemal", "epithelial", "stromal"]


def _adjacency(coords: np.ndarray, spacing_mult: float = 1.5):
    """Symmetric spot-adjacency (pairs within spacing_mult x median nearest-neighbour distance)."""
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    d, _ = tree.query(coords, k=2)
    spacing = float(np.median(d[:, 1])) if len(coords) > 1 else 1.0
    pairs = tree.query_pairs(r=spacing_mult * spacing, output_type="ndarray")
    return pairs, spacing


def morans_i(values: np.ndarray, pairs: np.ndarray) -> float:
    """Global Moran's I of `values` over the adjacency `pairs` (binary weights). +1 clustered,
    0 random, -1 checkerboard. Returns nan if degenerate."""
    x = np.asarray(values, float)
    n = len(x)
    if n < 3 or len(pairs) == 0:
        return float("nan")
    xc = x - x.mean()
    denom = np.sum(xc ** 2)
    if denom == 0:
        return float("nan")
    num = np.sum(xc[pairs[:, 0]] * xc[pairs[:, 1]])           # each undirected pair once
    W = len(pairs)                                            # sum of weights over ordered pairs / 2
    return float((n / (2.0 * W)) * (2.0 * num) / denom)


def detect_nodules(mask: np.ndarray, pairs: np.ndarray, min_size: int = 3) -> list[int]:
    """Connected-component sizes of the masked (e.g. blastemal-dominant) spots over `pairs`."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    remap = {g: i for i, g in enumerate(idx)}
    sub = [(remap[a], remap[b]) for a, b in pairs if a in remap and b in remap]
    m = len(idx)
    if sub:
        r, c = zip(*sub)
        g = coo_matrix((np.ones(len(sub)), (r, c)), shape=(m, m))
        g = g + g.T
    else:
        g = coo_matrix((m, m))
    _, lab = connected_components(g, directed=False)
    sizes = np.bincount(lab)
    return sorted((int(s) for s in sizes if s >= min_size), reverse=True)


def architecture_metrics(df: pd.DataFrame) -> dict:
    """df needs x_um, y_um, deconv_blastemal, dominant_state. Returns macro-shape descriptors."""
    coords = df[["x_um", "y_um"]].to_numpy(float)
    pairs, spacing = _adjacency(coords)
    n = len(df)
    blast = df["deconv_blastemal"].to_numpy(float)
    dom = df["dominant_state"].to_numpy()
    nod = detect_nodules(blast >= 0.5, pairs)
    n_blast = int((blast >= 0.5).sum())
    interface = (float(np.mean(dom[pairs[:, 0]] != dom[pairs[:, 1]])) if len(pairs) else float("nan"))
    return {
        "n_spots": n,
        "spot_spacing_um": round(spacing, 1),
        "morans_I_blastemal": round(morans_i(blast, pairs), 4),   # >0 nodular/segregated, ~0 diffuse
        "n_nodules": len(nod),
        "largest_nodule_spots": nod[0] if nod else 0,
        "largest_nodule_frac": round(nod[0] / n_blast, 4) if (nod and n_blast) else 0.0,
        "median_nodule_spots": int(np.median(nod)) if nod else 0,
        "blastemal_dominant_frac": round(n_blast / n, 4) if n else 0.0,
        "interface_fraction": round(interface, 4),
        # nodular (few large blastemal masses) vs diffuse (many small / high interface)
        "nodularity_index": round((nod[0] / n_blast) * max(morans_i(blast, pairs), 0.0), 4)
                             if (nod and n_blast) else 0.0,
    }


def main() -> None:
    cfg = load_config()
    sig = pd.read_parquet(resolve_path(cfg, "data/processed/spot_signatures.parquet"))
    sig = sig[sig["in_tissue"] == 1].copy()
    spot_um = float(cfg["phase_c"].get("spot_diameter_um", 55.0))
    spatial_root = Path(cfg["paths"]["phase_b"]["spatial_root"])
    if not spatial_root.is_absolute():
        spatial_root = resolve_path(cfg, str(spatial_root))

    # locate each library's dir (sample/SCPCL*_spatial) to get coordinates
    lib_dirs = {p.name.replace("_spatial", ""): p
                for p in spatial_root.glob("SCPCS*/SCPCL*_spatial")}
    rows = []
    for lib_id, g in sig.groupby("library_id"):
        ld = lib_dirs.get(lib_id)
        if ld is None:
            continue
        coords = load_spot_coords_um(ld, spot_um)
        m = g.merge(coords[["x_um", "y_um"]], left_on="barcode", right_index=True, how="inner")
        if len(m) < 20:
            continue
        met = architecture_metrics(m)
        met.update({"sample_id": g["sample_id"].iloc[0], "library_id": lib_id})
        rows.append(met)

    if not rows:
        raise SystemExit("[tissue_architecture] no libraries with coordinates on disk")
    per_lib = pd.DataFrame(rows)
    # per-sample aggregate (spot-weighted mean of the numeric descriptors)
    num = per_lib.select_dtypes("number").columns
    agg = (per_lib.groupby("sample_id")
           .apply(lambda d: pd.Series({c: np.average(d[c], weights=d["n_spots"]) for c in num}),
                  include_groups=False)
           .reset_index())
    out = ensure_dir(resolve_path(cfg, "results/spatial"))
    per_lib.to_csv(out / "tissue_architecture.csv", index=False)
    agg.round(4).to_csv(out / "tissue_architecture_per_sample.csv", index=False)
    print(f"[ok] tissue architecture: {len(per_lib)} libraries, {agg.shape[0]} samples "
          f"-> {out/'tissue_architecture.csv'}")
    cols = ["sample_id", "n_nodules", "largest_nodule_frac", "morans_I_blastemal",
            "interface_fraction", "nodularity_index"]
    print(agg[[c for c in cols if c in agg]].round(3).head(12).to_string(index=False))


if __name__ == "__main__":
    main()
