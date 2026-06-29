#!/usr/bin/env python3
"""Validate H&E spot fractions against independent marker-gene Visium deconvolution."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr

from spatial_utils import CELL_STATES, discover_libraries, load_visium_library, marker_fractions, spot_id_for
from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging


def marker_map_from_features(features: list[dict]) -> dict[str, list[str]]:
    by_id = {f["id"]: f for f in features}
    return {
        "blastemal": by_id["blastemal_program"]["genes_positive"],
        "epithelial": by_id["epithelial_program"]["genes_positive"],
        "stromal": by_id["stromal_program"]["genes_positive"],
    }


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "marker_deconv")
    out_json = resolve_path(cfg, cfg["paths"]["phase_b"]["marker_deconv_json"])
    ensure_dir(out_json.parent)

    if out_json.exists() and not args.force:
        print(f"[skip] Marker deconv exists: {out_json}")
        return

    frac_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_fractions_csv"])
    manifest_path = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"]) / "tiles_manifest.json"
    spatial_root = resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"])

    if not frac_path.exists() or not manifest_path.exists():
        raise SystemExit("Run Phase B pipeline through 05_spot_fractions.py first")

    frac_df = pd.read_csv(frac_path)
    with open(manifest_path) as f:
        manifest = json.load(f)

    markers = marker_map_from_features(cfg["features"]["features"])
    by_lib: dict[str, list[dict]] = {}
    for entry in manifest:
        by_lib.setdefault(entry["library_id"], []).append(entry)

    marker_rows = []
    libs = {l["library_id"]: l for l in discover_libraries(spatial_root)}

    for lib_id, entries in by_lib.items():
        if lib_id not in libs:
            continue
        adata, _, _, _ = load_visium_library(libs[lib_id]["library_dir"])
        mfrac = marker_fractions(adata, markers)
        for entry in entries:
            barcode = entry.get("barcode")
            if barcode not in mfrac.index and "_" in entry.get("spot_id", ""):
                # Visium barcodes use hyphens; spot_id sanitizes to underscores
                barcode = barcode.replace("_", "-") if barcode else None
            if barcode not in mfrac.index:
                continue
            row = {"spot_id": entry["spot_id"]}
            for st in CELL_STATES:
                row[f"marker_{st}"] = float(mfrac.loc[barcode, f"marker_{st}"])
            marker_rows.append(row)

    if not marker_rows:
        raise SystemExit("No marker fractions computed — check manifest barcodes")

    marker_df = pd.DataFrame(marker_rows)
    merged = frac_df.merge(marker_df, on="spot_id", how="inner")

    correlations = {}
    for st in CELL_STATES:
        r, p = pearsonr(merged[f"frac_{st}"], merged[f"marker_{st}"])
        rho, _ = spearmanr(merged[f"frac_{st}"], merged[f"marker_{st}"])
        correlations[st] = {
            "pearson_r": float(r),
            "pearson_p": float(p),
            "spearman_rho": float(rho),
        }

    dom_he = merged[[f"frac_{s}" for s in CELL_STATES]].idxmax(axis=1).str.replace("frac_", "")
    dom_marker = merged[[f"marker_{s}" for s in CELL_STATES]].idxmax(axis=1).str.replace("marker_", "")
    agreement = float((dom_he == dom_marker).mean())

    out = {
        "validation_method": "Independent marker-gene softmax (positive genes only per compartment)",
        "n_spots": int(len(merged)),
        "dominant_state_agreement": agreement,
        "correlations": correlations,
        "seed": seed,
    }
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)

    print(f"[ok] Marker deconv validation -> {out_json} ({len(merged)} spots)")
    for st, v in correlations.items():
        print(f"     {st}: r={v['pearson_r']:.3f}")


if __name__ == "__main__":
    main()
