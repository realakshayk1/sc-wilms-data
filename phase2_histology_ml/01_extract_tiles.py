#!/usr/bin/env python3
"""FR-B1: Extract H&E tiles aligned to Visium spots; Macenko stain-normalize."""

from __future__ import annotations

import argparse
import json

import cv2
import numpy as np

from spatial_utils import (
    _read_subdiagnosis,
    crop_spot_tile,
    discover_libraries,
    load_visium_library,
    macenko_normalize,
    program_fractions,
    score_spot_programs,
    select_tissue_spots,
    spot_id_for,
    estimate_macenko_stain_matrix,
)
from utils import (
    ensure_dir,
    is_demo_mode,
    load_config,
    resolve_path,
    set_seed_logged,
    setup_logging,
)


def make_demo_tiles(out_dir, n_tiles: int = 12, tile_size: int = 128) -> list[dict]:
    tiles_meta = []
    for i in range(n_tiles):
        rng = np.random.default_rng(1000 + i)
        base = rng.integers(180, 230, (tile_size, tile_size, 3), dtype=np.uint8)
        for _ in range(rng.integers(8, 20)):
            cx, cy = rng.integers(10, tile_size - 10, size=2)
            r = int(rng.integers(4, 10))
            cv2.circle(base, (int(cx), int(cy)), r, (80, 60, 120), -1)
        spot_id = f"DEMO_SPOT_{i:03d}"
        path = out_dir / f"{spot_id}.png"
        cv2.imwrite(str(path), cv2.cvtColor(base, cv2.COLOR_RGB2BGR))
        tiles_meta.append(
            {
                "spot_id": spot_id,
                "image_path": str(path.relative_to(out_dir.parent.parent.parent)),
                "sample_id": "DEMO",
                "library_id": "DEMO",
                "barcode": spot_id,
                "array_row": i // 4,
                "array_col": i % 4,
                "dominant_state": ["blastemal", "epithelial", "stromal"][i % 3],
            }
        )
    return tiles_meta


def fit_reference_stain(cfg: dict) -> np.ndarray:
    pb = cfg["phase_b"]
    spatial_root = resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"])
    ref_id = pb.get("reference_library_id")
    libs = discover_libraries(spatial_root)
    ref = next((l for l in libs if l["library_id"] == ref_id), libs[0] if libs else None)
    if ref is None:
        raise FileNotFoundError("No Visium libraries found for stain reference")
    _, image, _, _ = load_visium_library(ref["library_dir"])
    return estimate_macenko_stain_matrix(image)


def extract_real_tiles(cfg: dict, force: bool) -> None:
    pb = cfg["phase_b"]
    seed = set_seed_logged(cfg["features"]["seed"], "extract_tiles")
    spatial_root = resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"])
    out_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    meta_path = out_dir / "tiles_manifest.json"
    ensure_dir(out_dir)

    if meta_path.exists() and not force:
        print(f"[skip] Manifest exists: {meta_path}")
        return

    libs = discover_libraries(
        spatial_root,
        pb.get("library_allowlist"),
        pb.get("max_libraries"),
    )
    if not libs:
        raise SystemExit(f"No Visium libraries under {spatial_root}")

    radius = int(pb["tile_patch_radius_hires_px"])
    target_stain = fit_reference_stain(cfg) if pb.get("stain_normalization") == "macenko" else None

    manifest: list[dict] = []
    signature_rows: list[dict] = []

    for lib in libs:
        lib_dir = lib["library_dir"]
        sample_id = lib["sample_id"]
        library_id = lib["library_id"]
        subdiagnosis = _read_subdiagnosis(spatial_root, sample_id)
        print(f"[load] {library_id} ({sample_id}, {subdiagnosis})")

        adata, image, positions, scalef = load_visium_library(lib_dir)
        scale = float(scalef["tissue_hires_scalef"])
        scores = score_spot_programs(adata, cfg["features"]["features"])
        fracs = program_fractions(scores, pb["signature_programs"])
        spots = select_tissue_spots(
            scores,
            min_umis=int(pb["min_umis_per_spot"]),
            max_spots=pb.get("max_spots_per_library"),
            seed=seed + hash(library_id) % 10_000,
        )

        for barcode in spots:
            sid = spot_id_for(library_id, barcode)
            pos = positions.loc[barcode]
            tile = crop_spot_tile(
                image,
                float(pos["pxl_row_in_fullres"]),
                float(pos["pxl_col_in_fullres"]),
                scale,
                radius,
            )
            if tile is None:
                continue
            if target_stain is not None:
                tile = macenko_normalize(tile, target_stain)

            out_path = out_dir / f"{sid}.png"
            cv2.imwrite(str(out_path), cv2.cvtColor(tile, cv2.COLOR_RGB2BGR))
            dom = str(fracs.loc[barcode, "dominant_state"])
            manifest.append(
                {
                    "spot_id": sid,
                    "image_path": str(out_path.relative_to(cfg["root"])),
                    "sample_id": sample_id,
                    "library_id": library_id,
                    "barcode": barcode,
                    "array_row": int(pos["array_row"]),
                    "array_col": int(pos["array_col"]),
                    "subdiagnosis": subdiagnosis,
                    "dominant_state": dom,
                    "total_counts": int(scores.loc[barcode, "total_counts"]),
                }
            )
            row = {"spot_id": sid, "sample_id": sample_id, "library_id": library_id, "barcode": barcode}
            row.update(scores.loc[barcode].to_dict())
            row.update(fracs.loc[barcode].to_dict())
            signature_rows.append(row)

        print(f"[ok] {library_id}: {len([m for m in manifest if m['library_id'] == library_id])} tiles")

    import pandas as pd

    pd.DataFrame(signature_rows).to_parquet(sig_path, index=False)
    with open(meta_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[ok] {len(manifest)} tiles -> {out_dir}")
    print(f"[ok] Spot signatures -> {sig_path}")


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true", help="Generate synthetic H&E tiles")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    args = parser.parse_args()
    demo = args.demo or is_demo_mode()

    cfg = load_config()
    out_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    ensure_dir(out_dir)
    meta_path = out_dir / "tiles_manifest.json"

    if demo:
        set_seed_logged(cfg["features"]["seed"], "extract_tiles")
        meta = make_demo_tiles(out_dir)
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)
        print(f"[ok] Demo tiles -> {out_dir} ({len(meta)} spots)")
        return

    extract_real_tiles(cfg, force=args.force)


if __name__ == "__main__":
    main()
