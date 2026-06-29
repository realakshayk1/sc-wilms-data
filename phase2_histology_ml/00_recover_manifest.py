#!/usr/bin/env python3
"""Rebuild tiles_manifest.json from spot_signatures when PNGs/masks exist but manifest was lost."""

from __future__ import annotations

import argparse
import json

import pandas as pd

from spatial_utils import _read_subdiagnosis
from utils import load_config, resolve_path, setup_logging


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    manifest_path = tiles_dir / "tiles_manifest.json"
    spatial_root = resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"])

    if manifest_path.exists() and not args.force:
        with open(manifest_path) as f:
            n = len(json.load(f))
        if n > 1000:
            print(f"[skip] Manifest has {n} entries: {manifest_path}")
            return

    if not sig_path.exists():
        raise SystemExit(f"Missing {sig_path} — run 01_extract_tiles.py first")

    sig = pd.read_parquet(sig_path)
    sub_cache: dict[str, str] = {}
    manifest = []
    for row in sig.itertuples(index=False):
        sid = str(row.spot_id)
        sample_id = str(row.sample_id)
        if sample_id not in sub_cache:
            sub_cache[sample_id] = _read_subdiagnosis(spatial_root, sample_id)
        img = tiles_dir / f"{sid}.png"
        if not img.exists():
            continue
        manifest.append(
            {
                "spot_id": sid,
                "image_path": str(img.relative_to(cfg["root"])),
                "sample_id": sample_id,
                "library_id": str(row.library_id),
                "barcode": str(row.barcode),
                "subdiagnosis": sub_cache[sample_id],
                "dominant_state": str(row.dominant_state),
                "total_counts": int(row.total_counts),
            }
        )

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[ok] Recovered manifest ({len(manifest)} spots) -> {manifest_path}")


if __name__ == "__main__":
    main()
