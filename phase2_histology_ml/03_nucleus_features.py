#!/usr/bin/env python3
"""FR-B3/B4: Per-nucleus morphology + weak labels from spot transcriptomic dominant state."""

from __future__ import annotations

import argparse
import json

import cv2
import numpy as np
import pandas as pd
from skimage.measure import regionprops

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

FEATURE_COLS = [
    "area",
    "eccentricity",
    "solidity",
    "major_axis_length",
    "texture_var",
    "hematoxylin_intensity",
    "neighbor_density",
]


def hematoxylin_intensity(rgb: np.ndarray, mask: np.ndarray) -> float:
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    return float(gray[mask].mean()) if mask.any() else 0.0


def neighbor_density(centroids: np.ndarray, idx: int, radius: float = 25.0) -> float:
    if len(centroids) <= 1:
        return 0.0
    dists = np.linalg.norm(centroids - centroids[idx], axis=1)
    return float(np.sum((dists > 0) & (dists < radius)))


def load_weak_labels(cfg: dict, manifest: list[dict]) -> dict[str, str]:
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    if sig_path.exists():
        sig = pd.read_parquet(sig_path).set_index("spot_id")
        return {sid: str(sig.loc[sid, "dominant_state"]) for sid in sig.index if sid in {m["spot_id"] for m in manifest}}
    return {m["spot_id"]: m.get("dominant_state", "blastemal") for m in manifest}


def extract_nucleus_features(
    rgb: np.ndarray,
    labels: np.ndarray,
    spot_id: str,
    weak_label: str,
    meta: dict,
) -> list[dict]:
    props = regionprops(labels, intensity_image=cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY))
    centroids = np.array([p.centroid for p in props]) if props else np.empty((0, 2))
    rows = []
    for i, p in enumerate(props):
        if p.area < 20:
            continue
        mask = labels == p.label
        rows.append(
            {
                "nucleus_id": f"{spot_id}_n{i}",
                "spot_id": spot_id,
                "sample_id": meta.get("sample_id", ""),
                "library_id": meta.get("library_id", ""),
                "subdiagnosis": meta.get("subdiagnosis", ""),
                "area": float(p.area),
                "eccentricity": float(p.eccentricity),
                "solidity": float(p.solidity),
                "major_axis_length": float(p.major_axis_length),
                "texture_var": float(np.var(rgb[mask])),
                "hematoxylin_intensity": hematoxylin_intensity(rgb, mask),
                "neighbor_density": neighbor_density(centroids, i),
                "weak_label": weak_label,
            }
        )
    return rows


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    set_seed_logged(cfg["features"]["seed"], "nucleus_features")

    nuclei_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["nuclei_dir"])
    mask_dir = nuclei_dir / "masks"
    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    out_parquet = resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"])
    ensure_dir(out_parquet.parent)

    if out_parquet.exists() and not args.force:
        print(f"[skip] Features exist: {out_parquet}")
        return

    with open(tiles_dir / "tiles_manifest.json") as f:
        manifest = json.load(f)
    weak_labels = load_weak_labels(cfg, manifest)
    meta_by_spot = {m["spot_id"]: m for m in manifest}

    all_rows: list[dict] = []
    for entry in manifest:
        spot_id = entry["spot_id"]
        mask_path = mask_dir / f"{spot_id}_labels.npy"
        img_path = cfg["root"] / entry["image_path"]
        if not mask_path.exists() or not img_path.exists():
            continue
        labels = np.load(mask_path)
        bgr = cv2.imread(str(img_path))
        rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
        all_rows.extend(
            extract_nucleus_features(
                rgb,
                labels,
                spot_id,
                weak_labels.get(spot_id, entry.get("dominant_state", "blastemal")),
                meta_by_spot.get(spot_id, entry),
            )
        )

    df = pd.DataFrame(all_rows)
    df.to_parquet(out_parquet, index=False)
    print(f"[ok] {len(df)} nuclei features -> {out_parquet}")


if __name__ == "__main__":
    main()
