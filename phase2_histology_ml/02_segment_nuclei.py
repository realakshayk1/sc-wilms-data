#!/usr/bin/env python3
"""FR-B2: Segment nuclei per tile (StarDist default; Cellpose/watershed fallback)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cv2
import numpy as np
from scipy import ndimage as ndi
from skimage import morphology, segmentation
from skimage.color import rgb2hed
from skimage.filters import threshold_otsu

from utils import ensure_dir, is_demo_mode, load_config, resolve_path, set_seed_logged, setup_logging

_STARDIST_MODEL = None
_CELLPOSE_MODEL = None


def get_stardist_model():
    global _STARDIST_MODEL
    if _STARDIST_MODEL is None:
        from stardist.models import StarDist2D

        _STARDIST_MODEL = StarDist2D.from_pretrained("2D_versatile_he")
    return _STARDIST_MODEL


def get_cellpose_model():
    global _CELLPOSE_MODEL
    if _CELLPOSE_MODEL is None:
        from cellpose import models

        # Nuclei model — appropriate for H&E spot tiles; load once per run
        _CELLPOSE_MODEL = models.CellposeModel(gpu=False, model_type="nuclei")
    return _CELLPOSE_MODEL


def segment_stardist(rgb: np.ndarray) -> np.ndarray:
    model = get_stardist_model()
    labels, _ = model.predict_instances(rgb, n_tiles=4, show_tile_progress=False)
    return labels.astype(np.int32)


def segment_cellpose(rgb: np.ndarray) -> np.ndarray:
    model = get_cellpose_model()
    masks, _, _, _ = model.eval(rgb, diameter=25, channel_axis=-1)
    return masks.astype(np.int32)


def segment_watershed(rgb: np.ndarray) -> np.ndarray:
    """Hematoxylin-channel watershed — standard classical nuclei segmentation."""
    h = rgb2hed(rgb)[:, :, 0]
    thresh = threshold_otsu(h)
    binary = h < thresh
    binary = morphology.opening(binary, morphology.disk(1))
    binary = morphology.remove_small_objects(binary, max_size=15)
    distance = ndi.distance_transform_edt(binary)
    local_max = morphology.local_maxima(distance)
    markers = morphology.label(local_max)
    if markers.max() == 0:
        markers, _ = cv2.connectedComponents(binary.astype(np.uint8))
        return markers.astype(np.int32)
    labels = segmentation.watershed(-distance, markers, mask=binary)
    return labels.astype(np.int32)


def segment_tile(rgb: np.ndarray, method: str = "auto") -> tuple[np.ndarray, str]:
    if method in ("auto", "stardist"):
        try:
            return segment_stardist(rgb), "stardist"
        except Exception as exc:
            if method == "stardist":
                raise
            print(f"[warn] StarDist unavailable ({exc}); trying Cellpose nuclei")
    if method in ("auto", "cellpose"):
        try:
            return segment_cellpose(rgb), "cellpose"
        except Exception as exc:
            if method == "cellpose":
                raise
            print(f"[warn] Cellpose failed ({exc}); using hematoxylin watershed")
    if method in ("auto", "watershed", "demo_threshold"):
        return segment_watershed(rgb), "watershed"
    raise ValueError(f"Unknown segmentation method: {method}")


def save_overlay(rgb: np.ndarray, labels: np.ndarray, out_path: Path) -> None:
    overlay = rgb.copy()
    for lab in np.unique(labels):
        if lab == 0:
            continue
        mask = labels == lab
        ys, xs = np.where(mask)
        cy, cx = int(ys.mean()), int(xs.mean())
        cv2.circle(overlay, (cx, cy), 3, (0, 255, 0), 1)
    cv2.imwrite(str(out_path), cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR))


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument(
        "--method",
        choices=["auto", "stardist", "cellpose", "watershed", "demo_threshold"],
        default="auto",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    demo = args.demo or is_demo_mode()

    cfg = load_config()
    pb = cfg.get("phase_b", {})
    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    nuclei_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["nuclei_dir"])
    overlay_dir = ensure_dir(nuclei_dir / "overlays")
    mask_dir = ensure_dir(nuclei_dir / "masks")
    log_path = nuclei_dir / "segmentation_log.json"

    if log_path.exists() and not args.force and not demo:
        print(f"[skip] Segmentation log exists: {log_path}")
        return

    set_seed_logged(cfg["features"]["seed"], "segment_nuclei")

    manifest_path = tiles_dir / "tiles_manifest.json"
    if not manifest_path.exists():
        raise SystemExit("Run 01_extract_tiles.py first")

    with open(manifest_path) as f:
        manifest = json.load(f)

    method = pb.get("segmentation_method", args.method)
    if demo and args.method == "auto":
        method = "watershed"

    seg_log = []
    for i, entry in enumerate(manifest):
        img_path = cfg["root"] / entry["image_path"] if "image_path" in entry else tiles_dir / f"{entry['spot_id']}.png"
        bgr = cv2.imread(str(img_path))
        if bgr is None:
            print(f"[warn] Could not read {img_path}")
            continue
        rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
        labels, used = segment_tile(rgb, method=method if method != "auto" else args.method)

        spot_id = entry["spot_id"]
        np.save(mask_dir / f"{spot_id}_labels.npy", labels)
        save_overlay(rgb, labels, overlay_dir / f"{spot_id}_overlay.png")
        n_nuclei = len(np.unique(labels)) - 1
        seg_log.append({"spot_id": spot_id, "method": used, "n_nuclei": n_nuclei})
        if (i + 1) % 50 == 0 or i == 0:
            print(f"[ok] {spot_id}: {n_nuclei} nuclei ({used}) [{i + 1}/{len(manifest)}]")

    with open(nuclei_dir / "segmentation_log.json", "w") as f:
        json.dump(seg_log, f, indent=2)
    print(f"[ok] Segmented {len(seg_log)} tiles -> {overlay_dir}")


if __name__ == "__main__":
    main()
