"""Shared utilities for Phase C (histology/omics -> PhysiCell ABM).

Config loading, per-library Visium coordinate loading with a pixel->micron affine
transform, and seeded RNG helpers. Keeps Phase C decoupled from Phase B while reusing
the same repo-root / paths.yaml convention.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

LOG = logging.getLogger("wilms_abm_phase_c")

COMPARTMENTS = ["blastemal", "epithelial", "stromal"]


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (Path.cwd().resolve(), here.parent, here.parent.parent):
        if (candidate / "config" / "paths.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate repo root (config/paths.yaml)")


def load_config() -> dict[str, Any]:
    root = repo_root()
    paths = yaml.safe_load((root / "config" / "paths.yaml").read_text())
    phase_c = yaml.safe_load((root / "config" / "phase_c.yaml").read_text())
    physicell = yaml.safe_load((root / "config" / "physicell.yaml").read_text())
    return {"root": root, "paths": paths, "phase_c": phase_c, "physicell": physicell}


def resolve_path(cfg: dict[str, Any], rel: str) -> Path:
    return cfg["root"] / rel


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def rng_for(seed: int, label: str = "phase_c") -> np.random.Generator:
    LOG.info("[seed] %s = %d", label, seed)
    return np.random.default_rng(seed)


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def safe_barcode(barcode: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(barcode))


def um_per_pixel(scalefactors: dict[str, Any], spot_diameter_um: float) -> float:
    """Microns per full-resolution pixel from the Visium scalefactors.

    `spot_diameter_fullres` is the capture-spot diameter in full-res pixels; dividing the
    known physical spot diameter by it yields the isotropic micron/pixel scale.
    """
    dia_px = float(scalefactors["spot_diameter_fullres"])
    if dia_px <= 0:
        raise ValueError(f"non-positive spot_diameter_fullres: {dia_px}")
    return spot_diameter_um / dia_px


def load_spot_coords_um(
    library_dir: Path | str,
    spot_diameter_um: float,
) -> pd.DataFrame:
    """Per-barcode full-res pixel coordinates converted to microns.

    Returns a frame indexed by barcode with columns x_um, y_um, in_tissue. The micron
    frame is anchored so the tumor's own coordinates start near the origin downstream
    (recentring happens at placement, against tissue spots only).
    """
    library_dir = Path(library_dir)
    spatial = library_dir / "spatial"
    positions = pd.read_csv(
        spatial / "tissue_positions_list.csv",
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col",
               "pxl_col_in_fullres", "pxl_row_in_fullres"],
    )
    scalefactors = json.loads((spatial / "scalefactors_json.json").read_text())
    mpp = um_per_pixel(scalefactors, spot_diameter_um)
    out = pd.DataFrame({
        "barcode": positions["barcode"].astype(str),
        "in_tissue": positions["in_tissue"].astype(int),
        # PhysiCell x is horizontal (columns), y vertical (rows)
        "x_um": positions["pxl_col_in_fullres"].to_numpy(float) * mpp,
        "y_um": positions["pxl_row_in_fullres"].to_numpy(float) * mpp,
    }).set_index("barcode")
    out.attrs["um_per_pixel"] = mpp
    return out


def discover_library_dirs(spatial_root: Path) -> dict[str, dict[str, str]]:
    """Map sample_id -> {library_id, library_dir} for the Visium libraries on disk."""
    libs: dict[str, dict[str, str]] = {}
    if not spatial_root.exists():
        return libs
    for sample_dir in sorted(spatial_root.glob("SCPCS*")):
        if not sample_dir.is_dir():
            continue
        for lib_dir in sorted(sample_dir.glob("SCPCL*_spatial")):
            libs[sample_dir.name] = {
                "library_id": lib_dir.name.replace("_spatial", ""),
                "library_dir": str(lib_dir),
            }
    return libs
