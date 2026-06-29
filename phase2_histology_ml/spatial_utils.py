"""Visium spatial I/O, Phase A program scoring, and H&E tile utilities."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import cv2
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse

CELL_STATES = ["blastemal", "epithelial", "stromal"]


def load_phase_b_config(root: Path) -> dict[str, Any]:
    with open(root / "config" / "phase_b.yaml") as f:
        return yaml.safe_load(f)


def discover_libraries(spatial_root: Path, allowlist: list[str] | None = None, max_libraries: int | None = None) -> list[dict[str, str]]:
    libs: list[dict[str, str]] = []
    if not spatial_root.exists():
        return libs
    for sample_dir in sorted(spatial_root.glob("SCPCS*")):
        if not sample_dir.is_dir():
            continue
        sample_id = sample_dir.name
        for lib_dir in sorted(sample_dir.glob("SCPCL*_spatial")):
            lib_id = lib_dir.name.replace("_spatial", "")
            if allowlist and lib_id not in allowlist and sample_id not in allowlist:
                continue
            libs.append(
                {
                    "sample_id": sample_id,
                    "library_id": lib_id,
                    "library_dir": str(lib_dir),
                }
            )
    if max_libraries is not None and len(libs) > max_libraries:
        libs = libs[:max_libraries]
    return libs


def _read_spatial_metadata(lib_dir: Path) -> dict[str, Any]:
    meta_path = lib_dir / f"{lib_dir.name.replace('_spatial', '')}_metadata.json"
    if meta_path.exists():
        with open(meta_path) as f:
            return json.load(f)
    return {}


def _read_subdiagnosis(spatial_root: Path, sample_id: str) -> str:
    meta_tsv = spatial_root / "spatial_metadata.tsv"
    if not meta_tsv.exists():
        return "unknown"
    df = pd.read_csv(meta_tsv, sep="\t")
    sub = df.loc[df["scpca_sample_id"] == sample_id, "subdiagnosis"]
    if sub.empty:
        return "unknown"
    val = str(sub.iloc[0]).lower()
    if "anaplas" in val:
        return "anaplastic"
    if "favor" in val:
        return "favorable"
    return val


def load_visium_library(lib_dir: Path | str) -> tuple[Any, np.ndarray, pd.DataFrame, dict[str, float]]:
    lib_dir = Path(lib_dir)
    mtx_dir = lib_dir / "filtered_feature_bc_matrix"
    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
    adata.var_names_make_unique()

    spatial_dir = lib_dir / "spatial"
    positions = pd.read_csv(
        spatial_dir / "tissue_positions_list.csv",
        header=None,
        names=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ],
    )
    positions = positions.set_index("barcode")
    adata.obs = adata.obs.join(positions, how="left")
    adata.obs["in_tissue"] = adata.obs["in_tissue"].fillna(0).astype(int)

    with open(spatial_dir / "scalefactors_json.json") as f:
        scalefactors = json.load(f)

    hires_path = spatial_dir / "tissue_hires_image.png"
    if not hires_path.exists():
        hires_path = spatial_dir / "tissue_hires_image.jpg"
    bgr = cv2.imread(str(hires_path))
    if bgr is None:
        raise FileNotFoundError(f"Could not read H&E image: {hires_path}")
    image = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)

    return adata, image, positions, scalefactors


def build_gene_lookup(var_names: np.ndarray, var_df: pd.DataFrame | None = None) -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = {}
    symbols = var_names
    if var_df is not None and "gene_symbols" in var_df.columns:
        symbols = var_df["gene_symbols"].astype(str).values
    for sym, gid in zip(symbols, var_names, strict=False):
        if sym and sym != "nan":
            lookup.setdefault(str(sym).upper(), []).append(str(gid))
    return lookup


def resolve_genes(gene_names: list[str], lookup: dict[str, list[str]], var_names: np.ndarray) -> list[str]:
    hits: list[str] = []
    index = set(var_names)
    for g in gene_names:
        g = str(g).upper()
        if g in lookup:
            hits.extend(lookup[g])
        elif g in index:
            hits.append(g)
    return list(dict.fromkeys(hits))


def score_feature_matrix(
    X: sparse.spmatrix | np.ndarray,
    var_names: np.ndarray,
    genes_pos: list[str],
    genes_neg: list[str],
    lookup: dict[str, list[str]],
) -> np.ndarray:
    pos = resolve_genes(genes_pos, lookup, var_names)
    neg = resolve_genes(genes_neg, lookup, var_names)
    var_index = {g: i for i, g in enumerate(var_names)}
    pos_idx = [var_index[g] for g in pos if g in var_index]
    neg_idx = [var_index[g] for g in neg if g in var_index]

    if sparse.issparse(X):
        lib = np.asarray(X.sum(axis=1)).ravel().astype(float)
    else:
        lib = X.sum(axis=1).astype(float)
    lib[lib == 0] = 1.0
    med = float(np.median(lib))

    def _col_sums(idxs: list[int]) -> np.ndarray:
        if not idxs:
            return np.zeros(X.shape[0])
        if sparse.issparse(X):
            return np.asarray(X[:, idxs].sum(axis=1)).ravel()
        return X[:, idxs].sum(axis=1)

    pos_score = _col_sums(pos_idx) / lib * med
    neg_score = _col_sums(neg_idx) / lib * med
    return np.log1p(pos_score) - np.log1p(neg_score)


def score_spot_programs(adata: Any, features_cfg: list[dict[str, Any]]) -> pd.DataFrame:
    lookup = build_gene_lookup(adata.var_names.to_numpy(), adata.var)
    X = adata.X
    rows = {}
    for feat in features_cfg:
        fid = feat["id"]
        rows[fid] = score_feature_matrix(
            X,
            adata.var_names.to_numpy(),
            feat.get("genes_positive", []),
            feat.get("genes_negative", []),
            lookup,
        )
    score_df = pd.DataFrame(rows, index=adata.obs_names)
    score_df["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    score_df["in_tissue"] = adata.obs["in_tissue"].to_numpy()
    return score_df


def program_fractions(
    score_df: pd.DataFrame,
    program_map: dict[str, str],
) -> pd.DataFrame:
    cols = [program_map[s] for s in CELL_STATES if s in program_map]
    mat = score_df[cols].to_numpy(dtype=float)
    mat = mat - mat.max(axis=1, keepdims=True)
    exp = np.exp(mat)
    frac = exp / exp.sum(axis=1, keepdims=True)
    out = pd.DataFrame(frac, index=score_df.index, columns=[f"deconv_{s}" for s in CELL_STATES])
    dominant_idx = np.argmax(frac, axis=1)
    out["dominant_state"] = [CELL_STATES[i] for i in dominant_idx]
    return out


def estimate_macenko_stain_matrix(rgb: np.ndarray, beta: float = 0.15, alpha: float = 1.0) -> np.ndarray:
    """Return 2x3 stain matrix (rows = H/E basis vectors in OD space)."""
    img = rgb.reshape(-1, 3).astype(np.float64)
    od = -np.log((img + 1.0) / 256.0)
    od = od[~np.any(od < beta, axis=1)]
    if len(od) < 100:
        od = -np.log((rgb.reshape(-1, 3).astype(np.float64) + 1.0) / 256.0)

    cov = np.cov(od.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    v1 = eigvecs[:, int(np.argmax(eigvals))]
    v2 = eigvecs[:, int(np.argsort(eigvals)[-2])]
    if v1[0] < v2[0]:
        v1, v2 = v2, v1

    proj = od @ np.column_stack([v1, v2])
    phi = np.arctan2(proj[:, 1], proj[:, 0])
    min_phi = np.percentile(phi, alpha)
    max_phi = np.percentile(phi, 100 - alpha)
    h_mask = (phi >= min_phi) & (phi <= max_phi)
    e_mask = ~h_mask

    stain = np.zeros((2, 3))
    if h_mask.any():
        stain[0] = np.median(od[h_mask], axis=0)
    if e_mask.any():
        stain[1] = np.median(od[e_mask], axis=0)
    stain = stain / np.maximum(np.linalg.norm(stain, axis=1, keepdims=True), 1e-6)
    return stain


def macenko_normalize(rgb: np.ndarray, target_stain: np.ndarray) -> np.ndarray:
    """Normalize tile RGB to a reference stain matrix (Macenko et al.)."""
    source_stain = estimate_macenko_stain_matrix(rgb)
    od = -np.log((rgb.reshape(-1, 3).astype(np.float64) + 1.0) / 256.0)
    conc = np.linalg.lstsq(source_stain.T, od.T, rcond=None)[0]
    od_norm = (target_stain.T @ conc).T
    norm = 256.0 * np.exp(-od_norm) - 1.0
    norm = np.clip(norm, 0, 255).reshape(rgb.shape).astype(np.uint8)
    return norm


def crop_spot_tile(
    image: np.ndarray,
    row_fullres: float,
    col_fullres: float,
    scalef: float,
    radius_px: int,
) -> np.ndarray | None:
    row = int(row_fullres * scalef)
    col = int(col_fullres * scalef)
    h, w = image.shape[:2]
    y0, y1 = max(0, row - radius_px), min(h, row + radius_px)
    x0, x1 = max(0, col - radius_px), min(w, col + radius_px)
    if y1 - y0 < 8 or x1 - x0 < 8:
        return None
    tile = image[y0:y1, x0:x1]
    return cv2.resize(tile, (2 * radius_px, 2 * radius_px), interpolation=cv2.INTER_LINEAR)


def select_tissue_spots(
    score_df: pd.DataFrame,
    min_umis: int,
    max_spots: int | None,
    seed: int,
) -> pd.Index:
    tissue = score_df.index[(score_df["in_tissue"] == 1) & (score_df["total_counts"] >= min_umis)]
    if max_spots is not None and len(tissue) > max_spots:
        rng = np.random.default_rng(seed)
        tissue = pd.Index(rng.choice(tissue, size=max_spots, replace=False))
    return tissue


def marker_fractions(adata: Any, marker_map: dict[str, list[str]]) -> pd.DataFrame:
    """Independent marker-gene softmax fractions (not Phase A program scores)."""
    lookup = build_gene_lookup(adata.var_names.to_numpy(), adata.var)
    X = adata.X
    if sparse.issparse(X):
        lib = np.asarray(X.sum(axis=1)).ravel().astype(float)
    else:
        lib = X.sum(axis=1).astype(float)
    lib[lib == 0] = 1.0
    med = float(np.median(lib))
    var_names = adata.var_names.to_numpy()
    var_index = {g: i for i, g in enumerate(var_names)}

    scores = {}
    for state, genes in marker_map.items():
        hits = resolve_genes(genes, lookup, var_names)
        idx = [var_index[g] for g in hits if g in var_index]
        if not idx:
            scores[state] = np.zeros(adata.n_obs)
            continue
        if sparse.issparse(X):
            raw = np.asarray(X[:, idx].sum(axis=1)).ravel() / lib * med
        else:
            raw = X[:, idx].sum(axis=1) / lib * med
        scores[state] = np.log1p(raw)

    mat = np.column_stack([scores[s] for s in CELL_STATES if s in scores])
    mat = mat - mat.max(axis=1, keepdims=True)
    exp = np.exp(mat)
    frac = exp / exp.sum(axis=1, keepdims=True)
    out = pd.DataFrame(frac, index=adata.obs_names, columns=[f"marker_{s}" for s in CELL_STATES])
    out["marker_dominant"] = [CELL_STATES[i] for i in np.argmax(frac, axis=1)]
    return out


def spot_id_for(library_id: str, barcode: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9]+", "_", barcode)
    return f"{library_id}__{safe}"
