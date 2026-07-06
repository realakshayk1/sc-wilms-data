#!/usr/bin/env python3
"""Stage 7: Phase C validation — spatial (now) + emergent (after cluster runs).

Two validators:

  1. SPATIAL (runs on CPU now). Compartment neighbourhood-enrichment z-scores on the real
     Visium (a squidpy-style permutation test, implemented in-house on a spatial kNN graph
     to avoid the squidpy dependency). This is the observed baseline the simulated tissue
     is later compared against.

  2. EMERGENT (needs the cluster runs). A PATIENT-LEVEL contrast of a per-tumor simulation
     QoI (e.g. final tumor size / invasion) between anaplastic vs favorable and relapse vs
     not: Mann-Whitney U across tumors, Cliff's delta effect size, BH-FDR over QoIs. The
     statistics are pure and unit-tested here; `load_sim_qoi` reads PhysiCell output once it
     exists.

Writes results/abm/observed_spatial_qoi.csv now; emergent_validation.csv when sim QoIs are
present.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import (  # noqa: E402
    COMPARTMENTS, discover_library_dirs, ensure_dir, load_config,
    load_spot_coords_um, resolve_path, setup_logging,
)


# --------------------------------------------------------------------------- spatial QoI
def neighbor_enrichment(coords: np.ndarray, labels: np.ndarray, categories: list[str],
                        k: int = 6, n_perm: int = 1000, seed: int = 0) -> pd.DataFrame:
    """Permutation z-score of category adjacency on a spatial kNN graph (squidpy-style).

    Positive z => the two compartments neighbour each other more than under label
    reshuffling; negative => spatial segregation.
    """
    n = len(labels)
    if n <= k:
        return pd.DataFrame(columns=["cat_a", "cat_b", "z", "observed", "n_spots"])
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1)      # +1: first neighbour is the point itself
    edges = {(i, j) for i, row in enumerate(idx[:, 1:]) for j in row}
    edges = np.array([(min(a, b), max(a, b)) for a, b in edges])
    edges = np.unique(edges, axis=0)

    cat_index = {c: i for i, c in enumerate(categories)}
    lab = np.array([cat_index.get(x, -1) for x in labels])
    C = len(categories)
    e0, e1 = edges[:, 0], edges[:, 1]

    def counts(l):
        a, b = l[e0], l[e1]
        ok = (a >= 0) & (b >= 0)
        lo = np.minimum(a[ok], b[ok])
        hi = np.maximum(a[ok], b[ok])
        flat = np.bincount(lo * C + hi, minlength=C * C)
        return flat.reshape(C, C).astype(float)

    obs = counts(lab)
    rng = np.random.default_rng(seed)
    null = np.stack([counts(rng.permutation(lab)) for _ in range(n_perm)])
    mu, sd = null.mean(0), null.std(0)
    z = np.divide(obs - mu, sd, out=np.zeros_like(obs), where=sd > 0)

    rows = []
    for i in range(C):
        for j in range(i, C):
            rows.append({"cat_a": categories[i], "cat_b": categories[j],
                         "z": round(float(z[i, j]), 3), "observed": int(obs[i, j]),
                         "n_spots": n})
    return pd.DataFrame(rows)


def observed_spatial_baseline(cfg) -> pd.DataFrame:
    sig = pd.read_parquet(
        resolve_path(cfg, "data/processed/spot_signatures.parquet"),
        columns=["sample_id", "barcode", "in_tissue", "dominant_state"])
    libs = discover_library_dirs(resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"]))
    dia = float(cfg["phase_c"]["spot"]["diameter_um"])
    seed = int(cfg["phase_c"]["seed"])
    out = []
    for sid, lib in libs.items():
        s = sig[(sig["sample_id"] == sid) & (sig["in_tissue"] == 1)].copy()
        if s.empty:
            continue
        coords = load_spot_coords_um(lib["library_dir"], dia)
        s = s.join(coords[["x_um", "y_um"]], on="barcode").dropna(subset=["x_um", "y_um"])
        if len(s) <= 6:
            continue
        df = neighbor_enrichment(s[["x_um", "y_um"]].to_numpy(),
                                 s["dominant_state"].to_numpy(), COMPARTMENTS, seed=seed)
        df.insert(0, "sample_id", sid)
        out.append(df)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# --------------------------------------------------------------------------- emergent QoI
def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def patient_level_emergent_test(qoi: pd.DataFrame, group_col: str,
                                value_cols: list[str]) -> pd.DataFrame:
    """Mann-Whitney U across tumors + Cliff's delta + BH-FDR over the QoIs tested.

    qoi has one row per tumor (the unit of inference is the patient, never the agent).
    group_col is binary (e.g. 1=anaplastic/relapse, 0=other).
    """
    from statsmodels.stats.multitest import multipletests
    rows = []
    g = qoi.dropna(subset=[group_col])
    pos, neg = g[g[group_col] == 1], g[g[group_col] == 0]
    for col in value_cols:
        a = pos[col].dropna().to_numpy()
        b = neg[col].dropna().to_numpy()
        if len(a) < 2 or len(b) < 2:
            rows.append({"qoi": col, "n_pos": len(a), "n_neg": len(b),
                         "u": np.nan, "p_value": np.nan, "cliffs_delta": np.nan})
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"qoi": col, "n_pos": len(a), "n_neg": len(b),
                     "median_pos": float(np.median(a)), "median_neg": float(np.median(b)),
                     "u": float(u), "p_value": float(p),
                     "cliffs_delta": round(cliffs_delta(a, b), 3)})
    res = pd.DataFrame(rows)
    ok = res["p_value"].notna()
    res["p_bh"] = np.nan
    if ok.any():
        res.loc[ok, "p_bh"] = multipletests(res.loc[ok, "p_value"], method="fdr_bh")[1]
    return res


def load_sim_qoi(cfg) -> pd.DataFrame | None:
    """Per-tumor simulation QoIs from PhysiCell output (cluster). Returns None until runs
    exist. Expected: results/abm/<sample_id>/output/qoi.csv written by 06_run_cohort."""
    out_dir = resolve_path(cfg, "results/abm")
    frames = []
    for q in out_dir.glob("SCPCS*/output/qoi.csv"):
        d = pd.read_csv(q)
        d["sample_id"] = q.parent.parent.name
        frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else None


def main() -> None:
    setup_logging()
    cfg = load_config()
    out_dir = ensure_dir(resolve_path(cfg, "results/abm"))

    # 1. observed spatial baseline (CPU, now)
    spatial = observed_spatial_baseline(cfg)
    if not spatial.empty:
        spatial.to_csv(out_dir / "observed_spatial_qoi.csv", index=False)
        seg = spatial[spatial.cat_a != spatial.cat_b].groupby(["cat_a", "cat_b"])["z"].mean()
        print(f"[ok] observed spatial QoI -> {out_dir/'observed_spatial_qoi.csv'} "
              f"({spatial.sample_id.nunique()} tumors)")
        print("[info] mean cross-compartment enrichment z (negative = segregation):")
        print(seg.round(2).to_string())

    # 2. emergent test (only if sim QoIs exist)
    sim = load_sim_qoi(cfg)
    if sim is None:
        print("[pending] no simulation QoIs yet — emergent test runs after the cluster cohort.")
        return
    pt = pd.read_csv(resolve_path(cfg, "results/mechanotypes/per_tumor_scores.csv"))
    pt["anaplastic"] = (pt["subdiagnosis"].str.lower() == "anaplastic").astype(int)
    merged = sim.merge(pt[["sample_id", "anaplastic", "relapse"]], on="sample_id", how="left")
    value_cols = [c for c in sim.columns if c != "sample_id"]
    for grp in ["anaplastic", "relapse"]:
        res = patient_level_emergent_test(merged, grp, value_cols)
        res.to_csv(out_dir / f"emergent_validation_{grp}.csv", index=False)
        print(f"[ok] emergent test ({grp}) -> emergent_validation_{grp}.csv")
        print(res.to_string(index=False))


if __name__ == "__main__":
    main()
