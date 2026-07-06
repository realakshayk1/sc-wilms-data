#!/usr/bin/env python3
"""Stage 7: Phase C validation — spatial (now) + emergent (after cluster runs).

Two validators:

  1. SPATIAL (runs on CPU now). Compartment architecture on the real Visium: a
     neighbourhood-enrichment z-score (permutation test on a kNN graph) and a co-occurrence
     curve (co-location enrichment vs distance). Both are implemented in-house (scipy) so
     the default install needs no extra dependency; with squidpy importable (backend auto /
     squidpy) it additionally emits squidpy's nhood z-scores and Ripley's L as a cross-check.
     This is the observed baseline the simulated tissue is later compared against.

  2. EMERGENT (needs the cluster runs). A PATIENT-LEVEL contrast of a per-tumor simulation
     QoI (e.g. final tumor size / invasion) between anaplastic vs favorable and relapse vs
     not: Mann-Whitney U across tumors, Cliff's delta effect size, BH-FDR over QoIs. The
     statistics are pure and unit-tested here; `load_sim_qoi` reads PhysiCell output once it
     exists.

Writes results/abm/observed_spatial_qoi.csv + observed_cooccurrence.csv now (+ squidpy
Ripley/nhood if available); emergent_validation_*.csv when sim QoIs are present.
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


def co_occurrence(coords: np.ndarray, labels: np.ndarray, categories: list[str],
                  radii: list[float]) -> pd.DataFrame:
    """Co-occurrence enrichment vs distance (squidpy `co_occurrence` semantics).

    For each distance ring and ordered pair (cond, exp): score =
    P(a neighbour in the ring is `exp` | centre is `cond`) / P(`exp` overall). >1 means the
    two compartments co-locate at that distance; <1 means they segregate. Efficient via
    cKDTree.count_neighbors over cumulative radii.
    """
    cat_index = {c: i for i, c in enumerate(categories)}
    lab = np.array([cat_index.get(x, -1) for x in labels])
    keep = lab >= 0
    coords, lab = coords[keep], lab[keep]
    C, N = len(categories), len(lab)
    R = np.asarray(radii, float)
    counts_i = np.array([(lab == i).sum() for i in range(C)])
    if N == 0 or (counts_i > 0).sum() < 2 or len(R) < 2:
        return pd.DataFrame(columns=["cond", "exp", "r_um", "cooccur", "n_spots"])
    global_frac = counts_i / N
    trees = [cKDTree(coords[lab == i]) if counts_i[i] else None for i in range(C)]

    cum = np.zeros((C, C, len(R)))
    for a in range(C):
        if trees[a] is None:
            continue
        for b in range(C):
            if trees[b] is None:
                continue
            cn = np.asarray(trees[a].count_neighbors(trees[b], R), float)
            if a == b:
                cn = cn - counts_i[a]                 # drop self-pairs (distance 0)
            cum[a, b] = cn
    ring = np.diff(cum, axis=2)                        # pairs per ring
    rmid = 0.5 * (R[1:] + R[:-1])
    rows = []
    for a in range(C):
        ring_a_total = ring[a].sum(axis=0)            # a's neighbours in each ring
        for b in range(C):
            frac = np.divide(ring[a, b], ring_a_total,
                             out=np.zeros_like(ring_a_total), where=ring_a_total > 0)
            score = frac / global_frac[b] if global_frac[b] > 0 else np.full_like(frac, np.nan)
            for k, r in enumerate(rmid):
                rows.append({"cond": categories[a], "exp": categories[b],
                             "r_um": round(float(r), 1), "cooccur": round(float(score[k]), 3),
                             "n_spots": N})
    return pd.DataFrame(rows)


def has_squidpy() -> bool:
    import importlib.util
    return importlib.util.find_spec("squidpy") is not None


def squidpy_extras(coords: np.ndarray, labels: np.ndarray, categories: list[str],
                   knn: int, n_perm: int, seed: int) -> dict[str, pd.DataFrame]:
    """Battle-tested spatial QoIs via squidpy: nhood_enrichment z-matrix and Ripley's L
    (multi-scale clustering). Co-occurrence is covered by our own curve. squidpy only."""
    import anndata as ad
    import pandas as _pd
    import squidpy as sq

    cats = [c for c in categories if c in set(labels)]
    adata = ad.AnnData(
        X=np.zeros((len(labels), 1), dtype="float32"),
        obs=_pd.DataFrame({"compartment": _pd.Categorical(labels, categories=cats)}),
    )
    adata.obsm["spatial"] = np.asarray(coords, float)
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=knn)
    out: dict[str, pd.DataFrame] = {}

    # n_jobs=1: squidpy's default spawns processes, which deadlocks under Windows spawn
    # semantics (re-imports the caller). Single-process is fine at this data size.
    sq.gr.nhood_enrichment(adata, cluster_key="compartment", n_perms=n_perm, seed=seed,
                           n_jobs=1, show_progress_bar=False)
    z = adata.uns["compartment_nhood_enrichment"]["zscore"]
    out["nhood_z"] = _pd.DataFrame(z, index=cats, columns=cats)

    sq.gr.ripley(adata, cluster_key="compartment", mode="L")
    out["ripley_L"] = adata.uns["compartment_ripley_L"]["L_stat"]
    return out


def iter_tumor_spatial(cfg):
    """Yield (sample_id, coords_um, dominant_state_labels) for each Visium tumor."""
    sig = pd.read_parquet(
        resolve_path(cfg, "data/processed/spot_signatures.parquet"),
        columns=["sample_id", "barcode", "in_tissue", "dominant_state"])
    libs = discover_library_dirs(resolve_path(cfg, cfg["paths"]["phase_b"]["spatial_root"]))
    dia = float(cfg["phase_c"]["spot"]["diameter_um"])
    knn = int(cfg["phase_c"]["spatial_qoi"]["knn"])
    for sid, lib in libs.items():
        s = sig[(sig["sample_id"] == sid) & (sig["in_tissue"] == 1)].copy()
        if s.empty:
            continue
        coords_all = load_spot_coords_um(lib["library_dir"], dia)
        s = s.join(coords_all[["x_um", "y_um"]], on="barcode").dropna(subset=["x_um", "y_um"])
        if len(s) <= knn:
            continue
        yield sid, s[["x_um", "y_um"]].to_numpy(), s["dominant_state"].to_numpy()


def observed_spatial_baseline(cfg) -> tuple[pd.DataFrame, pd.DataFrame]:
    seed = int(cfg["phase_c"]["seed"])
    sq_cfg = cfg["phase_c"]["spatial_qoi"]
    knn, n_perm = int(sq_cfg["knn"]), int(sq_cfg["n_perm"])
    radii = list(sq_cfg["cooccur_radii_um"])
    enr, occ = [], []
    for sid, xy, lab in iter_tumor_spatial(cfg):
        e = neighbor_enrichment(xy, lab, COMPARTMENTS, k=knn, n_perm=n_perm, seed=seed)
        e.insert(0, "sample_id", sid)
        enr.append(e)
        o = co_occurrence(xy, lab, COMPARTMENTS, radii)
        o.insert(0, "sample_id", sid)
        occ.append(o)
    enr_df = pd.concat(enr, ignore_index=True) if enr else pd.DataFrame()
    occ_df = pd.concat(occ, ignore_index=True) if occ else pd.DataFrame()
    return enr_df, occ_df


def observed_squidpy_baseline(cfg) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ripley's L and squidpy nhood z-scores per tumor (long format). squidpy only."""
    seed = int(cfg["phase_c"]["seed"])
    sq_cfg = cfg["phase_c"]["spatial_qoi"]
    knn, n_perm = int(sq_cfg["knn"]), int(sq_cfg["n_perm"])
    ripley, nhood = [], []
    for sid, xy, lab in iter_tumor_spatial(cfg):
        try:
            ex = squidpy_extras(xy, lab, COMPARTMENTS, knn, n_perm, seed)
        except Exception as e:  # pragma: no cover - squidpy runtime quirks
            print(f"[warn] squidpy failed for {sid}: {e}")
            continue
        rl = ex["ripley_L"].copy()
        rl.insert(0, "sample_id", sid)
        ripley.append(rl)
        z = ex["nhood_z"].reset_index().melt(id_vars="index", var_name="cat_b", value_name="z")
        z = z.rename(columns={"index": "cat_a"})
        z.insert(0, "sample_id", sid)
        nhood.append(z)
    return (pd.concat(ripley, ignore_index=True) if ripley else pd.DataFrame(),
            pd.concat(nhood, ignore_index=True) if nhood else pd.DataFrame())


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

    # 1. observed spatial baseline (CPU, now): neighbourhood enrichment + co-occurrence curve
    backend = cfg["phase_c"]["spatial_qoi"].get("backend", "auto")
    use_sq = backend == "squidpy" or (backend == "auto" and has_squidpy())
    enr, occ = observed_spatial_baseline(cfg)
    if not enr.empty:
        enr.to_csv(out_dir / "observed_spatial_qoi.csv", index=False)
        seg = enr[enr.cat_a != enr.cat_b].groupby(["cat_a", "cat_b"])["z"].mean()
        print(f"[ok] enrichment  -> {out_dir/'observed_spatial_qoi.csv'} "
              f"({enr.sample_id.nunique()} tumors)")
        print("[info] mean cross-compartment enrichment z (negative = segregation):")
        print(seg.round(2).to_string())
    if not occ.empty:
        occ.to_csv(out_dir / "observed_cooccurrence.csv", index=False)
        cross = (occ[occ.cond != occ.exp].groupby(["cond", "exp", "r_um"])["cooccur"]
                 .mean().reset_index())
        near = cross[cross.r_um <= 150].groupby(["cond", "exp"])["cooccur"].mean()
        print(f"[ok] co-occurrence -> {out_dir/'observed_cooccurrence.csv'} "
              f"(curve over {occ.r_um.nunique()} distances)")
        print("[info] mean near-range (<=150um) cross co-occurrence (>1 co-locate, <1 segregate):")
        print(near.round(2).to_string())
    print(f"[info] spatial backend: {'squidpy' if use_sq else 'inhouse'}"
          f"{'' if has_squidpy() else ' (squidpy not installed)'}")
    if use_sq:
        ripley, nhood_z = observed_squidpy_baseline(cfg)
        if not ripley.empty:
            ripley.to_csv(out_dir / "observed_ripley_L.csv", index=False)
            print(f"[ok] squidpy Ripley's L -> {out_dir/'observed_ripley_L.csv'}")
        if not nhood_z.empty:
            nhood_z.to_csv(out_dir / "observed_nhood_z_squidpy.csv", index=False)
            print(f"[ok] squidpy nhood z    -> {out_dir/'observed_nhood_z_squidpy.csv'}")

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
