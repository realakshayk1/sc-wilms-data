#!/usr/bin/env python3
"""Two research-driven upgrades to the regressive-tissue pilot:

(a) DE-SPECKLE (SpotSweeper-style spatially-aware QC, following Totty et al., Nat Methods
    2025): our low-UMI label conflates biological necrosis with technical dropout
    (dryspots/hangnails). SpotSweeper flags spots that are QC outliers RELATIVE TO THEIR
    LOCAL SPATIAL NEIGHBORHOOD (isolated bad spots = artifacts) while preserving coherent
    low-quality DOMAINS (real necrotic regions). We reimplement the core local-outlier
    z-score in Python (no Bioconductor dep) and split regressive spots into
    coherent-necrosis (keep) vs isolated-dropout (technical, remove).

(b) BROADEN (SIOP regression is more than necrosis): chemo-induced regression also includes
    FIBROSIS and XANTHOMATOUS/foamy-macrophage change, which are CELLULAR (high RNA) and so
    are invisible to a low-UMI necrosis label. We score fibrosis + macrophage programs per
    spot and define a broadened regression label with three subtypes, then report its
    composition, treatment contrast, and tumor-level H&E readout (vs the necrosis-only one).

Honest labels: this is a Python reimplementation of SpotSweeper's local-outlier idea, not the
package; signature thresholds are heuristic and unvalidated against pathology.
"""
from __future__ import annotations

import glob
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
SPATIAL_ROOT = ROOT / "data" / "raw" / "scpca_downloads" / "spaceranger" / "SCPCP000006_spatial"
SIG_CACHE = PROC / "spot_regression_signatures.parquet"
OUT_DIR = ROOT / "results" / "regressive_pilot"
SEED = 42

FIBROSIS = ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A3", "FN1", "DCN", "LUM", "POSTN",
            "SPARC", "THBS2", "ACTA2", "TAGLN", "MMP2", "TIMP1", "VIM"]
MACRO = ["CD68", "CD163", "MRC1", "APOE", "APOC1", "TREM2", "GPNMB", "MARCO", "FABP4",
         "PLIN2", "SPP1", "FTL", "FTH1", "LGALS3", "CD14", "CTSD"]


def load_pilot():
    p = Path(__file__).resolve().parent / "19_regressive_tissue_pilot.py"
    spec = importlib.util.spec_from_file_location("pilot", p); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


def sanitize(bc):
    import re
    return re.sub(r"[^A-Za-z0-9]+", "_", bc)


def compute_signatures():
    """Per-spot fibrosis + macrophage program scores (mean log1p CP10k), cached."""
    if SIG_CACHE.exists():
        return pd.read_parquet(SIG_CACHE)
    rows = []
    for i, ld in enumerate(sorted(glob.glob(str(SPATIAL_ROOT / "*" / "*_spatial")))):
        ld = Path(ld); lib = ld.name.replace("_spatial", "")
        mtx = ld / "filtered_feature_bc_matrix"
        if not mtx.exists():
            continue
        ad = sc.read_10x_mtx(mtx, var_names="gene_symbols", make_unique=True)
        counts = np.asarray(ad.X.sum(1)).ravel(); counts[counts == 0] = 1
        up = list(ad.var_names.str.upper())
        def score(genes):
            idx = [up.index(g) for g in genes if g in up]
            cp = ad.X[:, idx] / counts[:, None] * 1e4
            return np.asarray(np.log1p(cp).mean(1)).ravel()
        fib, mac = score(FIBROSIS), score(MACRO)
        for bc, f, m in zip(ad.obs_names, fib, mac):
            rows.append((f"{lib}__{sanitize(bc)}", float(f), float(m)))
        print(f"   sig {i+1}/100 {lib}", flush=True)
    df = pd.DataFrame(rows, columns=["spot_id", "fibrosis_score", "macrophage_score"])
    df.to_parquet(SIG_CACHE, index=False)
    return df


def load_coords():
    names = ["barcode", "in_tissue", "array_row", "array_col", "pr", "pc"]
    frames = []
    for f in glob.glob(str(SPATIAL_ROOT / "*" / "*_spatial" / "spatial" / "tissue_positions_list.csv")):
        lib = Path(f).parents[1].name.replace("_spatial", "")
        d = pd.read_csv(f, header=None, names=names); d["library_id"] = lib
        frames.append(d[["library_id", "barcode", "array_row", "array_col"]])
    return pd.concat(frames, ignore_index=True)


def spotsweeper_local_z(lab, metric="log_umi", k=18):
    """Local-outlier z of a QC metric vs spatial kNN, per library (SpotSweeper core idea)."""
    lab = lab.copy()
    lab["log_umi"] = np.log10(lab["total_counts"] + 1)
    z = np.full(len(lab), np.nan)
    for lib, d in lab.dropna(subset=["array_row"]).groupby("library_id"):
        xy = d[["array_row", "array_col"]].to_numpy(float)
        if len(xy) < k + 2:
            continue
        _, idx = cKDTree(xy).query(xy, k=k + 1)
        neigh = idx[:, 1:]
        v = d[metric].to_numpy()
        nmean = v[neigh].mean(1); nstd = v[neigh].std(1) + 1e-9
        z[d.index.values] = (v - nmean) / nstd
    lab["local_z"] = z
    return lab


def loto(X, y, kind="ridge", seed=SEED):
    n = len(y); pred = np.full(n, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        if kind == "rf":
            m = make_pipeline(StandardScaler(), PCA(min(10, X.shape[1], tr.sum()-1), random_state=seed),
                              RandomForestRegressor(300, min_samples_leaf=2, random_state=seed, n_jobs=-1))
        else:
            m = make_pipeline(StandardScaler(), Ridge(alpha=10.0))
        m.fit(X[tr], y[tr]); pred[i] = m.predict(X[i:i+1])[0]
    return pred


def he_readout(per_tumor, col):
    emb = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    tum = emb.groupby("sample_id")[ecols].mean().reset_index().merge(
        per_tumor[["sample_id", col]], on="sample_id", how="inner").dropna()
    X = tum[ecols].to_numpy(); y = tum[col].to_numpy()
    rng = np.random.default_rng(SEED)
    res = {}
    for kind in ["ridge", "rf"]:
        p = loto(X, y, kind)
        res[kind] = {"pearson": float(pearsonr(y, p)[0]), "spearman": float(spearmanr(y, p)[0])}
    p = loto(X, y, "ridge")
    boot = [pearsonr(y[i], p[i])[0] for i in (rng.integers(0, len(y), len(y)) for _ in range(2000))
            if np.std(y[i]) > 1e-12 and np.std(p[i]) > 1e-12]
    res["ridge_ci95"] = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    return res


def main():
    pilot = load_pilot()
    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    lab = pilot.label_regressive(sig, pilot.compute_spot_qc())          # has regressive, pct_mito, n_genes
    lab = lab.merge(compute_signatures(), on="spot_id", how="left")
    lab = lab.merge(load_coords(), on=["library_id", "barcode"], how="left")
    tx = pilot.load_treatment_map()

    # ---------- (a) de-speckle: coherent necrosis vs isolated technical dropout ----------
    print("[a] SpotSweeper-style local-outlier QC ...", flush=True)
    lab = spotsweeper_local_z(lab, "log_umi", k=18)
    # isolated dropout = strong LOCAL low outlier (worse than neighbours); coherent = not
    lab["is_local_dropout"] = (lab["local_z"] < -2.5).astype(int)
    reg = lab["regressive"] == 1
    lab["necrosis_coherent"] = (reg & (lab["is_local_dropout"] == 0)).astype(int)
    n_reg = int(reg.sum()); n_drop = int((reg & (lab["is_local_dropout"] == 1)).sum())
    print(f"    of {n_reg:,} regressive spots, {n_drop:,} ({100*n_drop/max(n_reg,1):.1f}%) are "
          f"isolated technical dropouts -> removed; {n_reg-n_drop:,} coherent necrosis kept", flush=True)

    # ---------- (b) broaden: fibrotic + xanthomatous regression (cellular) ----------
    print("[b] broadening label with fibrosis + macrophage programs ...", flush=True)
    parts = []
    for libname, d in lab.groupby("library_id"):
        d = d.copy()
        if len(d) < 50:
            parts.append(d); continue
        fib_hi = d["fibrosis_score"] >= d["fibrosis_score"].quantile(0.85)
        mac_hi = d["macrophage_score"] >= d["macrophage_score"].quantile(0.85)
        low_blast = d["blastemal_program"] < d["blastemal_program"].quantile(0.5)
        d["fibrotic"] = (fib_hi & low_blast & (d["necrosis_coherent"] == 0)).astype(int)
        d["xanthomatous"] = (mac_hi & (d["necrosis_coherent"] == 0) & (d["fibrotic"] == 0)).astype(int)
        parts.append(d)
    lab = pd.concat(parts)
    lab["regression_broad"] = ((lab["necrosis_coherent"] == 1) | (lab["fibrotic"] == 1) |
                               (lab["xanthomatous"] == 1)).astype(int)

    # ---------- composition + treatment ----------
    per_tumor = (lab.groupby("sample_id").agg(
        frac_necrosis=("necrosis_coherent", "mean"), frac_fibrotic=("fibrotic", "mean"),
        frac_xanthomatous=("xanthomatous", "mean"), frac_regression_broad=("regression_broad", "mean"),
        frac_regressive_orig=("regressive", "mean")).reset_index().merge(tx, on="sample_id", how="left"))
    per_tumor.to_csv(OUT_DIR / "per_tumor_regression_broad.csv", index=False)

    def contrast(col):
        pc = per_tumor.loc[per_tumor.treatment == "Resection post chemotherapy", col]
        up = per_tumor.loc[per_tumor.treatment == "Upfront resection", col]
        return float(pc.mean()), float(up.mean()), float(mannwhitneyu(pc, up, alternative="greater")[1])

    comp = {c: {"mean_frac": float(lab[c].mean())} for c in
            ["necrosis_coherent", "fibrotic", "xanthomatous", "regression_broad"]}
    print("    composition (mean % of tissue): "
          f"necrosis={comp['necrosis_coherent']['mean_frac']*100:.1f}  "
          f"fibrotic={comp['fibrotic']['mean_frac']*100:.1f}  "
          f"xanthomatous={comp['xanthomatous']['mean_frac']*100:.1f}  "
          f"BROAD={comp['regression_broad']['mean_frac']*100:.1f}", flush=True)
    tx_out = {}
    for col in ["frac_necrosis", "frac_fibrotic", "frac_xanthomatous", "frac_regression_broad", "frac_regressive_orig"]:
        pcm, upm, p = contrast(col); tx_out[col] = {"postchemo": pcm, "upfront": upm, "p_greater": p}
        print(f"    {col:26} post-chemo={pcm*100:4.1f}% vs upfront={upm*100:4.1f}%  p={p:.3f}", flush=True)

    # ---------- H&E readout: broad regression + necrosis-coherent ----------
    print("[c] tumor-level H&E readout on broadened + de-speckled targets ...", flush=True)
    he = {}
    for col in ["frac_regression_broad", "frac_necrosis"]:
        he[col] = he_readout(per_tumor, col)
        r = he[col]
        print(f"    H&E -> {col:24} ridge r={r['ridge']['pearson']:+.3f} (Spearman {r['ridge']['spearman']:+.3f}), "
              f"rf r={r['rf']['pearson']:+.3f}, CI{[round(x,2) for x in r['ridge_ci95']]}", flush=True)

    out = {"despeckle": {"n_regressive": n_reg, "n_isolated_dropout_removed": n_drop,
                         "pct_removed": float(100*n_drop/max(n_reg, 1))},
           "composition": comp, "treatment_contrast": tx_out, "he_readout": he, "seed": SEED}
    (OUT_DIR / "regression_broaden_despeckle.json").write_text(json.dumps(out, indent=2))
    print("\n[ok] wrote", OUT_DIR / "regression_broaden_despeckle.json", flush=True)


if __name__ == "__main__":
    main()
