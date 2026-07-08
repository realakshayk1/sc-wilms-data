#!/usr/bin/env python3
"""AUDIT pass — further robustness checks on the regressive tumor-level signal.

Addresses risks found auditing scripts 19-21:
  1. OUTLIER ROBUSTNESS : Spearman (rank) alongside Pearson — %regressive is right-skewed,
     so Pearson r=0.34 could be driven by a few high-necrosis tumors.
  2. PROPER PERMUTATION : re-run the WHOLE LOTO with shuffled labels (not the fixed-prediction
     shortcut used in 21) for an honest p-value.
  3. MODEL DEPENDENCE   : plain Ridge alongside RF+PCA — is r a modelling artifact?
  4. LESS-CIRCULAR TARGET: predict per-tumor mean mito% (independent necrosis QC, NOT the
     UMI-defined label) from H&E. If H&E tracks mito%, the necrosis signal is not just the
     label's own definition leaking through.
  6. LABEL COLLINEARITY : how redundant is "low-UMI AND depleted-viable"?

(5 ResNet-vs-Phikon ablation is a separate embedding job: 23_regressive_resnet_ablation.py.)
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
OUT_DIR = ROOT / "results" / "regressive_pilot"
SEED = 42


def load_pilot():
    p = Path(__file__).resolve().parent / "19_regressive_tissue_pilot.py"
    spec = importlib.util.spec_from_file_location("pilot", p); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


def make_model(kind, seed, nfeat, ntrain):
    if kind == "rf":
        ncomp = min(10, nfeat, ntrain - 1)
        return make_pipeline(StandardScaler(), PCA(ncomp, random_state=seed),
                             RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                                   random_state=seed, n_jobs=-1))
    return make_pipeline(StandardScaler(), Ridge(alpha=10.0))


def loto(X, y, kind="rf", seed=SEED):
    n = len(y); pred = np.full(n, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        m = make_model(kind, seed, X.shape[1], tr.sum())
        m.fit(X[tr], y[tr]); pred[i] = m.predict(X[i:i+1])[0]
    return pred


def main():
    pilot = load_pilot()
    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    qc = pilot.compute_spot_qc()
    lab = pilot.label_regressive(sig, qc)

    # per-tumor: regressive fraction + independent necrosis QC (mean mito%, mean n_genes)
    qcm = lab.dropna(subset=["pct_mito", "n_genes"])
    per_tumor = (lab.groupby("sample_id").agg(frac_regressive=("regressive", "mean")).reset_index()
                 .merge(qcm.groupby("sample_id").agg(mean_mito=("pct_mito", "mean"),
                                                     mean_ngenes=("n_genes", "mean")).reset_index(),
                        on="sample_id", how="left"))
    emb = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    tum = emb.groupby("sample_id")[ecols].mean().reset_index().merge(per_tumor, on="sample_id", how="inner")
    E = tum[ecols].to_numpy(); n = len(tum)
    rng = np.random.default_rng(SEED)
    out = {"n_tumors": int(n), "seed": SEED}

    # ---- 1,3,4: predict each target with RF and Ridge; Pearson + Spearman ----
    print("[1,3,4] targets x models (Pearson + Spearman) ...", flush=True)
    targets = {"frac_regressive": tum["frac_regressive"].to_numpy(),
               "mean_mito_QC": tum["mean_mito"].to_numpy(),
               "mean_ngenes_QC": tum["mean_ngenes"].to_numpy()}
    out["targets"] = {}
    preds_cache = {}
    for tname, yv in targets.items():
        out["targets"][tname] = {}
        for kind in ["rf", "ridge"]:
            p = loto(E, yv, kind=kind)
            preds_cache[(tname, kind)] = (yv, p)
            pr = float(pearsonr(yv, p)[0]); sr = float(spearmanr(yv, p)[0])
            out["targets"][tname][kind] = {"pearson_r": pr, "spearman_r": sr}
            print(f"    {tname:16} {kind:5}  Pearson={pr:+.3f}  Spearman={sr:+.3f}", flush=True)

    # ---- 2: proper permutation (re-run full LOTO) on frac_regressive, both models ----
    print("[2] proper label-permutation p (re-fit LOTO) ...", flush=True)
    out["proper_permutation"] = {}
    for kind in ["rf", "ridge"]:
        yv, p = preds_cache[("frac_regressive", kind)]
        r_obs = pearsonr(yv, p)[0]
        B = 300 if kind == "ridge" else 150
        null = []
        for _ in range(B):
            yp = rng.permutation(yv)
            pp = loto(E, yp, kind=kind)
            null.append(pearsonr(yp, pp)[0])
        null = np.array(null)
        pval = float((null >= r_obs).mean())
        out["proper_permutation"][kind] = {"r_obs": float(r_obs), "B": B, "perm_p_one_sided": pval,
                                           "null_mean": float(null.mean()), "null_p95": float(np.percentile(null, 95))}
        print(f"    {kind:5}  r={r_obs:+.3f}  proper-perm p={pval:.4f} (B={B}, null95={np.percentile(null,95):+.3f})", flush=True)

    # ---- bootstrap CI on the less-circular mito% target (Ridge) ----
    yv, p = preds_cache[("mean_mito_QC", "ridge")]
    boot = [pearsonr(yv[idx], p[idx])[0] for idx in (rng.integers(0, n, n) for _ in range(2000))
            if np.std(yv[idx]) > 1e-12 and np.std(p[idx]) > 1e-12]
    out["mito_target_ci95"] = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    print(f"    mean_mito_QC (ridge) bootstrap 95% CI [{out['mito_target_ci95'][0]:.3f}, {out['mito_target_ci95'][1]:.3f}]", flush=True)

    # ---- 6: label internal collinearity ----
    print("[6] label collinearity (low-UMI vs depleted-viable) ...", flush=True)
    s = sig[sig["in_tissue"] == 1].copy()
    s["viable_max"] = s[["blastemal_program", "epithelial_program", "stromal_program"]].max(1)
    lo, dep = [], []
    for _, d in s.groupby("library_id"):
        if len(d) < 50:
            continue
        ref = np.quantile(d["total_counts"], 0.75)
        lo.append((d["total_counts"] < 0.30 * max(ref, 1.0)).to_numpy())
        dep.append((d["viable_max"] < d["viable_max"].quantile(0.40)).to_numpy())
    lo = np.concatenate(lo); dep = np.concatenate(dep)
    # of low-UMI spots, what fraction also pass depleted-viable (redundancy)?
    frac_lo_also_dep = float(dep[lo].mean()) if lo.sum() else float("nan")
    jacc = float((lo & dep).sum() / (lo | dep).sum())
    out["label_collinearity"] = {"P(depleted | low_umi)": frac_lo_also_dep,
                                 "jaccard_low_umi_depleted": jacc,
                                 "frac_low_umi": float(lo.mean()), "frac_depleted": float(dep.mean())}
    print(f"    P(depleted | low-UMI)={frac_lo_also_dep:.2f}  Jaccard={jacc:.2f} "
          f"(low-UMI {lo.mean()*100:.1f}%, depleted {dep.mean()*100:.1f}%)", flush=True)

    (OUT_DIR / "regressive_audit.json").write_text(json.dumps(out, indent=2))
    print("\n[ok] wrote", OUT_DIR / "regressive_audit.json", flush=True)


if __name__ == "__main__":
    main()
