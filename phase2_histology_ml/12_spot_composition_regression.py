#!/usr/bin/env python3
"""B4: Spot-level composition regression — can aggregated H&E morphology predict the
per-spot transcriptomic compartment composition, held out across tumors (LOTO)?

This replaces the nucleus classifier (which collapsed under 85%-stromal weak labels).
It is a genuine CROSS-MODAL test: features are morphology only, target is the
transcriptomic composition, so there is no label/target circularity. The target is
z-scored per program across spots before softmax to remove the systematic stromal
offset that biased the weak labels.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

CELL_STATES = ["blastemal", "epithelial", "stromal"]
MORPH = ["area", "eccentricity", "solidity", "major_axis_length",
         "texture_var", "hematoxylin_intensity", "neighbor_density"]


def spot_morphology(nuc: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-nucleus morphology to per-spot features (mean, std, count)."""
    g = nuc.groupby("spot_id")
    agg = g[MORPH].agg(["mean", "std"])
    agg.columns = [f"{c}_{s}" for c, s in agg.columns]
    agg["n_nuclei"] = g.size()
    meta = g[["sample_id", "subdiagnosis"]].first()
    return agg.join(meta).reset_index()


def zscore_softmax_target(sig: pd.DataFrame) -> pd.DataFrame:
    """Bias-corrected per-spot composition: z-score each program across spots, softmax."""
    cols = [f"{s}_program" for s in CELL_STATES]
    z = sig[cols].to_numpy(float)
    z = (z - z.mean(0)) / (z.std(0) + 1e-9)
    z = z - z.max(1, keepdims=True)
    e = np.exp(z)
    frac = e / e.sum(1, keepdims=True)
    out = pd.DataFrame(frac, columns=[f"y_{s}" for s in CELL_STATES])
    out["spot_id"] = sig["spot_id"].values
    return out


def loto_regression(df: pd.DataFrame, feat_cols: list[str], seed: int,
                    max_train: int = 15000, shuffle_y=False, random_feat=False,
                    fold_limit: int | None = None, tag: str = "real"):
    samples = sorted(df["sample_id"].astype(str).unique())
    if fold_limit is not None:
        samples = samples[:fold_limit]
    rng = np.random.default_rng(seed)
    per_state = {s: [] for s in CELL_STATES}
    for k, ho in enumerate(samples):
        tr = df[df["sample_id"].astype(str) != ho]
        te = df[df["sample_id"].astype(str) == ho]
        if len(tr) < 200 or len(te) < 20:
            continue
        if len(tr) > max_train:
            tr = tr.sample(max_train, random_state=seed)
        Xtr, Xte = tr[feat_cols].to_numpy(), te[feat_cols].to_numpy()
        if random_feat:
            Xtr = rng.normal(size=Xtr.shape); Xte = rng.normal(size=Xte.shape)
        Ytr = tr[[f"y_{s}" for s in CELL_STATES]].to_numpy()
        if shuffle_y:
            Ytr = Ytr[rng.permutation(len(Ytr))]
        m = RandomForestRegressor(n_estimators=50, max_depth=16, min_samples_leaf=50,
                                  max_features="sqrt", n_jobs=-1, random_state=seed)
        m.fit(Xtr, Ytr)
        pred = m.predict(Xte)
        for i, s in enumerate(CELL_STATES):
            yt = te[f"y_{s}"].to_numpy()
            if yt.std() > 1e-6 and pred[:, i].std() > 1e-6:
                per_state[s].append(float(pearsonr(yt, pred[:, i])[0]))
        print(f"   [{tag}] fold {k+1}/{len(samples)} ({ho}) done", flush=True)
    return {s: (float(np.mean(v)) if v else float("nan"),
                float(np.std(v)) if v else float("nan"), len(v))
            for s, v in per_state.items()}


def main() -> None:
    setup_logging()
    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "spot_composition_regression")
    feats = resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"])
    sig_p = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    out_json = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"]).parent / "spot_composition_regression.json"
    ensure_dir(out_json.parent)

    nuc = pd.read_parquet(feats, columns=["spot_id", "sample_id", "subdiagnosis"] + MORPH)
    print(f"[load] {len(nuc):,} nuclei", flush=True)
    spot_feat = spot_morphology(nuc)
    del nuc
    print(f"[agg] {len(spot_feat):,} spots aggregated", flush=True)
    sig = pd.read_parquet(sig_p)
    sig = sig[sig["in_tissue"] == 1]
    target = zscore_softmax_target(sig)

    df = spot_feat.merge(target, on="spot_id", how="inner")
    df = df.dropna()
    feat_cols = [c for c in df.columns if any(c.startswith(m) for m in MORPH)] + ["n_nuclei"]
    print(f"[data] {len(df):,} spots, {df['sample_id'].nunique()} samples, {len(feat_cols)} features", flush=True)

    real = loto_regression(df, feat_cols, seed, tag="real")
    # negative controls on a subset of folds (cheaper; baseline only)
    shuf = loto_regression(df, feat_cols, seed, shuffle_y=True, fold_limit=6, tag="shuffled")
    rand = loto_regression(df, feat_cols, seed, random_feat=True, fold_limit=6, tag="random")

    out = {
        "n_spots": int(len(df)),
        "n_samples": int(df["sample_id"].nunique()),
        "n_features": len(feat_cols),
        "target": "z-scored-per-program softmax of transcriptomic compartment programs",
        "metric": "mean held-out Pearson r (LOTO) per compartment [mean, std, n_folds]",
        "real": real,
        "negative_control_shuffled_target": shuf,
        "negative_control_random_features": rand,
        "seed": seed,
    }
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)
    print("[ok] Spot composition regression (LOTO):")
    for s in CELL_STATES:
        print(f"   {s:11} real r={real[s][0]:.3f}±{real[s][1]:.3f}  "
              f"shuffled={shuf[s][0]:.3f}  random={rand[s][0]:.3f}  (n_folds={real[s][2]})")
    print(f"[ok] -> {out_json}")


if __name__ == "__main__":
    main()
