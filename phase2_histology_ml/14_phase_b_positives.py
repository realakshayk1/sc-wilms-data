#!/usr/bin/env python3
"""Phase B positives: turn the composition negative into real H&E signal.

The composition task (cross-modal regression) was the wrong instrument. H&E's native,
biologically-mandated signal is ANAPLASIA (nuclear atypia = the definition of unfavorable
histology). We test whether pathology-FM tile embeddings predict histology, held out
across tumors — at the spot level and tumor level — plus an extreme-nuclear-tail morphology
classifier. Embeddings are cached so classifiers iterate cheaply.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging


def _load13():
    p = Path(__file__).resolve().parent / "13_fm_embedding_regression.py"
    spec = importlib.util.spec_from_file_location("fmreg", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def compute_or_load_embeddings(cfg, encoder, n_per_sample, seed, device, batch=32):
    cache = resolve_path(cfg, cfg["paths"]["dirs"]["processed"]) / f"phikon_spot_embeddings_{encoder}.parquet"
    if cache.exists():
        df = pd.read_parquet(cache)
        if df["sample_id"].nunique() >= 40 and df.groupby("sample_id").size().median() >= n_per_sample - 5:
            print(f"[cache] {len(df):,} embeddings from {cache.name}", flush=True)
            return df
    fm = _load13()
    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    man = pd.DataFrame(json.loads((tiles_dir / "tiles_manifest.json").read_text()))
    man = (man.groupby("sample_id", group_keys=False)
              .apply(lambda d: d.sample(min(len(d), n_per_sample), random_state=seed)))
    # attach histology
    sig = pd.read_parquet(resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"]))
    embed, dim, used = fm.build_encoder(encoder, device)
    print(f"[encoder] {used} dim={dim}; embedding {len(man):,} tiles", flush=True)
    from PIL import Image
    root = Path(cfg["root"])
    E, ids, samps = [], [], []
    paths = man["image_path"].tolist(); sids = man["spot_id"].tolist(); ss = man["sample_id"].tolist()
    subs = man["subdiagnosis"].tolist() if "subdiagnosis" in man else [None] * len(man)
    sub_list = []
    for i in range(0, len(paths), batch):
        imgs = [Image.open(root / p).convert("RGB") for p in paths[i:i + batch]]
        E.append(embed(imgs))
        ids += sids[i:i + batch]; samps += ss[i:i + batch]; sub_list += subs[i:i + batch]
        if (i // batch) % 10 == 0:
            print(f"   embedded {i + len(imgs)}/{len(paths)}", flush=True)
    df = pd.DataFrame(np.vstack(E), columns=[f"e{i}" for i in range(dim)])
    df["spot_id"] = ids; df["sample_id"] = samps; df["subdiagnosis"] = sub_list
    ensure_dir(cache.parent); df.to_parquet(cache, index=False)
    print(f"[cache] wrote {cache}", flush=True)
    return df


def loto_auc(X, y, groups, model_fn, seed):
    """Leave-one-group-out predicted probabilities -> pooled AUC + per-fold."""
    logo = LeaveOneGroupOut()
    proba = np.full(len(y), np.nan)
    for tr, te in logo.split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        m = model_fn(); m.fit(X[tr], y[tr])
        proba[te] = m.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(proba)
    auc = roc_auc_score(y[ok], proba[ok]) if len(np.unique(y[ok])) == 2 else float("nan")
    return float(auc), int(ok.sum()), proba


def main():
    setup_logging()
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", default="phikon")
    ap.add_argument("--n-per-sample", type=int, default=60)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()
    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "phase_b_positives")
    import torch
    dev = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"

    emb = compute_or_load_embeddings(cfg, args.encoder, args.n_per_sample, seed, dev)
    emb = emb[emb["subdiagnosis"].isin(["favorable", "anaplastic"])].copy()
    emb["y"] = (emb["subdiagnosis"] == "anaplastic").astype(int)
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    results = {"encoder": args.encoder, "n_spots": int(len(emb)), "n_tumors": int(emb["sample_id"].nunique())}

    # Task A: spot-level histology, LOTO by tumor
    Xs = emb[ecols].to_numpy(); ys = emb["y"].to_numpy(); gs = emb["sample_id"].to_numpy()
    spot_model = lambda: make_pipeline(StandardScaler(), PCA(30, random_state=seed),
                                       LogisticRegression(max_iter=2000, C=0.5))
    auc_spot, n_spot, _ = loto_auc(Xs, ys, gs, spot_model, seed)
    results["spot_histology_auc"] = auc_spot
    print(f"[A] spot-level histology AUC (LOTO) = {auc_spot:.3f}  (n={n_spot})", flush=True)

    # Task B: tumor-level histology (mean embedding), leave-one-tumor-out
    tum = emb.groupby("sample_id").agg({**{c: "mean" for c in ecols}, "y": "first"}).reset_index()
    Xt = tum[ecols].to_numpy(); yt = tum["y"].to_numpy(); gt = tum["sample_id"].to_numpy()
    tum_model = lambda: make_pipeline(StandardScaler(), PCA(10, random_state=seed),
                                      LogisticRegression(max_iter=2000, C=1.0))
    auc_tum, n_tum, _ = loto_auc(Xt, yt, gt, tum_model, seed)
    results["tumor_histology_auc"] = auc_tum
    results["n_tumors_favorable"] = int((yt == 0).sum()); results["n_tumors_anaplastic"] = int((yt == 1).sum())
    print(f"[B] tumor-level histology AUC (LOO) = {auc_tum:.3f}  ({int((yt==0).sum())} fav / {int((yt==1).sum())} ana)", flush=True)

    # Task C: extreme-nuclear-tail morphology -> histology (tumor-level)
    nuc = pd.read_parquet(resolve_path(cfg, cfg["paths"]["phase_b"]["features_parquet"]),
                          columns=["sample_id", "subdiagnosis", "area", "major_axis_length",
                                   "hematoxylin_intensity", "eccentricity", "solidity"])
    nuc = nuc[nuc["area"] >= 20]
    def tail(s): return s.quantile([0.5, 0.9, 0.95, 0.99]).tolist()
    g = nuc.groupby("sample_id")
    rows = []
    for sid, d in g:
        feat = {"sample_id": sid, "y": int(d["subdiagnosis"].iloc[0] == "anaplastic")}
        for col in ["area", "major_axis_length", "hematoxylin_intensity", "eccentricity", "solidity"]:
            qs = d[col].quantile([0.5, 0.9, 0.99]).tolist()
            feat[f"{col}_p50"], feat[f"{col}_p90"], feat[f"{col}_p99"] = qs
            feat[f"{col}_cv"] = d[col].std() / (abs(d[col].mean()) + 1e-9)
        feat["giant_frac"] = float((d["area"] > 3 * d["area"].median()).mean())  # 3x-median = anaplasia proxy
        rows.append(feat)
    mdf = pd.DataFrame(rows)
    mdf = mdf[mdf["sample_id"].isin(emb["sample_id"].unique()) | True]
    mcols = [c for c in mdf.columns if c not in ("sample_id", "y")]
    Xm = mdf[mcols].to_numpy(); ym = mdf["y"].to_numpy(); gm = mdf["sample_id"].to_numpy()
    morph_model = lambda: RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                                 class_weight="balanced", random_state=seed, n_jobs=-1)
    auc_morph, n_m, _ = loto_auc(Xm, ym, gm, morph_model, seed)
    results["tumor_morphology_auc"] = auc_morph
    print(f"[C] tumor-level nuclear-morphology AUC (LOO) = {auc_morph:.3f}", flush=True)

    out = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"]).parent / "phase_b_positives.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"[ok] -> {out}")


if __name__ == "__main__":
    main()
