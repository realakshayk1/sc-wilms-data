#!/usr/bin/env python3
"""Close the open gap: run the tissue-quality confound control on the CONTINUOUS necrosis-QC
target (per-tumor mean mito%), not just the label fraction.

Question: is H&E -> mean-mito% (Ridge r=0.56) reading real necrosis morphology, or is it a
tissue-density/stain-quality path (H&E quality correlates with mito% which is itself a QC
metric)? Test: do quality scalars ALONE predict mito%, and does the embedding survive
partialling them out? Runs the control on both targets for completeness. Caches quality scalars.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import cv2
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
TILES = PROC / "he_tiles"
QUAL_CACHE = PROC / "spot_tile_quality.parquet"
OUT_DIR = ROOT / "results" / "regressive_pilot"
SEED = 42
QCOLS = ["brightness", "focus", "saturation", "hematoxylin"]


def load_pilot():
    p = Path(__file__).resolve().parent / "19_regressive_tissue_pilot.py"
    spec = importlib.util.spec_from_file_location("pilot", p); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


def tile_quality(spot_ids, manifest):
    if QUAL_CACHE.exists():
        q = pd.read_parquet(QUAL_CACHE)
        if set(spot_ids).issubset(set(q["spot_id"])):
            return q
    man = manifest.set_index("spot_id")
    rows = []
    for k, sid in enumerate(spot_ids):
        if sid not in man.index:
            continue
        img = cv2.imread(str(ROOT / man.loc[sid, "image_path"]))
        if img is None:
            continue
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
        rows.append((sid, gray.mean(), cv2.Laplacian(gray, cv2.CV_64F).var(),
                     hsv[:, :, 1].mean(), 255.0 - img[:, :, 0].mean()))
        if k % 2000 == 0:
            print(f"   quality {k}/{len(spot_ids)}", flush=True)
    q = pd.DataFrame(rows, columns=["spot_id"] + QCOLS)
    q.to_parquet(QUAL_CACHE, index=False)
    return q


def loto(X, y, kind="ridge", seed=SEED):
    n = len(y); pred = np.full(n, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        if kind == "rf":
            m = make_pipeline(StandardScaler(), PCA(min(10, X.shape[1], tr.sum()-1), random_state=seed),
                              RandomForestRegressor(n_estimators=300, min_samples_leaf=2, random_state=seed, n_jobs=-1))
        else:
            m = make_pipeline(StandardScaler(), Ridge(alpha=10.0))
        m.fit(X[tr], y[tr]); pred[i] = m.predict(X[i:i+1])[0]
    return pred


def partial_corr(y, a, ctrl):
    def resid(v):
        A = np.column_stack([np.ones_like(ctrl), ctrl]); b, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ b
    return float(pearsonr(resid(y), resid(a))[0])


def main():
    pilot = load_pilot()
    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    lab = pilot.label_regressive(sig, pilot.compute_spot_qc())
    qcm = lab.dropna(subset=["pct_mito"])
    per_tumor = (lab.groupby("sample_id").agg(frac_regressive=("regressive", "mean")).reset_index()
                 .merge(qcm.groupby("sample_id").agg(mean_mito=("pct_mito", "mean")).reset_index(),
                        on="sample_id", how="left"))
    emb = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    manifest = pd.DataFrame(json.loads((TILES / "tiles_manifest.json").read_text()))
    q = tile_quality(emb["spot_id"].tolist(), manifest).merge(emb[["spot_id", "sample_id"]], on="spot_id", how="left")
    tum_q = q.groupby("sample_id")[QCOLS].mean().reset_index()
    tum_e = emb.groupby("sample_id")[ecols].mean().reset_index()
    df = tum_e.merge(tum_q, on="sample_id").merge(per_tumor, on="sample_id").dropna()
    E = df[ecols].to_numpy(); Q = df[QCOLS].to_numpy()
    print(f"[data] {len(df)} tumors with embedding+quality+targets", flush=True)

    out = {"n_tumors": int(len(df)), "seed": SEED, "targets": {}}
    for tname in ["mean_mito", "frac_regressive"]:
        yv = df[tname].to_numpy()
        p_emb = loto(E, yv, "ridge"); p_qual = loto(Q, yv, "ridge")
        r_emb = float(pearsonr(yv, p_emb)[0]); r_qual = float(pearsonr(yv, p_qual)[0])
        pr = partial_corr(yv, p_emb, p_qual)
        # also which single quality scalar correlates most with the target (diagnostic)
        qcorr = {c: float(pearsonr(df[c].to_numpy(), yv)[0]) for c in QCOLS}
        out["targets"][tname] = {"r_embedding": r_emb, "r_quality_only": r_qual,
                                 "partial_r_emb_given_quality": pr,
                                 "quality_scalar_correlations": qcorr,
                                 "reading": ("embedding survives quality control" if pr > 0.2 and r_emb - r_qual > 0.05
                                             else "quality may explain it" if r_qual >= r_emb - 0.05 else "mixed")}
        print(f"  [{tname}] r_emb={r_emb:+.3f} r_qual={r_qual:+.3f} partial(emb|qual)={pr:+.3f} -> "
              f"{out['targets'][tname]['reading']}", flush=True)
        print(f"      quality-scalar corrs: " + ", ".join(f"{k}={v:+.2f}" for k, v in qcorr.items()), flush=True)

    (OUT_DIR / "mito_quality_control.json").write_text(json.dumps(out, indent=2))
    print("\n[ok] wrote", OUT_DIR / "mito_quality_control.json", flush=True)


if __name__ == "__main__":
    main()
