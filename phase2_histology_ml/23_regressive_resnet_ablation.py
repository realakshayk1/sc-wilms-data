#!/usr/bin/env python3
"""AUDIT 5: is the tumor-level regressive signal pathology-FM-specific, or generic texture?

Embed the SAME histology-cache spot tiles with a generic ImageNet ResNet50 (2048-d, no
pathology pretraining), take the per-tumor mean, and run the same LOTO -> %regressive test.
If ResNet matches Phikon-v2's r (~0.34), the signal is generic image structure; if ResNet
is clearly weaker, the pathology foundation model is contributing real histology.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from PIL import Image
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
TILES = PROC / "he_tiles"
OUT = PROC / "resnet_spot_embeddings.parquet"
OUT_DIR = ROOT / "results" / "regressive_pilot"
SEED = 42


def load_pilot():
    p = Path(__file__).resolve().parent / "19_regressive_tissue_pilot.py"
    spec = importlib.util.spec_from_file_location("pilot", p); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


def embed_resnet(spot_ids):
    import torchvision.transforms as T
    from torchvision.models import ResNet50_Weights, resnet50
    w = ResNet50_Weights.IMAGENET1K_V2
    net = resnet50(weights=w); net.fc = torch.nn.Identity(); net.eval()
    tf = w.transforms()
    man = pd.DataFrame(json.loads((TILES / "tiles_manifest.json").read_text())).set_index("spot_id")
    rows, embs = [], []
    batch, paths, sids = [], [], []
    for sid in spot_ids:
        if sid not in man.index:
            continue
        paths.append(ROOT / man.loc[sid, "image_path"]); sids.append(sid)

    @torch.no_grad()
    def flush(imgs):
        x = torch.stack([tf(im) for im in imgs])
        return net(x).cpu().numpy()

    buf, buf_ids = [], []
    for k, (p, sid) in enumerate(zip(paths, sids)):
        buf.append(Image.open(p).convert("RGB")); buf_ids.append(sid)
        if len(buf) == 32:
            embs.append(flush(buf)); rows.extend(buf_ids); buf, buf_ids = [], []
            if k % 1600 == 0:
                print(f"   resnet {k}/{len(paths)}", flush=True)
    if buf:
        embs.append(flush(buf)); rows.extend(buf_ids)
    E = np.vstack(embs)
    df = pd.DataFrame(E, columns=[f"e{i}" for i in range(E.shape[1])]); df["spot_id"] = rows
    return df


def loto_rf(X, y, seed=SEED):
    n = len(y); pred = np.full(n, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        m = make_pipeline(StandardScaler(), PCA(min(10, X.shape[1], tr.sum()-1), random_state=seed),
                          RandomForestRegressor(n_estimators=300, min_samples_leaf=2, random_state=seed, n_jobs=-1))
        m.fit(X[tr], y[tr]); pred[i] = m.predict(X[i:i+1])[0]
    return pred


def main():
    pilot = load_pilot()
    ph = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    if OUT.exists():
        rn = pd.read_parquet(OUT)
    else:
        rn = embed_resnet(ph["spot_id"].tolist())
        rn = rn.merge(ph[["spot_id", "sample_id"]], on="spot_id", how="left")
        rn.to_parquet(OUT, index=False)
        print(f"[cache] wrote {OUT} ({len(rn):,} embeddings)", flush=True)

    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    lab = pilot.label_regressive(sig, pilot.compute_spot_qc())
    per_tumor = lab.groupby("sample_id").agg(frac_regressive=("regressive", "mean")).reset_index()
    rcols = [c for c in rn.columns if c.startswith("e") and c[1:].isdigit()]
    tum = rn.groupby("sample_id")[rcols].mean().reset_index().merge(per_tumor, on="sample_id", how="inner")
    X = tum[rcols].to_numpy(); y = tum["frac_regressive"].to_numpy()
    p = loto_rf(X, y)
    res = {"encoder": "resnet50_imagenet", "n_tumors": int(len(tum)), "dim": len(rcols),
           "pearson_r": float(pearsonr(y, p)[0]), "spearman_r": float(spearmanr(y, p)[0]),
           "phikon_v2_reference_r": 0.337}
    (OUT_DIR / "regressive_resnet_ablation.json").write_text(json.dumps(res, indent=2))
    print(f"[ok] ResNet50 -> %regressive Pearson r={res['pearson_r']:.3f} "
          f"(Phikon-v2 was 0.337)  Spearman={res['spearman_r']:.3f}", flush=True)


if __name__ == "__main__":
    main()
