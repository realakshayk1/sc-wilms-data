#!/usr/bin/env python3
"""B-1: Scale + stronger encoder + attention-MIL for tumor-level histology.

The current positive (script 14) mean-pools 60 Phikon-v1 spot embeddings per tumor
and runs logistic regression -> tumor AUC 0.724. Three levers, additively:

  1. SCALE   : embed up to --n-per-sample (default 200) spots/tumor, not 60.
  2. ENCODER : phikon-v2 (ViT-L, 1024-d) instead of phikon-v1 (ViT-B, 768-d).
  3. MIL     : gated attention-MIL (Ilse et al. 2018) over spot embeddings so
               anaplastic spots get up-weighted, vs. a flat mean-pool that dilutes
               the diagnostic minority of atypical spots.

Held-out across tumors (leave-one-tumor-out). Reports DeLong CI + label-permutation
p for every model, and a paired DeLong test of the best new model vs the v1 mean-pool
baseline (V-1 hardening). Honest by construction: if MIL does not beat mean-pool, the
JSON shows it.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from utils import load_config, resolve_path, set_seed_logged, setup_logging
from phase_b_stats import delong_auc_ci, permutation_auc_p, delong_paired_test


def _load14():
    p = Path(__file__).resolve().parent / "14_phase_b_positives.py"
    spec = importlib.util.spec_from_file_location("pbpos", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def resumable_embed(cfg, encoder, n_per_sample, seed, device, batch=32, flush_every=8):
    """Embed n_per_sample tiles/tumor, CHECKPOINTING to append-only chunk files every
    flush_every batches so a process restart resumes instead of restarting (the
    embedding is the long pole and this env tears the process down periodically).
    Append-only -> each flush costs O(chunk), not O(all-so-far). Returns the full
    embedding DataFrame; writes the final cache used by script 14."""
    proc_dir = resolve_path(cfg, cfg["paths"]["dirs"]["processed"])
    final = proc_dir / f"phikon_spot_embeddings_{encoder}.parquet"
    pdir = proc_dir / f"phikon_spot_embeddings_{encoder}_partial"
    if final.exists():
        df = pd.read_parquet(final)
        if df["sample_id"].nunique() >= 40 and df.groupby("sample_id").size().median() >= n_per_sample - 5:
            print(f"[cache] {len(df):,} embeddings from {final.name}", flush=True)
            return df

    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    man = pd.DataFrame(json.loads((tiles_dir / "tiles_manifest.json").read_text()))
    man = (man.groupby("sample_id", group_keys=False)
              .apply(lambda d: d.sample(min(len(d), n_per_sample), random_state=seed)))
    man = man.reset_index(drop=True)
    if "subdiagnosis" not in man:
        man["subdiagnosis"] = None

    pdir.mkdir(parents=True, exist_ok=True)
    chunks = sorted(pdir.glob("chunk_*.parquet"))
    done = set()
    if chunks:
        prev = pd.concat([pd.read_parquet(c) for c in chunks], ignore_index=True)
        done = set(prev["spot_id"])
        print(f"[resume] {len(done):,} tiles already embedded across {len(chunks)} chunks", flush=True)
    next_idx = len(chunks)
    todo = man[~man["spot_id"].isin(done)]
    print(f"[embed] {len(todo):,} of {len(man):,} tiles remaining ({encoder})", flush=True)

    if len(todo):
        embed_fn, dim, used = _load14()._load13().build_encoder(encoder, device)
        from PIL import Image
        root = Path(cfg["root"])
        paths = todo["image_path"].tolist(); sids = todo["spot_id"].tolist()
        samps = todo["sample_id"].tolist(); subs = todo["subdiagnosis"].tolist()
        buf_E, buf_meta, n_batches = [], [], 0
        for i in range(0, len(paths), batch):
            imgs = [Image.open(root / p).convert("RGB") for p in paths[i:i + batch]]
            buf_E.append(embed_fn(imgs))
            for j in range(len(imgs)):
                buf_meta.append((sids[i + j], samps[i + j], subs[i + j]))
            n_batches += 1
            if n_batches % flush_every == 0:
                next_idx = _flush_chunk(buf_E, buf_meta, dim, pdir, next_idx)
                buf_E, buf_meta = [], []
                print(f"   embedded {i + len(imgs)}/{len(paths)} (checkpointed)", flush=True)
        if buf_E:
            next_idx = _flush_chunk(buf_E, buf_meta, dim, pdir, next_idx)

    df = pd.concat([pd.read_parquet(c) for c in sorted(pdir.glob("chunk_*.parquet"))], ignore_index=True)
    df.to_parquet(final, index=False)
    for c in pdir.glob("chunk_*.parquet"):
        c.unlink()
    pdir.rmdir()
    print(f"[cache] wrote {final} ({len(df):,} embeddings)", flush=True)
    return df


def _flush_chunk(buf_E, buf_meta, dim, pdir, next_idx):
    chunk = pd.DataFrame(np.vstack(buf_E), columns=[f"e{i}" for i in range(dim)])
    chunk["spot_id"] = [m[0] for m in buf_meta]
    chunk["sample_id"] = [m[1] for m in buf_meta]
    chunk["subdiagnosis"] = [m[2] for m in buf_meta]
    tmp = pdir / f".chunk_{next_idx:04d}.tmp.parquet"
    chunk.to_parquet(tmp, index=False)
    tmp.rename(pdir / f"chunk_{next_idx:04d}.parquet")   # atomic publish
    return next_idx + 1


# --------------------------- attention-MIL ---------------------------
class GatedAttentionMIL(nn.Module):
    """Tiny gated-attention MIL head. Kept deliberately small (one instance layer,
    one attention head) because there are only ~40 training bags per LOTO fold."""

    def __init__(self, in_dim: int, hidden: int = 128, att_dim: int = 64, dropout: float = 0.4):
        super().__init__()
        self.inst = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.att_V = nn.Linear(hidden, att_dim)
        self.att_U = nn.Linear(hidden, att_dim)
        self.att_w = nn.Linear(att_dim, 1)
        self.clf = nn.Linear(hidden, 1)

    def forward(self, bag: torch.Tensor):
        h = self.inst(bag)                                  # [K, hidden]
        a = self.att_w(torch.tanh(self.att_V(h)) * torch.sigmoid(self.att_U(h)))  # [K,1]
        a = torch.softmax(a, dim=0)
        z = torch.sum(a * h, dim=0, keepdim=True)           # [1, hidden]
        return self.clf(z).squeeze(), a.squeeze(-1)


def train_mil(bags, labels, in_dim, seed, epochs=120, lr=1e-3, wd=1e-3):
    torch.manual_seed(seed)
    model = GatedAttentionMIL(in_dim)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    pos_w = torch.tensor([(len(labels) - sum(labels)) / max(1, sum(labels))], dtype=torch.float32)
    lossf = nn.BCEWithLogitsLoss(pos_weight=pos_w)
    model.train()
    order = list(range(len(bags)))
    rng = np.random.default_rng(seed)
    for _ in range(epochs):
        rng.shuffle(order)
        for i in order:
            opt.zero_grad()
            logit, _ = model(bags[i])
            loss = lossf(logit.unsqueeze(0), torch.tensor([labels[i]], dtype=torch.float32))
            loss.backward()
            opt.step()
    return model


def loto_mil(emb, ecols, seed, epochs):
    """Leave-one-tumor-out attention-MIL. Returns (proba[len=n_tumors], y, tumor_ids)."""
    tumors = emb["sample_id"].unique()
    y_by_t = emb.groupby("sample_id")["y"].first()
    proba = np.full(len(tumors), np.nan)
    y = np.array([int(y_by_t[t]) for t in tumors])
    for ti, held in enumerate(tumors):
        tr_ids = [t for t in tumors if t != held]
        if len(np.unique([int(y_by_t[t]) for t in tr_ids])) < 2:
            continue
        scaler = StandardScaler().fit(emb.loc[emb.sample_id.isin(tr_ids), ecols].to_numpy())

        def bag_of(t):
            X = scaler.transform(emb.loc[emb.sample_id == t, ecols].to_numpy())
            return torch.tensor(X, dtype=torch.float32)

        bags = [bag_of(t) for t in tr_ids]
        labels = [int(y_by_t[t]) for t in tr_ids]
        model = train_mil(bags, labels, len(ecols), seed + ti, epochs=epochs)
        model.eval()
        with torch.no_grad():
            logit, _ = model(bag_of(held))
            proba[ti] = torch.sigmoid(logit).item()
    return proba, y, list(tumors)


def loto_meanpool(emb, ecols, seed, C=1.0, n_pca=10):
    """Mean-pool tumor embedding -> logistic (the script-14 baseline), LOTO."""
    from sklearn.decomposition import PCA
    from sklearn.pipeline import make_pipeline
    tum = emb.groupby("sample_id").agg({**{c: "mean" for c in ecols}, "y": "first"}).reset_index()
    X = tum[ecols].to_numpy(); y = tum["y"].to_numpy(); g = tum["sample_id"].to_numpy()
    proba = np.full(len(y), np.nan)
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        if len(np.unique(y[tr])) < 2:
            continue
        m = make_pipeline(StandardScaler(), PCA(min(n_pca, tr.sum() - 1), random_state=seed),
                          LogisticRegression(max_iter=2000, C=C))
        m.fit(X[tr], y[tr])
        proba[i] = m.predict_proba(X[i:i + 1])[:, 1][0]
    return proba, y, list(g)


def evaluate(name, proba, y, n_perm, seed):
    ok = ~np.isnan(proba)
    ci = delong_auc_ci(y[ok], proba[ok])
    perm = permutation_auc_p(y[ok], proba[ok], n_perm=n_perm, seed=seed)
    print(f"[{name}] AUC={ci['auc']:.3f}  95% CI [{ci['ci_low']:.3f}, {ci['ci_high']:.3f}]  "
          f"perm p={perm['p_value']:.4f}  (n={int(ok.sum())})", flush=True)
    return {"auc": ci["auc"], "ci_low": ci["ci_low"], "ci_high": ci["ci_high"],
            "delong_se": ci["se"], "perm_p": perm["p_value"], "n": int(ok.sum())}


def main():
    setup_logging()
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", default="phikon-v2", choices=["phikon", "phikon-v2"])
    ap.add_argument("--n-per-sample", type=int, default=200)
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "phase_b_mil")
    dev = "cpu" if (args.device == "cuda" and not torch.cuda.is_available()) else args.device

    # --- embed (resumable checkpointed cache); scale spots/tumor ---
    emb = resumable_embed(cfg, args.encoder, args.n_per_sample, seed, dev)
    emb = emb[emb["subdiagnosis"].isin(["favorable", "anaplastic"])].copy()
    emb["y"] = (emb["subdiagnosis"] == "anaplastic").astype(int)
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    print(f"[data] {len(emb):,} spots, {emb['sample_id'].nunique()} tumors, "
          f"{len(ecols)}-d {args.encoder}, median {int(emb.groupby('sample_id').size().median())} spots/tumor",
          flush=True)

    results = {"encoder": args.encoder, "embedding_dim": len(ecols),
               "n_spots": int(len(emb)), "n_tumors": int(emb["sample_id"].nunique()),
               "n_per_sample": args.n_per_sample,
               "n_favorable": int((emb.groupby("sample_id")["y"].first() == 0).sum()),
               "n_anaplastic": int((emb.groupby("sample_id")["y"].first() == 1).sum()),
               "baseline_v1_meanpool_auc_ref": 0.724, "seed": seed, "models": {}}

    # --- models, all LOTO over the same tumors ---
    p_mp, y_mp, t_mp = loto_meanpool(emb, ecols, seed)
    results["models"]["meanpool_logistic"] = evaluate(f"{args.encoder} mean-pool", p_mp, y_mp, args.n_perm, seed)

    p_mil, y_mil, t_mil = loto_mil(emb, ecols, seed, args.epochs)
    results["models"]["attention_mil"] = evaluate(f"{args.encoder} attn-MIL", p_mil, y_mil, args.n_perm, seed)

    # --- paired DeLong: MIL vs mean-pool on the shared tumor set ---
    assert t_mp == t_mil, "tumor order mismatch between models"
    ok = (~np.isnan(p_mil)) & (~np.isnan(p_mp))
    paired = delong_paired_test(y_mil[ok], p_mil[ok], p_mp[ok])
    results["paired_delong_mil_vs_meanpool"] = paired
    print(f"[paired DeLong] MIL {paired['auc_a']:.3f} vs mean-pool {paired['auc_b']:.3f}  "
          f"delta={paired['delta']:+.3f}  p={paired['p_value']:.3f}", flush=True)

    out = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"]).parent / f"phase_b_mil_{args.encoder}.json"
    out.write_text(json.dumps(results, indent=2))
    # also persist held-out predictions for downstream paired tests / ensembling
    pd.DataFrame({"sample_id": t_mil, "y": y_mil, "p_meanpool": p_mp, "p_mil": p_mil}).to_csv(
        out.with_suffix(".predictions.csv"), index=False)
    print(f"[ok] -> {out}", flush=True)


if __name__ == "__main__":
    main()
