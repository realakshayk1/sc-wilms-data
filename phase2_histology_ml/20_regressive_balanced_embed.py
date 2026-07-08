#!/usr/bin/env python3
"""L2: regressive-balanced Phikon-v2 tile embedding.

The existing phikon_spot_embeddings cache subsampled ~200 spots/tumor for the HISTOLOGY
task, so regressive spots (~6-9%) are under-represented and the per-spot regressive readout
is starved of positives. Here we select a BALANCED set per tumor — up to N_PER regressive
and N_PER viable spots, using the sharpened label from 19_regressive_tissue_pilot.py — and
embed those tiles with Phikon-v2. Checkpointed (append-only chunks) so it resumes.

Output: data/processed/regressive_balanced_embeddings_phikon-v2.parquet
Then re-run 19_regressive_tissue_pilot.py to pick it up (L2 per-spot line).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
TILES = PROC / "he_tiles"
OUT = PROC / "regressive_balanced_embeddings_phikon-v2.parquet"
PDIR = PROC / "regressive_balanced_partial"
SEED = 42
N_PER = 40           # up to this many regressive + this many viable per tumor
BATCH = 32
FLUSH_EVERY = 6


def load_label():
    """Import label_regressive machinery from script 19 and build the sharpened label."""
    p = Path(__file__).resolve().parent / "19_regressive_tissue_pilot.py"
    spec = importlib.util.spec_from_file_location("pilot", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    qc = m.compute_spot_qc()
    return m.label_regressive(sig, qc)


def build_encoder():
    import os
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    from transformers import AutoImageProcessor, AutoModel
    proc = AutoImageProcessor.from_pretrained("owkin/phikon-v2")
    model = AutoModel.from_pretrained("owkin/phikon-v2").eval()

    @torch.no_grad()
    def embed(imgs):
        px = proc(images=imgs, return_tensors="pt")["pixel_values"]
        return model(pixel_values=px).last_hidden_state[:, 0].cpu().numpy()
    return embed, 1024


def select_balanced(lab: pd.DataFrame, manifest: pd.DataFrame, seed=SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    keep = []
    for sid, d in lab.groupby("sample_id"):
        for cls in (1, 0):
            pool = d[d["regressive"] == cls]["spot_id"].tolist()
            if not pool:
                continue
            take = rng.choice(pool, size=min(len(pool), N_PER), replace=False)
            keep.extend(take.tolist())
    sel = manifest[manifest["spot_id"].isin(set(keep))].copy()
    lab_map = lab.set_index("spot_id")["regressive"].to_dict()
    sel["regressive"] = sel["spot_id"].map(lab_map)
    return sel


def main():
    PDIR.mkdir(parents=True, exist_ok=True)
    lab = load_label()
    manifest = pd.DataFrame(json.loads((TILES / "tiles_manifest.json").read_text()))
    sel = select_balanced(lab, manifest)
    print(f"[select] {len(sel):,} tiles ({int(sel.regressive.sum())} regressive) "
          f"across {sel.sample_id.nunique()} tumors", flush=True)

    done = set()
    chunks = sorted(PDIR.glob("chunk_*.parquet"))
    if chunks:
        done = set(pd.concat([pd.read_parquet(c, columns=["spot_id"]) for c in chunks])["spot_id"])
        print(f"[resume] {len(done):,} already embedded", flush=True)
    todo = sel[~sel["spot_id"].isin(done)].reset_index(drop=True)
    print(f"[embed] {len(todo):,} tiles remaining", flush=True)

    if len(todo):
        embed, dim = build_encoder()
        next_idx = len(chunks)
        buf_E, buf_m, nb = [], [], 0
        paths = todo["image_path"].tolist(); sids = todo["spot_id"].tolist()
        samps = todo["sample_id"].tolist(); subs = todo.get("subdiagnosis", pd.Series([""]*len(todo))).tolist()
        for i in range(0, len(paths), BATCH):
            imgs = [Image.open(ROOT / p).convert("RGB") for p in paths[i:i+BATCH]]
            buf_E.append(embed(imgs))
            for j in range(len(imgs)):
                buf_m.append((sids[i+j], samps[i+j], subs[i+j]))
            nb += 1
            if nb % FLUSH_EVERY == 0:
                ch = pd.DataFrame(np.vstack(buf_E), columns=[f"e{k}" for k in range(dim)])
                ch["spot_id"] = [m[0] for m in buf_m]; ch["sample_id"] = [m[1] for m in buf_m]
                ch["subdiagnosis"] = [m[2] for m in buf_m]
                tmp = PDIR / f".c{next_idx:04d}.tmp"; ch.to_parquet(tmp, index=False)
                tmp.rename(PDIR / f"chunk_{next_idx:04d}.parquet"); next_idx += 1
                buf_E, buf_m = [], []
                print(f"   embedded {i+len(imgs)}/{len(paths)}", flush=True)
        if buf_E:
            ch = pd.DataFrame(np.vstack(buf_E), columns=[f"e{k}" for k in range(dim)])
            ch["spot_id"] = [m[0] for m in buf_m]; ch["sample_id"] = [m[1] for m in buf_m]
            ch["subdiagnosis"] = [m[2] for m in buf_m]
            tmp = PDIR / f".c{next_idx:04d}.tmp"; ch.to_parquet(tmp, index=False)
            tmp.rename(PDIR / f"chunk_{next_idx:04d}.parquet")

    df = pd.concat([pd.read_parquet(c) for c in sorted(PDIR.glob("chunk_*.parquet"))], ignore_index=True)
    df.to_parquet(OUT, index=False)
    for c in PDIR.glob("chunk_*.parquet"):
        c.unlink()
    PDIR.rmdir()
    print(f"[ok] wrote {OUT} ({len(df):,} embeddings)", flush=True)


if __name__ == "__main__":
    main()
