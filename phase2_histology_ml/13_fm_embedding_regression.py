#!/usr/bin/env python3
"""B6: Foundation-model tile embeddings -> spot composition (LOTO).

Tests the *representation* hypothesis behind Phase B's hand-crafted-feature negative:
replace the 7 morphology scalars with a pretrained pathology-FM embedding of each spot
tile, then run the SAME leave-one-tumor-out composition regression (script 12) with the
same shuffled/random negative controls.

Encoders (ungated, PyTorch — no TensorFlow needed):
  * phikon     : owkin/phikon       (ViT-B, 768-d)  [default; fast on CPU]
  * phikon-v2  : owkin/phikon-v2    (ViT-L, 1024-d) [stronger, slower]
  * resnet50   : torchvision ImageNet (2048-d)      [fallback; needs pytorch.org]
XMAG (5x-native, the preferred match for Visium-hires) is wired as 'xmag' and used
automatically once its weights are public.

Pilot-friendly: --max-spots-per-sample subsamples tiles so CPU inference is minutes.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

CELL_STATES = ["blastemal", "epithelial", "stromal"]
HF_ENCODERS = {"phikon": "owkin/phikon", "phikon-v2": "owkin/phikon-v2",
               "xmag": "owkin/phikon"}  # xmag placeholder -> swap repo id when released


def _load_reg_helpers():
    """Import loto_regression + zscore_softmax_target from script 12 (numeric-prefixed)."""
    path = Path(__file__).resolve().parent / "12_spot_composition_regression.py"
    spec = importlib.util.spec_from_file_location("spot_reg", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_encoder(name: str, device: str):
    """Return (embed_fn(list[PIL.Image]) -> np.ndarray, dim, used_name)."""
    if name in HF_ENCODERS or name == "auto":
        try:
            from transformers import AutoImageProcessor, AutoModel
            repo = HF_ENCODERS.get(name, HF_ENCODERS["phikon"])
            proc = AutoImageProcessor.from_pretrained(repo)
            model = AutoModel.from_pretrained(repo).to(device).eval()

            @torch.no_grad()
            def embed(imgs):
                px = proc(images=imgs, return_tensors="pt")["pixel_values"].to(device)
                out = model(pixel_values=px).last_hidden_state[:, 0]  # CLS token
                return out.cpu().numpy()

            dim = model.config.hidden_size
            return embed, dim, (name if name != "auto" else "phikon")
        except Exception as exc:
            print(f"[warn] HF encoder '{name}' unavailable ({exc}); trying resnet50")
    # torchvision ImageNet fallback
    import torchvision.transforms as T
    from torchvision.models import ResNet50_Weights, resnet50
    w = ResNet50_Weights.IMAGENET1K_V2
    net = resnet50(weights=w); net.fc = torch.nn.Identity(); net = net.to(device).eval()
    tf = w.transforms()

    @torch.no_grad()
    def embed(imgs):
        batch = torch.stack([tf(im) for im in imgs]).to(device)
        return net(batch).cpu().numpy()

    return embed, 2048, "resnet50"


def main() -> None:
    setup_logging()
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", default="phikon",
                    choices=["phikon", "phikon-v2", "xmag", "resnet50", "auto"])
    ap.add_argument("--max-spots-per-sample", type=int, default=40)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "fm_embedding_regression")
    reg = _load_reg_helpers()

    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    manifest = json.loads((tiles_dir / "tiles_manifest.json").read_text())
    man = pd.DataFrame(manifest)
    # subsample spots per tumor for a CPU-feasible pilot
    rng = np.random.default_rng(seed)
    man = (man.groupby("sample_id", group_keys=False)
              .apply(lambda d: d.sample(min(len(d), args.max_spots_per_sample), random_state=seed)))
    print(f"[data] {len(man):,} tiles across {man['sample_id'].nunique()} tumors", flush=True)

    embed, dim, used = build_encoder(args.encoder, args.device)
    print(f"[encoder] {used} (dim={dim}) on {args.device}", flush=True)

    from PIL import Image
    root = Path(cfg["root"])
    rows, embs = [], []
    paths = man["image_path"].tolist()
    spot_ids = man["spot_id"].tolist()
    samples = man["sample_id"].tolist()
    for i in range(0, len(paths), args.batch_size):
        chunk = paths[i:i + args.batch_size]
        imgs = [Image.open(root / p).convert("RGB") for p in chunk]
        embs.append(embed(imgs))
        for j in range(len(chunk)):
            rows.append({"spot_id": spot_ids[i + j], "sample_id": samples[i + j]})
        if (i // args.batch_size) % 10 == 0:
            print(f"   embedded {i + len(chunk)}/{len(paths)}", flush=True)
    E = np.vstack(embs)
    emb_df = pd.DataFrame(E, columns=[f"e{i}" for i in range(dim)])
    emb_df["spot_id"] = [r["spot_id"] for r in rows]
    emb_df["sample_id"] = [r["sample_id"] for r in rows]

    sig = pd.read_parquet(resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"]))
    target = reg.zscore_softmax_target(sig)
    df = emb_df.merge(target, on="spot_id", how="inner").dropna()
    feat_cols = [c for c in df.columns if c.startswith("e")]
    print(f"[merge] {len(df):,} spots with embedding+target, {len(feat_cols)} dims", flush=True)

    real = reg.loto_regression(df, feat_cols, seed, max_train=20000, tag="fm-real")
    shuf = reg.loto_regression(df, feat_cols, seed, shuffle_y=True, fold_limit=6, tag="fm-shuf")
    rand = reg.loto_regression(df, feat_cols, seed, random_feat=True, fold_limit=6, tag="fm-rand")

    out = {
        "encoder": used, "embedding_dim": dim,
        "n_spots": int(len(df)), "n_samples": int(df["sample_id"].nunique()),
        "max_spots_per_sample": args.max_spots_per_sample,
        "metric": "mean held-out Pearson r (LOTO) per compartment [mean, std, n_folds]",
        "real": real, "negative_control_shuffled": shuf, "negative_control_random": rand,
        "seed": seed,
    }
    out_json = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"]).parent / f"fm_embedding_regression_{used}.json"
    out_json.write_text(json.dumps(out, indent=2))
    print(f"[ok] FM embedding regression ({used}):")
    for s in CELL_STATES:
        print(f"   {s:11} real r={real[s][0]:.3f}±{real[s][1]:.3f}  "
              f"shuf={shuf[s][0]:.3f}  rand={rand[s][0]:.3f}")
    print(f"[ok] -> {out_json}")


if __name__ == "__main__":
    main()
