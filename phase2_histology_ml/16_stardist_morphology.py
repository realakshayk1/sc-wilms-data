#!/usr/bin/env python3
"""B-3 + B-2: StarDist nuclear morphology, and morphology+embedding ensemble.

Watershed morphometry gave tumor-level histology AUC 0.39 (worse than chance) — the
segmentation, not the hypothesis, was the bottleneck (AGENTS.md: switching segmentation
tool must be recorded). StarDist '2D_versatile_he' is a learned H&E nuclei model. We:

  1. Subsample N spots/tumor (seed-aligned with the embedding run), segment with StarDist
     in-memory (no 260k-tile mask dump), extract per-nucleus morphology.
  2. Build per-tumor anaplasia-proxy features (giant-nucleus fraction at 3x-median, p99
     area, pleomorphism = area CV, hyperchromasia tail) and a LOTO RandomForest -> AUC.
  3. ENSEMBLE: concat tumor-mean phikon embedding + StarDist morphology -> LOTO logistic.
  4. DeLong CI + permutation p for each, vs the watershed-morphology 0.39 baseline.

Honest by construction: if StarDist morphology still underperforms the embedding, the JSON
shows it. CPU-heavy — run after the embedding job frees cores.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from utils import load_config, resolve_path, set_seed_logged, setup_logging
from phase_b_stats import delong_auc_ci, permutation_auc_p, delong_paired_test

MORPH_COLS = ["area", "major_axis_length", "eccentricity", "solidity",
              "hematoxylin_intensity", "texture_var"]


def _load(mod_file, name):
    p = Path(__file__).resolve().parent / mod_file
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def segment_and_featurize(cfg, n_per_sample, seed, max_tiles=None, flush_every=50, prob_thresh=0.4):
    """Subsample tiles, StarDist-segment, return per-nucleus feature DataFrame.
    CHECKPOINTS per-nucleus rows to append-only chunk parquets (this env tears the
    process down periodically; StarDist on CPU is slow). Resume skips done spot_ids.
    prob_thresh below StarDist's default 0.69 — these 96px Visium-HIRES tiles are
    low-res, so the default detects ~0 nuclei; lowering it recovers detections."""
    import cv2
    seg = _load("02_segment_nuclei.py", "segmod")
    feat = _load("03_nucleus_features.py", "featmod")

    proc_dir = resolve_path(cfg, cfg["paths"]["dirs"]["processed"])
    pdir = proc_dir / f"nucleus_features_stardist_{n_per_sample}_pt{int(prob_thresh*100)}_partial"
    pdir.mkdir(parents=True, exist_ok=True)

    tiles_dir = resolve_path(cfg, cfg["paths"]["phase_b"]["tiles_dir"])
    man = pd.DataFrame(json.loads((tiles_dir / "tiles_manifest.json").read_text()))
    man = (man.groupby("sample_id", group_keys=False)
              .apply(lambda d: d.sample(min(len(d), n_per_sample), random_state=seed)))
    if max_tiles:
        man = man.iloc[:max_tiles]

    chunks = sorted(pdir.glob("chunk_*.parquet"))
    done = set()
    if chunks:
        done = set(pd.concat([pd.read_parquet(c, columns=["spot_id"]) for c in chunks])["spot_id"])
        print(f"[resume] {len(done):,} tiles segmented across {len(chunks)} chunks", flush=True)
    next_idx = len(chunks)
    todo = man[~man["spot_id"].isin(done)]
    print(f"[data] StarDist: {len(todo):,} of {len(man):,} tiles remaining / "
          f"{man['sample_id'].nunique()} tumors", flush=True)

    model = seg.get_stardist_model()  # 2D_versatile_he
    root = Path(cfg["root"])
    buf, n_since = [], 0
    for k, (_, e) in enumerate(todo.iterrows()):
        bgr = cv2.imread(str(root / e["image_path"]))
        if bgr is None:
            continue
        rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
        labels, _ = model.predict_instances(rgb, prob_thresh=prob_thresh, nms_thresh=0.3,
                                            show_tile_progress=False)
        meta = {"sample_id": e["sample_id"], "library_id": e.get("library_id", ""),
                "subdiagnosis": e.get("subdiagnosis", "")}
        rows = feat.extract_nucleus_features(rgb, labels.astype(np.int32),
                                             e["spot_id"], str(e.get("subdiagnosis", "")), meta)
        if not rows:  # keep a spot_id marker so resume counts the tile as done
            rows = [{"spot_id": e["spot_id"], "sample_id": e["sample_id"],
                     "subdiagnosis": str(e.get("subdiagnosis", "")), "area": np.nan}]
        buf.extend(rows); n_since += 1
        if n_since >= flush_every:
            tmp = pdir / f".chunk_{next_idx:04d}.tmp.parquet"
            pd.DataFrame(buf).to_parquet(tmp, index=False)
            tmp.rename(pdir / f"chunk_{next_idx:04d}.parquet")
            next_idx += 1; buf, n_since = [], 0
            print(f"   segmented {k + 1}/{len(todo)} (checkpointed)", flush=True)
    if buf:
        tmp = pdir / f".chunk_{next_idx:04d}.tmp.parquet"
        pd.DataFrame(buf).to_parquet(tmp, index=False)
        tmp.rename(pdir / f"chunk_{next_idx:04d}.parquet")

    df = pd.concat([pd.read_parquet(c) for c in sorted(pdir.glob("chunk_*.parquet"))], ignore_index=True)
    df = df[df["area"].notna()]  # drop empty-tile markers
    print(f"[seg] {len(df):,} nuclei (median {int(df.groupby('sample_id').size().median())}/tumor)", flush=True)
    return df, pdir


def tumor_morphology_features(nuc):
    """Per-tumor anaplasia-proxy morphology vector."""
    nuc = nuc[nuc["area"] >= 20]
    rows = []
    for sid, d in nuc.groupby("sample_id"):
        sub = str(d["subdiagnosis"].iloc[0])
        if sub not in ("favorable", "anaplastic"):
            continue  # only the histology contrast tumors carry a valid label
        feat = {"sample_id": sid, "y": int(sub == "anaplastic")}
        for col in MORPH_COLS:
            qs = d[col].quantile([0.5, 0.9, 0.99]).tolist()
            feat[f"{col}_p50"], feat[f"{col}_p90"], feat[f"{col}_p99"] = qs
            feat[f"{col}_cv"] = d[col].std() / (abs(d[col].mean()) + 1e-9)
        feat["giant_frac"] = float((d["area"] > 3 * d["area"].median()).mean())
        feat["n_nuclei"] = len(d)
        rows.append(feat)
    return pd.DataFrame(rows)


def loto_auc(X, y, groups, model_fn, seed):
    from sklearn.model_selection import LeaveOneGroupOut
    proba = np.full(len(y), np.nan)
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        m = model_fn(); m.fit(X[tr], y[tr]); proba[te] = m.predict_proba(X[te])[:, 1]
    return proba


def evaluate(name, y, proba, n_perm, seed):
    ok = ~np.isnan(proba)
    ci = delong_auc_ci(y[ok], proba[ok])
    perm = permutation_auc_p(y[ok], proba[ok], n_perm=n_perm, seed=seed)
    print(f"[{name}] AUC={ci['auc']:.3f}  95% CI [{ci['ci_low']:.3f}, {ci['ci_high']:.3f}]  "
          f"perm p={perm['p_value']:.4f}", flush=True)
    return {"auc": ci["auc"], "ci_low": ci["ci_low"], "ci_high": ci["ci_high"],
            "perm_p": perm["p_value"], "n": int(ok.sum())}


def main():
    setup_logging()
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-sample", type=int, default=200)
    ap.add_argument("--encoder", default="phikon-v2")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--max-tiles", type=int, default=None)
    ap.add_argument("--prob-thresh", type=float, default=0.4)
    args = ap.parse_args()
    cfg = load_config()
    seed = set_seed_logged(cfg["features"]["seed"], "stardist_morphology")
    proc_dir = resolve_path(cfg, cfg["paths"]["dirs"]["processed"])

    # 1. StarDist features (resumable chunked cache, keyed by prob_thresh)
    cache = proc_dir / f"nucleus_features_stardist_{args.n_per_sample}_pt{int(args.prob_thresh*100)}.parquet"
    if cache.exists():
        nuc = pd.read_parquet(cache); print(f"[cache] {len(nuc):,} StarDist nuclei", flush=True)
    else:
        nuc, pdir = segment_and_featurize(cfg, args.n_per_sample, seed, args.max_tiles,
                                          prob_thresh=args.prob_thresh)
        nuc.to_parquet(cache, index=False)
        for c in pdir.glob("chunk_*.parquet"):
            c.unlink()
        pdir.rmdir()
        print(f"[cache] wrote {cache}", flush=True)

    mdf = tumor_morphology_features(nuc)  # already restricted to favorable/anaplastic tumors
    mcols = [c for c in mdf.columns if c not in ("sample_id", "y")]
    med_nuc = int(nuc.groupby("sample_id").size().median())

    results = {"encoder": args.encoder, "n_per_sample": args.n_per_sample,
               "prob_thresh": args.prob_thresh, "n_nuclei": int(len(nuc)),
               "median_nuclei_per_tumor": med_nuc,
               "n_tumors": int(mdf["sample_id"].nunique()),
               "baseline_watershed_morph_auc_ref": 0.393, "seed": seed, "models": {}}

    # 2. StarDist morphology classifier (median-impute: tumors with few nuclei -> NaN feats)
    Xm = mdf[mcols].to_numpy(); ym = mdf["y"].to_numpy().astype(int); gm = mdf["sample_id"].to_numpy()
    rf = lambda: make_pipeline(SimpleImputer(strategy="median"),
                               RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                                      class_weight="balanced", random_state=seed, n_jobs=-1))
    p_morph = loto_auc(Xm, ym, gm, rf, seed)
    results["models"]["stardist_morphology_rf"] = evaluate("StarDist morphology", ym, p_morph, args.n_perm, seed)

    # 3. Ensemble with tumor-mean phikon-v2 embedding
    emb_cache = proc_dir / f"phikon_spot_embeddings_{args.encoder}.parquet"
    if emb_cache.exists():
        emb = pd.read_parquet(emb_cache)
        emb = emb[emb["subdiagnosis"].isin(["favorable", "anaplastic"])]
        ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
        tum = emb.groupby("sample_id")[ecols].mean().reset_index()
        ens = mdf.merge(tum, on="sample_id", how="inner")
        Xe = ens[mcols + ecols].to_numpy(); ye = ens["y"].to_numpy().astype(int); ge = ens["sample_id"].to_numpy()
        ens_model = lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                          LogisticRegression(max_iter=3000, C=0.3))
        p_ens = loto_auc(Xe, ye, ge, ens_model, seed)
        results["models"]["ensemble_morph_plus_embedding"] = evaluate("Ensemble morph+emb", ye, p_ens, args.n_perm, seed)
        # embedding-only on the SAME tumors, for a paired DeLong of the ensemble gain
        p_emb = loto_auc(ens[ecols].to_numpy(), ye, ge,
                         lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                               LogisticRegression(max_iter=3000, C=0.3)), seed)
        ok = (~np.isnan(p_ens)) & (~np.isnan(p_emb))
        results["paired_delong_ensemble_vs_embedding"] = delong_paired_test(ye[ok], p_ens[ok], p_emb[ok])
    else:
        print(f"[warn] no embedding cache at {emb_cache}; skipping ensemble", flush=True)

    out = resolve_path(cfg, cfg["paths"]["phase_b"]["classifier_pkl"]).parent / "stardist_morphology.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"[ok] -> {out}", flush=True)


if __name__ == "__main__":
    main()
