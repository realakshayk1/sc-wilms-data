#!/usr/bin/env python3
"""RIGOR pass on the regressive-tissue pilot — everything answerable with local data.

Addresses the honest caveats from the pilot review:
  A. UNCERTAINTY   : bootstrap 95% CI + permutation p on the tumor-level r (was point est).
  B. QUALITY CONFOUND (the important one): could H&E->%regressive be reading tissue/stain
     quality, not necrosis? Compute per-tile quality scalars (brightness, focus, saturation,
     hematoxylin) and test whether the embedding still predicts %regressive BEYOND them.
  C. VALIDATOR DE-CONFOUND: mito% was offered as "independent", but genes/spot is collinear
     with UMI. Test whether mito% separates regressive vs viable AT MATCHED UMI depth.
  D. SPATIAL ARTIFACT: is regressive clustering just tissue-edge/coverage gradient? Test edge
     enrichment and recompute clustering z on INTERIOR spots only.
  E. THRESHOLD SENSITIVITY: sweep the label thresholds; report how ratio / r / mito-effect move.
  F. DeLong CIs on the classification AUCs.

Not fixable locally (stated, not run): pathologist ground truth; per-spot resolution ceiling;
post-chemo power. Reuses cached artifacts from scripts 19/20.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
SPATIAL_ROOT = ROOT / "data" / "raw" / "scpca_downloads" / "spaceranger" / "SCPCP000006_spatial"
TILES = PROC / "he_tiles"
OUT_DIR = ROOT / "results" / "regressive_pilot"
FIG_DIR = ROOT / "results" / "figures"
SEED = 42
rng = np.random.default_rng(SEED)


def load_pilot():
    p = Path(__file__).resolve().parent / "19_regressive_tissue_pilot.py"
    spec = importlib.util.spec_from_file_location("pilot", p); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


def loto_pred(X, y, seed=SEED, randfeat=False):
    """Leave-one-tumor-out RF regression; returns held-out predictions."""
    n = len(y); pred = np.full(n, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        Xtr = rng.normal(size=(tr.sum(), X.shape[1])) if randfeat else X[tr]
        Xte = rng.normal(size=(1, X.shape[1])) if randfeat else X[i:i+1]
        ncomp = min(10, X.shape[1], tr.sum() - 1)
        m = make_pipeline(StandardScaler(), PCA(ncomp, random_state=seed),
                          RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                                random_state=seed, n_jobs=-1))
        m.fit(Xtr, y[tr]); pred[i] = m.predict(Xte)[0]
    return pred


def tile_quality(spot_ids, manifest):
    """Per-tile H&E quality scalars: brightness, focus(Laplacian var), saturation, hematoxylin."""
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
        # hematoxylin ~ blue-purple darkness: inverted grayscale weighted to blue channel
        hema = 255.0 - img[:, :, 0].mean()  # B channel is dark where hematoxylin-rich
        rows.append((sid, gray.mean(), cv2.Laplacian(gray, cv2.CV_64F).var(),
                     hsv[:, :, 1].mean(), hema))
        if k % 2000 == 0:
            print(f"   quality {k}/{len(spot_ids)}", flush=True)
    return pd.DataFrame(rows, columns=["spot_id", "brightness", "focus", "saturation", "hematoxylin"])


def partial_corr(y, a, control):
    """Partial correlation of y and a, controlling for `control` (all 1-D)."""
    def resid(v):
        A = np.column_stack([np.ones_like(control), control])
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta
    ry, ra = resid(y), resid(a)
    return float(pearsonr(ry, ra)[0])


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pilot = load_pilot()
    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    qc = pilot.compute_spot_qc()
    lab = pilot.label_regressive(sig, qc)
    tx = pilot.load_treatment_map(); lab = lab.merge(tx, on="sample_id", how="left")

    per_tumor = (lab.groupby("sample_id")
                    .agg(frac_regressive=("regressive", "mean"), treatment=("treatment", "first"))
                    .reset_index())
    emb = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    tum_emb = emb.groupby("sample_id")[ecols].mean().reset_index()
    df = tum_emb.merge(per_tumor, on="sample_id", how="inner").dropna(subset=["frac_regressive"])
    E = df[ecols].to_numpy(); y = df["frac_regressive"].to_numpy(); n = len(df)
    print(f"[data] {n} tumors with mean-embedding + regressive fraction", flush=True)

    out = {"n_tumors": int(n), "seed": SEED}

    # ---------- A. bootstrap CI + permutation p on tumor-level r ----------
    print("[A] bootstrap CI + permutation p on tumor-level r ...", flush=True)
    p_emb = loto_pred(E, y)
    r_emb = float(pearsonr(y, p_emb)[0])
    boot = []
    for _ in range(2000):
        idx = rng.integers(0, n, n)
        if np.std(y[idx]) < 1e-9 or np.std(p_emb[idx]) < 1e-9:
            continue
        boot.append(pearsonr(y[idx], p_emb[idx])[0])
    boot = np.array(boot)
    perm = np.array([pearsonr(rng.permutation(y), p_emb)[0] for _ in range(5000)])
    out["A_tumor_level_r"] = {
        "r": r_emb, "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "perm_p": float((np.abs(perm) >= abs(r_emb)).mean())}
    print(f"    r={r_emb:.3f}  95% CI [{out['A_tumor_level_r']['ci95'][0]:.3f}, "
          f"{out['A_tumor_level_r']['ci95'][1]:.3f}]  perm p={out['A_tumor_level_r']['perm_p']:.4f}", flush=True)

    # ---------- B. tissue-quality confound control ----------
    print("[B] tissue-quality confound control ...", flush=True)
    manifest = pd.DataFrame(json.loads((TILES / "tiles_manifest.json").read_text()))
    q = tile_quality(emb["spot_id"].tolist(), manifest).merge(
        emb[["spot_id", "sample_id"]], on="spot_id", how="left")
    qcols = ["brightness", "focus", "saturation", "hematoxylin"]
    tum_q = q.groupby("sample_id")[qcols].mean().reset_index()
    dq = df.merge(tum_q, on="sample_id", how="inner")
    Q = dq[qcols].to_numpy(); yq = dq["frac_regressive"].to_numpy()
    Eq = dq[ecols].to_numpy()
    p_qual = loto_pred(Q, yq)
    p_emb2 = loto_pred(Eq, yq)
    r_qual = float(pearsonr(yq, p_qual)[0]); r_emb_q = float(pearsonr(yq, p_emb2)[0])
    # does the embedding prediction explain %regressive BEYOND the quality prediction?
    pr = partial_corr(yq, p_emb2, p_qual)
    out["B_quality_confound"] = {
        "r_embedding_only": r_emb_q, "r_quality_only": r_qual,
        "partial_r_embedding_given_quality": pr,
        "reading": ("embedding adds beyond quality" if pr > 0.15 else
                    "embedding largely explained by quality" if r_qual >= r_emb_q - 0.05 else "mixed")}
    print(f"    r_emb={r_emb_q:.3f}  r_quality={r_qual:.3f}  partial_r(emb|quality)={pr:.3f} "
          f"-> {out['B_quality_confound']['reading']}", flush=True)

    # ---------- C. mito% at matched UMI depth ----------
    print("[C] mito% de-confounded from UMI depth ...", flush=True)
    m = lab.dropna(subset=["pct_mito", "total_counts"]).copy()
    m["logumi"] = np.log10(m["total_counts"] + 1)
    pr_mito = partial_corr(m["regressive"].to_numpy(float), m["pct_mito"].to_numpy(),
                           m["logumi"].to_numpy())
    # within UMI-decile medians
    m["umi_bin"] = pd.qcut(m["logumi"], 10, labels=False, duplicates="drop")
    strat = []
    for b, d in m.groupby("umi_bin"):
        r_, v_ = d[d.regressive == 1]["pct_mito"], d[d.regressive == 0]["pct_mito"]
        if len(r_) >= 10 and len(v_) >= 10:
            strat.append({"umi_decile": int(b), "mito_regressive": float(r_.median()),
                          "mito_viable": float(v_.median()), "n_reg": int(len(r_))})
    frac_bins_higher = float(np.mean([s["mito_regressive"] > s["mito_viable"] for s in strat])) if strat else float("nan")
    out["C_mito_matched_depth"] = {"partial_r_regressive_mito_given_umi": pr_mito,
                                   "frac_umi_deciles_mito_higher_in_regressive": frac_bins_higher,
                                   "n_deciles_tested": len(strat)}
    print(f"    partial r(regressive, mito | logUMI)={pr_mito:.3f}; "
          f"mito higher in regressive in {frac_bins_higher*100:.0f}% of UMI deciles", flush=True)

    # ---------- D. spatial edge artifact ----------
    print("[D] spatial edge-artifact test ...", flush=True)
    names = ["barcode", "in_tissue", "array_row", "array_col", "pr", "pc"]
    frames = []
    import glob
    for f in glob.glob(str(SPATIAL_ROOT / "*" / "*_spatial" / "spatial" / "tissue_positions_list.csv")):
        libid = Path(f).parents[1].name.replace("_spatial", "")
        d = pd.read_csv(f, header=None, names=names); d["library_id"] = libid
        frames.append(d[["library_id", "barcode", "array_row", "array_col"]])
    lab2 = lab.merge(pd.concat(frames, ignore_index=True), on=["library_id", "barcode"], how="left")
    edge_frac_reg, edge_frac_via, z_all, z_int = [], [], [], []
    for libid, d in lab2.dropna(subset=["array_row"]).groupby("library_id"):
        xy = d[["array_row", "array_col"]].to_numpy(float)
        if len(xy) < 60:
            continue
        _, idx = cKDTree(xy).query(xy, k=7)
        neigh = idx[:, 1:]
        # neighbor distance ~1 for adjacent Visium spots; edge = few close neighbors
        d0 = cKDTree(xy)
        close = np.array([len(d0.query_ball_point(p, r=2.0)) - 1 for p in xy])
        is_edge = close < 5
        reg = d["regressive"].to_numpy().astype(bool)
        if reg.sum() >= 5:
            edge_frac_reg.append(is_edge[reg].mean()); edge_frac_via.append(is_edge[~reg].mean())
        def clz(mask):
            sub = xy[mask]; r = reg[mask]
            if r.sum() < 5 or (~r).sum() < 5 or len(sub) < 20:
                return np.nan
            _, ii = cKDTree(sub).query(sub, k=min(7, len(sub)))
            ii = ii[:, 1:]
            same = lambda l: np.mean(l[ii] == l[:, None])
            null = np.array([same(rng.permutation(r)) for _ in range(100)])
            return (same(r) - null.mean()) / (null.std() + 1e-9)
        z_all.append(clz(np.ones(len(xy), bool))); z_int.append(clz(~is_edge))
    out["D_spatial_edge"] = {
        "edge_frac_regressive": float(np.nanmean(edge_frac_reg)),
        "edge_frac_viable": float(np.nanmean(edge_frac_via)),
        "clustering_z_all_median": float(np.nanmedian(z_all)),
        "clustering_z_interior_only_median": float(np.nanmedian(z_int))}
    print(f"    edge fraction: regressive={out['D_spatial_edge']['edge_frac_regressive']:.2f} vs "
          f"viable={out['D_spatial_edge']['edge_frac_viable']:.2f}; "
          f"clustering z all={out['D_spatial_edge']['clustering_z_all_median']:.1f} -> "
          f"interior-only={out['D_spatial_edge']['clustering_z_interior_only_median']:.1f}", flush=True)

    # ---------- E. threshold sensitivity sweep ----------
    print("[E] label-threshold sensitivity sweep ...", flush=True)
    sweep = []
    for cut in [0.20, 0.25, 0.30, 0.35, 0.40]:
        for qd in [0.30, 0.40, 0.50]:
            orig_cut, orig_q = pilot.REL_UMI_CUT, pilot.VIABLE_DEPLETE_Q
            pilot.REL_UMI_CUT, pilot.VIABLE_DEPLETE_Q = cut, qd
            L = pilot.label_regressive(sig, qc)
            pilot.REL_UMI_CUT, pilot.VIABLE_DEPLETE_Q = orig_cut, orig_q
            pt = L.groupby("sample_id")["regressive"].mean().rename("frac").reset_index()
            dd = tum_emb.merge(pt, on="sample_id", how="inner").dropna()
            yy = dd["frac"].to_numpy()
            rr = float(pearsonr(yy, loto_pred(dd[ecols].to_numpy(), yy))[0]) if yy.std() > 1e-6 else float("nan")
            mm = partial_corr(L["regressive"].to_numpy(float),
                              L["pct_mito"].fillna(L["pct_mito"].median()).to_numpy(),
                              np.log10(L["total_counts"].to_numpy() + 1))
            sweep.append({"umi_cut": cut, "deplete_q": qd, "pct_regressive": float(L["regressive"].mean()),
                          "tumor_r": rr, "mito_partial_r": float(mm)})
            print(f"    cut={cut} q={qd}: %reg={L['regressive'].mean()*100:4.1f}  r={rr:+.3f}  mito_pr={mm:+.3f}", flush=True)
    out["E_threshold_sweep"] = sweep
    out["E_summary"] = {"tumor_r_range": [float(min(s["tumor_r"] for s in sweep)),
                                          float(max(s["tumor_r"] for s in sweep))],
                        "pct_regressive_range": [float(min(s["pct_regressive"] for s in sweep)),
                                                 float(max(s["pct_regressive"] for s in sweep))]}

    # ---------- F. DeLong CI on high/low tumor AUC ----------
    print("[F] DeLong CI on high/low tumor AUC ...", flush=True)
    stats = load_pilot()  # reuse
    from phase_b_stats import delong_auc_ci
    ybin = (y > np.median(y)).astype(int)
    dl = delong_auc_ci(ybin, p_emb)
    out["F_highlow_auc"] = {"auc": float(dl["auc"]), "ci95": [float(dl["ci_low"]), float(dl["ci_high"])]}
    print(f"    high/low AUC={dl['auc']:.3f}  DeLong 95% CI [{dl['ci_low']:.3f}, {dl['ci_high']:.3f}]", flush=True)

    (OUT_DIR / "regressive_rigor.json").write_text(json.dumps(out, indent=2))
    print("\n[ok] wrote", OUT_DIR / "regressive_rigor.json", flush=True)


if __name__ == "__main__":
    main()
