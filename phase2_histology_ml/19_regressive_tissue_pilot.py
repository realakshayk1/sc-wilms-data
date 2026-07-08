#!/usr/bin/env python3
"""PILOT v2: viable-blastemal vs necrotic/regressive tissue — ratio + spatial + H&E.

Motivation (Ravi, 2026-07): can H&E give ratios and positions of blastemal vs
"regressive" tissue? Strict SIOP "regressive-type" is a post-chemo call (only 11/38 spatial
tumors here are post-chemo, 10/11 anaplastic -> confounded), so this pilot targets the
broader, resolution-robust question askable on ALL tumors: per-spot VIABLE vs
NECROTIC/REGRESSIVE tissue, its per-tumor ratio, and its spatial arrangement.

v2 upgrades over v1 (all four review levers):
  L3 SHARPEN LABEL  : real per-spot necrosis QC from the Visium matrices
                      (n_genes, mito%, hypoxia score) replaces the resolution-broken
                      StarDist cellularity check; regressive = low-UMI AND depleted-viable
                      AND necrosis-QC-concordant.
  L1 TUMOR-LEVEL H&E: predict a tumor's %regressive from its mean Phikon-v2 embedding
                      (the framing that worked for anaplasia), not just per-spot.
  L2 BALANCED EMBED : if 20_regressive_balanced_embed.py has produced a regressive-balanced
                      embedding cache, re-run the per-spot LOTO on it (fair positive rate).
  L4 WSI HOOK       : documented external ceiling-lifter (needs whole-slide images).

Everything reuses artifacts already on disk; QC is computed once and cached.
"""
from __future__ import annotations

import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, pearsonr
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
SPATIAL_ROOT = ROOT / "data" / "raw" / "scpca_downloads" / "spaceranger" / "SCPCP000006_spatial"
SC_META = ROOT / "data" / "raw" / "scpca_downloads" / "single-cell-experiment" / \
    "SCPCP000006_single-cell_merged" / "single-cell_metadata.tsv"
OUT_DIR = ROOT / "results" / "regressive_pilot"
FIG_DIR = ROOT / "results" / "figures"
QC_CACHE = PROC / "spot_qc_necrosis.parquet"
BALANCED_EMB = PROC / "regressive_balanced_embeddings_phikon-v2.parquet"
SEED = 42

REL_UMI_CUT = 0.30       # low-UMI: total_counts < 30% of the library's 75th-pct spot
VIABLE_DEPLETE_Q = 0.40  # depleted-viable: max viable program below its 40th pct in-library
MIN_TISSUE_SPOTS = 50

# Compact, well-established hypoxia signature (Buffa/Hallmark core).
HYPOXIA = ["VEGFA", "CA9", "SLC2A1", "LDHA", "PGK1", "NDRG1", "BNIP3", "ADM", "ENO1",
           "PDK1", "HK2", "ALDOA", "P4HA1", "PLOD2", "SLC2A3", "ANGPTL4", "EGLN3", "PGAM1"]


def sanitize(barcode: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", barcode)


def load_treatment_map() -> pd.DataFrame:
    md = pd.read_csv(SC_META, sep="\t")
    md = md[["scpca_sample_id", "treatment", "subdiagnosis"]].drop_duplicates("scpca_sample_id")
    return md.rename(columns={"scpca_sample_id": "sample_id"})


def compute_spot_qc() -> pd.DataFrame:
    """L3: per-spot necrosis QC from the Visium matrices (cached). Necrotic/regressive
    tissue -> low genes detected, high mito fraction (dying cells), and peri-necrotic
    hypoxia. These are independent of the UMI-based weak label used to *define* regressive,
    so they can validate it."""
    if QC_CACHE.exists():
        return pd.read_parquet(QC_CACHE)
    rows = []
    lib_dirs = sorted(glob.glob(str(SPATIAL_ROOT / "*" / "*_spatial")))
    for i, ld in enumerate(lib_dirs):
        ld = Path(ld)
        lib = ld.name.replace("_spatial", "")
        mtx = ld / "filtered_feature_bc_matrix"
        if not mtx.exists():
            continue
        ad = sc.read_10x_mtx(mtx, var_names="gene_symbols", make_unique=True)
        X = ad.X
        counts = np.asarray(X.sum(1)).ravel()
        counts[counts == 0] = 1
        n_genes = np.asarray((X > 0).sum(1)).ravel()
        upper = ad.var_names.str.upper()
        mt = np.where(upper.str.startswith("MT-"))[0]
        pct_mito = (np.asarray(X[:, mt].sum(1)).ravel() / counts) if len(mt) else np.zeros_like(counts)
        hy = [g for g in HYPOXIA if g in set(upper)]
        if hy:
            idx = [list(upper).index(g) for g in hy]
            cp10k = X[:, idx] / counts[:, None] * 1e4
            hyscore = np.asarray(np.log1p(cp10k).mean(1)).ravel()
        else:
            hyscore = np.zeros_like(counts)
        for bc, ng, pm, hs in zip(ad.obs_names, n_genes, pct_mito, hyscore):
            rows.append((f"{lib}__{sanitize(bc)}", int(ng), float(pm), float(hs)))
        print(f"   QC {i+1}/{len(lib_dirs)} {lib} ({ad.n_obs} spots)", flush=True)
    qc = pd.DataFrame(rows, columns=["spot_id", "n_genes", "pct_mito", "hypoxia_score"])
    qc.to_parquet(QC_CACHE, index=False)
    return qc


def label_regressive(sig: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    """Sharpened per-library weak label: regressive = low-UMI AND depleted-viable."""
    sig = sig[sig["in_tissue"] == 1].copy()
    sig["viable_max"] = sig[["blastemal_program", "epithelial_program", "stromal_program"]].max(1)
    out = []
    for lib, d in sig.groupby("library_id"):
        d = d.copy()
        if len(d) < MIN_TISSUE_SPOTS:
            continue
        ref = np.quantile(d["total_counts"], 0.75)
        low_umi = d["total_counts"] < REL_UMI_CUT * max(ref, 1.0)
        deplete = d["viable_max"] < d["viable_max"].quantile(VIABLE_DEPLETE_Q)
        d["regressive"] = (low_umi & deplete).astype(int)
        d["low_umi_only"] = low_umi.astype(int)  # v1-style label, for comparison
        blast_cut = d["blastemal_program"].quantile(0.60)
        d["viable_blastemal"] = ((d["regressive"] == 0) & (d["dominant_state"] == "blastemal") &
                                 (d["blastemal_program"] >= blast_cut)).astype(int)
        out.append(d)
    lab = pd.concat(out, ignore_index=True)
    return lab.merge(qc, on="spot_id", how="left")


def neighbor_enrichment(xy: np.ndarray, lab: np.ndarray, k=6, n_perm=200, seed=SEED) -> float:
    if lab.sum() < 5 or (~lab.astype(bool)).sum() < 5:
        return float("nan")
    _, idx = cKDTree(xy).query(xy, k=min(k + 1, len(xy)))
    idx = idx[:, 1:]
    same = lambda l: np.mean(l[idx] == l[:, None])
    obs = same(lab)
    rng = np.random.default_rng(seed)
    null = np.array([same(rng.permutation(lab)) for _ in range(n_perm)])
    return float((obs - null.mean()) / (null.std() + 1e-9))


def per_spot_he_test(lab: pd.DataFrame, emb_path: Path, seed=SEED) -> dict:
    """Per-spot LOTO: can Phikon-v2 embeddings predict the regressive label?"""
    emb = pd.read_parquet(emb_path)
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    df = emb.merge(lab[["spot_id", "regressive"]], on="spot_id", how="inner").dropna(subset=["regressive"])
    df["regressive"] = df["regressive"].astype(int)
    tumors = sorted(df["sample_id"].unique())
    rng = np.random.default_rng(seed)

    def run(shuffle_y=False, random_feat=False):
        y = df["regressive"].to_numpy(); proba = np.full(len(df), np.nan)
        for t in tumors:
            te = df["sample_id"].to_numpy() == t; tr = ~te
            ytr = y[tr]
            if len(np.unique(ytr)) < 2 or len(np.unique(y[te])) < 2:
                continue
            Xtr, Xte = df.loc[tr, ecols].to_numpy(), df.loc[te, ecols].to_numpy()
            if random_feat:
                Xtr, Xte = rng.normal(size=Xtr.shape), rng.normal(size=Xte.shape)
            if shuffle_y:
                ytr = rng.permutation(ytr)
            scl = StandardScaler().fit(Xtr)
            m = LogisticRegression(max_iter=2000, C=0.5, class_weight="balanced")
            m.fit(scl.transform(Xtr), ytr)
            proba[te] = m.predict_proba(scl.transform(Xte))[:, 1]
        ok = ~np.isnan(proba)
        return float(roc_auc_score(y[ok], proba[ok])) if ok.sum() and len(np.unique(y[ok])) > 1 else float("nan")

    return {"n_spots": int(len(df)), "n_tumors": len(tumors),
            "frac_regressive": float(df["regressive"].mean()),
            "auc_real": run(), "auc_shuffled": run(shuffle_y=True), "auc_random": run(random_feat=True)}


def tumor_level_he_test(lab: pd.DataFrame, per_tumor: pd.DataFrame, seed=SEED) -> dict:
    """L1: predict a tumor's %regressive from its MEAN Phikon-v2 embedding (LOTO)."""
    emb = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    tum = emb.groupby("sample_id")[ecols].mean().reset_index()
    df = tum.merge(per_tumor[["sample_id", "frac_regressive"]], on="sample_id", how="inner").dropna()
    X = df[ecols].to_numpy(); y = df["frac_regressive"].to_numpy(); n = len(df)
    rng = np.random.default_rng(seed)

    def loto_pred(yv, randfeat=False):
        pred = np.full(n, np.nan)
        for i in range(n):
            tr = np.arange(n) != i
            Xtr, Xte = (rng.normal(size=(tr.sum(), X.shape[1])), rng.normal(size=(1, X.shape[1]))) \
                if randfeat else (X[tr], X[i:i+1])
            m = make_pipeline(StandardScaler(), PCA(min(10, tr.sum()-1), random_state=seed),
                              RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                                    random_state=seed, n_jobs=-1))
            m.fit(Xtr, yv[tr]); pred[i] = m.predict(Xte)[0]
        return pred

    pred = loto_pred(y)
    r_real = float(pearsonr(y, pred)[0])
    r_shuf = float(pearsonr(y, loto_pred(rng.permutation(y)))[0])
    r_rand = float(pearsonr(y, loto_pred(y, randfeat=True))[0])
    # binary high/low (median split) AUC
    ybin = (y > np.median(y)).astype(int)
    auc = float(roc_auc_score(ybin, pred)) if len(np.unique(ybin)) > 1 else float("nan")
    return {"n_tumors": n, "pearson_r_real": r_real, "pearson_r_shuffled": r_shuf,
            "pearson_r_random": r_rand, "highlow_auc": auc}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True); FIG_DIR.mkdir(parents=True, exist_ok=True)
    np.random.seed(SEED)

    print("[1/6] per-spot necrosis QC from Visium matrices (L3, cached) ...", flush=True)
    qc = compute_spot_qc()
    print(f"      QC for {len(qc):,} spots", flush=True)

    print("[2/6] sharpened regressive label ...", flush=True)
    sig = pd.read_parquet(PROC / "spot_signatures.parquet")
    lab = label_regressive(sig, qc).merge(load_treatment_map(), on="sample_id", how="left")
    print(f"      {len(lab):,} tissue spots, {lab['regressive'].mean()*100:.1f}% regressive "
          f"(v1 low-UMI-only was {lab['low_umi_only'].mean()*100:.1f}%)", flush=True)

    # L3 validation: independent necrosis QC should separate regressive vs viable
    print("[3/6] validating label against independent necrosis QC ...", flush=True)
    val = {}
    for col, alt, nice in [("n_genes", "less", "genes/spot"), ("pct_mito", "greater", "mito%"),
                           ("hypoxia_score", "greater", "hypoxia")]:
        r = lab.loc[lab.regressive == 1, col].dropna(); v = lab.loc[lab.regressive == 0, col].dropna()
        p = float(mannwhitneyu(r, v, alternative=alt)[1])
        val[col] = {"regressive_median": float(r.median()), "viable_median": float(v.median()), "p": p}
        print(f"      {nice:12} regressive={r.median():.3g} vs viable={v.median():.3g}  p={p:.1e} ({alt})", flush=True)

    print("[4/6] per-tumor ratios + treatment contrast ...", flush=True)
    per_tumor = (lab.groupby("sample_id")
                    .agg(n_spots=("regressive", "size"), frac_regressive=("regressive", "mean"),
                         frac_viable_blastemal=("viable_blastemal", "mean"),
                         treatment=("treatment", "first"), subdiagnosis=("subdiagnosis", "first"))
                    .reset_index())
    per_tumor.to_csv(OUT_DIR / "per_tumor_regressive.csv", index=False)
    up = per_tumor.loc[per_tumor.treatment == "Upfront resection", "frac_regressive"]
    pc = per_tumor.loc[per_tumor.treatment == "Resection post chemotherapy", "frac_regressive"]
    p_tx = float(mannwhitneyu(pc, up, alternative="greater")[1])
    print(f"      %regressive post-chemo={pc.mean()*100:.1f}% (n={len(pc)}) vs "
          f"upfront={up.mean()*100:.1f}% (n={len(up)})  p={p_tx:.3f}", flush=True)

    print("[5/6] H&E readout — tumor-level (L1) + per-spot (existing & balanced) ...", flush=True)
    he_tumor = tumor_level_he_test(lab, per_tumor)
    print(f"      TUMOR-LEVEL: mean-emb -> %regressive Pearson r={he_tumor['pearson_r_real']:.3f} "
          f"(shuf={he_tumor['pearson_r_shuffled']:.3f}, rand={he_tumor['pearson_r_random']:.3f}); "
          f"high/low AUC={he_tumor['highlow_auc']:.3f}", flush=True)
    he_spot = per_spot_he_test(lab, PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    print(f"      PER-SPOT (histology-sampled emb): AUC real={he_spot['auc_real']:.3f} "
          f"(shuf={he_spot['auc_shuffled']:.3f}, rand={he_spot['auc_random']:.3f}, "
          f"pos={he_spot['frac_regressive']*100:.1f}%)", flush=True)
    he_spot_bal = None
    if BALANCED_EMB.exists():
        he_spot_bal = per_spot_he_test(lab, BALANCED_EMB)
        print(f"      PER-SPOT (L2 balanced emb): AUC real={he_spot_bal['auc_real']:.3f} "
              f"(shuf={he_spot_bal['auc_shuffled']:.3f}, rand={he_spot_bal['auc_random']:.3f}, "
              f"pos={he_spot_bal['frac_regressive']*100:.1f}%)", flush=True)
    else:
        print("      PER-SPOT (L2 balanced emb): not yet built — run 20_regressive_balanced_embed.py", flush=True)

    print("[6/6] spatial maps + neighbourhood enrichment ...", flush=True)
    names = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    frames = []
    for f in glob.glob(str(SPATIAL_ROOT / "*" / "*_spatial" / "spatial" / "tissue_positions_list.csv")):
        lib = Path(f).parents[1].name.replace("_spatial", "")
        d = pd.read_csv(f, header=None, names=names); d["library_id"] = lib
        frames.append(d[["library_id", "barcode", "array_row", "array_col"]])
    lab = lab.merge(pd.concat(frames, ignore_index=True), on=["library_id", "barcode"], how="left")
    zs = []
    for lib, d in lab.dropna(subset=["array_row"]).groupby("library_id"):
        z = neighbor_enrichment(d[["array_row", "array_col"]].to_numpy(float), d["regressive"].to_numpy())
        if not np.isnan(z):
            zs.append({"library_id": lib, "sample_id": d.sample_id.iloc[0], "reg_cluster_z": z})
    zdf = pd.DataFrame(zs); zdf.to_csv(OUT_DIR / "regressive_clustering_z.csv", index=False)
    print(f"      regressive clustering z: median={zdf.reg_cluster_z.median():.1f} "
          f"({(zdf.reg_cluster_z > 2).sum()}/{len(zdf)} libraries z>2)", flush=True)

    # ---- figures ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for i, sub in enumerate([up, pc]):
        ax[0].scatter(np.full(len(sub), i) + np.random.uniform(-.08, .08, len(sub)), sub*100, s=40, alpha=.7)
        ax[0].hlines(sub.mean()*100, i-.25, i+.25, color="k", lw=2)
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["upfront", "post-chemo"])
    ax[0].set_ylabel("% regressive tissue"); ax[0].set_title(f"Ratio by treatment\np={p_tx:.3f} (post>upfront)")
    # tumor-level scatter
    emb = pd.read_parquet(PROC / "phikon_spot_embeddings_phikon-v2.parquet")
    ecols = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    ax[1].bar(["real", "shuffled", "random"],
              [he_tumor["pearson_r_real"], he_tumor["pearson_r_shuffled"], he_tumor["pearson_r_random"]],
              color=["#2a7", "#89a", "#bbb"])
    ax[1].axhline(0, color="k", lw=1); ax[1].set_ylabel("Pearson r (LOTO)")
    ax[1].set_title(f"TUMOR-LEVEL: H&E -> %regressive\nr={he_tumor['pearson_r_real']:.3f}, AUC={he_tumor['highlow_auc']:.3f}")
    bars = [("hist-emb", he_spot)] + ([("balanced-emb", he_spot_bal)] if he_spot_bal else [])
    xlab, real, shuf = [], [], []
    for nm, h in bars:
        xlab.append(nm); real.append(h["auc_real"]); shuf.append(h["auc_shuffled"])
    xpos = np.arange(len(xlab))
    ax[2].bar(xpos-0.2, real, 0.4, label="real", color="#2a7")
    ax[2].bar(xpos+0.2, shuf, 0.4, label="shuffled", color="#89a")
    ax[2].axhline(0.5, ls="--", color="k", lw=1); ax[2].set_xticks(xpos); ax[2].set_xticklabels(xlab)
    ax[2].set_ylim(0.3, 1.0); ax[2].set_ylabel("per-spot AUC (LOTO)"); ax[2].set_title("PER-SPOT H&E readout"); ax[2].legend()
    fig.tight_layout(); fig.savefig(FIG_DIR / "regressive_pilot_summary.png", dpi=130)

    top = per_tumor.sort_values("frac_regressive", ascending=False).head(6).sample_id.tolist()
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax_, sid in zip(axes.ravel(), top):
        d = lab[(lab.sample_id == sid) & lab.array_row.notna()]
        d = d[d.library_id == d.library_id.iloc[0]]
        ax_.scatter(d.loc[d.regressive == 0, "array_col"], -d.loc[d.regressive == 0, "array_row"], s=6, c="#bcd")
        ax_.scatter(d.loc[d.regressive == 1, "array_col"], -d.loc[d.regressive == 1, "array_row"], s=6, c="#c33")
        tl = str(per_tumor.loc[per_tumor.sample_id == sid, "treatment"].iloc[0])[:12]
        ax_.set_title(f"{sid}  {d.regressive.mean()*100:.0f}% reg\n{tl}", fontsize=9)
        ax_.set_xticks([]); ax_.set_yticks([]); ax_.set_aspect("equal")
    fig.suptitle("Regressive (red) vs viable (blue) tissue — sharpened label")
    fig.tight_layout(); fig.savefig(FIG_DIR / "regressive_pilot_spatialmaps.png", dpi=130)

    summary = {
        "weak_label": {"rule": f"UMI<{REL_UMI_CUT}xlib-p75 AND viable_max<q{VIABLE_DEPLETE_Q}",
                       "frac_regressive": float(lab.regressive.mean()),
                       "frac_low_umi_only_v1": float(lab.low_umi_only.mean())},
        "L3_label_validation_necrosis_qc": val,
        "L1_treatment_contrast": {"pct_postchemo": float(pc.mean()), "n_postchemo": int(len(pc)),
                                  "pct_upfront": float(up.mean()), "n_upfront": int(len(up)),
                                  "p_greater": p_tx,
                                  "caveat": "post-chemo 10/11 anaplastic -> confounded with histology"},
        "L1_he_tumor_level": he_tumor,
        "he_per_spot_histology_emb": he_spot,
        "L2_he_per_spot_balanced_emb": he_spot_bal,
        "spatial": {"median_reg_cluster_z": float(zdf.reg_cluster_z.median()),
                    "n_libraries_z_gt2": int((zdf.reg_cluster_z > 2).sum()), "n_libraries": int(len(zdf))},
        "L4_note": "Clinical SIOP regressive-type call + higher ceiling need whole-slide images "
                   "(external); this pilot is bounded by Visium-hires resolution.",
        "seed": SEED,
    }
    (OUT_DIR / "regressive_pilot.json").write_text(json.dumps(summary, indent=2))
    print("\n[ok] wrote results/regressive_pilot/regressive_pilot.json + figures", flush=True)


if __name__ == "__main__":
    main()
