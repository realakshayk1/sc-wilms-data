#!/usr/bin/env python3
"""Cross-modal concordance: snRNA compartment proportions vs Visium fractions."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr

from utils import ensure_dir, load_config, resolve_path, set_seed_logged, setup_logging

CELL_STATES = ["blastemal", "epithelial", "stromal"]


def snrna_sample_fractions(scores_rds: Path) -> pd.DataFrame:
    """Aggregate snRNA compartment proportions per sample via R."""
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        out_csv = Path(tmp.name)
    with tempfile.NamedTemporaryFile(suffix=".R", delete=False, mode="w", encoding="utf-8") as rf:
        r_path = str(scores_rds).replace("\\", "/")
        o_path = str(out_csv).replace("\\", "/")
        rf.write(f"""
dat <- readRDS("{r_path}")
meta <- dat$meta
rows <- list()
for (sid in unique(as.character(meta$sample_id))) {{
  sub <- meta[meta$sample_id == sid, , drop = FALSE]
  hist <- unique(sub$histology[!is.na(sub$histology)])
  row <- data.frame(
    sample_id = sid,
    histology = if (length(hist) == 1) hist[1] else NA_character_,
    stringsAsFactors = FALSE
  )
  for (st in c("blastemal", "epithelial", "stromal")) {{
    row[[paste0("snrna_frac_", st)]] <- mean(sub$cell_state == st, na.rm = TRUE)
  }}
  rows[[length(rows) + 1]] <- row
}}
write.csv(do.call(rbind, rows), "{o_path}", row.names = FALSE)
""")
        r_file = Path(rf.name)
    root = Path(__file__).resolve().parents[1]
    rscript = shutil.which("Rscript")
    if rscript:
        subprocess.run([rscript, str(r_file)], check=True)
    else:
        subprocess.run([str(root / "scripts" / "rscript.bat"), str(r_file)], check=True)
    r_file.unlink(missing_ok=True)
    df = pd.read_csv(out_csv)
    out_csv.unlink(missing_ok=True)
    return df


def main() -> None:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cfg = load_config()
    set_seed_logged(cfg["features"]["seed"], "snrna_spatial_concordance")
    out_csv = resolve_path(cfg, cfg["paths"]["phase_b"]["snrna_concordance_csv"])
    ensure_dir(out_csv.parent)

    if out_csv.exists() and not args.force:
        print(f"[skip] Concordance exists: {out_csv}")
        return

    scores_rds = resolve_path(cfg, cfg["paths"]["phase_a"]["scores_rds"])
    sig_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_signatures_parquet"])
    frac_path = resolve_path(cfg, cfg["paths"]["phase_b"]["spot_fractions_csv"])

    snrna_df = snrna_sample_fractions(scores_rds)
    sig = pd.read_parquet(sig_path)

    visium_rna = sig.groupby("sample_id")[[f"deconv_{s}" for s in CELL_STATES]].mean().reset_index()
    visium_rna.columns = ["sample_id"] + [f"visium_rna_frac_{s}" for s in CELL_STATES]

    merged = snrna_df.merge(visium_rna, on="sample_id", how="inner")

    if frac_path.exists():
        frac = pd.read_csv(frac_path)
        visium_he = frac.groupby("sample_id")[[f"frac_{s}" for s in CELL_STATES]].mean().reset_index()
        visium_he.columns = ["sample_id"] + [f"visium_he_frac_{s}" for s in CELL_STATES]
        merged = merged.merge(visium_he, on="sample_id", how="left")

    corr_rows = []
    for st in CELL_STATES:
        x = merged[f"snrna_frac_{st}"]
        y = merged[f"visium_rna_frac_{st}"]
        mask = x.notna() & y.notna()
        if mask.sum() < 3:
            continue
        r, p = pearsonr(x[mask], y[mask])
        rho, _ = spearmanr(x[mask], y[mask])
        row = {
            "compartment": st,
            "comparison": "snrna_vs_visium_rna",
            "pearson_r": float(r),
            "pearson_p": float(p),
            "spearman_rho": float(rho),
            "n_samples": int(mask.sum()),
        }
        if f"visium_he_frac_{st}" in merged.columns:
            xh = merged[f"visium_he_frac_{st}"]
            mask_h = x.notna() & xh.notna()
            if mask_h.sum() >= 3:
                rh, ph = pearsonr(x[mask_h], xh[mask_h])
                row["snrna_vs_visium_he_pearson_r"] = float(rh)
                row["snrna_vs_visium_he_pearson_p"] = float(ph)
        corr_rows.append(row)

    merged.to_csv(out_csv, index=False)
    pd.DataFrame(corr_rows).to_csv(
        out_csv.with_name("snrna_spatial_concordance_correlations.csv"), index=False
    )
    print(f"[ok] snRNA–spatial concordance -> {out_csv} ({len(merged)} samples)")


if __name__ == "__main__":
    main()
