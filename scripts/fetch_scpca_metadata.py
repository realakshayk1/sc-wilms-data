#!/usr/bin/env python3
"""Fetch SCPCP000006 metadata via public ScPCA API (no auth required)."""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

API = "https://api.scpca.alexslemonade.org/v1"
PROJECT_ID = "SCPCP000006"


def get_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "sc-wilms-data/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def flatten_sample(s: dict) -> dict:
    row = {k: v for k, v in s.items() if k not in ("computed_files", "multiplexed_with", "project")}
    if isinstance(row.get("additional_metadata"), dict):
        for k, v in row.pop("additional_metadata").items():
            row.setdefault(k, v)
    # List fields → pipe-separated for CSV
    for k in ("seq_units", "technologies", "modalities"):
        if isinstance(row.get(k), list):
            row[k] = "|".join(row[k])
    return row


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "data" / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Fetching {API}/projects/{PROJECT_ID} ...")
    proj = get_json(f"{API}/projects/{PROJECT_ID}")

    print(f"Title: {proj.get('title')}")
    print(f"Modalities: {proj.get('modalities')}")
    print(f"Has spatial: {proj.get('has_spatial_data')}")
    print(f"Downloadable samples: {proj.get('downloadable_sample_count')}")

    dc = proj.get("diagnoses_counts") or {}
    if dc:
        print("\n--- diagnosis counts ---")
        for k, v in sorted(dc.items(), key=lambda x: -x[1]):
            print(f"  {k}: {v}")

    import pandas as pd

    rows = [flatten_sample(s) for s in proj.get("samples") or []]
    df = pd.DataFrame(rows)
    all_csv = out_dir / f"{PROJECT_ID}_samples.csv"
    df.to_csv(all_csv, index=False)
    print(f"\n[ok] {len(df)} samples -> {all_csv}")

    if "subdiagnosis" in df.columns:
        print("\n--- histology (subdiagnosis) ---")
        print(df["subdiagnosis"].value_counts().to_string())

    if "seq_units" in df.columns:
        print("\n--- seq_units ---")
        print(df["seq_units"].value_counts().to_string())
        nuc = df[df["seq_units"].str.contains("nucleus", case=False, na=False)]
        vis = df[df["seq_units"].str.contains("spot", case=False, na=False)]
        if len(nuc):
            p = out_dir / f"{PROJECT_ID}_nucleus_samples.csv"
            nuc.to_csv(p, index=False)
            print(f"[ok] {len(nuc)} samples with nucleus -> {p}")
        if len(vis):
            p = out_dir / f"{PROJECT_ID}_visium_samples.csv"
            vis.to_csv(p, index=False)
            print(f"[ok] {len(vis)} samples with spot/Visium -> {p}")

    summary = {
        k: proj[k]
        for k in (
            "scpca_id", "title", "modalities", "has_spatial_data",
            "has_bulk_rna_seq", "downloadable_sample_count", "diagnoses",
        )
        if k in proj
    }
    meta_json = out_dir / f"{PROJECT_ID}_project_summary.json"
    meta_json.write_text(json.dumps(summary, indent=2))
    print(f"[ok] summary -> {meta_json}")


if __name__ == "__main__":
    main()
