#!/usr/bin/env python3
"""WS3: pull TARGET-WT OPEN-ACCESS data from GDC (no dbGaP) for the survival lever-validation.

Downloads, via the GDC REST API:
  - clinical (OS: vital_status, days_to_death, days_to_last_follow_up, age, stage) -> clinical.tsv
  - bulk RNA-seq STAR-counts (125 cases)  -> expr_counts.tsv   (genes x case, `unstranded`)
  - open Masked Somatic Mutation MAFs (38 cases) -> mutations.tsv (case, mutated non-silent genes)
All open-access; the raw dbGaP sequencing is NOT touched. Output: data/raw/target_wt/.
Then 22_target_wt_survival.R builds surv_df (OS) + binary levers and runs Cox.

Usage: python phase1_mechanotypes/target_wt_pull.py
"""
from __future__ import annotations

import io
import json
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "raw" / "target_wt"
API = "https://api.gdc.cancer.gov"


def _get(path, params):
    return json.loads(urllib.request.urlopen(f"{API}/{path}?" + urllib.parse.urlencode(params),
                                             timeout=120).read())


def files_for(data_type, workflow=None):
    conts = [{"op": "in", "content": {"field": "cases.project.project_id", "value": ["TARGET-WT"]}},
             {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
             {"op": "in", "content": {"field": "files.data_type", "value": [data_type]}}]
    if workflow:
        conts.append({"op": "in", "content": {"field": "files.analysis.workflow_type", "value": [workflow]}})
    d = _get("files", {"filters": json.dumps({"op": "and", "content": conts}),
                       "fields": "file_id,file_name,cases.submitter_id", "format": "json", "size": "3000"})
    return d["data"]["hits"]


def download_bulk(file_ids, chunk=40):
    """POST file_ids to /data in chunks; yield (member_name, bytes) for every file in the tarballs."""
    for i in range(0, len(file_ids), chunk):
        ids = file_ids[i:i + chunk]
        req = urllib.request.Request(f"{API}/data", method="POST",
                                     data=json.dumps({"ids": ids}).encode(),
                                     headers={"Content-Type": "application/json"})
        raw = urllib.request.urlopen(req, timeout=600).read()
        try:
            tf = tarfile.open(fileobj=io.BytesIO(raw), mode="r:*")
        except tarfile.ReadError:
            # a single-file request returns the file itself, not a tar
            yield ids[0], raw
            continue
        for m in tf.getmembers():
            if m.isfile():
                yield m.name, tf.extractfile(m).read()
        print(f"  downloaded {min(i + chunk, len(file_ids))}/{len(file_ids)} files", flush=True)


def pull_clinical():
    filt = {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TARGET-WT"]}}
    d = _get("cases", {"filters": json.dumps(filt), "format": "json", "size": "2000",
                       "fields": ("submitter_id,demographic.vital_status,demographic.days_to_death,"
                                  "diagnoses.days_to_last_follow_up,diagnoses.age_at_diagnosis,"
                                  "diagnoses.ajcc_pathologic_stage,diagnoses.tumor_stage")})
    rows = []
    for h in d["data"]["hits"]:
        dem = h.get("demographic") or {}
        dg = (h.get("diagnoses") or [{}])[0]
        rows.append({"case": h["submitter_id"], "vital_status": dem.get("vital_status"),
                     "days_to_death": dem.get("days_to_death"),
                     "days_to_last_follow_up": dg.get("days_to_last_follow_up"),
                     "age_days": dg.get("age_at_diagnosis"),
                     "stage": dg.get("ajcc_pathologic_stage") or dg.get("tumor_stage")})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "clinical.tsv", sep="\t", index=False)
    print(f"[clinical] {len(df)} cases -> clinical.tsv ({(df.vital_status=='Dead').sum()} deaths)", flush=True)


def pull_expression():
    hits = files_for("Gene Expression Quantification", "STAR - Counts")
    fid2case = {h["file_id"]: (h.get("cases") or [{}])[0].get("submitter_id") for h in hits}
    name2case = {}  # match tar member basename -> case via file_name
    fname2case = {h["file_name"]: (h.get("cases") or [{}])[0].get("submitter_id") for h in hits}
    cols = {}
    print(f"[expr] downloading {len(hits)} STAR-counts files ...", flush=True)
    for member, data in download_bulk([h["file_id"] for h in hits]):
        base = member.split("/")[-1]
        case = fname2case.get(base) or fid2case.get(member.split("/")[0])
        if case is None:
            continue
        t = pd.read_csv(io.BytesIO(data), sep="\t", comment="#")
        t = t[~t["gene_id"].astype(str).str.startswith("N_")]
        s = pd.Series(t["unstranded"].to_numpy(), index=t["gene_id"].astype(str).str.replace(r"\..*$", "", regex=True))
        cols[case] = cols.get(case, s).add(s, fill_value=0) if case in cols else s
    expr = pd.DataFrame(cols)
    expr.index.name = "gene_id"
    expr.to_csv(OUT / "expr_counts.tsv", sep="\t")
    print(f"[expr] {expr.shape[0]} genes x {expr.shape[1]} cases -> expr_counts.tsv", flush=True)


def pull_mutations():
    hits = files_for("Masked Somatic Mutation")
    print(f"[maf] downloading {len(hits)} MAF files ...", flush=True)
    rows = []
    SILENT = {"Silent", "Intron", "3'UTR", "5'UTR", "RNA", "IGR", "5'Flank", "3'Flank"}
    for member, data in download_bulk([h["file_id"] for h in hits]):
        try:
            m = pd.read_csv(io.BytesIO(data), sep="\t", comment="#", compression="gzip", low_memory=False)
        except Exception:
            try:
                m = pd.read_csv(io.BytesIO(data), sep="\t", comment="#", low_memory=False)
            except Exception:
                continue
        if "Hugo_Symbol" not in m or "Tumor_Sample_Barcode" not in m:
            continue
        m = m[~m["Variant_Classification"].isin(SILENT)]
        for bc, g in m.groupby("Tumor_Sample_Barcode"):
            case = "-".join(str(bc).split("-")[:3])   # TARGET-50-XXXX
            for gene in g["Hugo_Symbol"].dropna().unique():
                rows.append({"case": case, "gene": gene})
    mut = pd.DataFrame(rows).drop_duplicates()
    mut.to_csv(OUT / "mutations.tsv", sep="\t", index=False)
    print(f"[maf] {mut['case'].nunique() if len(mut) else 0} cases, {len(mut)} case-gene mutations -> mutations.tsv", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pull_clinical()
    pull_expression()
    pull_mutations()
    print("[ok] TARGET-WT open-access pull complete -> data/raw/target_wt/", flush=True)


if __name__ == "__main__":
    main()
