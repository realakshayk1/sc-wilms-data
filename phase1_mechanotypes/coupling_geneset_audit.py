#!/usr/bin/env python3
"""WS1 step 0 (audit finding A): prove the curated lever/axis gene panels are disjoint.

Reads config/levers.yaml and computes the pairwise shared-gene count + Jaccard for every program
pair. The DISJOINTNESS ASSERTION covers only COMPARED pairs — {levers | extrinsic_axes} — because
those are the ones that become network edges or transfer predictor<->response pairs. `typing_only`
programs and `definitional_maps` are excluded (they are allowed to overlap levers by design) but are
still reported for transparency.

Each program's gene set is (genes_positive | genes_negative): a gene shared even with opposite sign
induces dependence, so both poles count.

Output: results/couplings/geneset_overlap.csv  (every pair; asserted flag; shared genes)
Exit non-zero if any asserted pair overlaps — this is the gate WS1 depends on.

Usage:  python phase1_mechanotypes/coupling_geneset_audit.py
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
LEVERS_YAML = ROOT / "config" / "levers.yaml"
OUT = ROOT / "results" / "couplings" / "geneset_overlap.csv"


def gene_union(prog: dict) -> set[str]:
    return set(prog.get("genes_positive", []) or []) | set(prog.get("genes_negative", []) or [])


def main() -> int:
    cfg = yaml.safe_load(LEVERS_YAML.read_text())

    # collect programs with their role; definitional_maps have no genes -> skip
    programs: dict[str, dict] = {}
    roles: dict[str, str] = {}
    for role in ("levers", "extrinsic_axes", "typing_only"):
        for p in cfg.get(role, []) or []:
            programs[p["id"]] = gene_union(p)
            roles[p["id"]] = role

    # "compared" = the pairs the disjointness rule applies to: levers | extrinsic_axes
    compared = {pid for pid, r in roles.items() if r in ("levers", "extrinsic_axes")}

    rows = []
    n_violation = 0
    for a, b in combinations(sorted(programs), 2):
        ga, gb = programs[a], programs[b]
        shared = sorted(ga & gb)
        jacc = len(ga & gb) / len(ga | gb) if (ga | gb) else 0.0
        asserted = a in compared and b in compared
        violated = asserted and len(shared) > 0
        n_violation += int(violated)
        rows.append({
            "program_a": a, "role_a": roles[a],
            "program_b": b, "role_b": roles[b],
            "n_shared": len(shared), "jaccard": round(jacc, 4),
            "asserted_disjoint": asserted,
            "VIOLATION": violated,
            "shared_genes": ";".join(shared),
        })

    df = pd.DataFrame(rows).sort_values(
        ["VIOLATION", "asserted_disjoint", "n_shared"], ascending=[False, False, False]
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    n_lever = sum(r == "levers" for r in roles.values())
    n_axis = sum(r == "extrinsic_axes" for r in roles.values())
    print(f"[levers.yaml] {n_lever} levers, {n_axis} extrinsic axes, "
          f"{sum(r == 'typing_only' for r in roles.values())} typing-only programs")
    print(f"[audit] {len(compared)} compared programs -> "
          f"{len(list(combinations(sorted(compared), 2)))} pairs asserted disjoint")
    print(f"[out] {OUT.relative_to(ROOT)}")

    if n_violation:
        print(f"[FAIL] {n_violation} asserted pair(s) share genes:")
        for r in rows:
            if r["VIOLATION"]:
                print(f"   {r['program_a']}  &  {r['program_b']}: {r['shared_genes']}")
        return 1

    # show the informational (allowed) overlaps so the biology is visible
    info = df[(~df.asserted_disjoint) & (df.n_shared > 0)]
    if len(info):
        print(f"[ok] PASS - all asserted pairs disjoint. "
              f"{len(info)} allowed typing<->lever overlap(s) (informational):")
        for _, r in info.iterrows():
            print(f"   {r['program_a']}  &  {r['program_b']}: {r['shared_genes']}")
    else:
        print("[ok] PASS - all asserted pairs disjoint.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
