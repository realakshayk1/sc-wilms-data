#!/usr/bin/env python3
"""Stage 4: assemble a per-tumor PhysiCell model directory.

Ties Stages 1-3 together into a self-contained model per tumor:
  results/abm/<sample_id>/
    cells.csv                 (Stage 2)  initial agent positions + types
    rules.csv                 (Stage 3)  cell-behavior grammar
    PhysiCell_settings.xml    (here)     domain + oxygen substrate + cell defs + ruleset
    provenance.json           (here)     seeds, sources, versions

The XML domain is sized from the tumor's own agent bounding box; oxygen uses standard
BioFVM constants (Dirichlet at physioxia); cell definitions carry the per-tumor base rates
from results/abm/positives_to_physicell.yaml. The settings are scaffolding: validate
against the target PhysiCell (>=1.14.1, grammar-enabled) sample project on the cluster
before the first run (plan milestone M4 / PRD AC5).

Usage: python 04_build_model.py [--sample SCPCS000168 ...]
"""
from __future__ import annotations

import argparse
import json
import sys
import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from abm_utils import COMPARTMENTS, ensure_dir, load_config, resolve_path, setup_logging  # noqa: E402

# BioFVM oxygen defaults (um^2/min diffusion, 1/min decay, mmHg boundary, uptake 1/min)
O2_DIFFUSION, O2_DECAY, O2_DIRICHLET, O2_UPTAKE = 1.0e5, 0.1, 38.0, 0.1


def _sub(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, {k: str(v) for k, v in attrib.items()})
    if text is not None:
        el.text = str(text)
    return el


def _secretion_block(ph, rates, substrates):
    """Per-cell secretion/uptake for every substrate. Oxygen is consumed (gradient former);
    IGF2 uptake is per-tumor (IGF program); ECM is secreted by stromal cells (v1.1)."""
    sec = _sub(ph, "secretion")
    def one(name, secretion_rate, uptake_rate, target=1.0):
        s = _sub(sec, "substrate", name=name)
        _sub(s, "secretion_rate", secretion_rate, units="1/min")
        _sub(s, "secretion_target", target, units="substrate density")
        _sub(s, "uptake_rate", uptake_rate, units="1/min")
        _sub(s, "net_export_rate", 0.0, units="total substrate/min")
    one("oxygen", 0.0, O2_UPTAKE)
    if "IGF2" in substrates:
        one("IGF2", 0.0, rates.get("igf_uptake_rate", 0.001))
    if "ECM" in substrates:
        one("ECM", rates.get("ecm_secretion_rate", 0.0), 0.0)


def build_xml(sample_id, tumor, dom, sim, substrates=None) -> ET.Element:
    substrates = substrates or {}
    root = ET.Element("PhysiCell_settings", version="devel-metadata")

    d = _sub(root, "domain")
    _sub(d, "x_min", 0); _sub(d, "x_max", dom["x_max"])
    _sub(d, "y_min", 0); _sub(d, "y_max", dom["y_max"])
    _sub(d, "z_min", -10); _sub(d, "z_max", 10)
    _sub(d, "dx", 20); _sub(d, "dy", 20); _sub(d, "dz", 20)
    _sub(d, "use_2D", "true")

    ov = _sub(root, "overall")
    _sub(ov, "max_time", sim["max_time_min"], units="min")
    _sub(ov, "time_units", "min"); _sub(ov, "space_units", "micron")

    save = _sub(root, "save")
    _sub(save, "folder", "output")
    fd = _sub(save, "full_data"); _sub(fd, "interval", sim["save_interval_min"], units="min")
    sd = _sub(save, "SVG"); _sub(sd, "interval", sim["save_interval_min"], units="min")

    me = _sub(root, "microenvironment_setup")
    o2 = _sub(me, "variable", name="oxygen", units="mmHg", ID="0")
    pd_ = _sub(o2, "physical_parameter_set")
    _sub(pd_, "diffusion_coefficient", O2_DIFFUSION, units="micron^2/min")
    _sub(pd_, "decay_rate", O2_DECAY, units="1/min")
    _sub(o2, "initial_condition", O2_DIRICHLET, units="mmHg")
    bc = _sub(o2, "Dirichlet_boundary_condition", O2_DIRICHLET, units="mmHg", enabled="true")  # noqa: F841

    # v1.1 substrates (IGF2, ECM) from config, IDs after oxygen
    for j, (name, spec) in enumerate(substrates.items(), start=1):
        v = _sub(me, "variable", name=name, units="dimensionless", ID=str(j))
        ps = _sub(v, "physical_parameter_set")
        _sub(ps, "diffusion_coefficient", spec.get("diffusion", 0.0), units="micron^2/min")
        _sub(ps, "decay_rate", spec.get("decay", 0.0), units="1/min")
        _sub(v, "initial_condition", spec.get("initial", 0.0), units="dimensionless")
        _sub(v, "Dirichlet_boundary_condition", spec.get("boundary", 0.0),
             units="dimensionless", enabled="true" if spec.get("dirichlet") else "false")

    defs = _sub(root, "cell_definitions")
    for i, ct in enumerate(COMPARTMENTS):
        rates = tumor["cell_types"][ct]
        cd = _sub(defs, "cell_definition", name=ct, ID=str(i))
        ph = _sub(cd, "phenotype")
        cyc = _sub(ph, "cycle", model="live", code="6")
        rates_el = _sub(cyc, "phase_transition_rates", units="1/min")
        _sub(rates_el, "rate", rates["proliferation_rate"], start_index="0", end_index="0",
             fixed_duration="false")
        death = _sub(ph, "death")
        apo = _sub(death, "model", name="apoptosis", code="100")
        _sub(apo, "death_rate", rates["apoptosis_rate"], units="1/min")
        nec = _sub(death, "model", name="necrosis", code="101")
        _sub(nec, "death_rate", 0.0, units="1/min")
        mech = _sub(ph, "mechanics")
        _sub(mech, "cell_cell_adhesion_strength", rates.get("adhesion_strength", 0.4),
             units="micron/min")
        _sub(mech, "cell_cell_repulsion_strength", 10.0, units="micron/min")
        # motility: per-tumor migration speed (EMT-scaled in 17_positives_to_abm.py)
        mot = _sub(ph, "motility")
        _sub(mot, "speed", rates.get("migration_speed", 0.3), units="micron/min")
        _sub(mot, "persistence_time", 1.0, units="min")
        _sub(mot, "migration_bias", 0.5, units="dimensionless")
        mot_opt = _sub(mot, "options")
        _sub(mot_opt, "enabled", "true")
        _sub(mot_opt, "use_2D", "true")
        _secretion_block(ph, rates, substrates)

    rules = _sub(root, "cell_rules")
    rs = _sub(rules, "rulesets")
    r = _sub(rs, "ruleset", protocol="CBHG", version="3.0", format="csv", enabled="true")
    _sub(r, "folder", "."); _sub(r, "filename", "rules.csv")

    ic = _sub(root, "initial_conditions")
    cp = _sub(ic, "cell_positions", type="csv", enabled="true")
    _sub(cp, "folder", "."); _sub(cp, "filename", "cells.csv")

    opt = _sub(root, "options")
    _sub(opt, "virtual_wall_at_domain_edge", "true")

    ur = _sub(root, "user_parameters")
    _sub(ur, "sample_id", sample_id)
    _sub(ur, "high_grade_regime", str(tumor.get("high_grade_regime", False)).lower())
    return root


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=None)
    args = ap.parse_args()

    setup_logging()
    cfg = load_config()
    sim = cfg["phase_c"]["simulation"]
    substrates = cfg["phase_c"].get("substrates", {})
    margin = float(cfg["phase_c"]["domain"]["margin_um"])
    abm = yaml.safe_load(resolve_path(cfg, "results/abm/positives_to_physicell.yaml").read_text())
    out_dir = resolve_path(cfg, "results/abm")

    samples = args.sample or list(abm.get("tumors", {}))
    built = []
    for sid in samples:
        tumor = abm["tumors"].get(sid)
        d = out_dir / sid
        cells_csv = d / "cells.csv"
        if tumor is None or not cells_csv.exists() or not (d / "rules.csv").exists():
            print(f"[skip] {sid}: missing tumor params / cells.csv / rules.csv")
            continue
        cells = pd.read_csv(cells_csv)
        dom = {"x_max": round(cells["x"].max() + margin, 1),
               "y_max": round(cells["y"].max() + margin, 1)}
        root = build_xml(sid, tumor, dom, sim, substrates)
        xml = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
        (d / "PhysiCell_settings.xml").write_text(xml)
        prov = {"sample_id": sid, "seed": cfg["phase_c"]["seed"],
                "n_agents": int(len(cells)), "domain_um": dom,
                "high_grade_regime": bool(tumor.get("high_grade_regime", False)),
                "deconvolution_backend": cfg["phase_c"]["deconvolution"]["backend"],
                "density_source": cfg["phase_c"]["density"]["source"],
                "sources": {
                    "params": "results/abm/positives_to_physicell.yaml",
                    "coords": "Visium tissue_positions_list.csv",
                    "density": "nucleus_features_stardist_80_pt40.parquet"},
                "physicell_target": ">=1.14.1 (grammar-enabled)"}
        (d / "provenance.json").write_text(json.dumps(prov, indent=2))
        built.append(sid)
        print(f"[ok] {sid}: model dir -> {d}")

    if built:
        (out_dir / "model_manifest.json").write_text(
            json.dumps({"tumors": built, "n": len(built)}, indent=2))
        print(f"[ok] {len(built)} model dirs; manifest -> {out_dir/'model_manifest.json'}")


if __name__ == "__main__":
    main()
