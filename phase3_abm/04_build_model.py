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


def _secretion_block(ph, rates, substrates, o2_uptake=O2_UPTAKE):
    """Per-cell secretion/uptake for every substrate. Oxygen is consumed (gradient former);
    IGF2 uptake is per-tumor (IGF program); ECM is secreted by stromal cells (v1.1). Necrotic
    tissue passes zero rates (dead cells neither consume nor secrete)."""
    sec = _sub(ph, "secretion")
    def one(name, secretion_rate, uptake_rate, target=1.0):
        s = _sub(sec, "substrate", name=name)
        _sub(s, "secretion_rate", secretion_rate, units="1/min")
        _sub(s, "secretion_target", target, units="substrate density")
        _sub(s, "uptake_rate", uptake_rate, units="1/min")
        _sub(s, "net_export_rate", 0.0, units="total substrate/min")
    one("oxygen", 0.0, o2_uptake)
    if "IGF2" in substrates:
        one("IGF2", 0.0, rates.get("igf_uptake_rate", 0.001))
    if "ECM" in substrates:
        one("ECM", rates.get("ecm_secretion_rate", 0.0), 0.0)


def _write_cell_def(defs, name, idx, rates, substrates, o2_uptake=O2_UPTAKE, motile=True):
    """One complete <cell_definition> — mirrors the PhysiCell template phenotype schema so the
    type registers and cells.csv can reference it by name. Cycle is the Live model (code 5):
    a single phase whose 0->0 transition rate is the (omics-determined) proliferation rate."""
    cd = _sub(defs, "cell_definition", name=name, ID=str(idx))
    ph = _sub(cd, "phenotype")

    cyc = _sub(ph, "cycle", code="5", name="Live")
    _sub(_sub(cyc, "phase_transition_rates", units="1/min"), "rate",
         rates.get("proliferation_rate", 0.0), start_index="0", end_index="0", fixed_duration="false")

    death = _sub(ph, "death")
    apo = _sub(death, "model", code="100", name="apoptosis")
    _sub(apo, "death_rate", rates.get("apoptosis_rate", 0.0), units="1/min")
    _sub(_sub(apo, "phase_durations", units="min"), "duration", 516, index="0", fixed_duration="true")
    ap = _sub(apo, "parameters")
    for tag, val in [("unlysed_fluid_change_rate", 0.05), ("lysed_fluid_change_rate", 0),
                     ("cytoplasmic_biomass_change_rate", 1.66667e-02),
                     ("nuclear_biomass_change_rate", 5.83333e-03), ("calcification_rate", 0)]:
        _sub(ap, tag, val, units="1/min")
    _sub(ap, "relative_rupture_volume", 2.0, units="dimensionless")
    nec = _sub(death, "model", code="101", name="necrosis")
    _sub(nec, "death_rate", 0.0, units="1/min")
    npd = _sub(nec, "phase_durations", units="min")
    _sub(npd, "duration", 0, index="0", fixed_duration="true")
    _sub(npd, "duration", 86400, index="1", fixed_duration="true")
    npar = _sub(nec, "parameters")
    for tag, val in [("unlysed_fluid_change_rate", 1.11667e-2), ("lysed_fluid_change_rate", 8.33333e-4),
                     ("cytoplasmic_biomass_change_rate", 5.33333e-5),
                     ("nuclear_biomass_change_rate", 2.16667e-3), ("calcification_rate", 0)]:
        _sub(npar, tag, val, units="1/min")
    _sub(npar, "relative_rupture_volume", 2.0, units="dimensionless")

    vol = _sub(ph, "volume")
    for tag, val, u in [("total", 2494, "micron^3"), ("fluid_fraction", 0.75, "dimensionless"),
                        ("nuclear", 540, "micron^3"), ("fluid_change_rate", 0.05, "1/min"),
                        ("cytoplasmic_biomass_change_rate", 0.0045, "1/min"),
                        ("nuclear_biomass_change_rate", 0.0055, "1/min"),
                        ("calcified_fraction", 0, "dimensionless"), ("calcification_rate", 0, "1/min"),
                        ("relative_rupture_volume", 2.0, "dimensionless")]:
        _sub(vol, tag, val, units=u)

    mech = _sub(ph, "mechanics")
    _sub(mech, "cell_cell_adhesion_strength", rates.get("adhesion_strength", 0.4), units="micron/min")
    _sub(mech, "cell_cell_repulsion_strength", 10.0, units="micron/min")
    _sub(mech, "relative_maximum_adhesion_distance", 1.25, units="dimensionless")
    _sub(_sub(mech, "cell_adhesion_affinities"), "cell_adhesion_affinity", 1, name="default")
    mopt = _sub(mech, "options")
    _sub(mopt, "set_relative_equilibrium_distance", 1.8, enabled="false", units="dimensionless")
    _sub(mopt, "set_absolute_equilibrium_distance", 15.12, enabled="false", units="micron")
    _sub(mech, "attachment_elastic_constant", 0.01, units="1/min")
    _sub(mech, "attachment_rate", 0.0, units="1/min")
    _sub(mech, "detachment_rate", 0.0, units="1/min")
    _sub(mech, "maximum_number_of_attachments", 12)

    mot = _sub(ph, "motility")
    _sub(mot, "speed", rates.get("migration_speed", 0.3) if motile else 0.0, units="micron/min")
    _sub(mot, "persistence_time", 1.0, units="min")
    _sub(mot, "migration_bias", 0.5 if motile else 0.0, units="dimensionless")
    mo = _sub(mot, "options")
    _sub(mo, "enabled", "true" if motile else "false"); _sub(mo, "use_2D", "true")
    chem = _sub(mo, "chemotaxis")
    _sub(chem, "enabled", "false"); _sub(chem, "substrate", "oxygen"); _sub(chem, "direction", 1)
    ach = _sub(mo, "advanced_chemotaxis")
    _sub(ach, "enabled", "false"); _sub(ach, "normalize_each_gradient", "false")
    _sub(_sub(ach, "chemotactic_sensitivities"), "chemotactic_sensitivity", 0.0, substrate="oxygen")

    _secretion_block(ph, rates, substrates, o2_uptake=o2_uptake)

    ci = _sub(ph, "cell_interactions")
    _sub(ci, "apoptotic_phagocytosis_rate", 0, units="1/min")
    _sub(ci, "necrotic_phagocytosis_rate", 0, units="1/min")
    _sub(ci, "other_dead_phagocytosis_rate", 0, units="1/min")
    _sub(_sub(ci, "live_phagocytosis_rates"), "phagocytosis_rate", 0, name="default", units="1/min")
    _sub(_sub(ci, "attack_rates"), "attack_rate", 0, name="default", units="1/min")
    _sub(ci, "attack_damage_rate", 1, units="1/min")
    _sub(ci, "attack_duration", 0.1, units="min")
    _sub(_sub(ci, "fusion_rates"), "fusion_rate", 0, name="default", units="1/min")
    ct = _sub(ph, "cell_transformations")
    _sub(_sub(ct, "transformation_rates"), "transformation_rate", 0, name="default", units="1/min")
    cint = _sub(ph, "cell_integrity")
    _sub(cint, "damage_rate", 0.0, units="1/min"); _sub(cint, "damage_repair_rate", 0.0, units="1/min")

    _sub(_sub(cd, "custom_data"), "sample", 1.0, conserved="false", units="dimensionless")
    return cd


def build_xml(sample_id, tumor, dom, sim, substrates=None, include_necrotic=False) -> ET.Element:
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
    _sub(ov, "dt_diffusion", 0.01, units="min")       # PhysiCell RNG/solver setup needs these
    _sub(ov, "dt_mechanics", 0.1, units="min")
    _sub(ov, "dt_phenotype", 6, units="min")

    par = _sub(root, "parallel")                      # required: setup_rng() sizes per-thread RNG
    _sub(par, "omp_num_threads", 4)

    save = _sub(root, "save")
    _sub(save, "folder", "output")
    fd = _sub(save, "full_data")
    _sub(fd, "interval", sim["save_interval_min"], units="min"); _sub(fd, "enable", "true")
    sd = _sub(save, "SVG")
    _sub(sd, "interval", sim["save_interval_min"], units="min"); _sub(sd, "enable", "true")

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
    # PhysiCell requires a 'default' cell_definition as ID 0 (initialize_cell_definitions_from_
    # _pugixml sets cell_defaults from it); our compartments follow and are referenced by NAME
    # from cells.csv. Without 'default' the named types don't register and every agent is skipped.
    _write_cell_def(defs, "default", 0, {"proliferation_rate": 0.0, "apoptosis_rate": 0.0}, substrates)
    idx = 1
    for ct in COMPARTMENTS:
        _write_cell_def(defs, ct, idx, tumor["cell_types"][ct], substrates)
        idx += 1
    if include_necrotic:
        # inert necrotic tissue (seeded at high-mito spots): no cycle/death/motility/uptake.
        _write_cell_def(defs, "necrotic", idx,
                        {"adhesion_strength": 0.1, "igf_uptake_rate": 0.0, "ecm_secretion_rate": 0.0},
                        substrates, o2_uptake=0.0, motile=False)

    rules = _sub(root, "cell_rules")
    rs = _sub(rules, "rulesets")
    r = _sub(rs, "ruleset", protocol="CBHG", version="3.0", format="csv", enabled="true")
    _sub(r, "folder", "."); _sub(r, "filename", "rules.csv")

    ic = _sub(root, "initial_conditions")
    cp = _sub(ic, "cell_positions", type="csv", enabled="true")
    _sub(cp, "folder", "."); _sub(cp, "filename", "cells.csv")

    opt = _sub(root, "options")
    _sub(opt, "virtual_wall_at_domain_edge", "true")
    _sub(opt, "random_seed", 0)                        # inside <options>; 06_run_cohort seds per replicate

    ur = _sub(root, "user_parameters")
    # PhysiCell sample projects' setup_tissue() reads number_of_cells for demo placement;
    # 0 => place nothing procedurally and load the real agents from cells.csv instead.
    _sub(ur, "number_of_cells", 0, type="int", units="none",
         description="demo placement count; 0 = load initial cells from cells.csv only")
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
    include_necrotic = bool(cfg["phase_c"].get("necrotic", {}).get("enabled"))
    margin = float(cfg["phase_c"]["domain"]["margin_um"])
    abm = yaml.safe_load(resolve_path(cfg, "results/abm/positives_to_physicell.yaml").read_text())
    out_dir = resolve_path(cfg, "results/abm")

    # model units: (run_id, sample_id). Patch mode -> one per (tumor, patch) from the
    # patch manifest; whole-slide -> one per tumor. Rules/params are per tumor (shared by
    # a tumor's patches); each run dir gets its own copy of rules.csv for PhysiCell.
    pm_path = out_dir / "patch_manifest.csv"
    if pm_path.exists():
        pm = pd.read_csv(pm_path)
        if args.sample:
            pm = pm[pm["sample_id"].isin(args.sample)]
        units = list(zip(pm["run_id"], pm["sample_id"]))
    else:
        samples = args.sample or list(abm.get("tumors", {}))
        units = [(s, s) for s in samples]

    built = []
    for run_id, sid in units:
        tumor = abm["tumors"].get(sid)
        run_d = ensure_dir(out_dir / run_id)
        cells_csv = run_d / "cells.csv"
        rules_src = out_dir / sid / "rules.csv"          # written per tumor by Stage 3
        if tumor is None or not cells_csv.exists() or not rules_src.exists():
            print(f"[skip] {run_id}: missing tumor params / cells.csv / rules.csv")
            continue
        if run_d != rules_src.parent:                    # patch dir needs its own rules.csv
            (run_d / "rules.csv").write_bytes(rules_src.read_bytes())   # preserve \n endings
        cells = pd.read_csv(cells_csv)
        dom = {"x_max": round(cells["x"].max() + margin, 1),
               "y_max": round(cells["y"].max() + margin, 1)}
        root = build_xml(run_id, tumor, dom, sim, substrates, include_necrotic)
        xml = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
        (run_d / "PhysiCell_settings.xml").write_bytes(
            xml.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8"))    # \n for Linux
        prov = {"run_id": run_id, "sample_id": sid, "seed": cfg["phase_c"]["seed"],
                "n_agents": int(len(cells)), "domain_um": dom,
                "high_grade_regime": bool(tumor.get("high_grade_regime", False)),
                "deconvolution_backend": cfg["phase_c"]["deconvolution"]["backend"],
                "density_source": cfg["phase_c"]["density"]["source"],
                "patch_mode": pm_path.exists(),
                "sources": {
                    "params": "results/abm/positives_to_physicell.yaml",
                    "coords": "Visium tissue_positions_list.csv",
                    "density": "nucleus_features_stardist_80_pt40.parquet"},
                "physicell_target": ">=1.14.1 (grammar-enabled)"}
        (run_d / "provenance.json").write_text(json.dumps(prov, indent=2))
        built.append(run_id)
        print(f"[ok] {run_id}: model dir -> {run_d}")

    if built:
        (out_dir / "model_manifest.json").write_text(
            json.dumps({"runs": built, "n": len(built)}, indent=2))
        (out_dir / "model_manifest.txt").write_text("\n".join(built) + "\n")   # 06_run_cohort
        print(f"[ok] {len(built)} model dirs; manifest -> model_manifest.{{json,txt}}")


if __name__ == "__main__":
    main()
