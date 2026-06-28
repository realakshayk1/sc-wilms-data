"""Integration smoke test for Phase B demo pipeline."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
PHASE2 = ROOT / "phase2_histology_ml"
ENV = {**dict(__import__("os").environ), "WILMS_DEMO": "1"}


@pytest.fixture(scope="module")
def demo_pipeline_outputs(tmp_path_factory):
    """Run demo Phase B scripts into a temp processed tree."""
    # Use repo paths but allow re-run via demo flag
    scripts = [
        "01_extract_tiles.py",
        "02_segment_nuclei.py",
        "03_nucleus_features.py",
        "04_train_classifier.py",
        "05_spot_fractions.py",
        "06_map_to_physicell.py",
    ]
    for script in scripts:
        subprocess.run(
            [PY, str(PHASE2 / script), "--demo"],
            cwd=ROOT,
            env=ENV,
            check=True,
        )
    return ROOT


def test_demo_tiles_created(demo_pipeline_outputs):
    tiles_dir = demo_pipeline_outputs / "data" / "processed" / "he_tiles"
    manifest = tiles_dir / "tiles_manifest.json"
    assert manifest.exists()
    with open(manifest) as f:
        meta = json.load(f)
    assert len(meta) >= 1


def test_classifier_metrics_exist(demo_pipeline_outputs):
    metrics = demo_pipeline_outputs / "results" / "classifier" / "classifier_metrics.json"
    assert metrics.exists()
    with open(metrics) as f:
        m = json.load(f)
    assert "balanced_accuracy" in m
    assert "per_class_balanced_accuracy" in m


def test_abm_stub_complete(demo_pipeline_outputs):
    summary = demo_pipeline_outputs / "results" / "abm" / "run_001" / "simulation_summary.json"
    assert summary.exists()
