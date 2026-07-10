#!/usr/bin/env bash
# Stage 6: run the per-tumor PhysiCell cohort on a SLURM cluster.
#
# This is the ONE compute-gated stage. It assumes a grammar-enabled PhysiCell (>=1.14.1)
# has been built and its binary is on $PHYSICELL_BIN. Each array task takes one tumor and
# one replicate (distinct random seed), runs it in an isolated copy of that tumor's model
# dir, and writes a per-tumor QoI CSV that Stage 7 (07_validate.py) consumes.
#
# Usage (from repo root):
#   PHYSICELL_BIN=/path/to/PhysiCell/project \
#   sbatch --array=0-$(( $(wc -l < results/abm/model_manifest.txt) * REPLICATES - 1 )) \
#          phase3_abm/06_run_cohort.sh
#
# Prereqs produced on CPU by Stages 1-5:
#   results/abm/<sample_id>/{cells.csv,rules.csv,PhysiCell_settings.xml}
#   results/abm/<sample_id>/uq/sweep_manifest.csv        (optional UQ)
#   results/abm/model_manifest.txt                        (one sample_id per line)
#
# Patch-scale models (~1e3 agents, ~2e3 voxels) are small: 4 cores / 2 h / 8 GB is ample.
# NOTE: array size = n_runs * REPLICATES may exceed the scheduler MaxArraySize (often 1001);
# chunk with several `sbatch --array=A-B` submissions if so.
#SBATCH --job-name=wilms_abm
#SBATCH --account=mcb200052p
#SBATCH --partition=RM-shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=results/abm/logs/%A_%a.out
set -euo pipefail

: "${PHYSICELL_BIN:?set PHYSICELL_BIN to the compiled PhysiCell project binary}"
# compute nodes are fresh shells: load the gcc whose libstdc++ the binary needs (Bridges-2)
module load gcc/13.3.1-p20240614 2>/dev/null || module load gcc 2>/dev/null || true
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"     # PhysiCell is OpenMP; match the cores
REPLICATES="${REPLICATES:-10}"
ROOT="$(git rev-parse --show-toplevel)"
MANIFEST="$ROOT/results/abm/model_manifest.txt"

mapfile -t SAMPLES < "$MANIFEST"
N="${#SAMPLES[@]}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
SAMPLE="${SAMPLES[$(( IDX / REPLICATES ))]}"
REP="$(( IDX % REPLICATES ))"

SRC="$ROOT/results/abm/$SAMPLE"
WORK="$ROOT/results/abm/$SAMPLE/replicates/rep_${REP}"
mkdir -p "$WORK/output"
cp "$SRC/cells.csv" "$SRC/rules.csv" "$SRC/PhysiCell_settings.xml" "$WORK/"

# distinct RNG seed per replicate (reproducible from IDX)
SEED="$(( 1000 + IDX ))"
sed -i "s#<random_seed>.*</random_seed>#<random_seed>${SEED}</random_seed>#" \
    "$WORK/PhysiCell_settings.xml" || true

echo "[run] sample=$SAMPLE rep=$REP seed=$SEED"
cd "$WORK"
"$PHYSICELL_BIN" PhysiCell_settings.xml

# --- extract per-replicate QoI from PhysiCell output (final-frame tumor burden etc.) ---
python "$ROOT/phase3_abm/qoi_extract.py" --run-dir "$WORK" --out "$WORK/qoi.csv"

# aggregate replicates -> results/abm/<sample_id>/output/qoi.csv is done post-array by:
#   python phase3_abm/qoi_extract.py --aggregate results/abm/$SAMPLE
