#!/usr/bin/env bash
# Stage 4.5: CALIBRATION / build-check — run BEFORE the cohort array (06_run_cohort.sh).
#
# Times one median-sized and one largest patch (single replicate each) on the compiled
# PhysiCell binary, then sanity-checks the result and extrapolates cohort cost. This is the
# go/no-go gate: it confirms (a) the binary consumes our generated cells.csv + rules.csv +
# PhysiCell_settings.xml, (b) dynamics are sane (no agent explosion / collapse), and (c) the
# real per-run wallclock, so the ~1,110-run cohort is launched on measured numbers, not an
# estimate. Costs ~2 short runs of SU.
#
# Usage (from repo root, on Bridges-2):
#   PHYSICELL_BIN=/path/to/PhysiCell/project bash phase3_abm/calibrate.sh
#
# Reads results/abm/patch_manifest.csv + model dirs (Stages 2-4). Writes
# results/abm/calibration_report.txt.
set -euo pipefail

: "${PHYSICELL_BIN:?set PHYSICELL_BIN to the compiled PhysiCell project binary}"
# runtime libstdc++ must match the gcc the binary was built with (Bridges-2). Harmless if
# already loaded in the shell (its LD_LIBRARY_PATH is inherited); guarded so it never aborts.
module load gcc/13.3.1-p20240614 2>/dev/null || module load gcc 2>/dev/null || true
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
ROOT="$(git rev-parse --show-toplevel)"
ABM="$ROOT/results/abm"
MANIFEST="$ABM/patch_manifest.csv"
REPORT="$ABM/calibration_report.txt"
[ -f "$MANIFEST" ] || { echo "run 02-04 first (patch_manifest.csv missing)"; exit 1; }

# --- pick the median and max run by agent count (from the patch manifest) ---------------
read -r MED_RUN MED_N MAX_RUN MAX_N < <(python - "$MANIFEST" <<'PY'
import sys, pandas as pd
d = pd.read_csv(sys.argv[1]).sort_values("n_cells").reset_index(drop=True)
med = d.iloc[len(d)//2]; mx = d.iloc[-1]
print(med.run_id, int(med.n_cells), mx.run_id, int(mx.n_cells))
PY
)
echo "[calib] median run $MED_RUN (${MED_N} agents) | max run $MAX_RUN (${MAX_N} agents)"

: > "$REPORT"
COHORT_RUNS="$(wc -l < "$ABM/model_manifest.txt")"
REPLICATES="${REPLICATES:-10}"

run_one () {
  local run_id="$1" n0="$2" label="$3"
  local src="$ABM/$run_id" work="$ABM/$run_id/calib"
  [ -d "$src" ] || { echo "[skip] $run_id: model dir missing"; return; }
  rm -rf "$work"; mkdir -p "$work/output"
  cp "$src/cells.csv" "$src/rules.csv" "$src/PhysiCell_settings.xml" "$work/"
  echo "[calib] running $label $run_id ..."
  local t0 t1 secs
  t0="$(date +%s)"
  ( cd "$work" && "$PHYSICELL_BIN" PhysiCell_settings.xml >run.log 2>&1 ) \
    || { echo "  FAILED — see $work/run.log"; echo "$label $run_id: FAILED" >>"$REPORT"; return; }
  t1="$(date +%s)"; secs=$(( t1 - t0 ))
  # final agent count + sanity via the QoI extractor (pcdl)
  python "$ROOT/phase3_abm/qoi_extract.py" --run-dir "$work" --out "$work/qoi.csv" || true
  python - "$work/qoi.csv" "$n0" "$secs" "$label" "$run_id" "$COHORT_RUNS" "$REPLICATES" \
           "$OMP_NUM_THREADS" >>"$REPORT" <<'PY'
import sys, pandas as pd
qoi, n0, secs, label, run_id, runs, reps, cores = sys.argv[1:9]
n0, secs, runs, reps, cores = int(n0), int(secs), int(runs), int(reps), int(cores)
try:
    q = pd.read_csv(qoi).iloc[0].to_dict()
    nf = q.get("final_total_cells", float("nan")); fold = q.get("fold_growth", float("nan"))
except Exception:
    nf, fold = float("nan"), float("nan")
verdict = "OK"
if fold == fold:                              # not NaN
    if fold > 50:  verdict = "EXPLOSION (rates too fast?)"
    elif fold < 0.05: verdict = "COLLAPSE (all death?)"
core_h = secs/3600 * runs * reps * cores
print(f"{label} {run_id}: {secs}s wall | agents {n0}->{int(nf) if nf==nf else '?'} "
      f"(fold {fold:.2f}) | verdict {verdict}")
print(f"    cohort extrapolation @ this run: {runs} runs x {reps} reps x {cores} cores "
      f"~= {core_h:,.0f} core-hours (SU)")
PY
  echo "  done in ${secs}s"
}

run_one "$MED_RUN" "$MED_N" "median"
run_one "$MAX_RUN" "$MAX_N" "max"

echo; echo "===== calibration report ====="; cat "$REPORT"
echo
echo "GO/NO-GO: proceed to 06_run_cohort.sh only if both verdicts are OK and the max run"
echo "finished well within the 2 h SLURM cap. Report saved -> $REPORT"
