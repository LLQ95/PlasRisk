#!/bin/bash
# =============================================================================
# rerun_with_new_weights.sh
# Re-run all analyses that depend on risk weights after pipdb_17 determines
# the data-driven final weights.
# =============================================================================
set -euo pipefail

RES_DIR="${1:-results}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

cd "$SCRIPT_DIR"

echo "============================================================"
echo " PlasRisk: re-running analyses with data-driven final weights"
echo " Results directory: $RES_DIR"
echo "============================================================"

run_step() {
  local script="$1"
  local desc="$2"
  echo ""
  echo ">>> [$desc] Running $script ..."
  if Rscript "$script" "$RES_DIR"; then
    echo "<<< [$desc] Completed successfully."
  else
    echo "!!! [$desc] FAILED (exit code $?). Continuing with next step ..."
  fi
}

run_step "pipdb_17_weight_optimization.R" "Step 1/6: Weight optimization"
run_step "pipdb_15_risk_weight_table.R" "Step 2/6: Per-replicon risk table"
run_step "pipdb_16_risk_10dimensions.R" "Step 3/6: 10-dim comparison"
run_step "pipdb_14_network_risk_hotspot.R" "Step 4/6: Network/risk validation"
run_step "pipdb_18_external_validation_benchmark.R" "Step 5/6: External validation"
run_step "pipdb_19_case_study.R" "Step 6/6: Case study"

echo ""
echo "============================================================"
echo " All done. Updated tables in $RES_DIR/tables/"
echo "============================================================"
