#!/usr/bin/env bash
# pipdb_07_beast_dating.sh — Run BEAST2 dating for each replicon
set -euo pipefail
EVO_DIR="${1:-results/evolution}"
THREADS="${2:-4}"
CHAIN="${3:-100000000}"

if ! command -v beast >/dev/null 2>&1; then
  echo "ERROR: beast not found. conda install -c bioconda beast2=2.7"; exit 1
fi

for xml in "$EVO_DIR"/*/tree/*_beast.xml; do
  [ -f "$xml" ] || continue
  rep_dir="$(dirname "$xml")"
  log_file="${xml%.xml}.log"
  trees_file="${xml%.xml}.trees"
  if [ -f "$trees_file" ]; then
    echo "[skip] $(basename "$xml"): trees file exists"; continue
  fi
  echo "=== Running BEAST2: $(basename "$xml") ==="
  (cd "$rep_dir" && beast -threads "$THREADS" -overwrite "$(basename "$xml")" 2>&1 | tail -5)
done

echo ""
echo "=== Post-processing with Tracer / LogCombiner / TreeAnnotator ==="
echo "1. Check ESS > 200 in Tracer"
echo "2. LogCombiner -burnin 10% to combine runs"
echo "3. TreeAnnotator -heights ca -burnin 10% for MCC tree"
