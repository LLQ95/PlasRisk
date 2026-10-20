#!/usr/bin/env bash
# 15_cef_resume_missing.sh
#
# P0-5: resume the four cefprozil samples missing from the 68-sample
# collection: P1E0, P1E7, P1E90 (all three P1 timepoints) and P10E0
# (console reported DONE but the log stops at a bad-gzip aria2c error
# and the summary excludes it).
#
# Strategy:
#   1. build a 4-sample sub-manifest from the full chunk manifest;
#   2. verify every existing chunk with gzip -t; delete corrupt chunks
#      and their aria2 control files so aria2c resumes cleanly;
#   3. remove the partial assembly/geNomad/PlasRisk outputs for the
#      four targets so they are rebuilt end to end (other 68 samples
#      are untouched);
#   4. run the resume-safe 07 pipeline on the sub-manifest;
#   5. rerun 07 on the FULL manifest into the same outdir (all samples
#      are skipped, the collector rebuilds the full summary), then
#      rerun 08_paired_analysis.py.
#
# Usage:
#   bash 15_cef_resume_missing.sh \
#        cefprozil_tpa_chunks_manifest.tsv results/cefprozil 16 8 \
#        [path/to/07_cefprozil_tpa_benchmark.sh] [path/to/08_paired_analysis.py]

set -uo pipefail

FULL_MANIFEST="${1:-}"
OUTDIR="${2:-}"
THREADS="${3:-16}"
JOBS="${4:-8}"
SCRIPT_07="${5:-$(dirname "$0")/../07_cefprozil_tpa_benchmark.sh}"
SCRIPT_08="${6:-$(dirname "$0")/../08_paired_analysis.py}"
TARGETS="P1E0 P1E7 P1E90 P10E0"

if [[ -z "$FULL_MANIFEST" || ! -s "$FULL_MANIFEST" || -z "$OUTDIR" ]]; then
  echo "usage: bash 15_cef_resume_missing.sh <full_manifest.tsv> <outdir> [threads] [jobs] [07.sh] [08.py]"
  exit 1
fi
for s in "$SCRIPT_07" "$SCRIPT_08"; do
  [[ -s "$s" ]] || { echo "[error] missing $s"; exit 1; }
done
SCRIPT_07="$(cd "$(dirname "$SCRIPT_07")" && pwd)/$(basename "$SCRIPT_07")"
SCRIPT_08="$(cd "$(dirname "$SCRIPT_08")" && pwd)/$(basename "$SCRIPT_08")"
BASEDIR="$(cd "$(dirname "$OUTDIR")" 2>/dev/null && pwd)/$(basename "$OUTDIR")"
mkdir -p "$BASEDIR/00_manifest"
SUB_MANIFEST="$BASEDIR/00_manifest/missing4.tsv"

# 1. sub-manifest
head -1 "$FULL_MANIFEST" > "$SUB_MANIFEST"
for sid in $TARGETS; do
  awk -F'\t' -v s="$sid" 'NR==1 || $1==s' "$FULL_MANIFEST" | tail -n +2 >> "$SUB_MANIFEST"
done
echo "[resume] sub-manifest rows: $(($(wc -l < "$SUB_MANIFEST") - 1))"

# 2. diagnostics + cleanup for the four targets
for sid in $TARGETS; do
  echo "----- $sid -----"
  cdir="$BASEDIR/chunks/$sid"
  if [[ -d "$cdir" ]]; then
    bad=0
    while IFS= read -r f; do
      if ! gzip -t "$f" 2>/dev/null; then
        echo "  corrupt chunk: $(basename "$f") -> removing"
        rm -f "$f" "$f.aria2"
        bad=$((bad+1))
      fi
    done < <(find "$cdir" -type f -name '*.fna.gz' 2>/dev/null)
    find "$cdir" -type f -name '*.fna.gz.aria2' -delete 2>/dev/null
    echo "  corrupt chunks removed: $bad"
  else
    echo "  no chunks directory yet"
  fi
  # state of downstream outputs
  for rel in "assemblies/$sid.fna" "genomad/$sid" "plasrisk/$sid"; do
    p="$BASEDIR/$rel"
    if [[ -e "$p" ]]; then
      echo "  exists: $rel"
    fi
  done
  pr="$BASEDIR/plasrisk/$sid/meta_sample_summary.tsv"
  if [[ -s "$pr" ]]; then
    echo "  WARNING: meta_sample_summary exists but sample was missing from"
    echo "           the published summary; rebuilding anyway"
  fi
  # force rebuild from (verified/new) chunks
  rm -f "$BASEDIR/assemblies/$sid.fna"
  rm -rf "$BASEDIR/genomad/$sid" "$BASEDIR/plasrisk/$sid"
  rm -f "$BASEDIR/logs/$sid.tpa.log"
done

# 3. resume the four samples (07 skips every chunk that is already valid)
bash "$SCRIPT_07" "$SUB_MANIFEST" "$BASEDIR" "$THREADS" "$JOBS"
rc1=$?

# 4. rerun on the full manifest: no processing (all skipped), full
#    collector rebuild of benchmark_sample_summary.tsv
bash "$SCRIPT_07" "$FULL_MANIFEST" "$BASEDIR" "$THREADS" "$JOBS"
rc2=$?

# 5. paired analysis
python3 "$SCRIPT_08" "$BASEDIR/benchmark_sample_summary.tsv" "$BASEDIR/paired"
rc3=$?

echo "[resume] exit codes: 07-sub=$rc1 07-full=$rc2 08=$rc3"
echo "[resume] failures (if any):"
cat "$BASEDIR/failed_samples.tsv" 2>/dev/null || true
echo "[resume] sample rows now: $(($(wc -l < "$BASEDIR/benchmark_sample_summary.tsv") - 1))"

# quick per-target confirmation
python3 - "$BASEDIR" <<'PYEOF'
import sys, os, pandas as pd
base = sys.argv[1]
p = os.path.join(base, "benchmark_sample_summary.tsv")
d = pd.read_csv(p, sep="\t")
for sid in ["P1E0", "P1E7", "P1E90", "P10E0"]:
    row = d[d["sample_id"] == sid]
    if len(row):
        r = row.iloc[0]
        print("%s: n_predicted_bins=%s CPRI=%s"
              % (sid, r.get("n_predicted_bins"), r.get("pred_CPRI_mean")))
    else:
        print("%s: STILL MISSING" % sid)
PYEOF
