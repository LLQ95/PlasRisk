#!/usr/bin/env bash
# 04_preassembled_benchmark.sh
#
# Run the geNomad -> PlasRisk meta pipeline on PRE-ASSEMBLED metagenome
# contigs. This skips fastp, megahit and the coverage step entirely, so it is
# the fast route when deposited assemblies (ENA ERZ / MGnify TPA / author
# FASTA) already exist. Build manifests with make_assembly_manifest.py.
#
# Manifest TSV (no header, four columns):
#     id <TAB> habitat_group <TAB> pair_id <TAB> url
#   url is either an https://...fasta.gz / .fna link or a local FASTA path.
#
# Usage:
#   bash 04_preassembled_benchmark.sh assemblies.tsv results/preassembled 16 4
#
# Output layout (parallel to 02_public_metag_benchmark.sh):
#   outdir/assemblies/<id>.fna
#   outdir/genomad/<id>/...
#   outdir/plasrisk/<id>/meta_sample_summary.tsv
#   outdir/benchmark_sample_summary.tsv   <-- fed directly to 03_gradient_stats.py
#
# Note: coverage is not available for deposited assemblies, so plasrisk meta
# runs without --coverage. CPRI then uses per-contig risk scores and the
# plasmid prevalence statistic only (no abundance weighting).

set -uo pipefail

MANIFEST="${1:-}"
OUTDIR="${2:-}"
THREADS="${3:-16}"
JOBS="${4:-4}"

if [[ -z "$MANIFEST" || -z "$OUTDIR" ]]; then
  echo "usage: bash 04_preassembled_benchmark.sh <manifest.tsv> <outdir> [threads=16] [jobs=4]"
  exit 1
fi
if [[ ! -s "$MANIFEST" ]]; then
  echo "[error] manifest not found or empty: $MANIFEST"
  exit 1
fi

BASEDIR="$(cd "$(dirname "$OUTDIR")" 2>/dev/null && pwd)/$(basename "$OUTDIR")"
mkdir -p "$BASEDIR/assemblies" "$BASEDIR/genomad" "$BASEDIR/plasrisk" "$BASEDIR/logs"
LOG_DIR="$BASEDIR/logs"
FAIL_LOG="$BASEDIR/failed_samples.tsv"
: > "${FAIL_LOG}.tmp"

command -v python3 >/dev/null 2>&1 || { echo "[error] python3 not in PATH"; exit 1; }
command -v genomad >/dev/null 2>&1   || { echo "[error] genomad not in PATH"; exit 1; }
command -v plasrisk >/dev/null 2>&1  || { echo "[error] plasrisk not in PATH"; exit 1; }
if command -v aria2c >/dev/null 2>&1; then
  DL_BACKEND="aria2c"
elif command -v curl >/dev/null 2>&1; then
  DL_BACKEND="curl"
else
  echo "[error] need aria2c or curl for downloads"
  exit 1
fi
GENOMAD_DB="${GENOMAD_DB:-}"
if [[ -z "$GENOMAD_DB" ]]; then
  echo "[error] GENOMAD_DB not set; export it to the geNomad database directory"
  exit 1
fi

GE_VER="$(genomad --version 2>/dev/null | grep -oE '^[0-9]+\.[0-9]+' | head -1 || true)"
GE_ARG=""
if [[ -n "$GE_VER" ]]; then
  case "$GE_VER" in
    1.12*|1.13*|1.14*|1.15*|1.16*|1.17*|1.18*|1.19*|1.2*|1.[2-9][0-9]*) GE_ARG="--genomad-version v1.12" ;;
    *) GE_ARG="" ;;
  esac
fi
echo "[asm] backend: ${DL_BACKEND} | jobs: ${JOBS} | threads/job: ${THREADS} | gen: ${GE_VER:-?}"

download_one() {
  local id="$1" url="$2" asm_gz="$3" asm_fa="$4"
  if [[ -s "$asm_fa" ]]; then return 0; fi
  if [[ "$url" == http* ]]; then
    rm -f "$asm_gz"
    if [[ "$DL_BACKEND" == "aria2c" ]]; then
      aria2c -x 8 -s 8 -k 1M --file-allocation=none --console-log-level=warn \
        --summary-interval=0 --retry-wait=5 -m 0 -d "$(dirname "$asm_gz")" \
        -o "$(basename "$asm_gz")" "$url"
    else
      curl -L --retry 4 --retry-delay 5 -C - -o "$asm_gz" "$url"
    fi
  else
    # local file: link/copy, decompress below if it is gzipped
    if [[ ! -s "$url" ]]; then echo "local file not found: $url"; return 1; fi
    if [[ "$url" == *.gz ]]; then
      ln -sf "$(realpath "$url")" "$asm_gz" 2>/dev/null || cp "$url" "$asm_gz"
    else
      cp "$url" "$asm_fa"
      return 0
    fi
  fi
  gzip -t "$asm_gz" 2>/dev/null || { echo "bad gzip: $asm_gz"; return 1; }
  gunzip -c "$asm_gz" > "$asm_fa" && rm -f "$asm_gz"
  grep -c "^>" "$asm_fa" >/dev/null 2>&1 || { echo "no sequences in $asm_fa"; return 1; }
}

process_one() {
  local id="$1" habitat="$2" url="$3"
  local lg="$LOG_DIR/${id}.asm.log"
  local asm_fa="$BASEDIR/assemblies/${id}.fna"
  local gd="$BASEDIR/genomad/$id"
  local pr="$BASEDIR/plasrisk/$id"
  local stem; stem="$(basename "$asm_fa" | sed 's/\.[^.]*$//')"
  local gsum="$gd/${id}/${stem}_summary/${stem}_plasmid_summary.tsv"

  {
    echo "[asm] $id start"
    download_one "$id" "$url" "$BASEDIR/assemblies/${id}.fna.gz" "$asm_fa" || exit 11
    echo "[asm] $id contigs ready"

    if [[ ! -s "$gsum" ]]; then
      echo "  [genomad] running"
      genomad end-to-end --cleanup --threads "$THREADS" "$asm_fa" "$gd" "$GENOMAD_DB" \
        || exit 12
    fi
    if [[ ! -s "$gsum" ]]; then echo "genomad plasmid summary missing: $gsum"; exit 13; fi

    if [[ ! -s "$pr/meta_sample_summary.tsv" ]]; then
      echo "  [plasrisk] meta (no coverage)"
      plasrisk meta "$asm_fa" --genomad-bin "$gd/$id" $GE_ARG \
        -o "$pr" || exit 14
    fi
    if [[ ! -s "$pr/meta_sample_summary.tsv" ]]; then
      echo "plasrisk meta summary missing"; exit 15
    fi
    echo "[asm] $id DONE"
  } > "$lg" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[asm] $id FAILED (exit $rc)"
    echo -e "${id}\t${habitat}\texit_${rc}" >> "${FAIL_LOG}.tmp"
    tail -n 20 "$lg"
  else
    echo "[asm] $id DONE"
  fi
}

slot() {
  while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 2; done
}

while IFS=$'\t' read -r id habitat pair url; do
  [[ -z "$id" || "$id" == \#* ]] && continue
  slot
  process_one "$id" "$habitat" "$url" &
done < "$MANIFEST"
wait
mv "${FAIL_LOG}.tmp" "$FAIL_LOG"
echo "failed samples: $(wc -l < "$FAIL_LOG")"

# --------------------------- collect sample summaries ------------------------
python3 - "$MANIFEST" "$BASEDIR" <<'PYEOF'
import os, sys
import pandas as pd
manifest, base = sys.argv[1], sys.argv[2]
meta = pd.read_csv(manifest, sep="\t", header=None,
                   names=["accession", "habitat_group", "pair_id", "url"])
rows = []
for _, m in meta.iterrows():
    p = os.path.join(base, "plasrisk", m["accession"],
                     "meta_sample_summary.tsv")
    if not os.path.exists(p):
        continue
    d = pd.read_csv(p, sep="\t")
    d["accession"] = m["accession"]
    d["habitat_group"] = m["habitat_group"]
    d["pair_id"] = m["pair_id"]
    rows.append(d)
if rows:
    out = pd.concat(rows, ignore_index=True)
    cols = ["accession", "habitat_group", "pair_id"] + \
           [c for c in out.columns if c not in ("accession", "habitat_group", "pair_id")]
    out = out[cols]
    out.to_csv(os.path.join(base, "benchmark_sample_summary.tsv"),
               sep="\t", index=False)
    print("samples collected:", len(out))
    if "CPRI_mean" in out.columns:
        print(out.groupby("habitat_group")["CPRI_mean"].describe())
else:
    print("no completed samples found")
PYEOF
