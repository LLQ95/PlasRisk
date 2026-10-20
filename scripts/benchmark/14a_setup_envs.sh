#!/usr/bin/env bash
# 14a_setup_envs.sh
#
# One-time, manually run environment + database setup for the cross-tool
# comparison (script 14). Split out of 14_rfplasmid_plasmer.sh so that conda
# solving and the large database downloads can be run step by step.
#
# Usage:
#   bash 14a_setup_envs.sh [results/traveler_isolates]
#
# Optional overrides (export before running):
#   CHECKM_DB_PATH=/path/to/already/downloaded/checkm_data_2015_01_16
#   PLASMER_DB=/path/to/already/extracted/plasmer_db
#
# Creates:
#   conda env rfplasmid_cmp   (RFPlasmid 0.0.18; model DB built by --initialize)
#   conda env plasmer_cmp     (iskoldt/plasmer conda package, launcher on PATH)
#   <outdir>/cross_tool/db/checkm/           (CheckM1 marker DB, ~1.2 GB)
#   <outdir>/cross_tool/db/plasmer_db/       (Plasmer DB, ~14 GB xz, ~40 GB raw)
#
# Notes verified against the source:
#   * RFPlasmid 0.0.18 installs from the PyPI wheel shipped by bioconda.
#     Its wrapper only recognises --initialize (there is NO --prepareDB), and
#     --initialize prints "installed. Please try rfplasmid again" and exits 1
#     BY DESIGN after building both diamond databases. Exit 1 is not a failure;
#     readiness is judged by the presence of plasmiddb_cge.dmnd and
#     plasmiddb_total.dmnd inside the RFPlasmid package directory.
#   * RFPlasmid runs `checkm taxonomy_wf ... -t 16`, so the CheckM1 marker
#     database must be configured inside rfplasmid_cmp.
#   * The official Plasmer repository is https://github.com/nekokoe/Plasmer
#     (an older URL referenced in earlier drafts, tiagolbiotech/Plasmer,
#     returns 404). The conda package is published on the iskoldt channel.
#
# Safe to re-run; every step is idempotent.

set -uo pipefail

# English locale for all log output (timestamps, tool messages).
export LC_ALL=C

ISO_DIR="${1:-results/traveler_isolates}"
OUTDIR="$ISO_DIR/cross_tool"
DBROOT="$OUTDIR/db"
ENV_RFP="rfplasmid_cmp"
ENV_PLM="plasmer_cmp"

CHECKM_DB_URL="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
PLASMER_MAIN_URL="https://zenodo.org/api/records/7030675/files/plasmerMainDB.tar.xz/content"
PLASMER_KRAK_URL="https://zenodo.org/api/records/7030675/files/customizedKraken2DB.tar.xz/content"
PLASMER_MAIN_SHA1="0b08f5c30d60b137f54de6024ab7557031850db6"
PLASMER_KRAK_SHA1="b14efdd9232fd5f6d066716bd8e3e6ca80c9c0de"

mkdir -p "$DBROOT"
command -v conda >/dev/null 2>&1 || { echo "[error] conda not in PATH"; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "[setup] disk space under $OUTDIR"
df -h "$OUTDIR" || true

# ================================================================ RFPlasmid
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_RFP"; then
  echo "[setup] creating $ENV_RFP (heavy solve; be patient)"
  conda create -y -n "$ENV_RFP" -c conda-forge -c bioconda \
    rfplasmid python=3.10 || { echo "[error] rfplasmid env creation failed"; exit 20; }
else
  echo "[setup] $ENV_RFP already exists, skipping creation"
fi

conda activate "$ENV_RFP"

RFP_DIR="$(python3 -c 'import RFPlasmid, os; print(os.path.dirname(RFPlasmid.__file__))' 2>/dev/null)"
if [[ -z "$RFP_DIR" ]]; then
  echo "[error] cannot import RFPlasmid inside $ENV_RFP"
  exit 21
fi
echo "[setup] RFPlasmid package directory: $RFP_DIR"

if [[ -s "$RFP_DIR/plasmiddb_cge.dmnd" && -s "$RFP_DIR/plasmiddb_total.dmnd" ]]; then
  echo "[setup] RFPlasmid diamond databases already present, skipping initialize"
else
  echo "[setup] running 'rfplasmid --initialize' (it ALWAYS exits 1 after"
  echo "        printing 'installed. Please try rfplasmid again'; that is"
  echo "        expected and is not an error)"
  rfplasmid --initialize || true
fi
if [[ ! -s "$RFP_DIR/plasmiddb_cge.dmnd" || ! -s "$RFP_DIR/plasmiddb_total.dmnd" ]]; then
  echo "[error] RFPlasmid model databases missing after --initialize."
  echo "        Activate the env and inspect manually:"
  echo "        conda activate $ENV_RFP && rfplasmid --initialize; ls -lh $RFP_DIR"
  exit 21
fi
echo "[setup] RFPlasmid model databases ready:"
ls -lh "$RFP_DIR/plasmiddb_cge.dmnd" "$RFP_DIR/plasmiddb_total.dmnd"

# ------------------------------------------------- RFPlasmid pandas 3.x fix
# RFPlasmid 0.0.18 passes n positionally to Series.str.split(pat, n) and uses
# delim_whitespace=True in pd.read_csv; both are removed in pandas 3.x. With a
# fresh env the crash is:
#   TypeError: StringMethods.split() takes from 1 to 2 positional arguments
# right after "start blast plasmiddb". Patch the installed source (idempotent).
RFP_PY="$RFP_DIR/rfplasmid.py"
if [[ -f "$RFP_PY" ]] && grep -qE "str\.split\('_',( 1| n=1)\)\.str" "$RFP_PY"; then
  echo "[setup] patching rfplasmid.py for pandas >=2.2 compatibility"
  [[ -f "$RFP_PY.bak" ]] || cp "$RFP_PY" "$RFP_PY.bak"
  python3 - "$RFP_PY" <<'PYEOF'
import sys
p = sys.argv[1]
with open(p, "r", encoding="utf-8") as fh:
    txt = fh.read()
# both the original buggy form and the half-patched n=1 form keep the
# pre-0.23 ".str unpacking" idiom; expand=True is the modern equivalent.
n_split = txt.count("str.split('_', 1).str") + txt.count("str.split('_', n=1).str")
n_read = txt.count("delim_whitespace=True")
txt = txt.replace("str.split('_', 1).str", "str.split('_', n=1, expand=True)")
txt = txt.replace("str.split('_', n=1).str", "str.split('_', n=1, expand=True)")
txt = txt.replace("df1['contig'], df1['gene'] =", "df1[['contig', 'gene']] =")
txt = txt.replace("df['contig'], df['gene'] =", "df[['contig', 'gene']] =")
txt = txt.replace("delim_whitespace=True", "sep=r'\\s+'")
with open(p, "w", encoding="utf-8") as fh:
    fh.write(txt)
print("[setup] patched %d split lines, %d read_csv lines" % (n_split, n_read))
PYEOF
  python3 -m py_compile "$RFP_PY" || { echo "[error] rfplasmid.py patch broke syntax"; exit 28; }
  if grep -qE "str\.split\('_',( 1| n=1)\)\.str" "$RFP_PY"; then
    echo "[error] rfplasmid.py split patch not applied"; exit 28
  fi
  echo "[setup] RFPlasmid pandas compatibility patch applied"
else
  echo "[setup] RFPlasmid source already patched (or not found)"
fi

# ---------------------------------------------------------------- CheckM DB
if [[ -n "${CHECKM_DB_PATH:-}" ]]; then
  CHECKM_DIR="$CHECKM_DB_PATH"
else
  CHECKM_DIR="$DBROOT/checkm"
fi
if [[ -s "$CHECKM_DIR/taxon_marker_sets.tsv" ]]; then
  echo "[setup] CheckM marker DB found at $CHECKM_DIR"
else
  if [[ -n "${CHECKM_DB_PATH:-}" ]]; then
    echo "[error] CHECKM_DB_PATH=$CHECKM_DB_PATH has no taxon_marker_sets.tsv"
    exit 26
  fi
  mkdir -p "$CHECKM_DIR"
  echo "[setup] downloading CheckM marker DB (~1.2 GB)"
  curl -L --retry 3 -C - -o "$DBROOT/checkm_data_2015_01_16.tar.gz" \
    "$CHECKM_DB_URL" || {
    echo "[error] CheckM DB download failed. Download manually:"
    echo "  $CHECKM_DB_URL"
    echo "and extract into $CHECKM_DIR (it must contain taxon_marker_sets.tsv),"
    echo "then re-run this script."; exit 26; }
  echo "[setup] extracting CheckM DB"
  tar -xzf "$DBROOT/checkm_data_2015_01_16.tar.gz" -C "$CHECKM_DIR" \
    || { echo "[error] CheckM DB extraction failed"; exit 26; }
fi
# The official tar contains a top-level checkm_data_2015_01_16/ directory.
# Locate taxon_marker_sets.tsv (the file CheckM expects directly under dataRoot).
# Do not pipe to head: SIGPIPE under pipefail aborts the command substitution.
if [[ ! -s "$CHECKM_DIR/taxon_marker_sets.tsv" ]]; then
  nested="$(find "$CHECKM_DIR" -mindepth 2 -maxdepth 3 -name taxon_marker_sets.tsv -print -quit)"
  if [[ -n "$nested" ]]; then
    CHECKM_DIR="$(dirname "$nested")"
    echo "[setup] nested directory detected; using $CHECKM_DIR"
  fi
fi
[[ -s "$CHECKM_DIR/taxon_marker_sets.tsv" ]] || {
  echo "[error] taxon_marker_sets.tsv not found after extraction"; exit 26; }
# Absolute path is required by the env activation script below.
CHECKM_ABS="$(cd "$CHECKM_DIR" && pwd)"
echo "[setup] pointing CheckM data root to $CHECKM_ABS"
# setRoot with a path argument is non-interactive; never pipe 'yes' into it.
# With pipefail, 'yes | checkm ...' inherits yes' SIGPIPE status (141) even
# though CheckM exits 0; that caused the earlier false failure (exit 26).
set +e
checkm data setRoot "$CHECKM_ABS" < /dev/null > "$DBROOT/checkm_setroot.log" 2>&1
SR_RC=$?
set -uo pipefail
if [[ $SR_RC -ne 0 ]]; then
  echo "[warning] 'checkm data setRoot' returned $SR_RC; verifying functionally"
  echo "[warning] setRoot log: $DBROOT/checkm_setroot.log"
fi
# CHECKM_DATA_PATH is checked first in DBManager and bypasses DATA_CONFIG.
# Pin it on every 'conda activate rfplasmid_cmp' (RFPlasmid spawns checkm).
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/checkm_data_path.sh" <<EOF
export CHECKM_DATA_PATH="$CHECKM_ABS"
EOF
cat > "$CONDA_PREFIX/etc/conda/deactivate.d/checkm_data_path.sh" <<'EOF'
unset CHECKM_DATA_PATH
EOF
export CHECKM_DATA_PATH="$CHECKM_ABS"
echo "[setup] CHECKM_DATA_PATH pinned in env activation: $CHECKM_ABS"
# taxon_list reads the marker DB; it is the functional check.
checkm taxon_list > "$DBROOT/checkm_taxon_list.log" 2>&1 || {
  echo "[error] CheckM DB not usable (checkm taxon_list failed); see $DBROOT/checkm_taxon_list.log"
  exit 27; }
echo "[setup] CheckM DB verified (taxon_list OK)"

conda deactivate
echo "[setup] RFPlasmid + CheckM ready"

# ================================================================== Plasmer
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_PLM"; then
  echo "[setup] creating $ENV_PLM (pulls R, kmer-db, kraken2, infernal, ...)"
  conda create -y -n "$ENV_PLM" \
    -c iskoldt -c bioconda -c conda-forge -c defaults \
    plasmer || { echo "[error] plasmer env creation failed"; exit 22; }
else
  echo "[setup] $ENV_PLM already exists, skipping creation"
fi

conda activate "$ENV_PLM"
command -v Plasmer >/dev/null 2>&1 || {
  echo "[error] Plasmer launcher not on PATH inside $ENV_PLM"; exit 22; }
conda deactivate

PLASMER_DB="${PLASMER_DB:-$DBROOT/plasmer_db}"
mkdir -p "$PLASMER_DB"

# Optional institutional proxy: export DOWNLOAD_PROXY=http://host:port
# before running (otherwise the standard https_proxy/HTTPS_PROXY vars are
# honoured by both aria2c and curl).
if [[ -n "${DOWNLOAD_PROXY:-}" ]]; then
  export http_proxy="$DOWNLOAD_PROXY" https_proxy="$DOWNLOAD_PROXY"
fi

download_and_verify() {
  local url="$1" archive="$2" sha1_expected="$3"
  local dir base
  dir="$(dirname "$archive")"
  base="$(basename "$archive")"
  if [[ -f "$archive" ]]; then
    local sha1_have
    sha1_have="$(sha1sum "$archive" | awk '{print $1}')"
    if [[ "$sha1_have" == "$sha1_expected" ]]; then
      echo "[setup] $base already downloaded and checksum OK"
      return 0
    fi
    echo "[setup] partial/corrupt $base, resuming"
  fi
  if command -v aria2c >/dev/null 2>&1; then
    # Zenodo accepts range requests (206, accept-ranges: bytes). Parallel
    # connections aggregate throughput on slow international links and
    # resume automatically after EOF drops (curl 56).
    echo "[setup] aria2c downloading $base (16 connections, unlimited retries)"
    aria2c -c -x 16 -s 16 -k 1M --file-allocation=none \
      --max-tries=0 --retry-wait=5 --connect-timeout=30 --timeout=60 \
      --console-log-level=warn --summary-interval=15 \
      -d "$dir" -o "$base" "$url" || {
      echo "[error] aria2c failed: $url"; return 1; }
  else
    # --retry-all-errors covers curl 56 (unexpected EOF); -C - resumes.
    # Abort and retry when speed drops below 1 KB/s for 60 s.
    echo "[setup] curl downloading $base (resume, 10 retries)"
    curl -L --retry 10 --retry-all-errors --retry-delay 5 -C - \
      --speed-time 60 --speed-limit 1024 \
      -o "$archive" "$url" || {
      echo "[error] download failed: $url"; return 1; }
  fi
  local sha1_now
  sha1_now="$(sha1sum "$archive" | awk '{print $1}')"
  if [[ "$sha1_now" != "$sha1_expected" ]]; then
    echo "[error] sha1 mismatch for $base:"
    echo "        expected $sha1_expected"
    echo "        have     $sha1_now"
    return 1
  fi
  echo "[setup] $base sha1 OK"
}

if [[ -s "$PLASMER_DB/plasmerMainDB.list" && \
      -s "$PLASMER_DB/customizedKraken2DB.list" ]]; then
  echo "[setup] Plasmer DB already extracted at $PLASMER_DB"
else
  echo "[setup] Plasmer DB target: $PLASMER_DB (~14 GB archives, ~40 GB extracted)"
  df -h "$PLASMER_DB" || true
  download_and_verify "$PLASMER_MAIN_URL" \
    "$DBROOT/plasmerMainDB.tar.xz" "$PLASMER_MAIN_SHA1" || exit 24
  download_and_verify "$PLASMER_KRAK_URL" \
    "$DBROOT/customizedKraken2DB.tar.xz" "$PLASMER_KRAK_SHA1" || exit 24
  echo "[setup] extracting plasmerMainDB"
  tar -xJf "$DBROOT/plasmerMainDB.tar.xz" -C "$PLASMER_DB" || {
    echo "[error] plasmerMainDB extraction failed"; exit 24; }
  echo "[setup] extracting customizedKraken2DB"
  tar -xJf "$DBROOT/customizedKraken2DB.tar.xz" -C "$PLASMER_DB" || {
    echo "[error] customizedKraken2DB extraction failed"; exit 24; }
  if [[ ! -s "$PLASMER_DB/plasmerMainDB.list" ]]; then
    nested="$(find "$PLASMER_DB" -mindepth 2 -maxdepth 3 -name plasmerMainDB.list -print -quit)"
    [[ -n "$nested" ]] && PLASMER_DB="$(dirname "$nested")"
  fi
  if [[ -s "$PLASMER_DB/plasmerMainDB.list" && \
        -s "$PLASMER_DB/customizedKraken2DB.list" ]]; then
    echo "[setup] Plasmer DB extracted"
    rm -f "$DBROOT/plasmerMainDB.tar.xz" "$DBROOT/customizedKraken2DB.tar.xz"
  else
    echo "[error] Plasmer DB marker files missing after extraction."
    echo "        plasmerMainDB.list and customizedKraken2DB.list must be in $PLASMER_DB"
    exit 25
  fi
fi

echo "[setup] all environments and databases ready."
echo "[setup] next run (from the genomad env, which has pandas):"
echo "  conda activate genomad"
echo "  bash 14b_run_cross_tool.sh $ISO_DIR 16 4 1"
echo "[setup] recommended smoke test before 14b:"
echo "  conda activate $ENV_RFP"
echo "  mkdir -p /tmp/rfp_test/in && ln -sf \$PWD/$ISO_DIR/assemblies/<acc>.fna /tmp/rfp_test/in/<acc>.fasta"
echo "  rfplasmid --species Enterobacteriaceae --input /tmp/rfp_test/in --out /tmp/rfp_test/out --threads 8"
