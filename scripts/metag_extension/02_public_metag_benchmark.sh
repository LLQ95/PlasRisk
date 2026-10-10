#!/usr/bin/env bash
# =============================================================================
# 02_public_metag_benchmark.sh
# Tier A public-metagenome benchmark, resumable + parallel + multi-connection.
#   Stage 1 download : ERR/SRR -> ENA direct FASTQ. Multi-connection with
#                      aria2c (preferred), Aspera fasp (optional), or curl
#                      resume as fallback. SRR falls back to prefetch if ENA
#                      does not mirror the run.
#   Stage 2 process  : fastp -> MEGAHIT -> depth -> geNomad -> plasrisk meta
# Per-sample failures are isolated and logged; re-running resumes from the
# last completed step (markers and existing outputs are reused).
#
# Credentials (kept OUT of this repo): the NCBI API key is read from the
#   NCBI_API_KEY environment variable, or from a shell snippet at
#   $PLASRIS_CREDENTIALS (default ~/.config/plasrisk/credentials) containing
#     export NCBI_API_KEY=your_key
#   The key raises NCBI E-utilities/prefetch request limits; it does not raise
#   bulk byte throughput. ENA multi-connection (aria2c) is the real speedup.
#
# Usage:
#   02_public_metag_benchmark.sh ACCESSIONS WORKDIR THREADS \
#                                [MODE] [DL_JOBS] [PROC_JOBS] [CONNS]
#   MODE      = all | download | process           (default all)
#   DL_JOBS   = parallel samples downloading       (default 6)
#   PROC_JOBS = parallel samples in stage 2, each
#               using THREADS; keep PROC_JOBS*THREADS <= cores (default 1)
#   CONNS     = connections per file for aria2c    (default 8, max 16)
# Set ENABLE_ASPERA=1 to prefer Aspera when the `ascp` binary is available.
# accessions.tsv columns (no header): accession,habitat_group,pair_id
# =============================================================================
set -uo pipefail   # no -e: per-sample failures must not abort the whole run

ACC_FILE="${1:?accessions.tsv required}"
BASEDIR="${2:?output dir required}"
THREADS="${3:-16}"
MODE="${4:-all}"
DL_JOBS="${5:-6}"
PROC_JOBS="${6:-1}"
CONNS="${7:-8}"

# Load the NCBI API key from env or a local, non-repo credentials file.
CRED_FILE="${PLASRIS_CREDENTIALS:-$HOME/.config/plasrisk/credentials}"
if [ -z "${NCBI_API_KEY:-}" ] && [ -f "${CRED_FILE}" ]; then
  # shellcheck disable=SC1090
  . "${CRED_FILE}" >/dev/null 2>&1 || true
fi
export NCBI_API_KEY
if [ -n "${NCBI_API_KEY:-}" ]; then echo "[cred] NCBI API key loaded"; fi

mkdir -p "${BASEDIR}"/{reads,qc,asm,depth,genomad,plasrisk,logs}
FAIL_LOG="${BASEDIR}/failed_samples.tsv"
: > "${FAIL_LOG}.tmp"

have(){ command -v "$1" >/dev/null 2>&1; }

# --------------------------- dependency self-check ---------------------------
miss_dl=0; miss_pr=0
have curl || { echo "MISSING download binary: curl"; miss_dl=1; }
# Choose the download backend. aria2c gives multi-connection throughput;
# Aspera (fasp) is used only when explicitly enabled and installed.
DL_BACKEND="curl"
if have aria2c; then DL_BACKEND="aria2c"; fi
if [ "${ENABLE_ASPERA:-0}" = "1" ] && have ascp; then DL_BACKEND="aspera"; fi
echo "[dl] backend: ${DL_BACKEND} (${DL_JOBS} samples x ${CONNS} conns)"

# Locate the Aspera transfer public key (shipped with aspera-cli / Aspera
# Connect; it is a non-secret public key). Fall back if it cannot be found.
ASPERA_KEY="${ASPERA_SCP_KEY:-}"
if [ -z "${ASPERA_KEY}" ]; then
  for k in \
    "${CONDA_PREFIX:-}/etc/asperaweb_id_dsa.openssh" \
    "${CONDA_PREFIX:-}/opt/aspera-cli/etc/asperaweb_id_dsa.openssh" \
    "$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh" \
    "/usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh"; do
    [ -f "$k" ] && { ASPERA_KEY="$k"; break; }
  done
fi
if [ "${DL_BACKEND}" = "aspera" ] && [ -z "${ASPERA_KEY}" ]; then
  echo "[dl] ascp present but public key not found (set ASPERA_SCP_KEY); falling back"
  if have aria2c; then DL_BACKEND="aria2c"; else DL_BACKEND="curl"; fi
fi
[ -n "${ASPERA_KEY}" ] && export ASPERA_SCP_KEY="${ASPERA_KEY}"
echo "[dl] aspera key: ${ASPERA_KEY:-<none>}"
for b in fastp megahit bowtie2 bowtie2-build samtools \
         jgi_summarize_bam_contig_depths genomad plasrisk; do
  have "$b" || { echo "MISSING process binary: $b"; miss_pr=1; }
done
# SRR runs only need sra-tools if ENA does not mirror them (prefetch fallback).
if grep -Eq '^[[:space:]]*SRR' "${ACC_FILE}"; then
  for b in prefetch fasterq-dump; do
    have "$b" || echo "NOTE (SRR fallback): $b not found; ENA direct will be used"
  done
fi
if [ "${MODE}" != "process" ] && [ "${miss_dl}" -ne 0 ]; then
  echo "install missing download tools, then re-run"; exit 1
fi
if [ "${MODE}" != "download" ]; then
  [ "${miss_pr}" -eq 0 ] || { echo "install missing process tools, then re-run"; exit 1; }
  : "${GENOMAD_DB:?set GENOMAD_DB to the genomad_db folder}"
fi

# Configure sra-tools with the API key once, if a supported flag exists.
if [ -n "${NCBI_API_KEY:-}" ] && have vdb-config; then
  if vdb-config --help 2>&1 | grep -q -- '--api-key'; then
    vdb-config --api-key "${NCBI_API_KEY}" >/dev/null 2>&1 || true
  fi
fi

# Wait until fewer than $1 background jobs are running.
slot(){ while [ "$(jobs -rp | wc -l)" -ge "$1" ]; do wait -n || true; done; }

# ENA FASTQ directory for an ERR/SRR accession (handles the 000-999 sub-buckets
# used once the numeric part exceeds six digits). ENA mirrors NCBI SRR runs.
ena_dir () {
  local ACC="$1"
  local b6="${ACC:0:6}"
  local digits="${ACC:3}"
  local d="/vol1/fastq/${b6}"
  if [ "${#digits}" -gt 6 ]; then
    local n=$(( ${#digits} - 6 ))
    local suf="${digits: -$n}"
    local sub
    sub="$(printf '%03d' "$((10#${suf}))")"
    d="${d}/${sub}"
  fi
  echo "${d}/${ACC}"
}

# Fetch one ENA FASTQ member into a target path.
#   ena_fetch <https_url> <out_path>
ena_fetch () {
  local URL="$1" OUT="$2"
  case "${DL_BACKEND}" in
    aria2c)
      # -c resume, -x/-s multi-connection, -m 0 keep retrying, falloc fast.
      aria2c -c -x "${CONNS}" -s "${CONNS}" -k 1M \
        --file-allocation=falloc --retry-wait=5 -m 0 \
        -d "$(dirname "${OUT}")" -o "$(basename "${OUT}")" "${URL}" ;;
    aspera)
      # fasp path is the same /vol1 tree on fasp.sra.ebi.ac.uk:33001.
      # -k1 enables resume; download into the directory (keeps file name).
      local rp="${URL#https://ftp.sra.ebi.ac.uk}"
      local keyargs=()
      [ -n "${ASPERA_KEY}" ] && keyargs=(-i "${ASPERA_KEY}")
      ascp -QT -k1 "${keyargs[@]}" \
        -l "${ASPERA_RATE:-300m}" -P 33001 \
        "fasp@fasp.sra.ebi.ac.uk:${rp}" "$(dirname "${OUT}")/" ;;
    *)
      curl -L -C - --retry 5 --retry-delay 5 --fail -o "${OUT}" "${URL}" ;;
  esac
}

# ------------------------------- Stage 1 download ----------------------------
dl_one () {
  local ACC="$1" RD="$2"
  local LOG="${BASEDIR}/logs/${ACC}.dl.log"
  mkdir -p "${RD}"
  if [ -f "${RD}/dl_done" ]; then echo "[dl] ${ACC} cached"; return 0; fi
  echo "[dl] ${ACC} start (${DL_BACKEND})" > "${LOG}"
  local d host url1 url2
  d="$(ena_dir "${ACC}")"
  host="https://ftp.sra.ebi.ac.uk"
  url1="${host}${d}_1.fastq.gz"
  url2="${host}${d}_2.fastq.gz"
  # ENA mirrors both ERR and SRR. Probe whether the member exists.
  ena_exists () {
    curl -L -s -o /dev/null --fail -r 0-0 "$1"
  }
  if ena_exists "${url1}"; then
    if ! ena_fetch "${url1}" "${RD}/${ACC}_1.fastq.gz" >>"${LOG}" 2>&1; then
      echo "[dl] ${ACC} FAILED (_1)"; return 1
    fi
    if ena_exists "${url2}"; then
      ena_fetch "${url2}" "${RD}/${ACC}_2.fastq.gz" >>"${LOG}" 2>&1 \
        || echo "[dl] ${ACC} no _2 (single-end?)" >>"${LOG}"
    else
      echo "[dl] ${ACC} no _2 on ENA (single-end?)" >>"${LOG}"
    fi
  elif [[ "${ACC}" == SRR* ]] && have prefetch; then
    # ENA does not mirror this run: fall back to NCBI prefetch.
    echo "[dl] ${ACC} ENA missing, prefetch fallback" >>"${LOG}"
    if ! prefetch "${ACC}" -O "${RD}" >>"${LOG}" 2>&1; then
      echo "[dl] ${ACC} FAILED (prefetch)"; return 1
    fi
    local sra
    sra="$(find "${RD}" -name '*.sra' 2>/dev/null | head -1)"
    if ! fasterq-dump "${sra:-$ACC}" -O "${RD}" -t "${RD}/tmp" \
          -e "${THREADS}" --progress >>"${LOG}" 2>&1; then
      echo "[dl] ${ACC} FAILED (fasterq)"; return 1
    fi
    have pigz && pigz -p "${THREADS}" "${RD}"/*.fastq \
              || gzip "${RD}"/*.fastq
  else
    echo "[dl] ${ACC} FAILED (not on ENA, no prefetch)"; return 1
  fi
  if ls "${RD}"/*_1.fastq.gz >/dev/null 2>&1; then
    touch "${RD}/dl_done"; echo "[dl] ${ACC} DONE"
  else
    echo "[dl] ${ACC} FAILED (no reads)"; return 1
  fi
}

if [ "${MODE}" != "process" ]; then
  echo "=== Stage 1: download (${DL_JOBS} parallel, resumable) ==="
  while IFS=$'\t' read -r ACC HAB PAIR; do
    [ -z "${ACC:-}" ] && continue
    case "${ACC}" in \#*) continue ;; esac
    slot "${DL_JOBS}"
    dl_one "${ACC}" "${BASEDIR}/reads/${ACC}" &
  done < "${ACC_FILE}"
  wait
fi

# ------------------------------- Stage 2 process -----------------------------
process_one () {
  local ACC="$1"
  local LOG="${BASEDIR}/logs/${ACC}.proc.log"
  echo "[proc] ${ACC} start"
  local RD="${BASEDIR}/reads/${ACC}"
  local R1 R2
  R1=$(ls "${RD}"/*_1.fastq.gz 2>/dev/null | head -1 || true)
  R2=$(ls "${RD}"/*_2.fastq.gz 2>/dev/null | head -1 || true)
  if [ -z "${R1}" ] || [ -z "${R2}" ]; then
    echo -e "${ACC}\tno paired reads" >> "${FAIL_LOG}.tmp"
    echo "[proc] ${ACC} skipped (no paired reads)"; return 0
  fi
  # Heavy body in an isolated subshell: any non-zero status is caught once.
  (
    set -e
    fastp -i "${R1}" -I "${R2}" \
      -o "${BASEDIR}/qc/${ACC}_1.fq.gz" -O "${BASEDIR}/qc/${ACC}_2.fq.gz" \
      -h "${BASEDIR}/qc/${ACC}.html" -j "${BASEDIR}/qc/${ACC}.json" \
      -w "${THREADS}"

    if [ ! -f "${BASEDIR}/asm/${ACC}/final.contigs.fa" ]; then
      megahit -1 "${BASEDIR}/qc/${ACC}_1.fq.gz" \
              -2 "${BASEDIR}/qc/${ACC}_2.fq.gz" \
              -o "${BASEDIR}/asm/${ACC}" -t "${THREADS}" --presets meta-large
    fi
    ASM="${BASEDIR}/asm/${ACC}/final.contigs.fa"

    bowtie2-build -q "${ASM}" "${BASEDIR}/depth/${ACC}_idx"
    bowtie2 -x "${BASEDIR}/depth/${ACC}_idx" \
      -1 "${BASEDIR}/qc/${ACC}_1.fq.gz" -2 "${BASEDIR}/qc/${ACC}_2.fq.gz" \
      -p "${THREADS}" --very-fast -S "${BASEDIR}/depth/${ACC}.sam"
    samtools view -bS "${BASEDIR}/depth/${ACC}.sam" \
      | samtools sort -@ "${THREADS}" -o "${BASEDIR}/depth/${ACC}.bam"
    samtools index "${BASEDIR}/depth/${ACC}.bam"
    jgi_summarize_bam_contig_depths \
      --outputDepth "${BASEDIR}/depth/${ACC}_depth.tsv" \
      "${BASEDIR}/depth/${ACC}.bam"

    if ! find "${BASEDIR}/genomad/${ACC}" \
          -name '*_plasmid_summary.tsv' 2>/dev/null | grep -q .; then
      genomad end-to-end --cleanup --threads "${THREADS}" \
        "${ASM}" "${BASEDIR}/genomad/${ACC}" "${GENOMAD_DB}"
    fi

    plasrisk meta "${ASM}" \
      --genomad-bin "${BASEDIR}/genomad/${ACC}" \
      --coverage "${BASEDIR}/depth/${ACC}_depth.tsv" \
      --coverage-column "${ACC}" \
      -o "${BASEDIR}/plasrisk/${ACC}" \
      --min-contig 1000 --complete-length 10000 --threads "${THREADS}"
  ) >> "${LOG}" 2>&1

  if [ $? -eq 0 ]; then
    echo "[proc] ${ACC} DONE"
  else
    echo -e "${ACC}\tprocessing error (see logs/${ACC}.proc.log)" >> "${FAIL_LOG}.tmp"
    echo "[proc] ${ACC} FAILED"
  fi
}

if [ "${MODE}" != "download" ]; then
  echo "=== Stage 2: process (${PROC_JOBS} parallel x ${THREADS} threads) ==="
  while IFS=$'\t' read -r ACC HAB PAIR; do
    [ -z "${ACC:-}" ] && continue
    case "${ACC}" in \#*) continue ;; esac
    slot "${PROC_JOBS}"
    process_one "${ACC}" &
  done < "${ACC_FILE}"
  wait
fi

mv "${FAIL_LOG}.tmp" "${FAIL_LOG}"
echo "failed samples: $(wc -l < "${FAIL_LOG}")"

# --------------------------- collect sample summaries ------------------------
python3 - "${ACC_FILE}" "${BASEDIR}" <<'PYEOF'
import os, sys
import pandas as pd
acc_file, base = sys.argv[1], sys.argv[2]
meta = pd.read_csv(acc_file, sep="\t", header=None,
                   names=["accession", "habitat_group", "pair_id"])
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
    out.to_csv(os.path.join(base, "benchmark_sample_summary.tsv"),
               sep="\t", index=False)
    print("samples collected:", len(out))
    print(out.groupby("habitat_group")["CPRI_mean"].describe())
else:
    print("no completed samples found")
PYEOF
