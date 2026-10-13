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
DL_JOBS="${5:-4}"
PROC_JOBS="${6:-1}"
CONNS="${7:-4}"

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

# Wait until fewer than $1 background jobs are running. Polling is used
# instead of `wait -n` for compatibility with bash < 4.3 (no wait -n).
slot(){ while [ "$(jobs -rp | wc -l)" -ge "$1" ]; do sleep 3; done; }

# Fetch one ENA FASTQ member.
#   ena_fetch <https_url> <aspera_remote> <out_path>
ena_fetch () {
  local HTTPS_URL="$1" ASP_REMOTE="$2" OUT="$3"
  case "${DL_BACKEND}" in
    aspera)
      # -k1 enables resume; download into the directory (keeps file name).
      local keyargs=()
      [ -n "${ASPERA_KEY}" ] && keyargs=(-i "${ASPERA_KEY}")
      ascp -QT -k1 "${keyargs[@]}" \
        -l "${ASPERA_RATE:-300m}" -P 33001 \
        "${ASP_REMOTE}" "$(dirname "${OUT}")/" ;;
    aria2c)
      # -c resume, -x/-s multi-connection, -m 0 keep retrying, falloc fast.
      aria2c -c -x "${CONNS}" -s "${CONNS}" -k 1M \
        --file-allocation=falloc --retry-wait=5 -m 0 \
        -d "$(dirname "${OUT}")" -o "$(basename "${OUT}")" "${HTTPS_URL}" ;;
    *)
      curl -L -C - --retry 5 --retry-delay 5 --fail -o "${OUT}" "${HTTPS_URL}" ;;
  esac
}

# Query the ENA portal filereport API and return the raw TSV only when it is a
# genuine table containing this run. Under high concurrency the EBI edge can
# answer with an HTML/JS challenge page instead of TSV; that is a transient
# block, so it is detected and retried with backoff rather than being parsed.
ena_api_tsv () {
  local ACC="$1" url out t
  url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}"
  url="${url}&result=read_run&fields=run_accession,fastq_ftp,fastq_aspera,fastq_md5"
  url="${url}&format=tsv&download=false"
  for ((t=1; t<=6; t++)); do
    out="$(curl -L -s --retry 3 --retry-delay 3 "${url}" || true)"
    if printf '%s' "${out}" | head -1 | grep -q 'run_accession' \
       && printf '%s' "${out}" | grep -q "^${ACC}[[:space:]]"; then
      printf '%s' "${out}"; return 0
    fi
    sleep $((t*5))
  done
  return 1
}

# True when the file is a valid gzip (a soft-blocked HTTP 200 returns an HTML
# body that is saved under the .fastq.gz name; gzip -t rejects it).
is_gzip(){ gzip -t "$1" 2>/dev/null; }

# ------------------------------- Stage 1 download ----------------------------
dl_one () {
  local ACC="$1" RD="$2"
  local LOG="${BASEDIR}/logs/${ACC}.dl.log"
  mkdir -p "${RD}"
  if [ -f "${RD}/dl_done" ]; then echo "[dl] ${ACC} cached"; return 0; fi
  echo "[dl] ${ACC} start (${DL_BACKEND})" > "${LOG}"

  # Authoritative paths come from ENA's filereport API (this also confirms the
  # run is present, so no separate range-probe against ftp.sra.ebi.ac.uk).
  local rep ftp_col asp_col
  rep="$(ena_api_tsv "${ACC}" | tail -n +2 || true)"
  ftp_col="$(printf '%s' "${rep}" | cut -f2)"
  asp_col="$(printf '%s' "${rep}" | cut -f3)"

  if [ -n "${ftp_col}" ]; then
    local oldifs="${IFS}"; IFS=';'
    # shellcheck disable=SC2206
    local ftps=(${ftp_col}); # shellcheck disable=SC2206
    local asps=(${asp_col}); IFS="${oldifs}"
    local i https_url asp_remote out base attempt ok
    for ((i=0; i<${#ftps[@]}; i++)); do
      base="$(basename "${ftps[$i]}")"
      https_url="https://${ftps[$i]}"
      if [ -n "${asps[$i]:-}" ]; then
        asp_remote="fasp@${asps[$i]}"
      else
        # derive the fasp remote from the FTP path (strip the FTP host).
        asp_remote="fasp@fasp.sra.ebi.ac.uk:${ftps[$i]#ftp.sra.ebi.ac.uk}"
      fi
      out="${RD}/${base}"
      if [ -s "${out}" ]; then
        if is_gzip "${out}"; then
          echo "[dl] ${ACC} ${base} exists" >>"${LOG}"; continue
        fi
        rm -f "${out}"   # truncated or an HTML block page saved under the name
      fi
      ok=0
      for attempt in 1 2 3 4; do
        if ena_fetch "${https_url}" "${asp_remote}" "${out}" >>"${LOG}" 2>&1 \
           && [ -s "${out}" ] && is_gzip "${out}"; then ok=1; break; fi
        rm -f "${out}"; sleep $((attempt*5))
      done
      if [ "${ok}" != 1 ]; then
        echo "[dl] ${ACC} FAILED (${base}; blocked/non-gzip after retries)"
        return 1
      fi
    done
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
    echo "[dl] ${ACC} FAILED (not confirmed on ENA after retries; re-run later)"
    return 1
  fi
  if ls "${RD}"/*.fastq.gz >/dev/null 2>&1; then
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
  local RD="${BASEDIR}/reads/${ACC}"
  local OUT="${BASEDIR}/plasrisk/${ACC}"
  # A finished sample is defined by its plasrisk summary, not by a loose marker.
  if [ -s "${OUT}/meta_sample_summary.tsv" ]; then
    echo "[proc] ${ACC} cached"; return 0
  fi
  echo "[proc] ${ACC} start"
  local R1 R2
  R1=$(ls "${RD}"/*_1.fastq.gz 2>/dev/null | head -1 || true)
  R2=$(ls "${RD}"/*_2.fastq.gz 2>/dev/null | head -1 || true)
  if [ -z "${R1}" ] || [ -z "${R2}" ]; then
    echo -e "${ACC}\tno paired reads" >> "${FAIL_LOG}.tmp"
    echo "[proc] ${ACC} skipped (no paired reads)"; return 0
  fi
  # An HTML block page or a truncated transfer can be saved under the .fastq.gz
  # name; fastp would fail instantly on it. Reject before any heavy work.
  if ! is_gzip "${R1}" || ! is_gzip "${R2}"; then
    echo -e "${ACC}\tcorrupt/non-gzip reads; re-run download" >> "${FAIL_LOG}.tmp"
    rm -f "${RD}/dl_done"
    echo "[proc] ${ACC} skipped (corrupt reads; re-download)"; return 0
  fi
  local ASM="${BASEDIR}/asm/${ACC}/final.contigs.fa"
  # Heavy body in an isolated subshell; each step logs a labelled marker so the
  # failing stage is obvious in logs/<ACC>.proc.log.
  (
    set -e
    echo ">> step 1/6 fastp"
    if [ ! -s "${BASEDIR}/qc/${ACC}_2.fq.gz" ]; then
      fastp -i "${R1}" -I "${R2}" \
        -o "${BASEDIR}/qc/${ACC}_1.fq.gz" -O "${BASEDIR}/qc/${ACC}_2.fq.gz" \
        -h "${BASEDIR}/qc/${ACC}.html" -j "${BASEDIR}/qc/${ACC}.json" \
        -w "${THREADS}"
    fi

    echo ">> step 2/6 megahit"
    if [ ! -f "${ASM}" ]; then
      [ -d "${BASEDIR}/asm/${ACC}" ] && rm -rf "${BASEDIR}/asm/${ACC}"
      megahit -1 "${BASEDIR}/qc/${ACC}_1.fq.gz" \
              -2 "${BASEDIR}/qc/${ACC}_2.fq.gz" \
              -o "${BASEDIR}/asm/${ACC}" -t "${THREADS}" --presets meta-large
    fi

    echo ">> step 3/6 depth (bowtie2/samtools)"
    if [ ! -s "${BASEDIR}/depth/${ACC}_depth.tsv" ]; then
      [ ! -f "${BASEDIR}/depth/${ACC}_idx.1.bt2" ] && \
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
      rm -f "${BASEDIR}/depth/${ACC}.sam"
    fi

    echo ">> step 4/6 geNomad"
    if ! find "${BASEDIR}/genomad/${ACC}" \
          -name '*_plasmid_summary.tsv' 2>/dev/null | grep -q .; then
      genomad end-to-end --cleanup --threads "${THREADS}" \
        "${ASM}" "${BASEDIR}/genomad/${ACC}" "${GENOMAD_DB}"
    fi

    echo ">> step 5/6 plasrisk meta"
    plasrisk meta "${ASM}" \
      --genomad-bin "${BASEDIR}/genomad/${ACC}" \
      --coverage "${BASEDIR}/depth/${ACC}_depth.tsv" \
      --coverage-column "${ACC}" \
      -o "${OUT}" \
      --min-contig 1000 --complete-length 10000 --threads "${THREADS}"

    echo ">> done"
  ) >> "${LOG}" 2>&1

  if [ $? -eq 0 ] && [ -s "${OUT}/meta_sample_summary.tsv" ]; then
    echo "[proc] ${ACC} DONE"
  else
    echo -e "${ACC}\tprocessing error (see logs/${ACC}.proc.log)" >> "${FAIL_LOG}.tmp"
    echo "[proc] ${ACC} FAILED -- last log lines:"
    tail -n 20 "${LOG}" 2>/dev/null | sed 's/^/    /'
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
