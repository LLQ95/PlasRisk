#!/usr/bin/env bash
# =============================================================================
# 02_public_metag_benchmark.sh
# Tier A public-metagenome benchmark, resumable + parallel.
#   Stage 1 download : ERR -> ENA direct FASTQ (curl, resume, parallel)
#                      SRR -> prefetch (resumable) + fasterq-dump + pigz
#   Stage 2 process  : fastp -> MEGAHIT -> depth -> geNomad -> plasrisk meta
# Per-sample failures are isolated and logged; re-running resumes from the
# last completed step (markers and existing outputs are reused).
#
# Usage:
#   02_public_metag_benchmark.sh ACCESSIONS WORKDIR THREADS \
#                                [MODE] [DL_JOBS] [PROC_JOBS]
#   MODE      = all | download | process           (default all)
#   DL_JOBS   = parallel downloads                 (default 4)
#   PROC_JOBS = parallel samples in stage 2, each
#               using THREADS; keep PROC_JOBS*THREADS <= cores (default 1)
# accessions.tsv columns (no header): accession,habitat_group,pair_id
# =============================================================================
set -uo pipefail   # no -e: per-sample failures must not abort the whole run

ACC_FILE="${1:?accessions.tsv required}"
BASEDIR="${2:?output dir required}"
THREADS="${3:-16}"
MODE="${4:-all}"
DL_JOBS="${5:-4}"
PROC_JOBS="${6:-1}"

mkdir -p "${BASEDIR}"/{reads,qc,asm,depth,genomad,plasrisk,logs}
FAIL_LOG="${BASEDIR}/failed_samples.tsv"
: > "${FAIL_LOG}.tmp"

have(){ command -v "$1" >/dev/null 2>&1; }

# --------------------------- dependency self-check ---------------------------
miss_dl=0; miss_pr=0
have curl || { echo "MISSING download binary: curl"; miss_dl=1; }
for b in fastp megahit bowtie2 bowtie2-build samtools \
         jgi_summarize_bam_contig_depths genomad plasrisk; do
  have "$b" || { echo "MISSING process binary: $b"; miss_pr=1; }
done
# SRR runs additionally need sra-tools + pigz (ERR use ENA direct FASTQ).
if grep -Eq '^[[:space:]]*SRR' "${ACC_FILE}"; then
  for b in prefetch fasterq-dump pigz; do
    have "$b" || { echo "MISSING download binary (SRR): $b"; miss_dl=1; }
  done
fi
if [ "${MODE}" != "process" ] && [ "${miss_dl}" -ne 0 ]; then
  echo "install missing download tools, then re-run"; exit 1
fi
if [ "${MODE}" != "download" ]; then
  [ "${miss_pr}" -eq 0 ] || { echo "install missing process tools, then re-run"; exit 1; }
  : "${GENOMAD_DB:?set GENOMAD_DB to the genomad_db folder}"
fi

# Wait until fewer than $1 background jobs are running.
slot(){ while [ "$(jobs -rp | wc -l)" -ge "$1" ]; do wait -n || true; done; }

# ENA FASTQ directory for an ERR accession (handles the 000-999 sub-buckets
# used once the numeric part exceeds six digits).
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

# ------------------------------- Stage 1 download ----------------------------
dl_one () {
  local ACC="$1" RD="$2"
  local LOG="${BASEDIR}/logs/${ACC}.dl.log"
  mkdir -p "${RD}"
  if [ -f "${RD}/dl_done" ]; then echo "[dl] ${ACC} cached"; return 0; fi
  echo "[dl] ${ACC} start" > "${LOG}"
  if [[ "${ACC}" == ERR* ]]; then
    local d host
    d="$(ena_dir "${ACC}")"
    host="https://ftp.sra.ebi.ac.uk"
    # -C - resumes a partial transfer; --retry handles transient failures.
    if ! curl -L -C - --retry 5 --retry-delay 5 --fail \
          -o "${RD}/${ACC}_1.fastq.gz" "${host}${d}_1.fastq.gz" >>"${LOG}" 2>&1; then
      echo "[dl] ${ACC} FAILED (_1)"; return 1
    fi
    curl -L -C - --retry 5 --retry-delay 5 --fail \
      -o "${RD}/${ACC}_2.fastq.gz" "${host}${d}_2.fastq.gz" >>"${LOG}" 2>&1 \
      || echo "[dl] ${ACC} no _2 (single-end?)" >>"${LOG}"
  else
    if ! prefetch "${ACC}" -O "${RD}" >>"${LOG}" 2>&1; then
      echo "[dl] ${ACC} FAILED (prefetch)"; return 1
    fi
    local sra
    sra="$(find "${RD}" -name '*.sra' 2>/dev/null | head -1)"
    if ! fasterq-dump "${sra:-$ACC}" -O "${RD}" -t "${RD}/tmp" \
          -e "${THREADS}" --progress >>"${LOG}" 2>&1; then
      echo "[dl] ${ACC} FAILED (fasterq)"; return 1
    fi
    pigz -p "${THREADS}" "${RD}"/*.fastq
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
