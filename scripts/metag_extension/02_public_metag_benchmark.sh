#!/usr/bin/env bash
# =============================================================================
# 02_public_metag_benchmark.sh
# Tier A: per-sample public metagenome pipeline
#   SRA download -> QC -> assembly -> depth -> geNomad -> plasrisk meta
#
# Usage:
#   ./02_public_metag_benchmark.sh accessions.tsv WORKDIR THREADS
# accessions.tsv columns (no header): accession,habitat_group,pair_id
#   pair_id is "-" for unpaired samples, identical string for matched pairs
#   (e.g. pre/post antibiotic, influent/effluent).
# =============================================================================
set -euo pipefail

ACC_FILE="${1:?accessions.tsv required}"
BASEDIR="${2:?output dir required}"
THREADS="${3:-16}"

mkdir -p "${BASEDIR}"/{reads,qc,asm,depth,genomad,plasrisk,logs}

run_one () {
  local ACC="$1" HAB="$2" PAIR="$3"
  local LOG="${BASEDIR}/logs/${ACC}.log"
  echo "[${ACC}] habitat=${HAB} pair=${PAIR}" | tee "${LOG}"
  local RD="${BASEDIR}/reads/${ACC}"
  mkdir -p "${RD}"

  # 1) download (fasterq-dump keeps temp inside the output dir)
  if [ ! -f "${RD}/done" ]; then
    fasterq-dump "${ACC}" -O "${RD}" -t "${RD}/tmp" -e "${THREADS}" \
      --progress >> "${LOG}" 2>&1
    pigz -p "${THREADS}" "${RD}"/*.fastq
    touch "${RD}/done"
  fi

  # 2) QC
  local R1 R2
  R1=$(ls "${RD}"/*_1.fastq.gz 2>/dev/null | head -1 || true)
  R2=$(ls "${RD}"/*_2.fastq.gz 2>/dev/null | head -1 || true)
  if [ -z "${R1}" ]; then
    echo "[${ACC}] single-end or download failure; skipped" | tee -a "${LOG}"
    return 0
  fi
  fastp -i "${R1}" -I "${R2}" \
    -o "${BASEDIR}/qc/${ACC}_1.fq.gz" \
    -O "${BASEDIR}/qc/${ACC}_2.fq.gz" \
    -h "${BASEDIR}/qc/${ACC}.html" -j "${BASEDIR}/qc/${ACC}.json" \
    -w "${THREADS}" >> "${LOG}" 2>&1

  # 3) assembly (MEGAHIT main analysis)
  if [ ! -f "${BASEDIR}/asm/${ACC}/final.contigs.fa" ]; then
    megahit -1 "${BASEDIR}/qc/${ACC}_1.fq.gz" \
            -2 "${BASEDIR}/qc/${ACC}_2.fq.gz" \
            -o "${BASEDIR}/asm/${ACC}" -t "${THREADS}" \
            --presets meta-large >> "${LOG}" 2>&1
  fi
  local ASM="${BASEDIR}/asm/${ACC}/final.contigs.fa"

  # 4) depth
  bowtie2-build -q "${ASM}" "${BASEDIR}/depth/${ACC}_idx" >> "${LOG}" 2>&1
  bowtie2 -x "${BASEDIR}/depth/${ACC}_idx" \
    -1 "${BASEDIR}/qc/${ACC}_1.fq.gz" \
    -2 "${BASEDIR}/qc/${ACC}_2.fq.gz" \
    -p "${THREADS}" --very-fast -S "${BASEDIR}/depth/${ACC}.sam" 2>> "${LOG}"
  samtools view -bS "${BASEDIR}/depth/${ACC}.sam" \
    | samtools sort -@ "${THREADS}" -o "${BASEDIR}/depth/${ACC}.bam"
  samtools index "${BASEDIR}/depth/${ACC}.bam"
  jgi_summarize_bam_contig_depths \
    --outputDepth "${BASEDIR}/depth/${ACC}_depth.tsv" \
    "${BASEDIR}/depth/${ACC}.bam"

  # 5) geNomad plasmids (positional order is INPUT OUTPUT_DIR DATABASE;
  #    --clean-up is a boolean flag).
  genomad classify --clean-up --threads "${THREADS}" \
    "${ASM}" "${BASEDIR}/genomad/${ACC}" "${GENOMAD_DB}" >> "${LOG}" 2>&1

  # 6) PlasRisk three layers
  plasrisk meta "${ASM}" \
    --genomad-bin "${BASEDIR}/genomad/${ACC}" \
    --coverage "${BASEDIR}/depth/${ACC}_depth.tsv" \
    --coverage-column "${ACC}" \
    -o "${BASEDIR}/plasrisk/${ACC}" \
    --min-contig 1000 --complete-length 10000 \
    --threads "${THREADS}" >> "${LOG}" 2>&1

  echo "[${ACC}] DONE" | tee -a "${LOG}"
}

while IFS=$'\t' read -r ACC HAB PAIR; do
  [ -z "${ACC}" ] && continue
  case "${ACC}" in \#*) continue ;; esac
  run_one "${ACC}" "${HAB}" "${PAIR:-"-"}"
done < "${ACC_FILE}"

# -----------------------------------------------------------------------------
# Collect sample summaries into one table with habitat / pair metadata
# -----------------------------------------------------------------------------
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
