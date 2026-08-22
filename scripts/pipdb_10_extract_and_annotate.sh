#!/usr/bin/env bash
# pipdb_10_extract_and_annotate.sh — Extract target contigs and annotate with prokka/Panaroo
set -euo pipefail
CONTIG_DIR="${1:?Usage: $0 <contigs_dir> <evolution_dir> [threads]}"
EVO_DIR="${2:?Usage: $0 <contigs_dir> <evolution_dir> [threads]}"
THREADS="${3:-32}"
WORK_DIR="${EVO_DIR}/_extracted"
mkdir -p "$WORK_DIR"

for cmd in seqkit awk sort parallel; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found"; exit 1; }
done

echo "=== Step 1: Collect target seq_ids ==="
ALL_IDS="${WORK_DIR}/all_target_ids.txt"
ID_MAP="${WORK_DIR}/id_to_replicon.tsv"
: > "$ID_MAP"
for td in "${EVO_DIR}"/*/tip_dates.tsv; do
  [ -f "$td" ] || continue
  rep=$(basename "$(dirname "$td")")
  awk -F'\t' -v r="$rep" 'NR>1 && $1!="" {print $1"\t"r}' "$td" >> "$ID_MAP"
done
cut -f1 "$ID_MAP" | sort -u > "$ALL_IDS"
echo "Target contigs: $(wc -l < "$ALL_IDS")"

echo "=== Step 2: Extract from part files ==="
PARTS=("${CONTIG_DIR}"/plasmid_contig.nucleotide.part*.fasta)
ALL_FA="${WORK_DIR}/all_selected.fasta"
if [ -s "$ALL_FA" ]; then echo "Already extracted"; else
  seqkit grep -j "$THREADS" -f "$ALL_IDS" "${PARTS[@]}" -o "$ALL_FA"
fi
echo "Extracted: $(grep -c '^>' "$ALL_FA") sequences"

echo "=== Step 3: Split by replicon ==="
extract_one() {
  local rep="$1"
  local rep_dir="${EVO_DIR}/${rep}"
  mkdir -p "${rep_dir}/fasta"
  awk -F'\t' -v r="$rep" '$2==r {print $1}' "$ID_MAP" | sort -u > "${rep_dir}/_ids.txt"
  seqkit grep -j 2 -f "${rep_dir}/_ids.txt" "$ALL_FA" -o "${rep_dir}/fasta/${rep}.fasta"
}
export -f extract_one; export EVO_DIR ID_MAP ALL_FA
cut -f2 "$ID_MAP" | sort -u | parallel -j "$((THREADS/2))" extract_one {}

echo "=== DONE ==="
