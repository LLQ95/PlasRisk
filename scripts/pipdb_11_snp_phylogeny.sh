#!/usr/bin/env bash
# pipdb_11_snp_phylogeny.sh — SNP-based phylogeny with snippy + Gubbins + IQ-TREE2
set -euo pipefail
EVO_DIR="${1:-results/evolution}"
REF_DIR="${2:-references}"
THREADS="${3:-32}"

declare -A REF_MAP=(
    ["IncFII"]="R100_AP000342.fasta"
    ["IncI1"]="R64_AP005147.fasta"
    ["IncX1"]="pOLA52_EU213072.fasta"
    ["IncX3"]="pEC14_35.fasta"
    ["IncN"]="R46_AY046276.fasta"
    ["IncHI2"]="R478_BX664015.fasta"
    ["IncA/C"]="pRMH760_KF976463.fasta"
    ["ColKP3"]="pMCR_1511_KX276659.fasta"
    ["IncR"]="pKPC_CAV1596.fasta"
)

run_snippy_one() {
  local rep="$1" ref="$2"
  local rep_dir="${EVO_DIR}/${rep}"
  local fa="${rep_dir}/fasta/${rep}.fasta"
  local snp_dir="${rep_dir}/tree/snippy"
  [ -s "$fa" ] || { echo "  [skip] $rep: no fasta"; return; }
  [ -f "$ref" ] || { echo "  [skip] $rep: reference not found"; return; }
  mkdir -p "$snp_dir"
  local map_dir="${snp_dir}/maps"; mkdir -p "$map_dir"
  local tmp="${snp_dir}/_split"; mkdir -p "$tmp"
  seqkit split -j 2 -i -O "$tmp" "$fa" >/dev/null 2>&1
  for sfa in "$tmp"/*.fasta; do
    [ -f "$sfa" ] || continue
    local sid=$(basename "$sfa" .fasta | tr '.' '_' | tr '/' '_')
    if [ ! -d "${map_dir}/${sid}" ]; then
      snippy --ctgs "$sfa" --ref "$ref" --outdir "${map_dir}/${sid}" --cpus 1 --force --quiet 2>/dev/null || true
    fi
  done
  rm -rf "$tmp"
  local aln="${snp_dir}/${rep}.full.aln"
  if [ ! -f "$aln" ]; then
    snippy-core --ref "$ref" --prefix "${snp_dir}/${rep}" "${map_dir}"/* 2>/dev/null
  fi
  if [ -f "$aln" ]; then
    snp-sites -c "$aln" > "${snp_dir}/${rep}.snps.fasta"
    local n_sites=$(grep -v '^>' "${snp_dir}/${rep}.snps.fasta" | head -1 | tr -d '\n' | wc -c)
    echo "  $rep: $n_sites core SNP sites"
    if [ "$n_sites" -gt 100 ]; then
      run_gubbins.py --prefix "${snp_dir}/${rep}_gubbins" --threads 4 "${snp_dir}/${rep}.snps.fasta" 2>/dev/null || true
      local tree_aln="${snp_dir}/${rep}.snps.fasta"
      [ -f "${snp_dir}/${rep}_gubbins.filtered_polymorphic_sites.fasta" ] && tree_aln="${snp_dir}/${rep}_gubbins.filtered_polymorphic_sites.fasta"
      iqtree2 -s "$tree_aln" -m GTR+F+R3 -B 1000 -T 4 --prefix "${snp_dir}/${rep}_ml" --quiet 2>/dev/null
      echo "  $rep: ML tree done"
    fi
  fi
}
export -f run_snippy_one; export EVO_DIR

echo "=== SNP-based phylogeny ==="
for rep in "${!REF_MAP[@]}"; do
  ref="${REF_DIR}/${REF_MAP[$rep]}"
  run_snippy_one "$rep" "$ref"
done
echo "DONE."
