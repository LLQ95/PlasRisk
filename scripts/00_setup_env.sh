#!/usr/bin/env bash
# =============================================================================
# 00_setup_env.sh — create conda env for PIPdb plasmid risk & evolution pipeline
# =============================================================================
set -euo pipefail
ENV=pipdb_risk

# ---- main analysis env (python + R) ----
conda create -y -n $ENV -c conda-forge -c bioconda \
  python=3.11 pandas numpy scipy scikit-learn matplotlib pyyaml biopython \
  shap tqdm \
  r-base=4.4 r-data.table r-ggplot2 r-dplyr r-tidyr r-patchwork r-ggrepel \
  r-viridis r-igraph r-ggraph r-ape r-phytools r-cowplot

# ---- phylogeny & dating env (separate to avoid conflicts) ----
conda create -y -n pipdb_phylo -c conda-forge -c bioconda \
  parsnp harvesttools prokka panaroo mafft iqtree2 \
  beast2=2.7 snp-sites snp-dists

# ---- NCBI download tools ----
conda activate $ENV
conda install -y -c conda-forge -c bioconda entrez-direct ncbi-datasets-cli seqkit

echo "Done. Activate with: conda activate $ENV"
