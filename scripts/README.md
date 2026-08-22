# PlasRisk Analysis Scripts

This directory contains all R, Python, and Shell scripts used to build and validate the PlasRisk ten-dimension plasmid risk assessment framework, as described in the manuscript.

## Pipeline Overview

```
PIPdb database (792,964 PSCs)
        │
        ▼
  pipdb_01_parse.py        ── Parse & integrate PIPdb metadata tables
        │
        ▼
  pipdb_02_risk_score.py   ── Initial 8-dimension risk scoring (Python prototype)
        │
        ▼
  pipdb_13_descriptive_epidemiology.R  ── Descriptive epidemiology (Figs 1–9)
        │
        ▼
  pipdb_14_network_risk_hotspot.R      ── Networks, risk quartiles, spatial hotspots (Figs 16–20)
        │
        ▼
  pipdb_15_risk_weight_table.R         ── Per-replicon risk table, RF validation (Figs 21–23)
        │
        ▼
  pipdb_16_risk_10dimensions.R         ── Add S_BM; 9-dim vs 10-dim comparison (Figs 25–29)
        │
        ▼
  pipdb_17_weight_optimization.R       ── Data-driven weight determination (Fig 24)
        │
        ▼
  pipdb_18_external_validation_benchmark.R  ── External NCBI validation + benchmark (Fig 25)
  pipdb_18_imbalanced_validation.R     ── Natural-prevalence holdout, calibration, DCA (Figs 27–28)
        │
        ▼
  pipdb_19_case_study.R                ── Case studies: pNDM-1, pHNSHP45, ColE1 (Fig 26)
```

## Script Descriptions

### Data Preparation

| Script | Language | Description |
|--------|----------|-------------|
| `00_setup_env.sh` | Bash | Create conda environments (`pipdb_risk`, `pipdb_phylo`) with all R/Python dependencies |
| `pipdb_01_parse.py` | Python | Parse PIPdb TSV tables, integrate integron annotations, classify mobility, sanitize years; outputs `psc_master.tsv` |
| `pipdb_02_risk_score.py` | Python | Initial 8-dimension composite risk score with AWaRe-weighted ARG hazard; outputs `psc_risk_scores.tsv` |
| `pipdb_03_evolution.py` | Python | Prepare tip-dated datasets for per-replicon Bayesian molecular clock dating |
| `pipdb_05_phylogeny.py` | Python | Phylogenetic analysis helpers (Parsnp/IQ-TREE integration) |
| `pipdb_06_beast_setup.py` | Python | Generate BEAST2 XML files for Bayesian dating |
| `pipdb_07_beast_dating.sh` | Bash | Run BEAST2 dating analyses |
| `pipdb_10_extract_and_annotate.sh` | Bash | Extract and annotate plasmid sequences with abricate |
| `pipdb_11_snp_phylogeny.sh` | Bash | SNP-based phylogeny with Parsnp/snp-sites/IQ-TREE |
| `pipdb_12_temporal_signal.R` | R | Temporal signal diagnostics (root-to-tip regression) |

### Core PlasRisk Analysis

| Script | Language | Description |
|--------|----------|-------------|
| `pipdb_13_descriptive_epidemiology.R` | R | Descriptive epidemiology: temporal trends, geographic maps, ARG burden, size evolution, high-risk ARGs, virulence, mobility, host range (Figs 1–15) |
| `pipdb_14_network_risk_hotspot.R` | R | ARG–replicon bipartite network, ARG–VF co-occurrence network, risk quartile validation, Getis-Ord Gi* spatial hotspots (Figs 16–20) |
| `pipdb_15_risk_weight_table.R` | R | Per-replicon 10-dimension risk table, RF feature importance, weight sensitivity, logistic regression (Figs 21–23) |
| `pipdb_16_risk_10dimensions.R` | R | Add S_BM (biocide/metal resistance), 9-dim vs 10-dim comparison, S_BM analysis (Figs 25–29) |
| `pipdb_17_weight_optimization.R` | R | **Data-driven weight determination**: RF-MDG, LASSO, entropy, grid search, LORO-CV; outputs `tab_weight_comparison.csv` with final consensus weights |
| `pipdb_18_external_validation_benchmark.R` | R | External validation on 40 NCBI plasmids + ROC/PR benchmark vs PIPdb ordinal, raw ARG count, single-component models |
| `pipdb_18_imbalanced_validation.R` | R | Natural-prevalence holdout (~12% high-risk), PR-AUC, decile calibration, decision-curve analysis, NCBI random-sample script |
| `pipdb_19_case_study.R` | R | Case studies: pNDM-1 (Grade A), pHNSHP45/mcr-1 (Grade B), ColE1-like (Grade D) |

### Utility

| Script | Language | Description |
|--------|----------|-------------|
| `pipdb_04_figures.R` | R | Initial figure generation (early prototype) |
| `pipdb_08_extended.R` | R | Extended analyses |
| `pipdb_09_ml_validation.py` | Python | Machine learning validation (Python prototype) |
| `rerun_with_new_weights.sh` | Bash | Re-run all weight-dependent analyses in correct order after `pipdb_17` determines final weights |

## Quick Start

```bash
# 1. Set up conda environment
bash 00_setup_env.sh
conda activate pipdb_risk

# 2. Parse PIPdb data (requires PIPdb TSV files in ../data/)
python pipdb_01_parse.py --config ../config.yaml

# 3. Run descriptive epidemiology
Rscript pipdb_13_descriptive_epidemiology.R results

# 4. Run network and risk analyses
Rscript pipdb_14_network_risk_hotspot.R results

# 5. Determine data-driven weights
Rscript pipdb_17_weight_optimization.R results

# 6. Re-run all weight-dependent analyses with final weights
bash rerun_with_new_weights.sh results

# 7. External validation and benchmarking
Rscript pipdb_18_external_validation_benchmark.R results
Rscript pipdb_18_imbalanced_validation.R results
```

## Output Structure

```
results/
├── tables/          # All CSV result tables (tab_*.csv)
├── figures/
│   └── descriptive/ # All PNG/PDF figures (fig*.png/pdf)
└── psc_master.tsv   # Integrated PIPdb master table
```

## Key Result Tables

| Table | Description |
|-------|-------------|
| `tab_weight_comparison.csv` | Expert/RF/LASSO/Entropy/Optimized/Final weights for all 10 dimensions |
| `tab_weight_auc_comparison.csv` | AUC for each weight scheme × 4 outcomes |
| `tab_benchmark_auc.csv` | PlasRisk vs PIPdb ordinal vs single-component AUCs |
| `tab_arg_vf_cooccurrence.csv` | Fisher's exact ARG–VF associations (OR, p_adj) |
| `tab_arg_replicon_network_metrics.csv` | Network modularity, betweenness centrality |
| `tab_risk_qgrade_metrics.csv` | Risk quartile characteristics |
| `tab_replicon_risk_10dimensions.csv` | Per-replicon 10-dimension scores and rankings |
| `tab_case_study_summary.csv` | Case study plasmid scores |
| `tab_SBM_weight_sensitivity.csv` | S_BM weight sensitivity analysis |
| `tab_weight_sensitivity.csv` | ±30% perturbation stability |

## Final Data-Driven Weights

| Dimension | Weight | Interpretation |
|-----------|--------|----------------|
| S_ARG | 0.2448 | ARG burden (incl. WHO high-priority bonus) |
| S_BM | 0.2112 | Biocide/metal resistance (co-selection) |
| S_MOB | 0.2041 | Mobility/conjugation machinery |
| S_SIZE | 0.1808 | Plasmid size (logistic sigmoid) |
| S_VF | 0.1096 | Virulence factor burden |
| S_HOST | 0.0282 | Host range |
| S_GROW | 0.0147 | Epidemic growth rate |
| S_REP | 0.0030 | Replicon type prior |
| S_HAB | 0.0022 | Habitat breadth |
| S_GEO | 0.0015 | Geographic spread |

## R Package Dependencies

Core packages (auto-installed by each script if missing):

```r
install.packages(c("data.table", "ggplot2", "pROC", "patchwork", "RColorBrewer", "scales"))
```

Packages needed only for specific scripts (auto-installed with fallback if unavailable):

| Package | Used by | Fallback |
|---------|---------|----------|
| `randomForest` | pipdb_15, pipdb_16, pipdb_17 | — |
| `glmnet` | pipdb_17 | — |
| `PRROC` | pipdb_18_imbalanced | pROC-based PR-AUC |
| `igraph`, `ggraph` | pipdb_14 | — |
| `scatterpie` | pipdb_13 (world map pies) | built-in pieGrob |
| `sf` | pipdb_13 (map loading) | `maps::map_data` |
| `ggrepel` | pipdb_13, pipdb_16 | `geom_text` |
| `viridis`, `cowplot` | optional | default ggplot scales |

## Reproducibility

- Random seed: `set.seed(42)` in all R scripts
- R version: 4.3.x / 4.4.x
- All session information captured via `sessionInfo()`
- Complete conda environment: `conda env export -n pipdb_risk > environment.yml`
- The PlasRisk Python package (independent of these analysis scripts) is available at https://github.com/LLQ95/PlasRisk
