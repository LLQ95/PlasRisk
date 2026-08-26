# PlasRisk Analysis Scripts

This directory contains all R, Python, and Shell scripts used to build and validate the PlasRisk ten-dimension plasmid risk assessment framework, as described in the manuscript.

## Pipeline Overview

```
PIPdb database (792,964 PSCs)
        |
        v
  pipdb_01_parse.py        -- Parse & integrate PIPdb metadata tables
        |
        v
  pipdb_02_risk_score.py   -- Initial 8-dimension risk scoring (Python prototype)
        |
        v
  pipdb_13_descriptive_epidemiology.R  -- Descriptive epidemiology
        |
        v
  pipdb_14_network_risk_hotspot.R      -- Networks, risk quartiles, spatial hotspots
        |
        v
  pipdb_15_risk_weight_table.R         -- Per-replicon risk table, RF validation
        |
        v
  pipdb_16_risk_10dimensions.R         -- Add S_BM; 9-dim vs 10-dim comparison
        |
        v
  pipdb_17_weight_optimization.R       -- Data-driven weight determination
        |
        v
  pipdb_18_external_validation_benchmark.R  -- External NCBI validation + benchmark
  pipdb_18_imbalanced_validation.R     -- Natural-prevalence holdout, calibration, DCA
        |
        v
  pipdb_19_case_study.R                -- Case studies: pNDM-1, pHNSHP45, ColE1
        |
        v
  fig2_*.R ... fig11_*.R                -- Final manuscript Figures 2-11
```

## Script Descriptions

### Data Preparation

| Script | Language | Description |
|--------|----------|-------------|
| `00_setup_env.sh` | Bash | Create conda environments (`pipdb_risk`, `pipdb_phylo`) with all R/Python dependencies |
| `pipdb_01_parse.py` | Python | Parse PIPdb TSV tables, integrate integron annotations, classify mobility, sanitize years; outputs `psc_master.tsv` |
| `pipdb_02_risk_score.py` | Python | Initial 8-dimension composite risk score with AWaRe-weighted ARG hazard; outputs `psc_risk_scores.tsv` |
| `pipdb_03_evolution.py` | Python | Prepare tip-dated datasets for per-replicon Bayesian molecular clock dating |
| `pipdb_05_phylogeny.py` | Python | Phylogeny analysis helpers (Parsnp/IQ-TREE integration) |
| `pipdb_06_beast_setup.py` | Python | Generate BEAST2 XML files for Bayesian dating |
| `pipdb_07_beast_dating.sh` | Bash | Run BEAST2 dating analyses |
| `pipdb_10_extract_and_annotate.sh` | Bash | Extract and annotate plasmid sequences with abricate |
| `pipdb_11_snp_phylogeny.sh` | Bash | SNP-based phylogeny with Parsnp/snp-sites/IQ-TREE |
| `pipdb_12_temporal_signal.R` | R | Temporal signal diagnostics (root-to-tip regression) |

### Core PlasRisk Analysis

| Script | Language | Description |
|--------|----------|-------------|
| `pipdb_13_descriptive_epidemiology.R` | R | Descriptive epidemiology: temporal trends, geographic maps, ARG burden, size evolution, high-risk ARGs, virulence, mobility, host range |
| `pipdb_14_network_risk_hotspot.R` | R | ARG-replicon bipartite network, ARG-VF co-occurrence network, risk quartile validation, Getis-Ord Gi* spatial hotspots |
| `pipdb_15_risk_weight_table.R` | R | Per-replicon 10-dimension risk table, RF feature importance, weight sensitivity, logistic regression |
| `pipdb_16_risk_10dimensions.R` | R | Add S_BM (biocide/metal resistance), 9-dim vs 10-dim comparison, S_BM analysis |
| `pipdb_17_weight_optimization.R` | R | **Data-driven weight determination**: RF-MDG, LASSO, entropy, grid search, LORO-CV; outputs `tab_weight_comparison.csv` with final consensus weights |
| `pipdb_18_external_validation_benchmark.R` | R | External validation on 367 independent NCBI plasmids + ROC/PR benchmark vs PIPdb ordinal, raw ARG count, single-component models |
| `pipdb_18_imbalanced_validation.R` | R | Natural-prevalence holdout (~3.1% high-risk), PR-AUC, decile calibration, decision-curve analysis |
| `pipdb_19_case_study.R` | R | Case studies: pNDM-1 (Grade A), pHNSHP45/mcr-1 (Grade B), ColE1-like (Grade D) |
| `pipdb_20_dimensionality_analysis.R` | R | All-subsets (1,023) evaluation, forward/backward selection, Pareto frontier, 5-dim lite model |
| `pipdb_21_overfitting_analysis.R` | R | Train-test gap, learning curves, bootstrap optimism, fold stability |
| `pipdb_22_bm_arg_coselection.R` | R | BMG-ARG co-occurrence network, CARD/BacMet cross-annotation sensitivity, MGE physical linkage |

### Manuscript figure scripts (Figures 2-11)

These scripts consume the `tab_*.csv` tables produced by the pipeline above and render the final main-text figures. Each is run as `Rscript <script> results` and writes PNG/PDF files to `results/figures/`.

| Script | Figure | Content |
|--------|--------|---------|
| `fig2_weight_validation.R` | Figure 2 | Consensus weight comparison across RF-MDG, LASSO, and grid search; weight robustness |
| `fig3_bm_arg_network.R` | Figure 3 | ARG family x BMG category co-occurrence bipartite network |
| `fig4_dimensionality.R` | Figure 4 | All-subsets performance curve and 5-dim lite vs 10-dim full comparison |
| `fig5_risk_stratification.R` | Figure 5 | A-E grade stratification, calibration, and decision-curve analysis |
| `fig6_highrisk_arg_replicon.R` | Figure 6 | High-risk ARG carriage across replicons, length-rate landscape, ARG-replicon network |
| `fig7_fusion_plasmids.R` | Figure 7 | MDR-VF fusion plasmid features and per-replicon enrichment |
| `fig8_dimension_epidemiology.R` | Figure 8 | Dimension-score heatmap, temporal trends, dimension correlation matrix, risk landscape |
| `fig9_global_outcomes_maps.R` | Figure 9 | Global country-level prevalence maps for the four outcomes |
| `fig10_grade_pie_map.R` | Figure 10 | World map with grade-distribution pies per country |
| `fig11_conjugative_replicon.R` | Figure 11 | Conjugative rate vs high-risk ARG/BMG/fusion cargo across replicons |

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

# 8. Render final manuscript figures
for f in fig2_*.R fig3_*.R fig4_*.R fig5_*.R fig6_*.R fig7_*.R fig8_*.R fig9_*.R fig10_*.R fig11_*.R; do
  Rscript "$f" results
done
```

## Output Structure

```
results/
|-- tables/          # All CSV result tables (tab_*.csv)
|-- figures/
|   `-- descriptive/ # All PNG/PDF figures (Figure*.png/pdf)
`-- psc_master.tsv   # Integrated PIPdb master table
```

## Key Result Tables

| Table | Description |
|-------|-------------|
| `tab_weight_comparison.csv` | Expert/RF/LASSO/Entropy/Optimized/Final weights for all 10 dimensions |
| `tab_weight_auc_comparison.csv` | AUC for each weight scheme x 4 outcomes |
| `tab_benchmark_auc.csv` | PlasRisk vs PIPdb ordinal vs single-component AUCs |
| `tab_arg_vf_cooccurrence.csv` | Fisher's exact ARG-VF associations (OR, p_adj) |
| `tab_arg_replicon_network_metrics.csv` | Network modularity, betweenness centrality |
| `tab_risk_qgrade_metrics.csv` | Risk quartile characteristics |
| `tab_replicon_risk_10dimensions.csv` | Per-replicon 10-dimension scores and rankings |
| `tab_psc_final_scores.csv` | Per-PSC final scores with four outcome labels (y_highrisk, y_fusion, y_conj, y_bm) |
| `tab_replicon_summary.csv` | Per-replicon n, ARG/BMG/conjugative rates, mean length, mean risk |
| `tab_case_study_summary.csv` | Case study plasmid scores |
| `tab_SBM_weight_sensitivity.csv` | S_BM weight sensitivity analysis |
| `tab_weight_sensitivity.csv` | +/-30% perturbation stability |

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
install.packages(c("data.table", "ggplot2", "pROC", "patchwork", "RColorBrewer", "scales", "ggrepel", "maps", "mapproj"))
```

Packages needed only for specific scripts (auto-installed with fallback if unavailable):

| Package | Used by | Fallback |
|---------|---------|----------|
| `randomForest` | pipdb_15, pipdb_16, pipdb_17 | - |
| `glmnet` | pipdb_17 | - |
| `PRROC` | pipdb_18_imbalanced | pROC-based PR-AUC |
| `igraph`, `ggraph` | pipdb_14, fig3 | - |
| `maps`, `mapproj` | fig9, fig10 (world maps) | required for map figures |
| `sf` | optional modern map loading | `maps::map_data` + `geom_polygon` |
| `scatterpie` | optional map pies | hand-computed trigonometric pie grobs (works on R 4.1.x) |
| `ggrepel` | fig6, fig8, fig11, pipdb_13, pipdb_16 | `geom_text` |
| `viridis`, `cowplot` | optional | default ggplot scales |

Note on compatibility: the `fig2`-`fig11` scripts are written in pure ASCII and
were tested on R 4.1.3. Greek letters in annotations use plotmath expressions
(`rho == 0.01` with `parse = TRUE`) rather than Unicode, because the R 4.1.x
PDF device cannot encode Unicode symbols. `scatterpie` (requires R >= 4.2) and
recent `sf`/`terra` stacks are deliberately avoided.

## Reproducibility

- Random seed: `set.seed(42)` in all R scripts
- R version: figure scripts tested on R 4.1.3; pipeline scripts on R 4.3.x / 4.4.x
- All session information captured via `sessionInfo()`
- Complete conda environment: `conda env export -n pipdb_risk > environment.yml`
- The PlasRisk Python package (independent of these analysis scripts) is available at https://github.com/LLQ95/PlasRisk
