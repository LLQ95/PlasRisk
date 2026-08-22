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

## Reproducibility

- Random seed: `set.seed(42)` in all R scripts
- R version: 4.3.x / 4.4.x
- Complete conda environment: `conda env export -n pipdb_risk > environment.yml`
- The PlasRisk Python package is available at https://github.com/LLQ95/PlasRisk
