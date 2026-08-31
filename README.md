# PlasRisk: Data-Driven, Multi-Dimensional Weighted Risk Scoring for Bacterial Plasmids

PlasRisk is a Python package that scores plasmid sequences on ten biologically grounded dimensions, then combines them with consensus weights learned three independent ways (random-forest Gini importance, LASSO logistic coefficients, and coordinate-descent grid search). It predicts four distinct, biologically meaningful outcomes: high-risk antimicrobial-resistance-gene (ARG) carriage, ARG–virulence-factor (VF) fusion, conjugative mobility, and biocide/metal-resistance-gene (BMG) carriage. The framework is designed for One-Health surveillance and source attribution across human, animal, food, and environmental isolates, so no single dimension or outcome can substitute for the others.

<p align="center">
  <img src="docs/graphical_abstract.svg" alt="PlasRisk graphical abstract" width="900"><br>
  <em>Graphical abstract — from FASTA input and ten biology-driven sub-scores, through three-method consensus weights and a calibrated composite score with full-10 / lite-5 modes, to multi-outcome validation and deployment. Editable vector source: <code>docs/graphical_abstract.svg</code>.</em>
</p>

---

## Table of Contents

- [Scoring models](#scoring-models)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Outputs](#outputs)
- [Consensus Weights](#consensus-weights)
- [Model Validation](#model-validation)
- [PIPdb Integration](#pipdb-integration)
- [Project Structure](#project-structure)
- [Reproducing the Paper Figures](#reproducing-the-paper-figures)
- [Development](#development)
- [Citation](#citation)
- [License](#license)

---

## Scoring models

Each plasmid sequence cluster (PSC) is scored on ten dimensions normalized to [0, 1]:

| Code | Dimension | Biological meaning |
|------|-----------|-------------------|
| `S_ARG` | ARG hazard | AWaRe-classification-weighted antimicrobial-resistance-gene hazard |
| `S_VF` | VF burden | Virulence-factor count, exotoxins, and effector-delivery systems |
| `S_MOB` | Mobility / MGE plasticity | Conjugative apparatus (T4CP, relaxase, oriT), integrons, IS density |
| `S_HOST` | Host pathogenicity | ESKAPE / WHO-critical host classes and human association |
| `S_REP` | Replicon prior | Risk-tiered Inc-group lookup table |
| `S_SIZE` | Plasmid size | Log-scaled sequence length |
| `S_BM` | Biocide / metal resistance | Metal-resistance operons (mer, ars, cop, sil, qac) |
| `S_GEO` | Geographic spread | Number of countries of isolation |
| `S_HAB` | Habitat breadth | Number of distinct isolation source categories |
| `S_GROW` | Epidemic growth | Annual collection growth rate from PIPdb time series |

Two deployable models are provided:

- **Full model (10 dimensions)** — the complete framework, including epidemiological context dimensions (`S_HOST`, `S_REP`, `S_GEO`, `S_HAB`, `S_GROW`) for One-Health source attribution.
- **Lite model (5 core dimensions)** — a parsimonious model using only `S_ARG`, `S_VF`, `S_MOB`, `S_SIZE`, and `S_BM`, selected from all 1,023 non-empty subsets. The lite model matches the full model within |ΔAUC| ≤ 0.001 on every outcome and reaches 98.2% grade agreement, and needs only standard annotation outputs (ARG, VF, mobility typing, length, biocide/metal genes).

---

## Installation

### From Bioconda (recommended)

```bash
# After adding the bioconda channels (see bioconda/meta.yaml for details)
conda install -c bioconda plasrisk
```

### From source

```bash
git clone https://github.com/LLQ95/PlasRisk.git
cd PlasRisk
pip install .
```

### Dependencies

- Python >= 3.9
- numpy, pandas, scikit-learn (>=1.3), scipy
- biopython

---

## Quick Start

### Score a single FASTA file

```bash
plasrisk run plasmid.fasta --output results.tsv
```

### Batch mode

```bash
plasrisk run *.fasta --output-dir results/ --model full
```

### Use the lite model

```bash
plasrisk run *.fasta --model lite --output lite_results.tsv
```

### Python API

```python
from plasrisk import PlasRiskScorer

scorer = PlasRiskScorer(model="full")
result = scorer.score("plasmid.fasta")

print(result.risk_score)       # continuous score in [0, 1]
print(result.grade)            # calibrated grade A–E
print(result.percentile)       # PIPdb reference percentile
print(result.sub_scores)       # the ten (or five) dimension scores
print(result.outcome_probs)    # probabilities for the four outcomes
```

---

## Outputs

For every input plasmid the tool reports:

- a continuous composite **risk score** in [0, 1];
- a calibrated **grade (A–E)** and a **PIPdb percentile**;
- the per-dimension **sub-scores**, so every contribution is transparent;
- predicted probabilities for the four independent outcomes (high-risk ARG, ARG–VF fusion, conjugative mobility, biocide/metal resistance);
- machine-readable TSV or JSON for downstream integration.

---

## Consensus Weights

Rather than assigning weights by expert judgement, PlasRisk learns them three independent ways and takes an equal-vote consensus:

1. **Random forest** — mean decrease in Gini impurity.
2. **LASSO logistic regression** — absolute sparse coefficients.
3. **Coordinate-descent grid search** — the weight vector maximizing mean cross-validated AUC.

### Full model (normalized)

| Dimension | Weight | Dimension | Weight |
|-----------|--------|-----------|--------|
| S_ARG | 0.245 | S_SIZE | 0.181 |
| S_VF | 0.110 | S_BM | 0.211 |
| S_MOB | 0.204 | S_GEO | 0.002 |
| S_HOST | 0.028 | S_HAB | 0.002 |
| S_REP | 0.003 | S_GROW | 0.015 |

### Lite model (normalized)

| Dimension | Weight |
|-----------|--------|
| S_ARG | 0.258 |
| S_VF | 0.115 |
| S_MOB | 0.215 |
| S_SIZE | 0.190 |
| S_BM | 0.222 |

The data-driven consensus departs from expert intuition in informative ways: biocide/metal resistance (`S_BM`) receives roughly three times the weight an expert panel would assign, while geographic spread, habitat breadth, and replicon prior contribute negligibly after the other dimensions are accounted for.

---

## Model Validation

### Internal validation (locked test set, triple 80/20 split, nested CV)

| Outcome | Full model | Lite model | PIPdb ordinal baseline |
|---------|-----------|------------|------------------------|
| High-risk ARG | 0.958 | 0.958 | 0.966 |
| ARG–VF fusion | 0.964 | 0.964 | 0.940 |
| Conjugative mobility | 0.856 | 0.854 | 0.829 |
| Biocide/metal resistance | 0.903 | 0.907 | 0.662 |
| **Mean AUC** | **0.920** | **0.921** | **0.849** |

### Robustness checks

- **Species-cluster bootstrap** (B = 1,000, holding out all 11,057 species clusters whole): confidence intervals remain narrow and all lower bounds stay well above the baselines.
- **Temporal validation**: training on pre-2020 PSCs (n = 722,350) and testing on post-2020 PSCs (n = 70,614) preserves the ranking, confirming the framework is not specific to one sampling period.
- **Leave-one-replicon-out × 40**: each Inc family is held out in turn; performance is stable across replicon backgrounds.
- **Firth penalized logistic regression** confirms every reported coefficient under rare-event/separation conditions, and variance-inflation factors stay below the collinearity threshold (max VIF ≈ 6.8 for the correlated ecological dimensions; all core dimensions < 1.5).
- **External validation** on 367 independent plasmids annotated with a separate toolchain: AUC 0.998 (95% CI 0.994–1.000), sensitivity 0.994, specificity 0.960.
- **Dimensionality ablation**: the five-dimension lite model is reached at the k = 5 plateau of all 1,023 non-empty subsets, after which adding dimensions yields no meaningful AUC gain.
- **Calibration**: reliability curves and Brier scores show the predicted probabilities remain calibrated after Platt/isotonic correction on the validation split.

---

## PIPdb Integration

PlasRisk is built on [PIPdb](https://www.pipdb.org), the curated plasmid genome database. The package ships with:

- reference percentiles and A–E grade thresholds calibrated on the PIPdb corpus;
- replicon and host lookup tables derived from PIPdb annotations;
- time-series growth statistics computed from PIPdb collection years;
- an update command to refresh reference tables when new PIPdb releases appear.

```bash
plasrisk update-reference --pipdb-version latest
```

---

## Project Structure

```
PlasRisk/
├── plasrisk/                  # Main package
│   ├── __init__.py
│   ├── scorer.py              # Core scoring pipeline
│   ├── dimensions/            # The ten dimension calculators
│   ├── models/                # Full and lite models, weight tables
│   ├── calibration/           # Grade thresholds and percentiles
│   ├── io.py                  # FASTA parsing, TSV/JSON output
│   └── cli.py                 # Command-line interface
├── scripts/                   # Paper analysis and figure scripts
│   ├── fig2_weight_validation.R
│   ├── fig3_arg_bmg_replicon_network.R
│   ├── ...                    # fig4–fig11
│   ├── pipdb_01_parse.py      # PIPdb build pipeline (01–22)
│   └── tables/                # Generated summary tables
├── tests/                     # Unit tests
├── bioconda/                  # Bioconda recipe
├── conda/                     # Environment files
├── docs/                      # Documentation and vector figures
├── pyproject.toml
└── README.md
```

---

## Reproducing the Paper Figures

The `scripts/` directory contains the R and Python scripts used to generate every manuscript figure and table. Figure scripts are numbered to match the paper:

```bash
# Example: regenerate the weight-validation figure
Rscript scripts/fig2_weight_validation.R

# Rebuild derived tables from a local PIPdb dump
bash scripts/00_setup_env.sh
python scripts/pipdb_01_parse.py
# ... continue through pipdb_22
```

See `scripts/README.md` for the full dependency list and run order.

---

## Development

```bash
# Create the development environment
conda env create -f conda/environment.yml
conda activate plasrisk-dev

# Install in editable mode with test dependencies
pip install -e ".[dev]"

# Run the test suite
pytest tests/ -v
```

Contributions are welcome. Please open an issue to discuss a proposed change before submitting a pull request, and add tests for any new behavior.

---

## Citation

If you use PlasRisk in published work, please cite the manuscript:

> Li L., et al. *PlasRisk: a data-driven, multi-dimensional weighted risk scoring framework for bacterial plasmids.* Journal details to be updated.

A BibTeX entry is provided in `CITATION.cff` once the article is assigned its DOI.

---

## License

PlasRisk is released under the MIT License. See [LICENSE](LICENSE) for details.
