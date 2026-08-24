# PlasRisk

**Multi-model risk assessment for bacterial plasmids from FASTA sequences**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)
[![Version](https://img.shields.io/badge/version-1.1.0-green.svg)](https://github.com/LLQ95/PlasRisk)

PlasRisk computes composite risk scores for bacterial plasmids. It accepts
plasmid FASTA sequences, automatically annotates them using
[abricate](https://github.com/tseemann/abricate), and supports **two scoring
models**:

1. **PlasRisk** (default) — a 10-dimension data-driven weighted model with a
   5-dimension lite option, optimized on 792,964 PIPdb plasmid sequence
   clusters (PSCs).
2. **PIPdb** — the original 8-item ordinal scoring system from the PIPdb paper
   (Zhu et al., *Nucleic Acids Res.*, 2025), provided for reproducibility and
   direct comparison.

---

## Scoring models

### PlasRisk (10-dimension weighted model)

```
Full model (10-dim):
S = 0.245*S_ARG + 0.110*S_VF + 0.204*S_MOB + 0.028*S_HOST
  + 0.003*S_REP + 0.181*S_SIZE + 0.211*S_BM
  + 0.002*S_GEO + 0.002*S_HAB + 0.015*S_GROW

Lite model (5-dim core, equivalent AUC):
S = 0.258*S_ARG + 0.115*S_VF + 0.215*S_MOB + 0.190*S_SIZE + 0.222*S_BM
```

Weights were derived by data-driven consensus (Random Forest mean decrease in
Gini, LASSO, and grid-search optimization) on 792,964 PIPdb PSCs; sum = 1.0.

Risk grades: **A** (Very High, S ≥ 0.60), **B** (High, ≥ 0.45),
**C** (Moderate, ≥ 0.30), **D** (Low, ≥ 0.15), **E** (Minimal, < 0.15).

### PIPdb (original 8-item ordinal model)

Reproduces the scoring system from PIPdb Table 1:

```
Combined risk index = Round(
    (Pathogenic_phylum + Pathogenic_species + Habitats + ARGs
     + VFGs + 2 × WHO_ARGs + ISs + Annual_average_growth_rate) / 8 + 0.6
)
```

Each item is binned into an ordinal score of 1–5; WHO-priority ARGs are
double-weighted. The combined index ranges from 1 (Minimal) to 5 (Very High).
Select with `--model pipdb` or `get_scorer('pipdb')`.

| Item | Score bins (1 → 5) |
|------|--------------------|
| Pathogenic phyla | 1, 2, 3, 4, ≥5 |
| Pathogenic species | [1,2), [2,4), [4,6), [6,8), ≥8 |
| Habitats | [1,2), [2,4), [4,6), [6,8), ≥8 |
| ARGs | 0, 1, [2,5), [5,10), ≥10 |
| VFGs | 0, 1, [2,4), [4,6), ≥6 |
| WHO ARGs (×2) | 0, 1, 2, ≥3 |
| Insertion sequences | [1,2), [2,5), [5,15), [15,30), ≥30 |
| Annual growth rate | [0,0.01), [0.01,0.05), [0.05,0.1), [0.1,0.2), ≥0.2 |

---

## Installation

### Option 1: conda (recommended)

```bash
conda create -n plasrisk -c bioconda -c conda-forge plasrisk
conda activate plasrisk
```

### Option 2: pip + manual abricate

```bash
pip install plasrisk
conda install -c bioconda abricate
```

### Option 3: from source

```bash
git clone https://github.com/LLQ95/PlasRisk.git
cd PlasRisk
pip install .
conda install -c bioconda abricate blast
```

Requires **Python 3.9+**.

---

## Quick start

```bash
# Score a single plasmid (default: PlasRisk 10-dim full model)
plasrisk plasmid.fasta

# Use the original PIPdb ordinal model
plasrisk --model pipdb plasmid.fasta

# Lite mode: 5-dimension core (ARG + VF + MOB + SIZE + BM)
plasrisk --mode lite *.fasta

# Score all FASTA files in a directory
plasrisk /path/to/plasmids/

# Specify output directory
plasrisk -o results *.fasta

# Use specific abricate databases
plasrisk --db card,vfdb,plasmidfinder,bacmet,isfinder plasmid.fasta

# Sequence-only mode (no abricate needed; uses replicon empirical priors)
plasrisk --no-abricate contigs.fasta

# JSON output
plasrisk --json -o results plasmid.fasta
```

### Model comparison

| | PlasRisk (full) | PlasRisk (lite) | PIPdb (ordinal) |
|---|---|---|---|
| Type | Continuous weighted | Continuous weighted | Ordinal bins |
| Dimensions | 10 | 5 | 8 items |
| Score range | [0, 1] | [0, 1] | 1–5 (integer) |
| Weights | Data-driven consensus | Renormalized subset | Equal (+WHO ×2) |
| WHO AWaRe weighting | Yes | Yes | WHO count only |
| Housekeeping gene exclusion | Yes | Yes | No |
| Replicon empirical fallback | Yes | Yes | N/A |
| Required annotations | ARG + VF + mobility + replicon + BacMet + metadata | ARG + VF + mobility + length + BacMet | ARG + VF + IS + metadata |
| Use case | Comprehensive One Health surveillance | Rapid screening | Reproducing PIPdb results |

### Full vs. Lite mode

Dimensionality analysis (all-subsets evaluation of 1,023 subsets with 5-fold
CV) showed that performance plateaus at k=5: adding S_VF as the 5th dimension
raises mean CV AUC from 0.918 to 0.920, and the remaining 5 context dimensions
contribute <0.1% additional AUC. The 10-dim full model is retained for
comprehensive surveillance because it is Pareto-optimal (non-dominated across
all 4 outcomes) and provides epidemiological context.

### Output files

| File | Description |
|------|-------------|
| `plasrisk_results.tsv` | Per-sequence scores: all components, S_total, S_norm, grade, gene lists |
| `plasrisk_summary.tsv` | Per-file summary: counts, grade distribution, mean/max scores |
| `plasrisk_results.json` | JSON format (with `--json`) |

When using `--model pipdb`, the output additionally includes all 8 ordinal
sub-scores (`score_phylum`, `score_species`, `score_habitats`, `score_args`,
`score_vfgs`, `score_who_args`, `score_iss`, `score_growth`) and the
`combined_risk_index` (1–5).

---

## Python API

```python
from plasrisk import get_scorer, PlasmidFeatures, annotate_fasta, load_replicon_lookup

lookup = load_replicon_lookup()
result = annotate_fasta("plasmid.fasta", lookup=lookup)

# PlasRisk 10-dim model (default)
scorer = get_scorer("plasrisk")
df = scorer.score_dataframe(result.features)

# PlasRisk lite 5-dim model
scorer_lite = get_scorer("plasrisk", mode="lite")
df_lite = scorer_lite.score_dataframe(result.features)

# Original PIPdb ordinal model
pipdb = get_scorer("pipdb")
df_pipdb = pipdb.score_dataframe(result.features)

# Direct class access also works
from plasrisk import PlasRiskScorer, PIPdbScorer
```

---

## The risk dimensions (PlasRisk)

| Component | Weight (full/lite) | What it measures |
|-----------|--------------------|------------------|
| **S_ARG** | 0.245 / 0.258 | ARG count with WHO AWaRe hazard weighting, high-risk genes (mcr, NDM, KPC, CTX-M, tetX), last-resort multipliers |
| **S_BM** | 0.211 / 0.222 | Biocide/metal resistance (mer, qac, ars/cop/sil) — co-selection potential; CARD fallback when BacMet unavailable |
| **S_MOB** | 0.204 / 0.215 | T4CP, relaxase, oriT, auxiliary transfer proteins, integrons, IS density |
| **S_SIZE** | 0.181 / 0.190 | Plasmid length (cargo capacity), sigmoid midpoint 30 kb |
| **S_VF** | 0.110 / 0.115 | VF count, exotoxins, secretion systems (T3SS/T4SS) |
| **S_HOST** | 0.028 / — | Host range breadth |
| **S_GROW** | 0.015 / — | Annual growth rate of the replicon |
| **S_REP** | 0.003 / — | Replicon backbone risk prior |
| **S_HAB** | 0.002 / — | Habitat breadth (human/animal/environment) |
| **S_GEO** | 0.002 / — | Geographic spread |

### Robustness features (v1.1.0)

- **Replicon empirical fallback**: When sequence-derived dimensions (S_MOB,
  S_HOST, S_GEO, S_HAB, S_GROW, S_REP) cannot be computed — e.g., abricate
  does not detect mobility genes or metadata is unavailable — replicon-specific
  median values from 792,964 PIPdb PSCs are used as priors. Multi-replicon
  plasmids (e.g., IncFII;IncFIA;IncR) take the maximum prior across replicons.
- **Housekeeping gene exclusion**: Chromosomal efflux pumps and porins
  (acrAB, tolC, mexAB-oprM, etc.) annotated by CARD are excluded from S_ARG
  to avoid inflating scores with non-transferred determinants.
- **Sigmoid S_SIZE**: Logistic function on log10(length) centered at 30 kb
  replaces the linear cap, better separating small mobilizable plasmids from
  large conjugative ones.
- **Additive S_MOB**: Mobility class base + integron bonus + IS density bonus
  (up to +0.30), with ISfinder/Tn database support.
- **CARD biocide/metal fallback**: When BacMet database is unavailable,
  S_BM detects biocide/metal genes from CARD annotations.

---

## Model validation

### Internal validation (792,964 PIPdb PSCs)

| Outcome | AUC (PlasRisk) | AUC (PIPdb) |
|---------|---------------|-------------|
| High-risk ARG carriage | 0.955 | 0.966 |
| MDR–VF fusion | 0.965 | 0.940 |
| Conjugation potential | 0.856 | — |
| Biocide/metal resistance | 0.903 | — |
| **Mean (4 outcomes)** | **0.920** | — |

- **Dimensionality analysis**: 5-dim lite achieves equivalent mean AUC to
  10-dim full (0.920 vs. 0.920); no overfitting (train-test gap < 0.003,
  bootstrap optimism < 0.005).
- **LORO-CV**: mean AUC = 0.962 across 40 replicons.
- **Weight sensitivity** (100 iterations, ±30%): Spearman ρ = 0.994, top-10
  overlap = 9.5/10.
- **Natural-prevalence holdout** (157,046 PSCs, 3.0% high-risk): ROC-AUC =
  0.968, PR-AUC = 0.413.
- **Temporal split validation**: train on pre-2020 PSCs, test on 2020+ PSCs;
  stable discrimination with no performance decay.

### External validation

40 independently curated NCBI plasmids (20 high-risk carrying blaNDM, blaKPC,
mcr-1, etc.; 20 low-risk cloning vectors and small cryptic plasmids), all
verified as true plasmid sequences (no chromosomal contamination). Performance
is reported both with and without sequence-derived annotations; dimensions that
cannot be annotated from a standalone FASTA are imputed from replicon empirical
medians.

---

## Running tests

```bash
python -m pytest tests/ -v          # PlasRisk core tests
python test_pipdb_model.py          # PIPdb model tests (8/8)
```

## Changelog

### v1.1.0

- Added `PIPdbScorer` class implementing the original PIPdb 8-item ordinal
  model (Table 1 formula), selectable via `--model pipdb` or
  `get_scorer('pipdb')`.
- Added `get_scorer()` factory function for model selection.
- Fixed S_MOB/S_HOST/S_GEO/S_HAB/S_GROW fallback to replicon empirical medians
  when abricate does not detect genes.
- Added multi-replicon support in prior lookup (max across replicons).
- Added housekeeping chromosomal gene exclusion from S_ARG (acrAB, tolC, etc.).
- Changed S_SIZE to sigmoid (logistic, midpoint 30 kb).
- Changed S_MOB to additive formula (class base + integron + IS density).
- Added CARD-based biocide/metal fallback when BacMet is unavailable.
- Added ISfinder and Tn transposon database support.
- Added `n_pathogenic_phylum`, `n_pathogenic_species` fields and `n_who_arg`
  property to `PlasmidFeatures`.
- Fixed data.table compatibility issues in external validation pipeline.

### v1.0.0

- Initial release: 10-dimension weighted model with 5-dim lite option.

## Citation

> [Authors]. PlasRisk: a multi-dimension weighted risk assessment framework for
> bacterial plasmids. *iMeta*, 2025. doi: [to be added]

Based on data from:
> Zhu Q, Chen Q, Lu X, et al. PIPdb: a comprehensive plasmid sequence resource
> for tracking the horizontal transfer of pathogenic factors and antimicrobial
> resistance genes. *Nucleic Acids Research*, 2025, 53(D1):D169-D178.
> doi:10.1093/nar/gkae952

## License

MIT License - see [LICENSE](LICENSE).
