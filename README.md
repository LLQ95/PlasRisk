# PlasRisk

**Multi-dimension weighted risk assessment for bacterial plasmids from FASTA sequences**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

PlasRisk computes a composite risk score for bacterial plasmids based on ten
biologically motivated dimensions: antimicrobial resistance gene (ARG) burden,
virulence factors (VFs), mobility/conjugation potential, host range, replicon
type, plasmid size, biocide/metal resistance (BMRG), geographic spread, habitat
breadth, and temporal growth rate. It accepts plasmid FASTA sequences and
automatically annotates them using [abricate](https://github.com/tseemann/abricate).

```
Full model (10-dim):
S = 0.245*S_ARG + 0.110*S_VF + 0.204*S_MOB + 0.028*S_HOST
  + 0.003*S_REP + 0.181*S_SIZE + 0.211*S_BM
  + 0.002*S_GEO + 0.002*S_HAB + 0.015*S_GROW

Lite model (5-dim core, equivalent AUC):
S = 0.258*S_ARG + 0.115*S_VF + 0.215*S_MOB + 0.190*S_SIZE + 0.222*S_BM

Weights derived by data-driven consensus (Random Forest MDG, LASSO, and
grid-search optimization) on 792,964 PIPdb PSCs; sum = 1.0.
```

Risk grades: **A** (Very High, S >= 0.60), **B** (High, >= 0.45),
**C** (Moderate, >= 0.30), **D** (Low, >= 0.15), **E** (Minimal, < 0.15).

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
# Score a single plasmid (10-dim full model)
plasrisk plasmid.fasta

# Lite mode: 5-dimension core (ARG + VF + MOB + SIZE + BM)
plasrisk --mode lite *.fasta

# Score all FASTA files in a directory
plasrisk /path/to/plasmids/

# Specify output directory
plasrisk -o results *.fasta

# Use specific abricate databases
plasrisk --db card,vfdb,plasmidfinder,bacmet plasmid.fasta

# Sequence-only mode (no abricate needed)
plasrisk --no-abricate contigs.fasta

# JSON output
plasrisk --json -o results plasmid.fasta
```

### Full vs. Lite mode

| | Full (10-dim) | Lite (5-dim) |
|---|---|---|
| Dimensions | S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE, S_BM, S_GEO, S_HAB, S_GROW | S_ARG, S_VF, S_MOB, S_SIZE, S_BM |
| Weights | 0.245, 0.110, 0.204, 0.028, 0.003, 0.181, 0.211, 0.002, 0.002, 0.015 | 0.258, 0.115, 0.215, 0.190, 0.222 |
| Mean AUC (4 outcomes) | 0.920 | 0.920 |
| Required annotations | ARG + VF + mobility + replicon + BacMet + metadata | ARG + VF + mobility + length + BacMet |
| Use case | Comprehensive One Health surveillance with epidemiological context | Rapid screening, resource-limited settings |

Dimensionality analysis (all-subsets evaluation of 1,023 subsets with 5-fold CV)
showed that performance plateaus at k=5: adding S_VF as the 5th dimension raises
mean CV AUC from 0.918 to 0.920, and the remaining 5 context dimensions contribute
<0.1% additional AUC. The 10-dim full model is retained for comprehensive
surveillance because it is Pareto-optimal (non-dominated across all 4 outcomes)
and provides epidemiological context.

### Output files

| File | Description |
|------|-------------|
| `plasrisk_results.tsv` | Per-sequence scores: all components, S_total, S_norm, grade, gene lists |
| `plasrisk_summary.tsv` | Per-file summary: counts, grade distribution, mean/max scores |
| `plasrisk_results.json` | JSON format (with `--json`) |

---

## Python API

```python
from plasrisk import PlasRiskScorer, PlasmidFeatures, annotate_fasta, load_replicon_lookup

lookup = load_replicon_lookup()
result = annotate_fasta("plasmid.fasta", lookup=lookup)

# Full 10-dim model
scorer = PlasRiskScorer(replicon_lookup=lookup)
df = scorer.score_dataframe(result.features)

# Lite 5-dim model
scorer_lite = PlasRiskScorer(replicon_lookup=lookup, mode="lite")
df_lite = scorer_lite.score_dataframe(result.features)
```

---

## The risk dimensions

| Component | Weight (full/lite) | What it measures |
|-----------|--------------------|------------------|
| **S_ARG** | 0.245 / 0.258 | ARG count, WHO-priority genes, high-risk genes (mcr, NDM, KPC, CTX-M, tetX) |
| **S_BM** | 0.211 / 0.222 | Biocide/metal resistance (mer, qac, ars/cop/sil) — co-selection potential |
| **S_MOB** | 0.204 / 0.215 | T4CP, relaxase, oriT, auxiliary transfer proteins |
| **S_SIZE** | 0.181 / 0.190 | Plasmid length (cargo capacity), sigmoid midpoint 30 kb |
| **S_VF** | 0.110 / 0.115 | VF count, exotoxins, secretion systems (T3SS/T4SS) |
| **S_HOST** | 0.028 / — | Host range breadth |
| **S_GROW** | 0.015 / — | Annual growth rate of the replicon |
| **S_REP** | 0.003 / — | Replicon backbone risk prior |
| **S_HAB** | 0.002 / — | Habitat breadth (human/animal/environment) |
| **S_GEO** | 0.002 / — | Geographic spread |

---

## Model validation

- **792,964 PSCs** from PIPdb (Zhu et al., *Nucleic Acids Res.*, 2025)
- **AUC**: 0.958 (high-risk ARG), 0.964 (MDR-VF fusion), 0.856 (conjugation), 0.903 (BMRG); mean 0.920
- **Dimensionality analysis**: 5-dim lite achieves equivalent mean AUC to 10-dim full (0.920 vs 0.920); no overfitting (train-test gap < 0.003, bootstrap optimism < 0.005)
- **LORO-CV**: mean AUC = 0.962 across 40 replicons
- **Weight sensitivity** (100 iterations, +/-30%): Spearman rho = 0.994, top-10 overlap = 9.5/10
- **External validation**: 40 independent NCBI plasmids (18/20 high-risk Grade A, 19/20 low-risk Grade D/E)
- **Natural-prevalence holdout** (157,046 PSCs, 3.0% high-risk): ROC-AUC = 0.968, PR-AUC = 0.413

---

## Running tests

```bash
python -m pytest tests/ -v
```

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
