# PlasRisk

**Data-driven weighted risk assessment for bacterial plasmids from FASTA sequences**
**(5-dimension FASTA-only lite core by default; 10-dimension full model available)**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

PlasRisk computes a composite risk score for bacterial plasmids based on
biologically motivated dimensions: antimicrobial resistance gene (ARG) burden,
virulence factors (VFs), mobility/conjugation potential, host range, replicon
type, plasmid size, biocide/metal resistance (BMRG), geographic spread, habitat
breadth, and temporal growth rate. It accepts plasmid FASTA sequences and
automatically annotates them using [abricate](https://github.com/tseemann/abricate).

The default output is the FASTA-only lite core:

```
S_lite = 0.258*S_ARG + 0.115*S_VF + 0.215*S_MOB + 0.190*S_SIZE + 0.222*S_BM
```

The full 10-dimension model (`--mode full`) is:

```
S_full = 0.245*S_ARG + 0.110*S_VF + 0.204*S_MOB + 0.028*S_HOST
       + 0.003*S_REP + 0.181*S_SIZE + 0.211*S_BM
       + 0.002*S_GEO + 0.002*S_HAB + 0.015*S_GROW

Weights derived by data-driven consensus (Random Forest MDG, LASSO, and
grid-search optimization) on 792,964 PIPdb PSCs; sum ≈ 1.0.
```

Risk grades: **A** (Very High, S >= 0.60), **B** (High, >= 0.45),
**C** (Moderate, >= 0.30), **D** (Low, >= 0.15), **E** (Minimal, < 0.15).

---

## Installation

### Option 1: conda (recommended)

```bash
# Create a dedicated environment
conda create -n plasrisk -c bioconda -c conda-forge plasrisk
conda activate plasrisk
```

This installs PlasRisk together with `abricate` and `blast` for full
annotation capability.

### Option 2: pip + manual abricate

```bash
pip install plasrisk

# Install abricate separately for annotation
conda install -c bioconda abricate
# or on Debian/Ubuntu: apt-get install abricate
```

### Option 3: from source

```bash
git clone https://github.com/LLQ95/PlasRisk.git
cd plasrisk
pip install .

# Install annotation dependencies
conda install -c bioconda abricate blast
```

### Set up abricate databases

After installing abricate, download the databases you need:

```bash
# Download/update all default databases
abricate-get_db --db card --force
abricate-get_db --db vfdb --force
abricate-get_db --db plasmidfinder --force
abricate-get_db --db resfinder --force
abricate-get_db --db ncbi --force

# Verify
abricate --list
```

For BacMet (biocide/metal resistance) database, see:
https://github.com/tseemann/abricate#making-your-own-database

---

## Quick start

```bash
# Score a single plasmid
plasrisk plasmid.fasta

# Score multiple plasmids
plasrisk *.fasta

# Score all FASTA files in a directory
plasrisk /path/to/plasmids/

# Specify output directory
plasrisk -o results *.fasta

# Use specific abricate databases
plasrisk --db card,vfdb,plasmidfinder,bacmet plasmid.fasta

# Sequence-only mode (no abricate needed; scores based on length + replicon lookup)
plasrisk --no-abricate contigs.fasta

# JSON output
plasrisk --json -o results plasmid.fasta

# Full 10-dimension model with epidemiological context (default is lite)
plasrisk --mode full *.fasta

# Verify that the package S_SIZE transform and the discovery-pipeline
# transform rank plasmids identically (writes rank_concordance.tsv)
plasrisk --rank-concordance *.fasta
```

### Lite (default) vs. Full mode

The default output is the 5-dimension FASTA-only lite core; no PIPdb-style
metadata are required. The full 10-dimension model is available with
`--mode full`.

| | Lite (5-dim, default) | Full (10-dim) |
|---|---|---|
| Dimensions | S_ARG, S_VF, S_MOB, S_SIZE, S_BM | S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE, S_BM, S_GEO, S_HAB, S_GROW |
| Weights | 0.258, 0.115, 0.215, 0.190, 0.222 | 0.245, 0.110, 0.204, 0.028, 0.003, 0.181, 0.211, 0.002, 0.002, 0.015 |
| Mean AUC (4 outcomes) | 0.920 | 0.920 |
| Grade agreement | 92.7% exact / 100% within one grade vs. full | reference |
| Required annotations | ARG + VF + mobility + length + BacMet | ARG + VF + mobility + replicon + BacMet + metadata |
| Use case | Rapid FASTA-only screening, resource-limited settings | Comprehensive risk assessment with epidemiological context |

### S_SIZE normalization and rank-concordance check

The released package computes the length component with a monotonic sigmoid
on log10 length, centered at 30 kb:

```
S_SIZE = 1 / (1 + exp(-4 * (log10(L) - log10(30,000))))
```

The discovery pipeline used the log-ratio `clip(log10(L)/log10(1,000,000), 0, 1)`
for weight fitting. Both transforms are strictly increasing in L, so plasmid
rankings are identical (Spearman rho = 1.000, Kendall tau = 1.000 over the
792,964 PIPdb PSCs); only the absolute scaling differs. The equivalence can be
reproduced on any input set with `plasrisk --rank-concordance` (writes
`rank_concordance.tsv`) or via the Python API
(`from plasrisk import size_rank_concordance`). See Supplementary Text S1.

### Example output

```
  Grade A:     3 (  3.2%) ######
  Grade B:    12 ( 12.9%) ##########################
  Grade C:    28 ( 30.1%) ############################################################
  Grade D:    35 ( 37.6%) ############################################################################
  Grade E:    15 ( 16.1%) ################################

  Top 10 highest-risk plasmids:
    pNDM-1_260kb                            S=0.712  grade A  8 ARG  IncX3 [blaNDM]
    pMCR-1_33kb                             S=0.581  grade B  4 ARG  IncX4 [mcr]
    pKPC-2_110kb                            S=0.534  grade B  6 ARG  IncFII(K) [blaKPC]
```

### Output files

| File | Description |
|------|-------------|
| `plasrisk_results.tsv` | Per-sequence scores: active components, S_total, S_norm, grade, gene lists |
| `rank_concordance.tsv` | S_SIZE rank-concordance check (with `--rank-concordance`) |
| `plasrisk_summary.tsv` | Per-file summary: counts, grade distribution, mean/max scores |
| `plasrisk_results.json` | JSON format (with `--json`) |

---

## Python API

```python
from plasrisk import PlasRiskScorer, PlasmidFeatures, annotate_fasta, load_replicon_lookup

# Option A: annotate a FASTA file directly (default = 5-dim FASTA-only lite)
lookup = load_replicon_lookup()
result = annotate_fasta("plasmid.fasta", lookup=lookup)
scorer = PlasRiskScorer(replicon_lookup=lookup)  # lite by default
df = scorer.score_dataframe(result.features)
print(df[["seq_id", "S_norm", "grade", "high_risk_genes"]])

# Full 10-dimension model with epidemiological context
scorer_full = PlasRiskScorer(replicon_lookup=lookup, mode="full")
df_full = scorer_full.score_dataframe(result.features)

# S_SIZE rank-concordance check (package sigmoid vs pipeline log-ratio)
from plasrisk import size_rank_concordance
print(size_rank_concordance([f.length_bp for f in result.features]))

# Option B: construct features manually
feat = PlasmidFeatures(
    seq_id="pExample",
    length_bp=85000,
    arg_names=["NDM-1", "CTX-M-15", "TEM-1"],
    vf_names=["aerobactin"],
    vf_categories=["Nutritional/Metabolic factor"],
    bm_gene_names=["merA", "qacEdelta1"],
    replicon="IncX3",
    has_t4cp=True,
    has_relaxase=True,
    has_oriT=True,
    has_auxiliary=True,
)
scores = scorer.score(feat)
print(f"S_norm = {scores['S_norm']:.3f}, grade = {scores['grade']}")
```

---

## The 10 risk dimensions

| Component | Weight | What it measures | Scoring basis |
|-----------|--------|------------------|---------------|
| **S_ARG** | 0.245 | ARG count, WHO-priority genes, high-risk genes (mcr, NDM, KPC, CTX-M, tetX, etc.) | Base + per-gene + high-risk bonuses |
| **S_BM** | 0.211 | Biocide/metal resistance (mer, qac, ars/cop/sil) — co-selection potential | Base + per-gene + family bonuses |
| **S_MOB** | 0.204 | T4CP, relaxase, oriT, auxiliary transfer proteins | Element-based additive score |
| **S_SIZE** | 0.181 | Plasmid length (cargo capacity) | Sigmoid: midpoint 30 kb |
| **S_VF** | 0.110 | VF count, exotoxins, secretion systems (T3SS/T4SS) | Base + per-gene + category bonuses |
| **S_HOST** | 0.028 | Number of host genera / replicon prior | Empirical host range or lookup |
| **S_GROW** | 0.015 | Annual growth rate of the replicon | PIPdb-derived lookup |
| **S_REP** | 0.003 | Replicon backbone risk (IncX3, IncN, ColKP3 high; ColpVC low) | PIPdb-derived lookup table |
| **S_HAB** | 0.002 | Habitat breadth (human/animal/environment) | PIPdb-derived lookup |
| **S_GEO** | 0.002 | Number of countries observed | PIPdb-derived lookup |

---

## Command-line options

```
plasrisk [options] <fasta1> [fasta2 ...]

positional arguments:
  FASTA                 FASTA file(s) or directory

options:
  -o, --output DIR      Output directory (default: ./plasrisk_output)
  -t, --threads N       Number of abricate threads (default: 4)
  --min-id FLOAT        Minimum abricate identity % (default: 75)
  --min-cov FLOAT       Minimum abricate coverage % (default: 50)
  --no-abricate         Skip abricate; sequence-only scoring
  --mode {lite,full}    Scoring mode: lite (5-dim FASTA-only core, default)
                        or full (10-dim)
  --rank-concordance    Verify S_SIZE rank concordance (sigmoid vs log-ratio)
  --db LIST             Comma-separated abricate databases (default: auto)
  --json                Also write JSON output
  -q, --quiet           Suppress progress messages
  -v, --version         Show version
  -h, --help            Show help
```

---

## Model validation

The PlasRisk model was developed and validated using 792,964 plasmid sequence
clusters from PIPdb (Zhu et al., *Nucleic Acids Res.*, 2025). Validation
included:

- **Quartile stratification**: Q1 (highest risk) plasmids had 65.6% ARG prevalence,
  12.2% high-risk ARG rate, 13.3% conjugative rate (vs. 0% in Q4).
- **Data-driven weights**: RF-MDG, LASSO, and grid-search optimization across four
  outcomes (high-risk ARG, MDR-VF fusion, conjugative capacity, BMRG carriage)
  converged on S_ARG (0.245), S_BM (0.211), S_MOB (0.204), and S_SIZE (0.181)
  as dominant predictors.
- **AUC validation**: Final weights achieved AUC 0.958 (high-risk ARG), 0.964
  (MDR-VF fusion), 0.856 (conjugation), 0.903 (BMRG); mean 0.920.
- **Dimensionality analysis**: All-subsets evaluation (1,023 subsets, 5-fold CV)
  showed that 5 core dimensions (ARG+VF+MOB+SIZE+BM) achieve equivalent mean AUC
  to the full 10-dim model (0.920 vs. 0.920), provided as `--mode lite`.
  Overfitting diagnostics (train-test gap, bootstrap optimism, learning curves)
  confirmed no excess optimism in either model.
- **Leave-one-replicon-out CV**: mean AUC = 0.962 across 40 replicons.
- **Weight perturbation sensitivity** (100 iterations, +/-30%): mean Spearman
  rho = 0.994, mean top-10 overlap = 9.2/10.
- **External validation**: 40 independent NCBI plasmids correctly classified
  (18/20 high-risk Grade A, 19/20 low-risk Grade D/E).

---

## Uploading to conda (bioconda)

To make PlasRisk installable via `conda install -c bioconda plasrisk`:

### Step 1: Upload to PyPI

```bash
# Install build tools
pip install build twine

# Build distributions
python -m build

# Upload to PyPI
twine upload dist/*
```

### Step 2: Fork and clone bioconda-recipes

```bash
git clone https://github.com/bioconda/bioconda-recipes.git
cd bioconda-recipes
```

### Step 3: Create the recipe

```bash
# Create recipe directory
mkdir -p recipes/plasrisk
```

Create `recipes/plasrisk/meta.yaml`:

```yaml
{% set version = "1.0.0" %}

package:
  name: plasrisk
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/p/plasrisk/plasrisk-{{ version }}.tar.gz
  sha256: <SHA256 from PyPI>

build:
  number: 0
  noarch: python
  entry_points:
    - plasrisk = plasrisk.cli:main
  script: "{{ PYTHON }} -m pip install . --no-deps --ignore-installed -vv"

requirements:
  host:
    - python >=3.9
    - pip
    - setuptools >=61.0
    - wheel
  run:
    - python >=3.9
    - pandas >=1.3
    - numpy >=1.20
    # abricate/blast optional; CLI falls back to --no-abricate mode

test:
  imports:
    - plasrisk
  commands:
    - plasrisk --help
    - plasrisk --version

about:
  home: https://github.com/LLQ95/PlasRisk
  license: MIT
  license_file: LICENSE
  summary: "Ten-dimension data-driven weighted risk assessment for bacterial plasmids"
```

> **Note:** The complete, ready-to-submit recipe is in the `bioconda/` directory.
> See `UPLOAD_GUIDE.md` for the full step-by-step release process.

### Step 4: Test locally

```bash
# Install bioconda-utils
conda install -c bioconda bioconda-utils

# Test the recipe
bioconda-utils build recipes/plasrisk --docker
```

### Step 5: Submit a pull request

```bash
git checkout -b plasrisk
git add recipes/plasrisk/
git commit -m "Add plasrisk recipe"
git push origin plasrisk
# Open PR at https://github.com/bioconda/bioconda-recipes
```

Once the PR is merged and CI passes, PlasRisk will be installable via:

```bash
conda install -c bioconda plasrisk
```

### Local conda build (without bioconda)

```bash
# Build from the conda/ directory in this repo
conda build conda/

# Install locally
conda install --use-local plasrisk
```

---

## Running tests

```bash
cd plasrisk_py
python -m pytest tests/ -v
# or
python tests/test_scoring.py
```

---

## Citation

If you use PlasRisk, please cite:

> [Authors]. PlasRisk: a ten-dimension weighted risk assessment framework for
> bacterial plasmids. *Journal*, 2025. doi: [to be added]

The model is based on data from:
> Zhu Q, Chen Q, Lu X, et al. PIPdb: a comprehensive plasmid sequence resource
> for tracking the horizontal transfer of pathogenic factors and antimicrobial
> resistance genes. *Nucleic Acids Research*, 2025, 53(D1):D169-D178.
> doi:10.1093/nar/gkae952

---

## License

MIT License - see [LICENSE](LICENSE) for details.
