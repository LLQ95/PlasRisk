# PlasRisk locked-analysis data release (v1.1.1)

This folder contains the reproducibility package for the PlasRisk manuscript
submitted to iMeta. It allows every number in the weight-derivation,
dimensionality, locked-test validation, grade-stratification and external
case-control challenge analyses to be reproduced from plain CSV files, without
access to the PIPdb source database or the annotation pipeline.

All analyses use a single, pre-specified split (`seed = 42`):
**train 475,779 / validation 158,591 / locked test 158,594 = 792,964 PSCs**.
The locked test set was evaluated exactly once.

## Files

Large tables are split into plain-CSV shards (GitHub's web API limit); each
shard carries the header row and can be concatenated as shown below. All files
in this folder are listed with SHA-256 hashes in `SHA256SUMS.txt`.

| File | Rows | Contents |
|---|---|---|
| `split/locked_split_ids.csv` | 792,964 | PSC `id` and split assignment (`train` / `validation` / `locked_test`) |
| `components/plasrisk_component_matrix_part01.csv` … `part08.csv` | 792,964 in total | Ten normalized component scores per PSC plus primary replicon (8 shards, ~99 k rows each) |
| `predictions/plasrisk_locked_test_predictions.csv` | 158,594 | Frozen full/lite scores, grades and four endpoint labels for the untouched locked test set |
| `predictions/plasrisk_validation_predictions.csv` | 158,591 | The same columns for the validation partition |
| `predictions/plasrisk_train_predictions_part01.csv`, `part02.csv` | 475,779 in total | The same columns for the training partition (2 shards) |
| `external_case_control_challenge_scores_367.csv` | 367 | Frozen scores for the external case-control challenge (167 high-risk vs 200 low-risk plasmids, 95%-ANI deduplicated) |
| `locked_weight_vectors.csv` | 10 | Expert, RF-MDG, LASSO, grid-search and consensus weights, and the five-dimension lite weights |

Reassemble the sharded tables:

```bash
# component matrix (792,964 rows)
head -n 1 components/plasrisk_component_matrix_part01.csv > component_matrix.csv
tail -n +2 -q components/plasrisk_component_matrix_part*.csv >> component_matrix.csv
# training predictions (475,779 rows)
head -n 1 predictions/plasrisk_train_predictions_part01.csv > train_predictions.csv
tail -n +2 -q predictions/plasrisk_train_predictions_part*.csv >> train_predictions.csv
```

## Column dictionary

### Component matrix

| Column | Meaning | Range |
|---|---|---|
| `id` | PIPdb plasmid segment cluster (PSC) identifier | text |
| `replicon_primary` | primary PlasmidFinder replicon assigned to the PSC | text / `Unknown` |
| `S_ARG` | critical-antibiotic-resistance-gene component | 0–1 |
| `S_VF` | virulence-factor component | 0–1 |
| `S_MOB` | mobility apparatus component (relaxase/T4CP/oriT) | 0–1 |
| `S_HOST` | host-range component | 0–1 |
| `S_REP` | replicon-type component | 0–1 |
| `S_SIZE` | plasmid-size component (sigmoid centred at 30 kb) | 0–1 |
| `S_BMG` | biocide/metal-resistance component | 0–1 |
| `S_GEO` | geographic-spread component | 0–1 |
| `S_HAB` | habitat-breadth component | 0–1 |
| `S_GROW` | isolate-record-growth component | 0–1 |

### Prediction files

`S_full` is the consensus ten-dimension score; `S_lite` is the FASTA-only
five-dimension score (S_ARG, S_VF, S_MOB, S_SIZE, S_BMG). Grades use fixed
thresholds: **A ≥ 0.60, B ≥ 0.45, C ≥ 0.30, D ≥ 0.15, E below 0.15**.
The four labels are `y_highrisk` (critical ARG carriage), `y_fusion`
(ARG–VF fusion), `y_conj` (complete conjugative apparatus) and `y_bm`
(biocide/metal resistance), all coded 0/1.

The external challenge file additionally contains accession, raw annotation
counts and the assigned grade. It is a deliberately extreme case-control
panel whose labels derive from annotation; it is **not** an independent
prospective validation cohort.

## Reproducing the scores

```python
import pandas as pd, glob

m = pd.concat([pd.read_csv(f) for f in
               sorted(glob.glob("components/plasrisk_component_matrix_part*.csv"))],
              ignore_index=True)
w = pd.read_csv("locked_weight_vectors.csv").set_index("component")
comps = ["S_ARG","S_VF","S_MOB","S_HOST","S_REP","S_SIZE",
         "S_BMG","S_GEO","S_HAB","S_GROW"]
lite5 = ["S_ARG","S_VF","S_MOB","S_SIZE","S_BMG"]

S_full = sum(m[c] * w.loc[c, "consensus"] for c in comps)
S_lite = sum(m[c] * w.loc[c, "lite"]      for c in lite5)
```

Consensus weights are the equal-vote consensus of random-forest
mean-decrease-Gini, LASSO and grid-AUC searches, each fitted across the four
reference endpoints on the training partition. The package implementation
(`plasrisk` v1.1.1, PyPI) produces identical scores directly from plasmid
FASTA files; see the main repository README for the command-line interface.

## Verification values (locked test, n = 158,594)

| Endpoint | ROC-AUC, full | ROC-AUC, lite |
|---|---|---|
| Critical ARG | 0.9549 | 0.9552 |
| ARG–VF fusion | 0.9668 | 0.9658 |
| Conjugative | 0.8462 | 0.8448 |
| Biocide/metal | 0.9095 | 0.9109 |
| Mean | 0.9194 | 0.9192 |

Grade counts on the locked test: A 1,136 (0.72%), B 5,835 (3.68%),
C 24,765 (15.62%), D 62,100 (39.16%), E 64,758 (40.83%). Grades A+B select
4.40% of PSCs and capture 51.1% of critical-ARG positives. Full/lite grade
agreement is 92.7% exact and 100% within one grade (Pearson 0.998,
Spearman 0.991, top-decile overlap 98.0%).

The external case-control challenge yields ROC-AUC 0.998 (95% CI
0.994–1.000); because the panel is annotation-labelled and extreme-case by
construction, it demonstrates score separation, not transportability to
unselected populations.

## Provenance and licensing

Component scores were computed with the PlasRisk annotation pipeline on the
PIPdb release (Zhu et al., *Nucleic Acids Research*, 2025): 1,009,571
pathogen genomes, 19,292,110 plasmid contigs, 792,964 PSCs at 95% ANI over
80% alignment coverage. No length filter was applied at the PSC level;
52,234 PSCs (6.59%) with mean length below 2 kb were retained and scored.
PSC identifiers remain the property of the PIPdb resource. These derived
score/split tables are released under the same license as the PlasRisk
repository (see LICENSE at the repository root).

## Corresponding Data Availability statement

"The PIPdb-derived component-score matrix (792,964 PSCs × 10 dimensions),
the seed-42 train/validation/locked-test split identifiers, all frozen
predictions and grades, the external case-control challenge scores, and the
weight vectors are available in the PlasRisk GitHub repository
(https://github.com/LLQ95/PlasRisk, `data/locked_v1.1.1/`) and archived at
Zenodo [DOI]. The source plasmid sequences are available from PIPdb
(https://nmdc.cn/pipdb; Zhu et al., *Nucleic Acids Research*, 2025). The
PlasRisk scoring package is available from PyPI (plasrisk v1.1.1) and
Bioconda, and the complete R/Python analysis code is provided in the same
GitHub repository."
