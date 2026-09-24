# Changelog

All notable changes to PlasRisk are documented in this file.

## [1.2.0] - 2026-09-24

### Changed
- The biocide/metal resistance component is renamed from S_BM to S_BMG,
  matching the endpoint name used throughout the manuscript and avoiding
  ambiguity. This applies to the Python API keys, CLI/TSV output columns,
  lookup tables, documentation, graphical abstract, and the locked data
  tables. Weight values are unchanged. This is a breaking change for code
  that reads the S_BM output column; use S_BMG instead.

## [1.1.2] - 2026-09-23

### Fixed
- Package weight constants now match the locked consensus weights released
  in data/locked_v1.1.1 and reported in manuscript Table 1 (S_ARG 0.237,
  S_BMG 0.222, S_MOB 0.189, S_SIZE 0.164, S_VF 0.126, S_HOST 0.030,
  S_HAB 0.019, S_GEO 0.005, S_REP 0.004, S_GROW 0.002). Package builds
  1.1.0 and 1.1.1 carried the pre-lock weights.
- Terms in all lite and full formulas are listed in descending weight
  order. Lite weights are 0.253 (S_ARG), 0.237 (S_BMG), 0.201 (S_MOB),
  0.175 (S_SIZE), and 0.135 (S_VF).

## [1.1.1] - 2026-09-21

### Changed
- The 5-dimension FASTA-only lite core (S_ARG, S_VF, S_MOB, S_SIZE, S_BMG)
  is now the default for both the CLI (`plasrisk`) and the Python API
  (`PlasRiskScorer()`, `get_scorer()`). The 10-dimension model remains
  available with `--mode full` / `mode="full"`.
- S_SIZE in the released package uses the monotonic sigmoid
  `1 / (1 + exp(-4 * (log10(L) - log10(30,000))))` (midpoint 30 kb).
  The discovery pipeline's log-ratio transform
  `clip(log10(L)/log10(1,000,000), 0, 1)` is retained as
  `PlasRiskScorer.score_size_logratio` for exact reproduction of the
  manuscript's weight-fitting analyses.

### Added
- `size_rank_concordance()` (Python API) and `plasrisk --rank-concordance`
  (CLI, writes `rank_concordance.tsv`) to verify that the sigmoid and
  log-ratio S_SIZE transforms rank plasmids identically. Over the
  792,964 PIPdb PSCs, Spearman rho = 1.000 and Kendall tau = 1.000;
  see Supplementary Text S1.
- Unit tests for the lite default, the log-ratio transform, and the
  rank-concordance check.

## [1.1.0] - 2026-08

### Added
- 5-dimension lite core with renormalized weights.
- `get_scorer()` factory and PIPdb 8-item ordinal model for comparison.

## [1.0.0] - 2026-08

### Added
- Initial release: 10-dimension data-driven weighted plasmid risk score,
  abricate-based annotation pipeline, CLI and Python API.
