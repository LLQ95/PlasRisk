# Changelog

All notable changes to PlasRisk are documented in this file.

## [1.5.0] - 2026-10-10

### Added
- Metagenomic three-layer scoring under the new `plasrisk meta` command
  (Python API: `plasrisk.metagenome`). Layer 1 scores fragmented plasmid
  contigs on the dimensions that remain identifiable after fragmentation:
  S_ARG/S_BMG/S_VF are always evaluated, S_MOB only from observed transfer
  apparatus or IS evidence (no isolate-derived priors are injected), and
  S_SIZE only for circularised headers or contigs at/above
  `--complete-length` (default 10,000 bp). Scores are renormalised over the
  assessable dimensions as S_cargo3, S_frag4 and S_lite5, with A-E grades
  emitted only when all five lite dimensions are assessable.
- Layer 2 aggregates contigs into predicted plasmid bins (generic
  `--bins` table, best-effort geNomad `--genomad-bin` parser for
  `*_plasmid_components.tsv` / `*_plasmid_summary.tsv`, or singleton
  bins for unbinned contigs): gene sets are unioned with deduplication,
  lengths are summed, mobility markers are OR-combined, and optional
  per-contig coverage (`--coverage`, including
  jgi_summarize_bam_contig_depths output) is length-weighted. Bins carry
  complete/probable/fragmented confidence labels.
- Layer 3 reports descriptive, abundance-aware community statistics
  (CPRI_mean/P90/P99, grade prevalence, coverage-weighted grade A/B and
  critical-ARG loads). CPRI is an unvalidated descriptive extension and is
  labelled as such in every sample summary.
- `plasrisk-fragment-sim` console script (`plasrisk.fragment_sim`) for the
  controlled fragmentation grid used to quantify score recovery versus
  intact plasmids (Spearman rho, grade exact/within-one agreement, median
  bias, critical-ARG recovery per fragment length).
- New test module `tests/test_metagenome.py` (29 tests total).

### Design notes
- The frozen consensus weights are never refit on metagenomic data.
- Fragment lengths are never passed through S_SIZE, avoiding systematic
  penalties for megaplasmids split by short-read assembly.

## [1.4.0] - 2026-09-30

### Added
- S_ARG can use the WHO AWaRe 2025 category table, selected with
  `--aware-version 2025` on the CLI or `aware_version="2025"` in
  `PlasRiskScorer` / `get_scorer`. Under the 2025 mapping, biocide/metal
  genes contribute 0 to S_ARG (they are scored under S_BMG) and the
  separate 1.5x last-resort multiplier is not applied because Reserve
  already carries the highest weight. The legacy mapping used for the
  frozen component matrix remains the default. The audit in Method S14
  shows that the 2025 mapping changes S_ARG by at least 0.05 for 4.9%
  of PSCs and moves the high-risk-ARG AUC on the locked test from
  0.9807 to 0.9789 (delta -0.0018); the frozen values are retained.

## [1.3.0] - 2026-09-24

### Added
- The original PIPdb 8-item ordinal index is now exposed as a third mode,
  `--mode ordinal` (grades 1–5), alongside the default continuous lite
  core (`--mode lite`) and the 10-dimension full model (`--mode full`).
  The Python factory accepts it as `get_scorer("plasrisk", mode="ordinal")`
  in addition to `get_scorer("pipdb")` and the `PIPdbScorer` class.

### Deprecated
- The separate `--model {plasrisk,pipdb}` flag is deprecated; `--model
  pipdb` continues to work and is equivalent to `--mode ordinal`. Mode
  selection should use `--mode {lite,full,ordinal}`.

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
