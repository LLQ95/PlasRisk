# Changelog

All notable changes to PlasRisk are documented in this file.

## [1.9.0] - 2026-10-26

### Added
- `armA` to the high-risk (critical) ARG family list. ArmA is a 16S rRNA
  methyltransferase that confers high-level resistance to nearly all
  clinically important aminoglycosides, matching the `rmt` family already
  in the list and the AWaRe 2025 Reserve mapping.

### Changed
- The critical-ARG label is now evaluated against the 14-family critical
  list on the full PIPdb catalogue, removing elements that carried only
  non-critical determinants (e.g. FosA, TEM, non-48-like OXA,
  AAC(6')-Ib-cr, DHA). Corrected critical prevalence is 2.52% of the
  catalogue; downstream discrimination and grade statistics are updated
  in the documentation.

## [1.8.0] - 2026-10-25

### Added
- System-level conjugative annotation with CONJScan (Cury et al., 2020)
  running in MacSyFinder v2 (Neron et al., 2023). The new
  `plasrisk.conjscan` module predicts proteins per replicon with Prodigal,
  runs the CONJScan/Plasmids models in ordered-replicon mode (gembase batch
  with a per-replicon topology file; linear by default, circular for closed
  plasmids) or content-only unordered mode, and reports complete, likely,
  putative conjugative (pConj) and putative mobilizable (pMob) systems with
  MPF types. The final mobility class is the union of this independent
  system-level evidence and the oriTfinder/oriTDB/MOB-typer-style gene
  evidence, so a system-level call can recover T4CP/relaxase calls missed
  by gene-name matching. The new `conjugative_putative` class receives an
  S_MOB base of 0.7, between conjugative_likely (0.85) and mobilizable
  (0.60).
- CLI flags `--no-conjscan`, `--conjscan-topology {linear,circular}` and
  `--conjscan-db-type {auto,unordered}`; `plasrisk meta` defaults to
  unordered content-only searches for fragmented metagenomic contigs and
  propagates the CONJScan class into bin-level and sample-level output.
- Per-sequence output columns `conjscan_class` and `conjscan_mpf`.
- When macsyfinder/prodigal are absent the tool degrades gracefully to
  gene-based evidence and emits an install warning; install with
  `pip install plasrisk[conj]` and `conda install -c bioconda prodigal`.

### Changed
- The locked training labels and frozen catalogue scores are unchanged;
  CONJScan operates on new inputs and provides the independent validation
  and union-label sensitivity analysis described in the supplementary
  materials.

## [1.7.3] - 2026-10-15

### Fixed
- `plasrisk meta` no longer runs silently when abricate is unavailable. When
  abricate was requested but no annotation database could be used (abricate
  missing, no installed databases, or BLAST dependency missing), annotation
  failures were previously swallowed and the run produced all-zero
  ARG/VF/BMG columns with no warning, leaving only length-driven scores.
  Annotation warnings are now printed to stderr, a prominent warning is
  emitted when no database is available or when databases return zero hits
  across more than 100 sequences, and the new `--strict-annotation` flag
  exits with code 3 instead of writing misleading results.
- Community sample summaries no longer treat every contig outside the
  caller's bin map as a plasmid bin. Bins are tagged with `bin_origin`
  (`predicted` for caller-supplied bins such as geNomad plasmid bins,
  `unbinned` for singleton contigs outside the bin map), and
  `aggregate_sample` additionally reports `n_predicted_bins`,
  `n_unbinned_bins`, and a full set of `pred_*` statistics computed over
  predicted plasmid bins only. These are the primary plasmid-level
  statistics; in whole-community assemblies the unbinned set is dominated
  by chromosome fragments and previously diluted plasmid prevalence by up
  to two orders of magnitude. Unprefixed columns retain the all-element
  (community-level) denominator.

## [1.7.2] - 2026-10-14

### Fixed
- `plasrisk meta` no longer raises `KeyError` when the coverage column passed
  via `--coverage-column` does not exactly match the jgi_summarize_bam_contig_depths
  header. jgi names the depth column after the BAM argument verbatim, which may
  be a full path or carry a `.bam` suffix. The column is now resolved by exact
  name, then by basename (directory and `.bam` stripped), then by the first
  jgi sample column, so a bare run accession works as expected.
- Benchmark script `02_public_metag_benchmark.sh` reads the real depth column
  from the jgi table header instead of assuming the bare accession, cleans an
  incomplete MEGAHIT output directory before reassembly, gzip-validates reads
  before fastp, skips completed samples and resumable per-step outputs, and
  prints the failing step and last log lines on a per-sample failure.

## [1.7.1] - 2026-10-11

### Changed
- Documentation states the generality explicitly: PlasRisk is not limited to
  curated plasmid sequences. A single FASTA may mix chromosomal and plasmid
  contigs, or complete and fragmented molecules from several replicons, and
  each contig is classified separately (replicon/circularity on long reads,
  geNomad on short reads) so that only the predicted plasmid sequences are
  extracted and scored and chromosomes are ignored.
- Releases are published non-interactively from the project's configured
  PyPI credential; no functional scoring change versus 1.7.0.

## [1.7.0] - 2026-10-11

### Added
- Unified `plasrisk auto` entry point (Python API:
  `plasrisk.extract.auto_main`, `characterize_assembly`,
  `choose_auto_strategy`). One command accepts a ready plasmid collection or
  a raw short-read WGS, long-read WGS, or metagenomic assembly and decides
  the workflow from the assembly itself: collections with no
  chromosome-scale molecule are scored directly; assemblies containing a
  chromosome are split by contig size/count into long-read vs short-read.
- Long-read non-plasmid assemblies are mined for plasmids with a new
  replicon classifier (`--classifier replicon`, Python API
  `classify_by_replicon`): contigs carrying a PlasmidFinder replicon /
  replication gene (abricate `plasmidfinder`, or a precomputed
  `--replicon-table`) are extracted, as are circularised molecules in the
  plasmid size range; chromosome-scale contigs without plasmid signal are
  excluded. If the replicon scan finds nothing and geNomad with a database
  is available, `auto` falls back to geNomad.
- Short-read non-plasmid assemblies are predicted with geNomad end-to-end
  and then routed to the isolate or metagenomic layers based on contig
  completeness; metagenomic input always uses the three-layer aggregation.
- Override switches `--input-type {auto,plasmid,genome,metagenome}`,
  `--read-type {auto,short,long}`, and `--classifier`. A per-sample
  `auto_status.tsv` records the detected read type, action, classifier,
  route, plasmid count, and any unscorable reason.
- FASTQ reads are not assembled in place. `auto` detects `.fastq/.fq(.gz)`
  input, exits with code 2, and prints the recommended assemblers
  (MEGAHIT/metaSPAdes for short reads; Flye/metaFlye for long reads).

### Changed
- geNomad invocation updated to the v1.7+ command line (verified on
  v1.12): `genomad end-to-end INPUT OUTPUT DATABASE --cleanup --threads N
  --min-score T`, replacing the removed `genomad classify`. Summary tables
  are located recursively so the nested `<stem>_summary/` layout is handled.

## [1.6.1] - 2026-10-10

### Fixed
- `load_genomad` no longer requires a `*_plasmid_components.tsv`, which
  geNomad v1.12 end-to-end does not produce. It now reads
  `*_plasmid_summary.tsv` (one row per plasmid, with topology and
  plasmid_score), falls back to `*_aggregated_classification.tsv` filtered
  to the plasmid class, and still supports the legacy components layout.
- The pilot grid and public metagenomic benchmark reuse a completed
  geNomad run by checking for `*_plasmid_summary.tsv`, so reruns do not
  recompute plasmid prediction or fail on existing output.

## [1.6.0] - 2026-09-23

### Added
- Raw-assembly prediction front end under the new `plasrisk run` command
  (Python API: `plasrisk.extract`). A short-read WGS, long-read, or
  multi-strain metagenomic assembly is length-filtered (>=1 kb), classified
  with PlasFlow or geNomad, and the predicted plasmid contigs are extracted
  to `predicted_plasmids.fasta`. Classifiers are invoked through subprocess
  when present and are not bundled; precomputed outputs can be supplied with
  `--plasflow-table`, `--genomad-table`, or a generic `--generic-table`.
- Automatic routing (`--route auto`): long/circular plasmid contigs go to
  the isolate scorer; short fragmented contigs and plasmid bins go to the
  metagenomic three-layer. Multiple classifiers are combined by majority
  vote (`--vote majority`, or `any` for high recall).
- Unscorable reporting. Assemblies with no detected plasmid are not forced
  into a score. They are reported with a reason, including a chromosome-only
  determination based on total/longest contig length and the classifier
  labels (`plasrisk_unscorable.tsv`). A per-contig `extraction_manifest.tsv`
  and a top-level `run_status.tsv` are always written.

### Fixed
- `plasrisk meta` previously called the annotation routine with arguments
  (`db_path`, `quiet`) that the routine does not accept, which raised a
  `TypeError` on the abricate path. It now uses the supported interface.

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
