# PlasRisk

PlasRisk is a data-driven risk assessment framework for bacterial
plasmids in One Health surveillance. It combines ten biologically
meaningful dimensions into a single, interpretable score and ranks
plasmids by their potential to carry and disseminate clinically relevant
antibiotic resistance and virulence. PlasRisk was developed on 792,964
plasmid segment clusters (PSCs) from the PIPdb resource and is
distributed as a command-line tool and a Python package.

![PlasRisk graphical abstract](graphical_abstract.png)

## Why PlasRisk

Plasmids are the main vectors of horizontal gene transfer, but only a
small subset pose an immediate public health risk. Existing databases
provide either raw sequence archives or single-feature annotations that
do not separate a transferable multidrug resistance plasmid from a
cryptic fragment. PlasRisk integrates resistance, virulence, mobility,
plasmid size, host range and geographic spread into one calibrated
score, with frozen weights and a transparent A-E grading scheme. It can
be run on a single plasmid FASTA, on fragmented metagenomic plasmid
contigs and bins, or directly on a raw short-read/long-read assembly
without an isolated plasmid sequence.

## Three scoring modes

PlasRisk offers three modes. The default lite core needs only a FASTA
file and is fully reproducible without external databases. The full
model adds five metadata-dependent dimensions. The original PIPdb
8-item ordinal index is retained for direct comparison.

The lite formula, with terms in descending weight order, is

S_lite = 0.253*S_ARG + 0.237*S_BMG + 0.201*S_MOB + 0.175*S_SIZE + 0.135*S_VF

The full formula is

S_full = 0.237*S_ARG + 0.222*S_BMG + 0.189*S_MOB + 0.164*S_SIZE + 0.126*S_VF + 0.030*S_HOST + 0.019*S_HAB + 0.005*S_GEO + 0.004*S_REP + 0.002*S_GROW

The dimensions are:

- S_ARG: critical antibiotic resistance genes, using the WHO AWaRe
  categories (an optional 2025 AWaRe table is available).
- S_BMG: biocide and metal resistance genes.
- S_VF: virulence factors and virulence categories.
- S_MOB: mobilisation potential, from relaxase, type IV coupling
  protein, oriT, auxiliary transfer genes and insertion sequences.
- S_SIZE: plasmid size through a monotonic sigmoid centred at 30 kb.
- S_HOST: host range across genera.
- S_HAB: number of distinct habitats (human, animal, environment).
- S_GEO: geographic spread across countries.
- S_REP: replicon prevalence.
- S_GROW: annual growth rate of the cluster.

## Installation

Install from PyPI:

```bash
pip install plasrisk
```

The lite core depends only on pandas. For sequence annotation, install
abricate separately:

```bash
conda install -c bioconda abricate
```

For the raw-assembly front end, install a plasmid classifier (either is
sufficient); these are optional and are not bundled:

```bash
# geNomad (recommended, maintained)
pip install genomad
# or PlasFlow
conda install -c bioconda plasflow
```

## Quick start

Score one plasmid with the default lite core:

```bash
plasrisk plasmid.fasta
```

Score a folder of plasmids and write results to a directory:

```bash
plasrisk -o results /path/to/plasmids/
```

Run the full ten-dimension model:

```bash
plasrisk --mode full plasmid.fasta
```

Run the original PIPdb ordinal index:

```bash
plasrisk --mode ordinal plasmid.fasta
```

Use the WHO AWaRe 2025 table for S_ARG:

```bash
plasrisk --aware-version 2025 plasmid.fasta
```

Verify that the package sigmoid S_SIZE and the pipeline log-ratio S_SIZE
rank the input sequences identically:

```bash
plasrisk --rank-concordance plasmid.fasta
```

## Scoring a raw assembly without an isolated plasmid

Most real FASTA files are short-read WGS, long-read, or multi-strain
metagenomic assemblies that do not contain a separately isolated
plasmid. The `plasrisk run` command predicts plasmid sequences first,
extracts them, and then scores them. It works with PlasFlow or geNomad:

```bash
# geNomad (provide its database)
plasrisk run assembly.fasta --classifier genomad --genomad-db /path/to/genomad_db -o run_out

# PlasFlow
plasrisk run assembly.fasta --classifier plasflow -o run_out
```

You can pass several assemblies or a directory:

```bash
plasrisk run /path/to/assemblies/ --classifier genomad --genomad-db /path/to/genomad_db -o run_out
```

If you have already run a classifier, supply its output table directly
and skip classification:

```bash
plasrisk run assembly.fasta --plasflow-table plasflow_predictions.tsv -o run_out
plasrisk run assembly.fasta --genomad-table genome_plasmid_summary.tsv -o run_out
plasrisk run assembly.fasta --generic-table labels.tsv -o run_out
```

For each sample, PlasRisk writes a per-contig `extraction_manifest.tsv`,
the extracted `predicted_plasmids.fasta`, and the score tables. Long or
circular plasmid contigs are routed to the isolate scorer, while short
fragmented contigs are routed to the metagenomic layers. A top-level
`run_status.tsv` summarises every sample.

Assemblies with no detected plasmid are not forced into a score. They
are reported as unscorable with a reason in `plasrisk_unscorable.tsv`.
When the assembly is chromosome-only, identified from the classifier
labels and the total or longest contig length, the result states that
PlasRisk scores plasmids rather than chromosomes. Small or ambiguous
assemblies are reported as too fragmented to classify, with a suggestion
to lower the classifier threshold or change the classifier.

## Metagenomic scoring

The `plasrisk meta` command scores fragmented plasmid contigs,
predicted plasmid bins and whole community samples in three layers:

```bash
plasrisk meta plasmid_contigs.fasta -o meta_out
```

With predicted bins and coverage, for example from geNomad and the JGI
coverage script:

```bash
plasrisk meta plasmid_contigs.fasta \
    --genomad-bin genomad_output \
    --coverage coverage.jgi.tsv \
    -o meta_out
```

Layer 1 reports per-contig component scores. S_ARG, S_BMG and S_VF are
always counted. S_MOB is assessed only from observed transfer apparatus
or insertion sequence evidence, and S_SIZE only for circularised
headers or contigs at or above the complete-length threshold. Missing
dimensions are reported rather than imputed, and the lite score is
renormalised over the assessable dimensions as S_cargo3, S_frag4 and
S_lite5. Layer 2 aggregates contigs into bins with unioned gene sets,
summed lengths and OR-combined mobility markers, with complete,
probable and fragmented confidence labels. Layer 3 reports descriptive,
abundance-aware community statistics (CPRI). CPRI is a descriptive
extension and is not trained or validated on metagenomic data.

## End-to-end workflow

A complete workflow from a raw assembly to PlasRisk scores uses seqkit,
a plasmid classifier and abricate. The individual steps are shown below;
`plasrisk run` performs the equivalent routing automatically.

Short-read WGS or metagenomic reads are assembled first, for example
with SPAdes or MEGAHIT, and short contigs are removed:

```bash
seqkit seq -m 1000 assembly.fasta > assembly.1k.fasta
```

Predict plasmids with geNomad:

```bash
genomad classify -t 8 --cleanup assembly.1k.fasta genomad_out /path/to/genomad_db
```

or with PlasFlow:

```bash
PlasFlow.py --input assembly.1k.fasta --output plasflow_out --threshold 0.7
```

Then run PlasRisk, which extracts the predicted plasmid contigs and
routes them by completeness:

```bash
plasrisk run assembly.fasta --classifier genomad --genomad-db /path/to/genomad_db -o run_out
```

For long-read assemblies, circular plasmid contigs are detected from
the assembly headers and topology fields and are routed directly to the
isolate scorer:

```bash
plasrisk run longread_assembly.fasta --classifier genomad --genomad-db /path/to/genomad_db -o run_out
```

## Fragmentation recovery

To quantify how scores recover after assembly fragmentation, PlasRisk
ships a controlled fragmentation grid:

```bash
plasrisk-fragment-sim --input intact_plasmids.fasta --out frag_out
```

It fragments intact plasmids across a grid of fragment lengths and
reports Spearman correlation, grade exact and within-one agreement,
median bias and critical-ARG recovery.

## Python API

```python
from plasrisk import annotate_fasta, get_scorer

ann = annotate_fasta("plasmid.fasta")
scorer = get_scorer(model="plasrisk", mode="lite")
df = scorer.score_dataframe(ann.features)
print(df[["seq_id", "S_norm", "grade"]])
```

Use the full model, the ordinal index, or the 2025 AWaRe table:

```python
scorer = get_scorer(model="plasrisk", mode="full")
scorer = get_scorer(model="plasrisk", mode="ordinal")
scorer = get_scorer(model="plasrisk", mode="lite", aware_version="2025")
```

Score a raw assembly or metagenomic contigs programmatically:

```python
from plasrisk import process_one, score_contigs, aggregate_bins, aggregate_sample

result = process_one("assembly.fasta", outdir="run_one",
                     classifier="genomad", genomad_db="/path/to/genomad_db")
print(result.status, result.reason)
```

## Output

The default command writes a per-sequence TSV with each component
score, the normalised total score, the rank and the A-E grade, together
with a per-file summary. Grades use prespecified thresholds: grade A,
S at or above 0.60; grade B, 0.45 to 0.60; grade C, 0.30 to 0.45;
grade D, 0.15 to 0.30; and grade E below 0.15. The thresholds are
priority categories rather than equal-width bins and are evaluated in
the Supplementary Materials against fixed-width, quantile and normal
binning schemes.

## Reproducible data

The locked data split, frozen weight vectors, frozen predictions and the
792,964 by ten component matrix are released under
[data/locked_v1.1.1](data/locked_v1.1.1), with SHA256 checksums. The WHO
AWaRe 2025 source table and derived mapping are under
[supplementary_data/aware2025](supplementary_data/aware2025).

## Citation

Please cite PlasRisk and the PIPdb resource. The manuscript citation
will be added on publication.

## License

PlasRisk is released under the MIT license.
