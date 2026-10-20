# Presubmission Server Scripts (P0 Checklist Pack)

This folder contains the scripts that run on the compute server
(`/db/student/metagenome/pipdb/contigs`, conda env `genomad`) for the
presubmission checklist item P0-1 and P0-3. Every script is resume-safe:
if a run is interrupted, rerun the same command and completed steps are
skipped automatically.

## Run order and resource estimates

| Script | Item | Purpose | Estimated time |
| --- | --- | --- | --- |
| `15_cef_resume_missing.sh` | P0-5 | Resume the four missing 07 cefprozil samples (P1 time points and P10E0), rerun the collection and the 08 paired analysis (completed; 72 samples) | 2-4 h |
| `12_longread_truth_plasmidome.sh` | P0-1 | Long-read plasmidome truth set for 24 samples (Kirstahler et al., PRJEB41171): ERZ long-read assemblies, matched Illumina reads re-assembled with MEGAHIT, minimap2 recovery funnel, geNomad/PlasRisk on both sides | Overnight (12-24 h, ~120 Gbp) |
| `13_ablation_genomad_only.py` | P0-3 | geNomad-only versus geNomad+PlasRisk ablation across the three cohorts; reads existing results, no reassembly | Minutes |
| `14a_setup_envs.sh` | P0-3 | One-time, idempotent setup of `rfplasmid_cmp` and `plasmer_cmp` environments and databases (CheckM DB, two Plasmer DBs, ~14 GB) | 0.5-3 h, mostly download |
| `14b_run_cross_tool.sh` | P0-3 | Run RFPlasmid (Enterobacteriaceae model) and Plasmer on the 189 isolates; four-tool comparison tables; RFP and Plasmer stages run concurrently | Hours to overnight |

## Commands

```bash
cd /db/student/metagenome/pipdb/contigs
conda activate genomad
export GENOMAD_DB=/path/to/genomad_db

# P0-5 (completed; rerun only if samples are still missing)
bash 15_cef_resume_missing.sh \
  cefprozil_tpa_chunks_manifest.tsv results/cefprozil 16 8 \
  ./07_cefprozil_tpa_benchmark.sh ./08_paired_analysis.py

# P0-1 (the manifest sits next to the script; use the ./ relative path)
bash 12_longread_truth_plasmidome.sh \
  ./assemblies_plasmidome.tsv results/plasmidome_truth 16 2

# P0-3 ablation (independent of scripts 12 and 14; requires 04/06/07 outputs)
python3 13_ablation_genomad_only.py \
  --isolates results/traveler_isolates \
  --sewage   results/preassembled \
  --cef      results/cefprozil \
  --out      results/ablation

# P0-3 cross-tool setup (run once; safe to interrupt and rerun)
bash 14a_setup_envs.sh results/traveler_isolates
# An existing rfplasmid_cmp is skipped. RFPlasmid --initialize prints
# "installed. Please try rfplasmid again" and exits 1 by design; readiness
# is judged by plasmiddb_cge.dmnd and plasmiddb_total.dmnd. 14a also
# patches RFPlasmid 0.0.18 for pandas >= 2.2 (str.split / delim_whitespace).
# Then run the comparison in the genomad environment
bash 14b_run_cross_tool.sh results/traveler_isolates 16 4 1
```

## Outputs to sync back to the local machine

1. `results/cefprozil/benchmark_sample_summary.tsv` and
   `results/cefprozil/paired/*` update all cefprozil numbers in the
   manuscript and SI. With 72 samples collected, the cohort text uses 72
   samples and 24 subjects; if samples were still missing, the 68-sample
   wording would be retained with the failure documented in Methods.

2. `results/plasmidome_truth/truth_summary.txt`, `truth_funnel.tsv`,
   `truth_recovery_contigs.tsv`, and `truth_long_grades.tsv` supply the
   numbers for the long-read truth-set statement in Methods (recovery
   funnel, length stratification, grade distribution of the truth
   molecules). Zero events for grade A or B are reported with the exact
   zero-event upper bound, not as a zero rate.

3. The tables under `results/ablation/` feed the SI cross-tool ablation
   section (Tables S70 onward).

4. `results/traveler_isolates/cross_tool/tool_comparison_*.tsv`
   supplements Table S53 (literature matrix) and updates the
   method-comparison text and the SI tables.

## Notes

* The geNomad version is matched automatically with `--genomad-version
  v1.12`, consistent with the prospective cohorts.

* The official Plasmer repository is `nekokoe/Plasmer` (the old
  `tiagolbiotech/Plasmer` URL no longer exists). The conda package is on
  the `iskoldt` channel (23.04.20); the databases come from Zenodo
  record 7030675 (plasmerMainDB plus customizedKraken2DB; the URLs and
  sha1 checksums are built into 14a). A fair comparison requires `-l 0`;
  success is judged by `predClass.tsv`, not by the exit code.

* Downloads from Zenodo use `aria2c` when available (16 connections,
  unlimited retries, resume); an institutional proxy can be supplied with
  `export DOWNLOAD_PROXY=http://host:port`.

* In RFPlasmid 0.0.18, `--threads` is parsed but not passed to CheckM or
  jellyfish (CheckM is hard-coded to 16), and diamond uses all cores by
  default, so RFP concurrency should not exceed 4. Each Plasmer task
  needs roughly 32 GB of memory, so the default Plasmer concurrency is 1.

* In script 12, the circular-molecule check combines a header regular
  expression with minimap2 self-alignment of the terminal overlap
  (at least 200 bp and identity at least 0.90). If the authors' figshare
  Table C is downloadable, their labels take precedence; the decision
  procedure is stated in the SI.

* Script 12 assembles the matched short reads with MEGAHIT rather than
  reproducing the authors' single-Nanopore-read workflow. The purpose is
  to test whether the frozen short-read pipeline recovers long-read truth
  molecules, not to reproduce the authors' assembly.
