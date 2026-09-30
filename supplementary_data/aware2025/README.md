# WHO AWaRe 2025 supplementary data

This folder contains the source data and audit outputs used for the S_ARG
AWaRe classification check reported in Supplementary Method S14 and Table S56
of the PlasRisk manuscript, so that reviewers can reproduce every number.

## Contents

| File | Description |
| --- | --- |
| `B09489-eng.xlsx` | Original WHO workbook (WHO AWaRe classification 2025, downloaded from the WHO website), uploaded alongside these files. |
| `AWaRe_classification_2025.csv` | Main sheet, all 272 listed antibiotics with Class, ATC code, Category and EML/EMLc notes. |
| `AWaRe_Access_2025.csv`, `AWaRe_Watch_2025.csv`, `AWaRe_Reserve_2025.csv` | Antibiotics listed under each AWaRe category (97 Access, 149 Watch, 34 Reserve). |
| `Not_recommended_2025.csv` | Agents not recommended for the listed indications. |
| `EML_24_2025.csv`, `EMLc_10_2025.csv` | Antibiotics on the 24th WHO Model List of Essential Medicines and 10th EML for Children. |
| `Introduction.csv` | WHO introductory notes. |
| `tab_aware_rule_discrepancies.csv` | Class-level comparison between the legacy regex mapping and the WHO AWaRe 2025 table. |
| `tab_aware_mapping_crosstab.csv` | PSC counts cross-tabulated between legacy and 2025 categories. |
| `sarg_aware_audit.py` | Script that recomputes S_ARG under both mappings and the audit statistics. |

## Reproducing the audit

Run `sarg_aware_audit.py` against the released PSC master table. It writes the
per-PSC audit table (`tab_sarg_aware_audit.csv`, 792,964 rows) and reports
that the 2025 mapping changes S_ARG by at least 0.05 for 4.9% of PSCs and moves
the high-risk-ARG AUC on the locked test from 0.9807 to 0.9789 (delta
-0.0018). The frozen legacy values remain the default; the 2025 mapping is
available in the package with `--aware-version 2025`.
