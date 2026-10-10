#!/usr/bin/env bash
# make_accessions.sh
# ---------------------------------------------------------------------------
# Expand curated BioProjects (projects.tsv) into run-level accessions for
# 02_public_metag_benchmark.sh. fasterq-dump needs run accessions (ERR/SRR),
# not project accessions (PRJEB/PRJNA).
#
# ENA projects (PRJEB/ERP) are expanded via the ENA Portal filereport API;
# NCBI projects (PRJNA/PRJ...) via the SRA runinfo endpoint. Only paired
# Illumina runs matching the requested platform are emitted.
#
# Input projects.tsv (header + tab columns):
#   project  habitat  pair  platform
# Output accessions.tsv (NO header, tab columns):
#   run_accession  habitat_group  pair_id
#
# Usage:
#   bash make_accessions.sh [projects.tsv] [accessions.tsv]
# ---------------------------------------------------------------------------
set -euo pipefail

IN="${1:-projects.tsv}"
OUT="${2:-accessions.tsv}"
python3 - "${IN}" "${OUT}" <<'PYEOF'
import csv, sys, urllib.request

in_tsv, out_tsv = sys.argv[1], sys.argv[2]


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=120) as fh:
        return fh.read().decode("utf-8", "replace")


def expand_ena(project, platform):
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s"
           "&result=read_run&fields=run_accession,instrument_platform,"
           "library_layout&format=tsv&download=false") % project
    txt = get(url)
    rows = []
    for r in csv.DictReader(txt.splitlines(), delimiter="\t"):
        if r.get("instrument_platform", "").upper() == platform \
                and r.get("library_layout", "").upper() == "PAIRED":
            rows.append(r["run_accession"])
    return rows


def expand_ncbi(project, platform):
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=%s" \
        % project
    txt = get(url)
    rows = []
    for r in csv.DictReader(txt.splitlines()):
        if r.get("Platform", "").upper() == platform \
                and (r.get("LibraryLayout", "").upper() == "PAIRED"):
            rows.append(r["Run"])
    return rows


seen, out_rows = set(), []
with open(in_tsv, encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        project = r["project"].strip()
        habitat = r["habitat"].strip()
        pair = r["pair"].strip()
        platform = r["platform"].strip().upper()
        if project.startswith(("PRJEB", "ERP", "PRJER")):
            runs = expand_ena(project, platform)
        else:
            runs = expand_ncbi(project, platform)
        print("  %-11s %-9s %d paired runs" % (project, habitat, len(runs)))
        for run in runs:
            if run not in seen:
                seen.add(run)
                out_rows.append((run, habitat, pair))

with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
    for run, habitat, pair in out_rows:
        fh.write("%s\t%s\t%s\n" % (run, habitat, pair))
print("wrote %s (%d runs)" % (out_tsv, len(out_rows)))
PYEOF

echo "Next: bash 02_public_metag_benchmark.sh ${OUT} <WORKDIR> <THREADS>"
