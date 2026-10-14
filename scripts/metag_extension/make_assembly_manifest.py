#!/usr/bin/env python3
"""Build a manifest of pre-assembled metagenome FASTA from ENA.

ENA (and MGnify "EMG produced TPA" assemblies) store assembled contigs as
SEQUENCE_ASSEMBLY analyses with ERZ accessions, one ``.fasta.gz`` / ``.fna.gz``
per sample. This script queries ENA for every assembly submitted under a
BioProject and writes a manifest that 04_preassembled_benchmark.sh consumes.

Example
-------
# MGnify TPA assembly of the global sewage AMR project (110 samples, total
# metagenome, metaSPAdes). The TPA assembly BioProject is PRJEB47494.
python3 make_assembly_manifest.py \
    --study PRJEB47494 --habitat wastewater --pair 0 \
    --out assemblies_wastewater.tsv

# Global sewage plasmidome (24 phi29-amplified plasmid-enriched assemblies).
python3 make_assembly_manifest.py \
    --study PRJEB41171 --habitat plasmidome --pair 0 \
    --out assemblies_plasmidome.tsv

The output TSV has four columns (no header):
    id <TAB> habitat_group <TAB> pair_id <TAB> url
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.parse
import urllib.request

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"


def _sanitize(text: str) -> str:
    text = (text or "").strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "sample"


def fetch_assemblies(study: str, limit: str = "0", retries: int = 3) -> list[dict]:
    params = {
        "result": "analysis",
        "query": f"study_accession={study}",
        "fields": "analysis_accession,analysis_type,sample_accession,"
                  "submitted_ftp,sample_alias",
        "format": "tsv",
        "limit": limit,
    }
    url = f"{ENA_API}?{urllib.parse.urlencode(params)}"
    raw = ""
    last_exc = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "plasrisk/meta"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                raw = resp.read().decode("utf-8", "replace")
            break
        except Exception as exc:  # noqa: BLE001 - transient ENA errors
            last_exc = exc
            if attempt < retries - 1:
                import time
                time.sleep(3 * (attempt + 1))
    if not raw:
        raise last_exc if last_exc else RuntimeError("empty ENA response")
    lines = [ln for ln in raw.splitlines() if ln.strip()]
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    rows = []
    for r in reader:
        if (r.get("analysis_type") or "").strip() != "SEQUENCE_ASSEMBLY":
            continue
        ftp = (r.get("submitted_ftp") or "").strip().split(";")[0].strip()
        if not ftp:
            continue
        if not ftp.startswith("http"):
            ftp = "https://" + ftp
        erz = (r.get("analysis_accession") or "").strip()
        alias = _sanitize(r.get("sample_alias") or erz)
        rows.append({
            "id": alias,
            "erz": erz,
            "sample": (r.get("sample_accession") or "").strip(),
            "url": ftp,
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True,
                    help="ENA BioProject that holds the SEQUENCE_ASSEMBLY "
                         "analyses (for MGnify TPA projects use the TPA "
                         "assembly BioProject, e.g. PRJEB47494).")
    ap.add_argument("--habitat", required=True,
                    help="habitat_group label written to the manifest, e.g. "
                         "wastewater, livestock, plasmidome.")
    ap.add_argument("--pair", default="0",
                    help="pair_id label written to the manifest.")
    ap.add_argument("--out", required=True, help="output manifest TSV path.")
    ap.add_argument("--id-prefix", default="",
                    help="optional prefix prepended to every sample id.")
    ap.add_argument("--limit", default="0",
                    help="ENA result limit (0 = all); use a finite number if "
                         "the unlimited query times out.")
    args = ap.parse_args()

    try:
        rows = fetch_assemblies(args.study, limit=args.limit)
    except Exception as exc:  # noqa: BLE001 - report the HTTP/parse error plainly
        sys.stderr.write(f"[manifest] ENA query failed for {args.study}: {exc}\n")
        return 1

    if not rows:
        sys.stderr.write(
            f"[manifest] no SEQUENCE_ASSEMBLY rows for {args.study}. "
            "The project has no deposited assemblies; assemble the reads "
            "instead.\n")
        return 2

    seen = set()
    n = 0
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for r in rows:
            sid = _sanitize(args.id_prefix + r["id"])
            if sid in seen:
                sid = f"{sid}_{r['erz']}"
            seen.add(sid)
            w.writerow([sid, args.habitat, args.pair, r["url"]])
            n += 1
    print(f"[manifest] {args.study}: {n} assemblies -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
