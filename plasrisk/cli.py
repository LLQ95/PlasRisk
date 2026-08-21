#!/usr/bin/env python3
"""
cli.py - Command-line interface for PlasRisk.

Usage:
    plasrisk [options] <fasta1> [fasta2 ...]
    plasrisk [options] /path/to/plasmids/
    plasrisk [options] *.fasta

Options:
    -o, --output DIR     Output directory (default: ./plasrisk_output)
    -t, --threads N      Number of abricate threads (default: 4)
    --min-id FLOAT       Minimum abricate identity % (default: 75)
    --min-cov FLOAT      Minimum abricate coverage % (default: 50)
    --no-abricate        Skip abricate; sequence-only scoring
    --db LIST            Comma-separated abricate databases to use
                         (default: auto-detect)
    --json               Also write JSON summary
    --tsv                Write per-sequence TSV (default)
    --summary            Write per-file summary TSV
    -q, --quiet          Suppress progress messages
    -h, --help           Show this help message
    -v, --version        Show version
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from typing import List

import pandas as pd

from . import __version__
from .annotate import (
    ABRICATE_DATABASES,
    abricate_available,
    abricate_databases,
    annotate_fasta,
)
from .lookup import load_replicon_lookup
from .scoring import PlasRiskScorer, RISK_WEIGHTS, WEIGHT_SUM


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="plasrisk",
        description="PlasRisk v%s - 10-dimension plasmid risk assessment" % __version__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  plasrisk plasmid.fasta
  plasrisk -o results *.fasta
  plasrisk /data/plasmids/
  plasrisk --db card,vfdb,plasmidfinder plasmid.fasta
  plasrisk --no-abricate contigs.fasta
        """,
    )
    parser.add_argument("inputs", nargs="+", metavar="FASTA",
                        help="FASTA file(s) or directory containing FASTA files")
    parser.add_argument("-o", "--output", default="plasrisk_output",
                        help="Output directory (default: ./plasrisk_output)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Number of threads (default: 4)")
    parser.add_argument("--min-id", type=float, default=75.0,
                        help="Minimum abricate identity %% (default: 75)")
    parser.add_argument("--min-cov", type=float, default=50.0,
                        help="Minimum abricate coverage %% (default: 50)")
    parser.add_argument("--no-abricate", action="store_true",
                        help="Skip abricate annotation; sequence-only scoring")
    parser.add_argument("--db", default=None,
                        help="Comma-separated abricate databases (default: auto)")
    parser.add_argument("--json", action="store_true",
                        help="Also write JSON output")
    parser.add_argument("--tsv", action="store_true", default=True,
                        help="Write per-sequence TSV (default: on)")
    parser.add_argument("--summary", action="store_true", default=True,
                        help="Write per-file summary TSV (default: on)")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress progress messages")
    parser.add_argument("-v", "--version", action="version",
                        version="PlasRisk %s" % __version__)
    return parser.parse_args(argv)


def collect_fasta_files(inputs: List[str]) -> List[str]:
    """Expand inputs (files, globs, directories) into a list of FASTA paths."""
    fasta_files = []
    fasta_exts = {".fasta", ".fa", ".fna", ".ffn", ".fsa", ".fas"}
    for inp in inputs:
        # Expand glob patterns
        matches = glob.glob(inp)
        if not matches and os.path.exists(inp):
            matches = [inp]
        for match in matches:
            if os.path.isdir(match):
                for root, _dirs, files in os.walk(match):
                    for f in files:
                        ext = os.path.splitext(f)[1].lower()
                        if ext in fasta_exts or ext in (".gz", ".bz2"):
                            fasta_files.append(os.path.join(root, f))
            elif os.path.isfile(match):
                fasta_files.append(match)
    # Deduplicate while preserving order
    seen = set()
    unique = []
    for f in fasta_files:
        ap = os.path.abspath(f)
        if ap not in seen:
            seen.add(ap)
            unique.append(f)
    return unique


def print_banner(args):
    print("=" * 68)
    print("  PlasRisk v%s - Plasmid Risk Assessment" % __version__)
    print("  10-dimension weighted risk model")
    print("=" * 68)
    if not args.quiet:
        print("  Weights: %s" % ", ".join(
            "%s=%.2f" % (k, v) for k, v in RISK_WEIGHTS.items()))
        print("  Weight sum: %.2f (normalized scores = S/%.2f)" % (WEIGHT_SUM, WEIGHT_SUM))
        if not args.no_abricate:
            if abricate_available():
                dbs = abricate_databases()
                print("  abricate: found (%d databases: %s)" % (
                    len(dbs), ", ".join(dbs[:6]) + ("..." if len(dbs) > 6 else "")))
            else:
                print("  abricate: NOT FOUND (install with: "
                      "conda install -c bioconda abricate)")
        else:
            print("  abricate: disabled (--no-abricate)")
        print("-" * 68)


def main(argv=None):
    args = parse_args(argv)
    print_banner(args)

    # Collect input files
    fasta_files = collect_fasta_files(args.inputs)
    if not fasta_files:
        print("ERROR: No FASTA files found.", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print("  Found %d FASTA file(s)" % len(fasta_files))

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Load replicon lookup
    lookup = load_replicon_lookup()
    scorer = PlasRiskScorer(replicon_lookup=lookup)

    # Determine abricate databases
    abricate_dbs = None
    if args.db:
        abricate_dbs = [d.strip() for d in args.db.split(",")]

    # Process each file
    all_results = []
    file_summaries = []
    total_seqs = 0
    t0 = time.time()

    for i, fasta_path in enumerate(fasta_files, 1):
        if not args.quiet:
            print("\n  [%d/%d] %s" % (i, len(fasta_files), os.path.basename(fasta_path)))

        try:
            ann = annotate_fasta(
                fasta_path,
                use_abricate=not args.no_abricate,
                abricate_dbs=abricate_dbs,
                min_id=args.min_id,
                min_cov=args.min_cov,
                threads=args.threads,
                lookup=lookup,
            )
        except Exception as exc:
            print("    ERROR: %s" % exc, file=sys.stderr)
            continue

        if ann.warnings and not args.quiet:
            for w in ann.warnings:
                print("    WARNING: %s" % w)

        if not ann.features:
            print("    No sequences processed.")
            continue

        # Score
        df = scorer.score_dataframe(ann.features)
        df.insert(0, "file", os.path.basename(fasta_path))
        all_results.append(df)
        total_seqs += len(df)

        # Per-file summary
        if len(df) > 0:
            file_summaries.append({
                "file": os.path.basename(fasta_path),
                "n_sequences": len(df),
                "n_ARG_positive": int((df["n_ARG"] > 0).sum()),
                "n_VF_positive": int((df["n_VF"] > 0).sum()),
                "n_BM_positive": int((df["n_BM"] > 0).sum()),
                "n_high_risk": int((df["high_risk_genes"] != "").sum()),
                "mean_S_norm": round(df["S_norm"].mean(), 4),
                "max_S_norm": round(df["S_norm"].max(), 4),
                "n_grade_A": int((df["grade"] == "A").sum()),
                "n_grade_B": int((df["grade"] == "B").sum()),
                "n_grade_C": int((df["grade"] == "C").sum()),
                "n_grade_D": int((df["grade"] == "D").sum()),
                "n_grade_E": int((df["grade"] == "E").sum()),
                "databases_used": ";".join(ann.databases_used),
            })

        if not args.quiet:
            top = df.iloc[0]
            print("    %d sequences scored | top: %s (S=%.3f, grade %s)" % (
                len(df), top["seq_id"], top["S_norm"], top["grade"]))

    if not all_results:
        print("\nERROR: No results generated.", file=sys.stderr)
        sys.exit(1)

    # Combine all results
    combined = pd.concat(all_results, ignore_index=True)

    # Write outputs
    tsv_path = os.path.join(args.output, "plasrisk_results.tsv")
    combined.to_csv(tsv_path, sep="\t", index=False, float_format="%.4f")
    if not args.quiet:
        print("\n  Results written to: %s" % tsv_path)

    if file_summaries:
        summary_df = pd.DataFrame(file_summaries)
        summary_path = os.path.join(args.output, "plasrisk_summary.tsv")
        summary_df.to_csv(summary_path, sep="\t", index=False, float_format="%.4f")
        if not args.quiet:
            print("  Summary written to: %s" % summary_path)

    if args.json:
        json_path = os.path.join(args.output, "plasrisk_results.json")
        records = combined.to_dict(orient="records")
        with open(json_path, "w") as jf:
            json.dump(records, jf, indent=2, default=str)
        if not args.quiet:
            print("  JSON written to: %s" % json_path)

    # Print final summary table
    elapsed = time.time() - t0
    print("\n" + "=" * 68)
    print("  SUMMARY")
    print("=" * 68)
    print("  Files processed:    %d" % len(fasta_files))
    print("  Sequences scored:   %d" % total_seqs)
    print("  Time elapsed:       %.1f sec" % elapsed)
    print("-" * 68)

    # Grade distribution
    grade_counts = combined["grade"].value_counts().sort_index()
    for grade in ["A", "B", "C", "D", "E"]:
        n = grade_counts.get(grade, 0)
        pct = 100 * n / total_seqs if total_seqs > 0 else 0
        bar = "#" * int(pct / 2)
        print("  Grade %s: %6d (%5.1f%%) %s" % (grade, n, pct, bar))
    print("-" * 68)

    # Top 10 highest risk
    if len(combined) > 0:
        print("\n  Top 10 highest-risk plasmids:")
        top10 = combined.nsmallest(10, "rank") if "rank" in combined.columns else combined.head(10)
        for _, row in top10.head(10).iterrows():
            hr = row.get("high_risk_genes", "")
            hr_str = " [%s]" % hr if hr else ""
            print("    %-40s S=%.3f  grade %s  %d ARG  %s%s" % (
                str(row["seq_id"])[:40],
                row["S_norm"], row["grade"],
                row["n_ARG"], row.get("replicon", "?"),
                hr_str,
            ))

    print("\nDone. Output directory: %s" % os.path.abspath(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
