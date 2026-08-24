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
    --model MODEL        Scoring model: 'plasrisk' (10-dim weighted, default)
                         or 'pipdb' (original PIPdb 8-item ordinal)
    --mode MODE          PlasRisk mode: 'full' (10-dim, default) or 'lite'
                         (5-dim core; ignored when --model pipdb)
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
from .scoring import (PlasRiskScorer, PIPdbScorer, RISK_WEIGHTS,
                      RISK_WEIGHTS_LITE, WEIGHT_SUM, get_scorer)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="plasrisk",
        description="PlasRisk v%s - plasmid risk assessment (PlasRisk or PIPdb model)" % __version__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  plasrisk plasmid.fasta
  plasrisk -o results *.fasta
  plasrisk /data/plasmids/
  plasrisk --db card,vfdb,plasmidfinder plasmid.fasta
  plasrisk --no-abricate contigs.fasta
  plasrisk --model pipdb plasmid.fasta
  plasrisk --model plasrisk --mode lite plasmid.fasta
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
    parser.add_argument("--model", choices=["plasrisk", "pipdb"],
                        default="plasrisk",
                        help="Scoring model: 'plasrisk' (10-dim weighted, "
                             "default) or 'pipdb' (original PIPdb 8-item ordinal)")
    parser.add_argument("--mode", choices=["full", "lite"], default="full",
                        help="PlasRisk mode: 'full' (10-dim, default) or 'lite' "
                             "(5-dim core; ignored when --model pipdb)")
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
    if args.model == "pipdb":
        model_str = "PIPdb original (8-item ordinal)"
    else:
        model_str = "LITE (5-dim core)" if args.mode == "lite" else "FULL (10-dim)"
    print("=" * 68)
    print("  PlasRisk v%s - Plasmid Risk Assessment" % __version__)
    print("  Model: %s" % model_str)
    print("=" * 68)
    if not args.quiet:
        if args.model == "plasrisk":
            w = RISK_WEIGHTS_LITE if args.mode == "lite" else RISK_WEIGHTS
            ws = sum(w.values())
            print("  Weights: %s" % ", ".join(
                "%s=%.3f" % (k, v) for k, v in w.items()))
            print("  Weight sum: %.4f" % ws)
        else:
            print("  Formula: Round((phylum+species+habitats+ARGs+VFGs")
            print("           +2*WHO_ARGs+ISs+growth)/8 + 0.6)")
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

    # Load replicon lookup and create scorer
    lookup = load_replicon_lookup()
    scorer = get_scorer(model=args.model, mode=args.mode,
                        replicon_lookup=lookup)

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
            is_pipdb = args.model == "pipdb"
            score_col = "risk_index_normalized" if is_pipdb else "S_norm"
            grade_col = "grade"
            if is_pipdb:
                grade_list = ["5", "4", "3", "2", "1"]
                n_bm = 0  # PIPdb model does not report BM separately
                n_hr = 0  # PIPdb model does not report high_risk_genes
            else:
                grade_list = ["A", "B", "C", "D", "E"]
                n_bm = int((df["n_BM"] > 0).sum()) if "n_BM" in df.columns else 0
                n_hr = int((df["high_risk_genes"] != "").sum()) if "high_risk_genes" in df.columns else 0
            file_summaries.append({
                "file": os.path.basename(fasta_path),
                "model": args.model,
                "n_sequences": len(df),
                "n_ARG_positive": int((df["n_ARG"] > 0).sum()),
                "n_VF_positive": int((df["n_VF"] > 0).sum()),
                "n_BM_positive": n_bm,
                "n_high_risk": n_hr,
                "mean_score": round(df[score_col].mean(), 4),
                "max_score": round(df[score_col].max(), 4),
                **{f"n_grade_{g}": int((df[grade_col] == g).sum())
                   for g in grade_list},
                "databases_used": ";".join(ann.databases_used),
            })

        if not args.quiet:
            top = df.iloc[0]
            score_col = "risk_index_normalized" if args.model == "pipdb" else "S_norm"
            print("    %d sequences scored | top: %s (score=%.3f, grade %s)" % (
                len(df), top["seq_id"], top[score_col], top["grade"]))

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
    is_pipdb = args.model == "pipdb"
    if is_pipdb:
        grade_order = ["5", "4", "3", "2", "1"]
        score_col = "risk_index_normalized"
    else:
        grade_order = ["A", "B", "C", "D", "E"]
        score_col = "S_norm"
    grade_counts = combined["grade"].value_counts()
    for grade in grade_order:
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
            if is_pipdb:
                extra = "index=%s" % row.get("combined_risk_index", "?")
            else:
                hr = row.get("high_risk_genes", "")
                extra = "[%s]" % hr if hr else ""
            print("    %-40s score=%.3f  grade %s  %d ARG  %s  %s" % (
                str(row["seq_id"])[:40],
                row[score_col], row["grade"],
                row["n_ARG"], row.get("replicon", "?"),
                extra,
            ))

    print("\nDone. Output directory: %s" % os.path.abspath(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
