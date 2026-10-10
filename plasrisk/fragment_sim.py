"""
fragment_sim.py - Controlled fragmentation grid for metagenomic validation.

Analogue of the simulation design used by parameter-free metagenomic
classifiers (e.g. SCiMS, Microbiome 2026): complete reference plasmids are
fragmented in silico across a grid of fragment lengths and random phases,
and three quantities are compared against the intact-plasmid score:

  S_bin_oracle    fragments of one plasmid are assigned to the same bin
                  (oracle binning, an upper bound on bin recovery),
  S_contig_mean   length-weighted mean of the best assessable per-contig
                  score (naive, no binning),
  S_whole         standard lite score on the intact plasmid.

The script reports Spearman rank correlation, grade exact / within-one
agreement, median score bias and critical-ARG recovery per fragment length.
It does not model assembly chimerism; pair it with ART/wgsim + metaSPAdes +
geNomad for the realistic layer (see the accompanying analysis plan).
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .annotate import AnnotationResult, annotate_fasta, parse_fasta
from .metagenome import aggregate_bin, score_contig, score_contigs
from .scoring import PlasRiskScorer, PlasmidFeatures

GRADE_ORDER = {"A": 4, "B": 3, "C": 2, "D": 1, "E": 0}
DEFAULT_GRID = (500, 1000, 2000, 5000, 10000, 20000)


def tile_fragments(seq: str, frag_len: int, rng: np.random.Generator
                   ) -> List[Tuple[int, int, str]]:
    """Tile a sequence into fragments with a random phase offset."""
    n = len(seq)
    if frag_len >= n:
        return [(0, n, seq.upper())]
    shift = int(rng.integers(0, frag_len))
    out = []
    start = shift
    if shift > 0:
        out.append((0, shift, seq[:shift].upper()))
    while start < n:
        end = min(start + frag_len, n)
        out.append((start, end, seq[start:end].upper()))
        start = end
    return out


def _best_contig_score(row) -> Optional[float]:
    for col in ("S_lite5", "S_frag4", "S_cargo3"):
        v = row.get(col)
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            return float(v)
    return None


def run_simulation(fasta_path: str,
                   frag_lengths: Tuple[int, ...] = DEFAULT_GRID,
                   n_replicates: int = 10,
                   min_fragment: int = 500,
                   complete_length: int = 10_000,
                   threads: int = 4,
                   use_abricate: bool = True,
                   db_path: Optional[str] = None,
                   aware_version: str = "legacy",
                   seed: int = 42,
                   quiet: bool = False
                   ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    scorer = PlasRiskScorer(replicon_lookup=None, mode="lite",
                            aware_version=aware_version)

    # ---- intact reference scores ----------------------------------------
    if use_abricate:
        whole_ann = annotate_fasta(fasta_path, db_path=db_path,
                                   threads=threads, quiet=True)
    else:
        recs = parse_fasta(fasta_path)
        whole_ann = AnnotationResult(
            file=fasta_path,
            features=[PlasmidFeatures(seq_id=i, length_bp=len(s))
                      for i, s in recs.items()])
    whole_df = scorer.score_dataframe(whole_ann.features)
    whole = {r["seq_id"]: r for _, r in whole_df.iterrows()}

    element_rows, summary_rows = [], []
    for frag_len in frag_lengths:
        for rep in range(n_replicates):
            # Build fragmented features (annotated in one batch when abricate
            # is available, else sequence-only).
            frag_records = []  # (plasmid_id, index, sequence)
            frag_index = []    # (plasmid_id, index, start, end)
            for pid, seq in parse_fasta(fasta_path).items():
                pieces = tile_fragments(seq, frag_len, rng)
                k = 0
                for (a, b, s) in pieces:
                    if len(s) < min_fragment:
                        continue
                    frag_records.append(("%s__f%05d" % (pid, k), s))
                    frag_index.append((pid, k, a, b))
                    k += 1

            if use_abricate:
                import tempfile
                with tempfile.TemporaryDirectory() as td:
                    fp = os.path.join(td, "fragments.fasta")
                    with open(fp, "w") as fh:
                        for fid, s in frag_records:
                            fh.write(">%s\n%s\n" % (fid, s))
                    ann = annotate_fasta(fp, db_path=db_path,
                                         threads=threads, quiet=True)
                feats = {f.seq_id: f for f in ann.features}
            else:
                feats = {fid: PlasmidFeatures(seq_id=fid, length_bp=len(s))
                         for fid, s in frag_records}

            by_plasmid: Dict[str, List[PlasmidFeatures]] = {}
            for fid, feat in feats.items():
                pid = fid.rsplit("__f", 1)[0]
                by_plasmid.setdefault(pid, []).append(feat)

            for pid, flist in by_plasmid.items():
                if pid not in whole:
                    continue
                w = whole[pid]
                bin_row = aggregate_bin(
                    pid, flist, scorer, circular=True,
                    plasmid_score=None, coverage=None)
                contig_ann = AnnotationResult(file=fasta_path, features=flist)
                cdf = score_contigs(contig_ann, scorer,
                                    complete_length=complete_length,
                                    min_length=min_fragment)
                scores = cdf.apply(_best_contig_score, axis=1).dropna()
                lens = cdf.loc[scores.index, "length_bp"].to_numpy(float)
                s_contig = (float(np.sum(scores.to_numpy(float) * lens)
                                  / np.sum(lens))
                            if len(scores) and lens.sum() > 0 else np.nan)
                hr_whole = str(w.get("high_risk_genes", "") or "")
                hr_bin = bin_row["high_risk_genes"]
                element_rows.append({
                    "plasmid_id": pid,
                    "frag_len_bp": frag_len,
                    "replicate": rep + 1,
                    "n_fragments": len(flist),
                    "length_whole_bp": int(w["length_bp"]),
                    "S_whole": round(float(w["S_norm"]), 4),
                    "grade_whole": w["grade"],
                    "S_bin_oracle": bin_row["S_lite"],
                    "grade_bin": bin_row["grade"],
                    "delta_bin": (round(bin_row["S_lite"]
                                        - float(w["S_norm"]), 4)
                                  if bin_row["S_lite"] is not None else np.nan),
                    "S_contig_mean": round(s_contig, 4)
                    if not np.isnan(s_contig) else np.nan,
                    "high_risk_recovered": int(
                        hr_whole != "" and hr_bin != ""),
                    "high_risk_present_whole": int(hr_whole != ""),
                })

    elem = pd.DataFrame(element_rows)
    if elem.empty:
        return elem, pd.DataFrame()

    for frag_len, grp in elem.groupby("frag_len_bp"):
        graded = grp[grp["grade_bin"] != ""]
        if len(graded) >= 3 and graded["grade_whole"].nunique() > 1:
            rho = float(graded["S_bin_oracle"].corr(
                graded["S_whole"], method="spearman"))
        else:
            rho = np.nan
        gw = graded["grade_whole"].map(GRADE_ORDER)
        gb = graded["grade_bin"].map(GRADE_ORDER)
        exact = float((gw == gb).mean()) if len(graded) else np.nan
        within1 = float((np.abs(gw - gb) <= 1).mean()) if len(graded) \
            else np.nan
        hr = grp[grp["high_risk_present_whole"] == 1]
        summary_rows.append({
            "frag_len_bp": frag_len,
            "n_comparisons": len(grp),
            "spearman_rho_bin_vs_whole": round(rho, 4)
            if not np.isnan(rho) else np.nan,
            "grade_exact_agreement": round(exact, 4)
            if not np.isnan(exact) else np.nan,
            "grade_within_one": round(within1, 4)
            if not np.isnan(within1) else np.nan,
            "median_delta_bin": round(
                float(grp["delta_bin"].median()), 4),
            "iqr_delta_bin": "%.4f-%.4f" % (
                float(grp["delta_bin"].quantile(0.25)),
                float(grp["delta_bin"].quantile(0.75))),
            "critical_ARG_recovery": round(
                float(hr["high_risk_recovered"].mean()), 4)
            if len(hr) else np.nan,
            "mean_n_fragments": round(float(grp["n_fragments"].mean()), 2),
        })
    summary = pd.DataFrame(summary_rows)
    if not quiet:
        print(summary.to_string(index=False))
    return elem, summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plasrisk-fragment-sim",
        description="In-silico fragmentation grid for PlasRisk metagenomic "
                    "validation (oracle-bin recovery vs intact plasmids).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("fasta", help="complete reference plasmid FASTA")
    p.add_argument("-o", "--output", default="fragment_sim",
                   help="output directory")
    p.add_argument("--frag-lengths", default=",".join(map(str, DEFAULT_GRID)),
                   help="comma-separated fragment lengths (bp)")
    p.add_argument("--replicates", type=int, default=10)
    p.add_argument("--min-fragment", type=int, default=500)
    p.add_argument("--complete-length", type=int, default=10000)
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--no-abricate", action="store_true",
                   help="sequence-only plumbing test (cargo dims = 0)")
    p.add_argument("--db", default=None)
    p.add_argument("--aware-version", choices=["legacy", "2025"],
                   default="legacy")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("-q", "--quiet", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    grid = tuple(int(x) for x in args.frag_lengths.split(",") if x.strip())
    elem, summary = run_simulation(
        args.fasta, frag_lengths=grid, n_replicates=args.replicates,
        min_fragment=args.min_fragment,
        complete_length=args.complete_length, threads=args.threads,
        use_abricate=not args.no_abricate, db_path=args.db,
        aware_version=args.aware_version, seed=args.seed,
        quiet=args.quiet)
    os.makedirs(args.output, exist_ok=True)
    ep = os.path.join(args.output, "fragment_sim_per_element.tsv")
    sp = os.path.join(args.output, "fragment_sim_summary.tsv")
    elem.to_csv(ep, sep="\t", index=False)
    summary.to_csv(sp, sep="\t", index=False)
    if not args.quiet:
        print("wrote %s\nwrote %s" % (ep, sp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
