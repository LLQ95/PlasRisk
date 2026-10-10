"""
fragment_sim.py - Controlled fragmentation recovery simulation.

Intact, circular reference plasmids are tiled into fragments of varying
length (with random phase offsets), and each fragment set is scored with
the metagenomic Layer-1 pipeline. An oracle bin (all fragments of one
source plasmid grouped together, circular=True) is compared with:

    S_whole       intact plasmid, standard lite score
    S_bin_oracle  oracle bin aggregation
    S_contig_mean length-weighted mean of the best assessable fragment
                  score (S_lite5 > S_frag4 > S_cargo3)

Outputs report, per fragment length: Spearman rank correlation with the
intact-plasmid scores, grade exact / within-one agreement, median score
delta (IQR), critical-ARG recovery, and mean fragments per plasmid.
"""

from __future__ import annotations

import os
import sys
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .annotate import (AnnotationResult, annotate_fasta, parse_fasta)
from .metagenome import (aggregate_bin, score_contigs, circular_map_from_fasta)
from .scoring import PlasRiskScorer, PlasmidFeatures


def tile_fragments(seq_id: str, sequence: str, frag_length: int,
                   rng: np.random.Generator) -> List[PlasmidFeatures]:
    """Tile one sequence into fragments with a random start phase."""
    L = len(sequence)
    if frag_length >= L:
        return [PlasmidFeatures(seq_id=seq_id, length_bp=L)]
    offset = int(rng.integers(0, frag_length))
    bounds = list(range(offset, L, frag_length))
    if not bounds or bounds[0] != 0:
        bounds = [0] + bounds
    if bounds[-1] != L:
        bounds = bounds + [L]
    out = []
    for i in range(len(bounds) - 1):
        a, b = bounds[i], bounds[i + 1]
        if b - a < 1:
            continue
        out.append(PlasmidFeatures(
            seq_id="%s_frag%d" % (seq_id, i + 1),
            length_bp=b - a))
    return out


def _best_contig_score(row: Dict) -> Optional[float]:
    for key in ("S_lite5", "S_frag4", "S_cargo3"):
        v = row.get(key)
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            return float(v)
    return None


def run_simulation(fasta_path: str,
                   frag_lengths: Sequence[int] = (500, 1000, 2000, 5000,
                                                  10000, 20000),
                   replicates: int = 10,
                   use_abricate: bool = True,
                   threads: int = 4,
                   seed: int = 42) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run the fragmentation grid; return per-element and summary tables."""
    rng = np.random.default_rng(seed)
    records = parse_fasta(fasta_path)
    scorer = PlasRiskScorer(mode="lite")

    # Intact-plasmid reference scores
    whole_rows = []
    for sid, seq in records.items():
        whole_rows.append({
            "seq_id": sid,
            "length_bp": len(seq),
            "S_whole": scorer.score_size(len(seq)) and None,  # placeholder
        })
    # Score intact plasmids through the standard annotation path
    if use_abricate:
        whole_ann = annotate_fasta(fasta_path, threads=threads, quiet=True)
    else:
        feats = [PlasmidFeatures(seq_id=sid, length_bp=len(seq))
                 for sid, seq in records.items()]
        whole_ann = AnnotationResult(file=fasta_path, features=feats)
    whole_df = score_contigs(whole_ann, scorer, min_length=1,
                             circular_map={sid: True for sid in records})
    whole_lookup = {}
    for _, r in whole_df.iterrows():
        whole_lookup[r["seq_id"]] = {
            "S": r["S_lite5"] if r["S_lite5"] is not None else
            _best_contig_score(r.to_dict()),
            "grade": r["grade"],
            "hr": r["high_risk_genes"],
        }

    per_element, summary = [], []
    for frag_length in frag_lengths:
        for rep in range(replicates):
            all_frag_feats: List[PlasmidFeatures] = []
            frag_to_source: Dict[str, str] = {}
            source_frags: Dict[str, List[str]] = {}
            for sid, seq in records.items():
                frags = tile_fragments(sid, seq, frag_length, rng)
                source_frags[sid] = [f.seq_id for f in frags]
                for f in frags:
                    frag_to_source[f.seq_id] = sid
                    all_frag_feats.append(f)

            if use_abricate:
                tmpdir = tempfile.mkdtemp(prefix="fragsim_")
                fpath = os.path.join(tmpdir, "fragments.fasta")
                with open(fpath, "w") as out:
                    for sid, seq in records.items():
                        for f in tile_fragments(sid, seq, frag_length, rng):
                            sub = seq  # re-tile below for sequence content
                # Rebuild fragments with actual sequence content for abricate
                with open(fpath, "w") as out:
                    for sid, seq in records.items():
                        L = len(seq)
                        if frag_length >= L:
                            out.write(">%s\n%s\n" % (sid, seq))
                            continue
                        offset = int(rng.integers(0, frag_length))
                        bounds = list(range(offset, L, frag_length))
                        if not bounds or bounds[0] != 0:
                            bounds = [0] + bounds
                        if bounds[-1] != L:
                            bounds = bounds + [L]
                        for i in range(len(bounds) - 1):
                            a, b = bounds[i], bounds[i + 1]
                            if b - a < 1:
                                continue
                            fid = "%s_frag%d" % (sid, i + 1)
                            out.write(">%s\n%s\n" % (fid, seq[a:b]))
                ann = annotate_fasta(fpath, threads=threads, quiet=True)
                contig_df = score_contigs(ann, scorer, min_length=1)
            else:
                ann = AnnotationResult(file="<fragments>",
                                       features=all_frag_feats)
                contig_df = score_contigs(ann, scorer, min_length=1)

            contig_lookup = {r["seq_id"]: r for _, r in
                             contig_df.iterrows()}

            # Oracle bins: all fragments of each source plasmid, circular
            oracle_rows = []
            for sid in records:
                members = [f for f in all_frag_feats
                           if frag_to_source.get(f.seq_id) == sid]
                if use_abricate:
                    members = [f for f in ann.features
                               if frag_to_source.get(f.seq_id) == sid]
                if not members:
                    continue
                row = aggregate_bin(sid, members, scorer,
                                    circular=True, plasmid_score=None,
                                    coverages=None)
                oracle_rows.append(row)

            oracle_lookup = {r["bin_id"]: r for r in oracle_rows}

            for sid in records:
                ref = whole_lookup.get(sid)
                if ref is None or ref["S"] is None:
                    continue
                frag_ids = [fid for fid in frag_to_source
                            if frag_to_source[fid] == sid]
                # length-weighted best assessable fragment score
                num = den = 0.0
                hr_recovered = False
                grades = []
                for fid in frag_ids:
                    cr = contig_lookup.get(fid)
                    if cr is None:
                        continue
                    bs = _best_contig_score(cr.to_dict())
                    L = int(cr["length_bp"])
                    if bs is not None:
                        num += bs * L
                        den += L
                    if cr["high_risk_genes"]:
                        hr_recovered = True
                    if cr["grade"]:
                        grades.append(cr["grade"])
                s_contig = num / den if den else np.nan
                s_oracle = oracle_lookup.get(sid, {}).get("S_lite")
                oracle_grade = oracle_lookup.get(sid, {}).get("grade", "")
                per_element.append({
                    "frag_length": frag_length,
                    "replicate": rep,
                    "seq_id": sid,
                    "length_bp": len(records[sid]),
                    "S_whole": ref["S"],
                    "grade_whole": ref["grade"],
                    "S_bin_oracle": s_oracle,
                    "grade_oracle": oracle_grade,
                    "S_contig_mean": s_contig,
                    "n_fragments": len(frag_ids),
                    "hr_in_whole": bool(ref["hr"]),
                    "hr_recovered": hr_recovered,
                })

        block = pd.DataFrame([e for e in per_element
                              if e["frag_length"] == frag_length])
        if block.empty:
            continue
        rho = block[["S_whole", "S_contig_mean"]].corr(
            method="spearman").iloc[0, 1]
        rho_oracle = block[["S_whole", "S_bin_oracle"]].corr(
            method="spearman").iloc[0, 1]
        grade_exact = np.mean(
            block["grade_whole"] == block["grade_oracle"])
        grade_within_one = np.mean([
            abs("ABCDE".index(a) - "ABCDE".index(b)) <= 1
            for a, b in zip(block["grade_whole"], block["grade_oracle"])
            if a and b]) if (block["grade_whole"] != "").any() else np.nan
        delta = block["S_contig_mean"] - block["S_whole"]
        hr_blocks = block[block["hr_in_whole"]]
        hr_recovery = hr_blocks["hr_recovered"].mean() if len(hr_blocks) else np.nan
        summary.append({
            "frag_length": frag_length,
            "n_elements": len(block),
            "spearman_contig": round(rho, 4),
            "spearman_oracle_bin": round(rho_oracle, 4),
            "grade_exact_oracle": round(grade_exact, 4),
            "grade_within_one_oracle": round(grade_within_one, 4)
            if not np.isnan(grade_within_one) else np.nan,
            "median_delta_contig": round(float(delta.median()), 4),
            "iqr_delta_contig": round(float(delta.quantile(0.75)
                                           - delta.quantile(0.25)), 4),
            "critical_ARG_recovery": round(hr_recovery, 4)
            if not np.isnan(hr_recovery) else np.nan,
            "mean_fragments": round(block["n_fragments"].mean(), 2),
        })

    return pd.DataFrame(per_element), pd.DataFrame(summary)


def build_parser():
    import argparse
    p = argparse.ArgumentParser(
        prog="plasrisk-fragment-sim",
        description="Fragmentation recovery simulation for PlasRisk "
                    "metagenomic scoring.")
    p.add_argument("fasta", help="intact reference plasmid FASTA")
    p.add_argument("-o", "--output", default="fragment_sim",
                   help="output directory")
    p.add_argument("--frag-lengths", default="500,1000,2000,5000,10000,20000",
                   help="comma-separated fragment lengths (bp)")
    p.add_argument("--replicates", type=int, default=10)
    p.add_argument("--no-abricate", action="store_true",
                   help="sequence-only pipeline (no gene annotation)")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=42)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    frag_lengths = [int(x) for x in args.frag_lengths.split(",") if x]
    per_element, summary = run_simulation(
        args.fasta, frag_lengths=frag_lengths,
        replicates=args.replicates, use_abricate=not args.no_abricate,
        threads=args.threads, seed=args.seed)
    os.makedirs(args.output, exist_ok=True)
    p1 = os.path.join(args.output, "fragment_sim_per_element.tsv")
    p2 = os.path.join(args.output, "fragment_sim_summary.tsv")
    per_element.to_csv(p1, sep="\t", index=False)
    summary.to_csv(p2, sep="\t", index=False)
    print("Saved %s and %s" % (p1, p2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
