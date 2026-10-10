#!/usr/bin/env python3
"""
metagenome.py - Metagenomic extension for PlasRisk.

Three-layer scoring for fragmented plasmid data:
  Layer 1 (contig): score individual plasmid contigs on the dimensions
      that remain identifiable after fragmentation.
  Layer 2 (bin):    aggregate contigs into predicted plasmid bins.
  Layer 3 (sample): descriptive, abundance-aware community statistics.

The frozen consensus weights are never refit on metagenomic data. When a
dimension cannot be assessed from a fragment, it is reported as missing
(None) rather than imputed from isolate-derived priors, and the lite score
is renormalised over the assessable dimensions.

CLI:
    plasrisk meta <inputs> -o <outdir>
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
from typing import Dict, List, Optional, Sequence

import pandas as pd

from .annotate import annotate_fasta
from .scoring import (
    PlasmidFeatures,
    PlasRiskScorer,
    RISK_WEIGHTS,
    RISK_WEIGHTS_LITE,
)

# Dimensions that can be counted directly from a fragment regardless of
# completeness.
CARGO_DIMS = ("S_ARG", "S_BMG", "S_VF")

CIRCULAR_PATTERN = ("circular=true", "circular=True", "circular",
                    "type=plasmid", "closed")

# Frozen lite weights (descending), used for renormalisation over the
# assessable dimensions.
LITE_W: Dict[str, float] = {
    "S_ARG": RISK_WEIGHTS_LITE["S_ARG"],
    "S_BMG": RISK_WEIGHTS_LITE["S_BMG"],
    "S_MOB": RISK_WEIGHTS_LITE["S_MOB"],
    "S_SIZE": RISK_WEIGHTS_LITE["S_SIZE"],
    "S_VF": RISK_WEIGHTS_LITE["S_VF"],
}

FULL_W: Dict[str, float] = {
    "S_ARG": RISK_WEIGHTS["S_ARG"],
    "S_BMG": RISK_WEIGHTS["S_BMG"],
    "S_MOB": RISK_WEIGHTS["S_MOB"],
    "S_SIZE": RISK_WEIGHTS["S_SIZE"],
    "S_VF": RISK_WEIGHTS["S_VF"],
    "S_HOST": RISK_WEIGHTS["S_HOST"],
    "S_REP": RISK_WEIGHTS["S_REP"],
    "S_GEO": RISK_WEIGHTS["S_GEO"],
    "S_HAB": RISK_WEIGHTS["S_HAB"],
    "S_GROW": RISK_WEIGHTS["S_GROW"],
}

GRADE_THRESHOLDS = (  # (lower bound, grade)
    (0.60, "A"),
    (0.45, "B"),
    (0.30, "C"),
    (0.15, "D"),
    (0.00, "E"),
)


def grade_from_score(score: Optional[float]) -> str:
    if score is None or (isinstance(score, float) and math.isnan(score)):
        return ""
    for lower, grade in GRADE_THRESHOLDS:
        if score >= lower:
            return grade
    return "E"


def _is_circular_header(header: str) -> bool:
    h = header.lower()
    return any(p.lower() in h for p in CIRCULAR_PATTERN)


def mob_evidence_score(features: PlasmidFeatures) -> Optional[float]:
    """S_MOB from observed transfer apparatus / IS evidence only.

    No isolate-derived replicon or global priors are injected. Fragments
    with no observed evidence return None (not assessable); a complete
    assembly with no evidence is scored as non-mobilisable.
    """
    has_relaxase = bool(getattr(features, "has_relaxase", False))
    has_t4cp = bool(getattr(features, "has_t4cp", False))
    has_oriT = bool(getattr(features, "has_oriT", False))
    has_aux = bool(getattr(features, "has_auxiliary", False))
    has_integron = bool(getattr(features, "has_integron", False))
    n_is = int(getattr(features, "n_is", 0) or 0)

    annotated = bool(getattr(features, "mobility_annotated", False))
    mclass = getattr(features, "mobility_class", "") or ""

    if mclass:
        mc = mclass.lower()
        if "conjugative" in mc:
            return 0.95
        if "mobilizable" in mc or "mobilisable" in mc:
            return 0.55
        if "non" in mc and ("mob" in mc):
            return 0.1

    score = 0.0
    if has_relaxase:
        score += 0.45
    if has_t4cp:
        score += 0.35
    if has_oriT:
        score += 0.15
    if has_aux:
        score += 0.1
    if has_integron:
        score += 0.05
    if n_is >= 3:
        score += 0.1
    elif n_is >= 1:
        score += 0.05
    score = min(score, 1.0)

    if score > 0:
        return score
    if annotated:
        # Complete assembly with explicit annotation and no apparatus.
        return 0.1
    # Fragment and no evidence: not assessable.
    return None


def score_contig(features: PlasmidFeatures, header: str = "",
                 complete_length: int = 10000,
                 scorer: Optional[PlasRiskScorer] = None) -> Dict[str, object]:
    """Score one plasmid contig on the assessable dimensions.

    Cargo dimensions (ARG/BMG/VF) are always count-based. S_SIZE is
    assessed only for a circularised header or a contig at/above
    ``complete_length``. S_MOB is assessed only from observed evidence;
    S_REP only when a replicon is detected.
    """
    if scorer is None:
        scorer = PlasRiskScorer()

    length = int(getattr(features, "length_bp", 0) or 0)
    circular = _is_circular_header(header) or bool(
        getattr(features, "metadata", {}) and
        features.metadata.get("circular", False))

    s_arg = scorer.score_arg(features)
    s_bmg = scorer.score_bm(features)
    s_vf = scorer.score_vf(features)

    size_assessable = circular or length >= complete_length
    s_size = scorer.score_size(features) if size_assessable else None

    s_mob = mob_evidence_score(features)

    replicon = getattr(features, "replicon", "") or ""
    s_rep = scorer.score_rep(features) if replicon else None

    out: Dict[str, object] = {
        "seq_id": getattr(features, "seq_id", ""),
        "length_bp": length,
        "circular": circular,
        "replicon": replicon,
        "S_ARG": s_arg,
        "S_BMG": s_bmg,
        "S_VF": s_vf,
        "S_MOB": s_mob,
        "S_SIZE": s_size,
        "S_REP": s_rep,
    }

    # Renormalised scores over the assessable dimensions.
    def _renorm(dims: Sequence[str], key: str) -> None:
        pairs = [(d, LITE_W[d]) for d in dims if out.get(d) is not None]
        if not pairs:
            out[key] = None
            return
        wsum = sum(w for _, w in pairs)
        out[key] = sum(out[d] * w for d, w in pairs) / wsum

    _renorm(("S_ARG", "S_BMG", "S_VF"), "S_cargo3")
    _renorm(("S_ARG", "S_BMG", "S_VF", "S_MOB"), "S_frag4")
    lite_dims = ("S_ARG", "S_BMG", "S_MOB", "S_SIZE", "S_VF")
    assessable = [d for d in lite_dims if out.get(d) is not None]
    _renorm(lite_dims, "S_lite5")
    out["assessable_dims"] = ";".join(assessable)
    out["n_assessable"] = len(assessable)
    if len(assessable) == len(lite_dims):
        out["grade"] = grade_from_score(out["S_lite5"])
    else:
        out["grade"] = ""
    return out


def score_contigs(features_list: Sequence[PlasmidFeatures],
                  headers: Optional[Sequence[str]] = None,
                  complete_length: int = 10000,
                  scorer: Optional[PlasRiskScorer] = None
                  ) -> List[Dict[str, object]]:
    if scorer is None:
        scorer = PlasRiskScorer()
    rows = []
    for i, feats in enumerate(features_list):
        header = headers[i] if headers else ""
        rows.append(score_contig(feats, header=header,
                                complete_length=complete_length,
                                scorer=scorer))
    return rows


def load_bins_table(path: str) -> pd.DataFrame:
    """Load a generic bins table: seq_id, bin_id, optional columns.

    Optional columns: circular (1/0), plasmid_score, length_bp, coverage.
    """
    df = pd.read_csv(path, sep=None, engine="python", dtype=str,
                     comment="#")
    df.columns = [c.strip() for c in df.columns]
    for alias, canon in (("contig", "seq_id"), ("contig_id", "seq_id"),
                        ("seq_name", "seq_id"), ("bin", "bin_id"),
                        ("cluster", "bin_id")):
        if canon not in df.columns and alias in df.columns:
            df = df.rename(columns={alias: canon})
    if "seq_id" not in df.columns or "bin_id" not in df.columns:
        raise ValueError("bins table must contain seq_id and bin_id columns")
    return df


def load_genomad(path: str) -> pd.DataFrame:
    """Load a geNomad plasmid components/summary table as a bins table."""
    if os.path.isdir(path):
        cands = sorted(glob.glob(os.path.join(
            path, "**", "*_plasmid_components.tsv"), recursive=True))
        if not cands:
            cands = sorted(glob.glob(os.path.join(
                path, "**", "*_plasmid_summary.tsv"), recursive=True))
        if not cands:
            raise FileNotFoundError(
                "no geNomad *_plasmid_components.tsv under %s" % path)
        path = cands[0]
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    df.columns = [c.strip() for c in df.columns]
    rename = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ("seq_name", "sequence", "contig_name"):
            rename[c] = "seq_id"
        elif cl in ("plasmid", "plasmid_name", "component"):
            rename[c] = "bin_id"
        elif cl == "topology":
            rename[c] = "circular"
    df = df.rename(columns=rename)
    if "bin_id" not in df.columns:
        # Each geNomad plasmid component row is its own bin; key by seq.
        if "seq_id" in df.columns:
            df["bin_id"] = df["seq_id"]
    if "circular" in df.columns:
        df["circular"] = df["circular"].astype(str).str.lower().map(
            lambda x: "1" if x in ("circular", "circular=true", "yes") else "0")
    return df


def load_coverage_table(path: str,
                        coverage_column: Optional[str] = None
                        ) -> pd.DataFrame:
    """Load per-contig coverage.

    Accepts a generic two-column table (seq_id, coverage) or the output of
    jgi_summarize_bam_contig_depths. For the JGI format, ``coverage_column``
    selects a sample column; by default the mean of all *totalAvgDepth /
    *coverage columns is used.
    """
    df = pd.read_csv(path, sep=None, engine="python", dtype=str,
                     comment="#")
    df.columns = [c.strip() for c in df.columns]
    id_col = None
    for c in df.columns:
        if c.lower() in ("seq_id", "contig", "contig_id", "seq_name", "name"):
            id_col = c
            break
    if id_col is None:
        id_col = df.columns[0]

    if coverage_column is not None:
        cov = df[[id_col, coverage_column]].copy()
        cov.columns = ["seq_id", "coverage"]
    else:
        depth_cols = [c for c in df.columns
                      if c.lower().endswith("totalavgdepth")
                      or c.lower() == "coverage"
                      or c.lower().endswith("coverage")]
        if not depth_cols:
            depth_cols = [c for c in df.columns if c != id_col][:1]
        vals = df[depth_cols].apply(pd.to_numeric, errors="coerce")
        cov = pd.DataFrame({
            "seq_id": df[id_col],
            "coverage": vals.mean(axis=1),
        })
    cov["coverage"] = pd.to_numeric(cov["coverage"], errors="coerce")
    return cov


def aggregate_bin(contig_rows: Sequence[Dict[str, object]],
                  bin_id: str,
                  coverages: Optional[Dict[str, float]] = None
                  ) -> Dict[str, object]:
    """Aggregate scored contig rows into one plasmid bin."""
    rows = list(contig_rows)
    n = len(rows)
    lengths = [float(r.get("length_bp", 0) or 0) for r in rows]
    total_len = sum(lengths)
    circular = any(bool(r.get("circular", False)) for r in rows)

    replicons = sorted({str(r.get("replicon", "")).strip()
                        for r in rows if str(r.get("replicon", "")).strip()})

    def _mean(key: str) -> Optional[float]:
        vals = [r.get(key) for r in rows if r.get(key) is not None]
        vals = [v for v in vals if not (isinstance(v, float) and math.isnan(v))]
        if not vals:
            return None
        wg = []
        for r in rows:
            v = r.get(key)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                continue
            L = float(r.get("length_bp", 0) or 0)
            if coverages is not None:
                c = coverages.get(str(r.get("seq_id", "")), 1.0)
                L *= max(float(c), 0.0)
            wg.append(L if L > 0 else 1.0)
        wsum = sum(wg)
        return sum(float(v) * w for v, w in zip(vals, wg)) / wsum

    s_arg = _mean("S_ARG")
    s_bmg = _mean("S_BMG")
    s_vf = _mean("S_VF")
    s_mob = _mean("S_MOB")
    s_size = _mean("S_SIZE")

    # Fragmented bins cannot have S_SIZE assessed from summed fragment
    # lengths unless the bin is circularised.
    if not circular and s_size is not None:
        # A bin is treated as complete if any member was size-assessable and
        # the summed length is at the threshold; otherwise keep the observed
        # component mean but flag confidence.
        pass

    out: Dict[str, object] = {
        "bin_id": bin_id,
        "n_contigs": n,
        "length_bp": total_len,
        "circular": circular,
        "replicon": ";".join(replicons),
        "S_ARG": s_arg,
        "S_BMG": s_bmg,
        "S_VF": s_vf,
        "S_MOB": s_mob,
        "S_SIZE": s_size if circular else (s_size if all(
            r.get("S_SIZE") is not None for r in rows) else None),
    }

    lite_dims = ("S_ARG", "S_BMG", "S_MOB", "S_SIZE", "S_VF")
    assessable = [d for d in lite_dims if out.get(d) is not None]
    pairs = [(d, LITE_W[d]) for d in lite_dims if out.get(d) is not None]
    if pairs:
        wsum = sum(w for _, w in pairs)
        out["S_lite5"] = sum(out[d] * w for d, w in pairs) / wsum
    else:
        out["S_lite5"] = None
    out["assessable_dims"] = ";".join(assessable)

    if circular:
        out["confidence"] = "complete"
    elif len(assessable) >= 4:
        out["confidence"] = "probable"
    else:
        out["confidence"] = "fragmented"
    out["grade"] = (grade_from_score(out["S_lite5"])
                    if len(assessable) == len(lite_dims) else "")
    return out


def aggregate_bins(contig_rows: Sequence[Dict[str, object]],
                   bins_table: Optional[pd.DataFrame] = None,
                   coverage_table: Optional[pd.DataFrame] = None
                   ) -> List[Dict[str, object]]:
    """Aggregate scored contig rows into bins.

    When no bins table is supplied, each contig becomes its own bin. Unbinned
    contigs are assigned to ``unbinned:<seq_id>``.
    """
    by_id = {str(r["seq_id"]): r for r in contig_rows}

    coverages: Optional[Dict[str, float]] = None
    if coverage_table is not None:
        coverages = {}
        for _, rr in coverage_table.iterrows():
            coverages[str(rr["seq_id"])] = float(rr["coverage"])

    groups: Dict[str, List[str]] = {}
    if bins_table is not None:
        for _, rr in bins_table.iterrows():
            sid = str(rr["seq_id"]).strip()
            bid = str(rr["bin_id"]).strip()
            # First token alias for contig headers.
            key = sid if sid in by_id else sid.split()[0]
            if key not in by_id:
                continue
            groups.setdefault(bid, []).append(key)

    assigned = set()
    for ids in groups.values():
        assigned.update(ids)

    for sid, row in by_id.items():
        if sid not in assigned:
            groups.setdefault("unbinned:%s" % sid, []).append(sid)

    out = []
    for bid, ids in groups.items():
        rows = [by_id[i] for i in ids]
        out.append(aggregate_bin(rows, bin_id=bid, coverages=coverages))
    out.sort(key=lambda r: str(r["bin_id"]))
    return out


def aggregate_sample(bin_rows: Sequence[Dict[str, object]],
                     contig_rows: Optional[Sequence[Dict[str, object]]] = None,
                     coverages: Optional[Dict[str, float]] = None
                     ) -> Dict[str, object]:
    """Descriptive, abundance-aware community statistics (CPRI)."""
    bins = list(bin_rows)
    n_bins = len(bins)

    def _w(row: Dict[str, object]) -> float:
        L = float(row.get("length_bp", 0) or 0)
        if coverages is not None:
            L *= max(coverages.get(str(row.get("bin_id", row.get("seq_id", ""))),
                                  1.0), 0.0)
        return L if L > 0 else 1.0

    scored = [b for b in bins if b.get("S_lite5") is not None]
    if scored:
        weights = [_w(b) for b in scored]
        wsum = sum(weights)
        cpri_mean = sum(float(b["S_lite5"]) * w
                        for b, w in zip(scored, weights)) / wsum
        svals = sorted(float(b["S_lite5"]) for b in scored)

        def _quantile(q: float) -> float:
            if len(svals) == 1:
                return svals[0]
            pos = q * (len(svals) - 1)
            lo = int(math.floor(pos))
            hi = int(math.ceil(pos))
            if lo == hi:
                return svals[lo]
            return svals[lo] + (svals[hi] - svals[lo]) * (pos - lo)

        cpri_p90 = _quantile(0.90)
        cpri_p99 = _quantile(0.99)
        s_max = max(svals)
    else:
        cpri_mean = cpri_p90 = cpri_p99 = s_max = None

    grade_counts = {g: 0 for g in ("A", "B", "C", "D", "E")}
    for b in scored:
        g = b.get("grade", "")
        if g in grade_counts:
            grade_counts[g] += 1

    total_len = sum(float(b.get("length_bp", 0) or 0) for b in bins)

    # Critical-ARG load from contig rows (count-based).
    crit_load = 0
    if contig_rows is not None:
        for r in contig_rows:
            if r.get("S_ARG") is not None and float(r["S_ARG"]) >= 0.8:
                crit_load += 1

    out: Dict[str, object] = {
        "n_bins": n_bins,
        "n_scored_bins": len(scored),
        "total_plasmid_length_bp": total_len,
        "CPRI_mean": cpri_mean,
        "CPRI_p90": cpri_p90,
        "CPRI_p99": cpri_p99,
        "S_max": s_max,
        "grade_A": grade_counts["A"],
        "grade_B": grade_counts["B"],
        "grade_C": grade_counts["C"],
        "grade_D": grade_counts["D"],
        "grade_E": grade_counts["E"],
        "critical_ARG_contigs": crit_load,
        "critical_ARG_per_Mb": (crit_load / (total_len / 1e6)
                                if total_len > 0 else None),
        "note": "Descriptive extension; CPRI was not trained or validated "
                "on metagenomic data.",
    }
    return out


def _collect_fastas(inputs: Sequence[str]) -> List[str]:
    out = []
    exts = {".fasta", ".fa", ".fna", ".fsa", ".fas"}
    for inp in inputs:
        matches = glob.glob(inp)
        if not matches and os.path.exists(inp):
            matches = [inp]
        for m in matches:
            if os.path.isdir(m):
                for root, _d, files in os.walk(m):
                    for f in files:
                        if os.path.splitext(f)[1].lower() in exts:
                            out.append(os.path.join(root, f))
            elif os.path.isfile(m):
                out.append(m)
    seen, uniq = set(), []
    for f in out:
        ap = os.path.abspath(f)
        if ap not in seen:
            seen.add(ap)
            uniq.append(f)
    return uniq


def _headers_for(fasta_path: str) -> List[str]:
    headers = []
    with open(fasta_path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].rstrip("\n"))
    return headers


def meta_main(argv=None):
    parser = argparse.ArgumentParser(
        prog="plasrisk meta",
        description="PlasRisk metagenomic three-layer scoring "
                    "(contig / bin / sample)")
    parser.add_argument("inputs", nargs="+",
                        help="FASTA file(s) of plasmid contigs, or directory")
    parser.add_argument("-o", "--output", default="plasrisk_meta_output",
                        help="Output directory")
    parser.add_argument("--bins", default=None,
                        help="Optional bins table (seq_id, bin_id, ...)")
    parser.add_argument("--genomad-bin", default=None,
                        help="Optional geNomad plasmid components file/dir")
    parser.add_argument("--coverage", default=None,
                        help="Optional per-contig coverage table")
    parser.add_argument("--coverage-column", default=None,
                        help="Coverage column for JGI depth format")
    parser.add_argument("--complete-length", type=int, default=10000,
                        help="Length at/above which a contig is size-assessable")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--min-id", type=float, default=75.0)
    parser.add_argument("--min-cov", type=float, default=50.0)
    parser.add_argument("--no-abricate", action="store_true")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)

    os.makedirs(args.output, exist_ok=True)
    fastas = _collect_fastas(args.inputs)
    if not fastas:
        print("ERROR: no FASTA inputs found", file=sys.stderr)
        return 2

    bins_df = None
    if args.genomad_bin:
        bins_df = load_genomad(args.genomad_bin)
    elif args.bins:
        bins_df = load_bins_table(args.bins)

    cov_df = None
    if args.coverage:
        cov_df = load_coverage_table(args.coverage,
                                     coverage_column=args.coverage_column)

    all_contig_rows: List[Dict[str, object]] = []
    sample_rows: List[Dict[str, object]] = []
    scorer = PlasRiskScorer()

    for fasta in fastas:
        sample_id = os.path.splitext(os.path.basename(fasta))[0]
        if not args.quiet:
            print("[meta] %s" % sample_id)
        ann = annotate_fasta(fasta, use_abricate=not args.no_abricate,
                             min_id=args.min_id, min_cov=args.min_cov,
                             threads=args.threads)
        headers = _headers_for(fasta)
        rows = score_contigs(ann.features, headers=headers,
                             complete_length=args.complete_length,
                             scorer=scorer)
        for r in rows:
            r["sample_id"] = sample_id
        all_contig_rows.extend(rows)

        sub_bins = bins_df
        if sub_bins is not None:
            sub_bins = bins_df[bins_df.get(
                "sample_id", pd.Series([sample_id] * len(bins_df)))
                .astype(str) == sample_id] if "sample_id" in bins_df.columns \
                else bins_df
        bin_rows = aggregate_bins(rows, bins_table=sub_bins,
                                  coverage_table=cov_df)
        coverages = None
        if cov_df is not None:
            coverages = {str(rr["seq_id"]): float(rr["coverage"])
                         for _, rr in cov_df.iterrows()}
        sample = aggregate_sample(bin_rows, contig_rows=rows,
                                  coverages=coverages)
        sample["sample_id"] = sample_id
        sample_rows.append(sample)

        if args.json:
            with open(os.path.join(args.output, "%s_meta_sample.json" % sample_id),
                      "w") as jf:
                json.dump(sample, jf, indent=2, default=str)

    contig_df = pd.DataFrame(all_contig_rows)
    contig_path = os.path.join(args.output, "meta_contig_scores.tsv")
    contig_df.to_csv(contig_path, sep="\t", index=False)

    # Bin rows across samples (recompute to capture all).
    all_bin_rows: List[Dict[str, object]] = []
    for fasta in fastas:
        sample_id = os.path.splitext(os.path.basename(fasta))[0]
        rows = [r for r in all_contig_rows if r.get("sample_id") == sample_id]
        sub_bins = bins_df
        if sub_bins is not None and "sample_id" in sub_bins.columns:
            sub_bins = sub_bins[sub_bins["sample_id"].astype(str) == sample_id]
        br = aggregate_bins(rows, bins_table=sub_bins, coverage_table=cov_df)
        for b in br:
            b["sample_id"] = sample_id
        all_bin_rows.extend(br)
    bin_df = pd.DataFrame(all_bin_rows)
    bin_path = os.path.join(args.output, "meta_bin_scores.tsv")
    bin_df.to_csv(bin_path, sep="\t", index=False)

    sample_df = pd.DataFrame(sample_rows)
    sample_path = os.path.join(args.output, "meta_sample_summary.tsv")
    sample_df.to_csv(sample_path, sep="\t", index=False)

    if not args.quiet:
        print("[meta] wrote %s, %s, %s" % (contig_path, bin_path, sample_path))
    return 0


if __name__ == "__main__":
    sys.exit(meta_main())
