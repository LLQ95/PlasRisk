"""
metagenome.py - Metagenomic extension of PlasRisk (three aggregation layers).

Layer 1  contig / element level
    Fragmented plasmid contigs from a metagenomic assembly are scored on the
    dimensions that remain identifiable after fragmentation. S_SIZE is not
    evaluated from a fragment length (it would systematically penalise
    megaplasmids that the assembler split), and isolate-derived mobility
    priors are never injected when the transfer apparatus was simply not
    covered. Scores are renormalised over the assessable dimensions:

        S_cargo3  = w_ARG, w_BMG, w_VF           (always assessable)
        S_frag4   = cargo3 + S_MOB               (mobility evidence present)
        S_lite5   = frag4  + S_SIZE              (circular / near-complete)

Layer 2  plasmid bin level
    Contigs assigned to the same predicted plasmid (geNomad components,
    SCAPP plasmids, MOB-suite bins, ...) are aggregated: gene sets are
    unioned (deduplicated), lengths are summed and mobility markers are
    OR-combined. Bins are then scored with the standard lite model, so bin
    scores are on the same scale as isolate plasmid scores. A confidence
    label (complete / probable / fragmented) is always carried.

Layer 3  sample / community level
    Descriptive, abundance-aware summary statistics per sample:
    coverage-length weighted mean score (CPRI_mean), P90/P99, grade
    prevalence, and coverage-weighted risk loads (grade A/B load,
    critical-ARG load). These indices were not trained or validated in the
    isolate-derived study and must be reported as an extension, not as
    validated predictions.

The frozen consensus weights are never refit.
"""

from __future__ import annotations

import math
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .annotate import AnnotationResult, annotate_fasta, parse_fasta
from .scoring import (
    GLOBAL_MEDIANS,
    LITE_DIMENSIONS,
    PlasRiskScorer,
    PlasmidFeatures,
    RISK_GRADES,
    RISK_WEIGHTS_LITE,
)

# Lite dimensions that can be observed directly from sequence evidence
CARGO_DIMS = ("S_ARG", "S_BMG", "S_VF")
CIRCULAR_PATTERN = re.compile(r"circular[=: ]*(true|yes|y|1)\b", re.I)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _grade(s_norm: float) -> Tuple[str, str]:
    for threshold, grade, label in RISK_GRADES:
        if s_norm >= threshold:
            return grade, label
    return "E", "Minimal"


def _renorm(component_values: Dict[str, Optional[float]],
            weights: Dict[str, float] = RISK_WEIGHTS_LITE
            ) -> Tuple[Optional[float], List[str]]:
    """Weighted mean over dimensions that are not None, renormalised so the
    weights of missing dimensions are redistributed (frozen weights)."""
    avail = {k: float(v) for k, v in component_values.items()
             if k in weights and v is not None and not (isinstance(v, float)
                                                        and math.isnan(v))}
    if not avail:
        return None, []
    wsum = sum(weights[k] for k in avail)
    score = sum(weights[k] * v for k, v in avail.items()) / wsum
    return float(score), list(avail.keys())


def is_circular_header(seq_id: str) -> bool:
    """Detect explicit circularisation flags in an assembly header."""
    return bool(CIRCULAR_PATTERN.search(str(seq_id)))


def circular_map_from_fasta(path: str) -> Dict[str, bool]:
    """Scan FASTA headers and map the first-word record id to a circular
    flag (parse_fasta/abricate keep only the first word as id)."""
    out: Dict[str, bool] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                rid = header.split()[0] if header else "seq"
                out[rid] = is_circular_header(header)
    return out


def mob_evidence_score(feat: PlasmidFeatures,
                       complete: bool) -> Optional[float]:
    """S_MOB from observed sequence evidence only.

    Mirrors PlasRiskScorer.score_mob but never falls back to isolate-derived
    replicon / global priors. Returns None when mobility cannot be assessed
    (fragment with no transfer apparatus or IS evidence); for a complete
    assembly, assembly-wide absence is informative and receives the
    non-mobilizable base of 0.1.
    """
    has_t4cp = bool(feat.has_t4cp)
    has_rel = bool(feat.has_relaxase)
    has_oriT = bool(feat.has_oriT)

    if feat.mobility_class:
        mc = feat.mobility_class.lower()
        if "complete" in mc:
            mob_base = 1.0
        elif "likely" in mc:
            mob_base = 0.85
        elif "mobilizable" in mc:
            mob_base = 0.6
        else:
            mob_base = 0.1
    elif has_t4cp and has_rel and has_oriT:
        mob_base = 1.0
    elif has_t4cp and has_rel:
        mob_base = 0.85
    elif has_rel or has_oriT:
        mob_base = 0.6
    elif feat.n_is > 0:
        mob_base = 0.1
    else:
        if complete:
            mob_base = 0.1
        else:
            return None

    mob_component = min(mob_base + (0.10 if feat.has_integron else 0.0), 1.0)
    length_kb = feat.length_bp / 1000.0 if feat.length_bp > 0 else 1.0
    is_density = feat.n_is / length_kb if length_kb > 0 else 0.0
    is_bonus = 0.30 * min(is_density, 2.0) / 2.0
    return float(np.clip(mob_component + is_bonus, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Layer 1: contig-level scoring
# ---------------------------------------------------------------------------

@dataclass
class ContigScore:
    seq_id: str
    length_bp: int
    circular: bool
    complete: bool
    replicon: str
    n_ARG: int
    n_VF: int
    n_BM: int
    high_risk_genes: str
    S_ARG: float
    S_BMG: float
    S_VF: float
    S_MOB: Optional[float]
    S_SIZE: Optional[float]
    S_REP: Optional[float]
    S_cargo3: Optional[float]
    S_frag4: Optional[float]
    S_lite5: Optional[float]
    grade: str
    grade_label: str
    assessable_dims: str


def score_contig(feat: PlasmidFeatures, scorer: PlasRiskScorer,
                 complete_length: int = 10_000,
                 circular_flag: bool = False) -> ContigScore:
    """Score one metagenomic contig with dimension completeness tracking."""
    circular = bool(circular_flag) or is_circular_header(feat.seq_id)
    complete = circular or feat.length_bp >= complete_length

    s_arg = scorer.score_arg(feat)
    s_bmg = scorer.score_bmg(feat)
    s_vf = scorer.score_vf(feat)
    s_mob = mob_evidence_score(feat, complete=complete)
    s_size = scorer.score_size(feat.length_bp) if complete else None
    s_rep = scorer.score_rep(feat) if feat.replicon else None

    cargo_vals = {"S_ARG": s_arg, "S_BMG": s_bmg, "S_VF": s_vf}
    s_cargo3, _ = _renorm(cargo_vals)

    s_frag4 = None
    s_lite5 = None
    lite5_dims = []
    if s_mob is not None:
        frag4_vals = dict(cargo_vals)
        frag4_vals["S_MOB"] = s_mob
        s_frag4, _ = _renorm(frag4_vals)
        if s_size is not None:
            lite5_vals = dict(frag4_vals)
            lite5_vals["S_SIZE"] = s_size
            s_lite5, lite5_dims = _renorm(lite5_vals)

    if s_lite5 is not None and set(lite5_dims) == set(LITE_DIMENSIONS):
        grade, grade_label = _grade(s_lite5)
    else:
        grade, grade_label = "", ""

    assessable = list(CARGO_DIMS)
    if s_mob is not None:
        assessable.append("S_MOB")
    if s_size is not None:
        assessable.append("S_SIZE")

    return ContigScore(
        seq_id=feat.seq_id,
        length_bp=feat.length_bp,
        circular=circular,
        complete=complete,
        replicon=feat.replicon or "",
        n_ARG=feat.n_arg,
        n_VF=feat.n_vf,
        n_BM=feat.n_bm,
        high_risk_genes=";".join(feat.high_risk_args) if feat.high_risk_args
        else "",
        S_ARG=round(s_arg, 4),
        S_BMG=round(s_bmg, 4),
        S_VF=round(s_vf, 4),
        S_MOB=round(s_mob, 4) if s_mob is not None else None,
        S_SIZE=round(s_size, 4) if s_size is not None else None,
        S_REP=round(s_rep, 4) if s_rep is not None else None,
        S_cargo3=round(s_cargo3, 4) if s_cargo3 is not None else None,
        S_frag4=round(s_frag4, 4) if s_frag4 is not None else None,
        S_lite5=round(s_lite5, 4) if s_lite5 is not None else None,
        grade=grade,
        grade_label=grade_label,
        assessable_dims=";".join(assessable),
    )


def score_contigs(annotation: AnnotationResult, scorer: PlasRiskScorer,
                  complete_length: int = 10_000,
                  min_length: int = 1_000,
                  circular_map: Optional[Dict[str, bool]] = None
                  ) -> pd.DataFrame:
    rows = []
    for feat in annotation.features:
        if feat.length_bp < min_length:
            continue
        cflag = bool((circular_map or {}).get(feat.seq_id, False))
        rows.append(score_contig(feat, scorer,
                                 complete_length=complete_length,
                                 circular_flag=cflag).__dict__)
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("S_cargo3", ascending=False).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Layer 2: plasmid bin aggregation
# ---------------------------------------------------------------------------

def _union(series_of_lists: List[List[str]]) -> List[str]:
    seen = set()
    out = []
    for items in series_of_lists:
        for x in items:
            if x not in seen:
                seen.add(x)
                out.append(x)
    return out


def aggregate_bin(bin_id: str, feats: List[PlasmidFeatures],
                  scorer: PlasRiskScorer, circular: Optional[bool],
                  plasmid_score: Optional[float],
                  coverages: Optional[Dict[str, float]],
                  complete_length: int = 10_000,
                  circ_contigs: Optional[set] = None) -> Dict:
    """Aggregate all contigs of one predicted plasmid into a bin and score it.

    Gene sets are unioned with deduplication, lengths are summed, mobility
    markers are OR-combined, and per-contig coverage is length-weighted.
    A fragmented bin (neither circularised nor at/above complete_length in
    total length) is not scored on S_SIZE; mobility is then assessed only
    from observed evidence.
    """
    circ_contigs = circ_contigs or set()
    any_circ = bool(circular) or any(f.seq_id in circ_contigs for f in feats)
    total_len = int(sum(f.length_bp for f in feats))
    size_assessable = any_circ or total_len >= complete_length

    arg_names = _union([f.arg_names for f in feats])
    vf_names = _union([f.vf_names for f in feats])
    bm_names = _union([f.bm_gene_names for f in feats])
    n_is = int(sum(f.n_is for f in feats))
    has_t4cp = any(f.has_t4cp for f in feats)
    has_rel = any(f.has_relaxase for f in feats)
    has_oriT = any(f.has_oriT for f in feats)
    has_aux = any(f.auxiliary for f in feats)
    has_integron = any(f.has_integron for f in feats)
    replicons = _union([[f.replicon] if f.replicon else [] for f in feats])
    replicon = replicons[0] if replicons else ""

    mob_classes = [f.mobility_class for f in feats if f.mobility_class]
    mobility_class = mob_classes[0] if mob_classes else (
        "non-mobilizable" if size_assessable else "unassessed")

    high_risk = _union([f.high_risk_args for f in feats])

    if coverages:
        num = sum(coverages.get(f.seq_id, 0.0) * f.length_bp for f in feats)
        den = sum(f.length_bp for f in feats)
        coverage = num / den if den else np.nan
    else:
        coverage = np.nan

    agg = PlasmidFeatures(
        seq_id=bin_id,
        length_bp=total_len,
        arg_names=arg_names,
        vf_names=vf_names,
        bm_gene_names=bm_names,
        replicon=replicon,
        has_t4cp=has_t4cp,
        has_relaxase=has_rel,
        has_oriT=has_oriT,
        auxiliary=has_aux,
        integron=has_integron,
        n_is=n_is,
        mobility_class=mobility_class,
        mobility_annotated=bool(mob_classes),
    )

    s_arg = scorer.score_arg(agg)
    s_bmg = scorer.score_bmg(agg)
    s_vf = scorer.score_vf(agg)
    s_mob = mob_evidence_score(agg, complete=size_assessable)
    s_size = scorer.score_size(total_len) if size_assessable else None
    s_rep = scorer.score_rep(agg) if replicon else None

    cargo_vals = {"S_ARG": s_arg, "S_BMG": s_bmg, "S_VF": s_vf}
    s_cargo3, _ = _renorm(cargo_vals)
    s_frag4, s_lite, grade, grade_label = None, None, "", ""
    dims = list(CARGO_DIMS)
    if s_mob is not None:
        f4 = dict(cargo_vals)
        f4["S_MOB"] = s_mob
        s_frag4, f4dims = _renorm(f4)
        dims.append("S_MOB")
        if s_size is not None:
            l5 = dict(f4)
            l5["S_SIZE"] = s_size
            s_lite, l5dims = _renorm(l5)
            dims.append("S_SIZE")
            if set(l5dims) == set(LITE_DIMENSIONS):
                grade, grade_label = _grade(s_lite)

    if any_circ:
        confidence = "complete"
    elif size_assessable:
        confidence = "probable"
    else:
        confidence = "fragmented"

    return {
        "rank": 0,
        "bin_id": bin_id,
        "n_contigs": len(feats),
        "total_length_bp": total_len,
        "total_length_kb": round(total_len / 1000.0, 2),
        "circular": any_circ,
        "confidence": confidence,
        "plasmid_score": plasmid_score,
        "coverage": round(coverage, 3) if not math.isnan(coverage) else np.nan,
        "replicon": replicon or "Unknown",
        "n_ARG": len(arg_names),
        "n_VF": len(vf_names),
        "n_BM": len(bm_names),
        "n_IS": n_is,
        "high_risk_genes": ";".join(high_risk),
        "mobility_class": mobility_class,
        "S_ARG": round(s_arg, 4),
        "S_BMG": round(s_bmg, 4),
        "S_VF": round(s_vf, 4),
        "S_SIZE": round(s_size, 4) if s_size is not None else None,
        "S_MOB": round(s_mob, 4) if s_mob is not None else None,
        "S_REP": round(s_rep, 4) if s_rep is not None else None,
        "S_cargo3": round(s_cargo3, 4) if s_cargo3 is not None else None,
        "S_frag4": round(s_frag4, 4) if s_frag4 is not None else None,
        "S_lite": round(s_lite, 4) if s_lite is not None else None,
        "grade": grade,
        "grade_label": grade_label,
        "assessable_dims": ";".join(dims),
    }


# ---------------------------------------------------------------------------
# Bin / coverage table parsers
# ---------------------------------------------------------------------------

def _find_column(columns: List[str], candidates: Tuple[str, ...]) -> Optional[str]:
    lower = {str(c).lower().strip(): c for c in columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    for c in columns:
        cl = str(c).lower().strip()
        if any(cand in cl for cand in candidates):
            return c
    return None


def load_bins_table(path: str) -> pd.DataFrame:
    """Generic contig-to-bin map.

    Required columns: seq_id, bin_id (alternative names accepted:
    contig/sequence/seq_name/contigName and bin/component/plasmid/bin_id).
    Optional columns: circular / topology, plasmid_score.
    """
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = [str(c).strip() for c in df.columns]
    seq_col = _find_column(list(df.columns),
                           ("seq_id", "seqname", "sequence", "contig",
                            "seq_name", "contigname", "name", "id"))
    bin_col = _find_column(list(df.columns),
                           ("bin_id", "bin", "component", "plasmid",
                            "plasmid_id", "cluster"))
    if seq_col is None or bin_col is None:
        raise ValueError(
            "bins table must contain seq_id and bin_id columns "
            "(found: %s)" % ", ".join(df.columns))
    out = pd.DataFrame({
        "seq_id": df[seq_col].astype(str),
        "bin_id": df[bin_col].astype(str),
    })
    circ_col = _find_column(list(df.columns),
                            ("circular", "topology", "circ"))
    if circ_col is not None:
        out["circular"] = df[circ_col].astype(str).str.lower().map(
            lambda x: x in ("1", "true", "yes", "y", "circular", "c"))
    else:
        out["circular"] = False
    score_col = _find_column(list(df.columns),
                             ("plasmid_score", "plasmid score", "score"))
    if score_col is not None:
        out["plasmid_score"] = pd.to_numeric(df[score_col],
                                             errors="coerce")
    else:
        out["plasmid_score"] = np.nan
    return out.drop_duplicates(["seq_id", "bin_id"])


def load_genomad(genomad_dir: str) -> pd.DataFrame:
    """Parse geNomad aggregate outputs into the generic bins table.

    Expects, anywhere under ``genomad_dir``:
      * ``*_plasmid_components.tsv`` with columns seq_name, components
        (components are ';'-separated and contigs within one component are
        comma-separated);
      * optionally ``*_plasmid_summary.tsv`` with per-plasmid length,
        topology and plasmid_score.
    """
    comp_files, summ_files = [], []
    for root, _dirs, files in os.walk(genomad_dir):
        for f in files:
            if f.endswith("_plasmid_components.tsv"):
                comp_files.append(os.path.join(root, f))
            elif f.endswith("_plasmid_summary.tsv"):
                summ_files.append(os.path.join(root, f))
    if not comp_files:
        raise FileNotFoundError(
            "no *_plasmid_components.tsv found under %s" % genomad_dir)

    rows = []
    for cf in comp_files:
        cdf = pd.read_csv(cf, sep="\t", dtype=str)
        seq_col = _find_column(list(cdf.columns), ("seq_name", "seq_id"))
        comp_col = _find_column(list(cdf.columns), ("components", "component"))
        if seq_col is None or comp_col is None:
            continue
        for _, r in cdf.iterrows():
            cell = str(r[comp_col])
            if cell in ("", "nan"):
                continue
            for i, group in enumerate(cell.split(";")):
                contigs = [x.strip() for x in group.split(",") if x.strip()]
                # geNomad names components "<seq_name>_component_<i>" when
                # the assembly holds several plasmids; stay compatible with
                # summary seq_name values where present.
                bin_id = "%s_component_%d" % (str(r[seq_col]), i + 1)
                for contig in contigs:
                    rows.append({"seq_id": contig, "bin_id": bin_id,
                                 "circular": False, "plasmid_score": np.nan})
    bins = pd.DataFrame(rows)

    if summ_files:
        sdf = pd.read_csv(summ_files[0], sep="\t", dtype=str)
        name_col = _find_column(list(sdf.columns), ("seq_name", "bin_id",
                                                    "component"))
        topo_col = _find_column(list(sdf.columns), ("topology", "circular"))
        score_col = _find_column(list(sdf.columns), ("plasmid_score", "score"))
        if name_col is not None:
            lookup = {}
            for _, r in sdf.iterrows():
                name = str(r[name_col])
                circ = False
                if topo_col is not None:
                    circ = str(r[topo_col]).lower() in ("circular", "1",
                                                        "true", "yes")
                score = (float(r[score_col]) if score_col is not None
                         and str(r[score_col]) not in ("", "nan")
                         else np.nan)
                lookup[name] = (circ, score)
            bins["circular"] = [lookup.get(b, (False, np.nan))[0]
                                for b in bins["bin_id"]]
            bins["plasmid_score"] = [lookup.get(b, (False, np.nan))[1]
                                     for b in bins["bin_id"]]
    return bins


def load_coverage_table(path: str,
                        coverage_column: Optional[str] = None) -> pd.DataFrame:
    """Per-contig coverage. Accepts generic TSV/CSV and jgi_summarize output
    (contigName + per-sample depth columns)."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = [str(c).strip() for c in df.columns]
    seq_col = _find_column(list(df.columns),
                           ("seq_id", "seqname", "sequence", "contig",
                            "seq_name", "contigname", "name"))
    if seq_col is None:
        raise ValueError("coverage table must contain a contig/seq_id column")
    if coverage_column is not None:
        if coverage_column not in df.columns:
            raise ValueError("coverage column '%s' not found (columns: %s)" % (
                coverage_column, ", ".join(df.columns)))
        return pd.DataFrame({
            "seq_id": df[seq_col].astype(str),
            "coverage": pd.to_numeric(df[coverage_column],
                                       errors="coerce"),
        })
    depth_candidates = ("coverage", "depth", "abundance", "meandepth",
                        "totalavgdepth", "avg_depth")
    depth_col = _find_column(list(df.columns), depth_candidates)
    if depth_col is None:
        numeric = [c for c in df.columns
                   if pd.to_numeric(df[c], errors="coerce").notna().any()
                   and c != seq_col]
        if not numeric:
            raise ValueError("no numeric coverage column found")
        depth_col = numeric[0]
    return pd.DataFrame({
        "seq_id": df[seq_col].astype(str),
        "coverage": pd.to_numeric(df[depth_col], errors="coerce"),
    })


def aggregate_bins(contig_df: pd.DataFrame,
                   features: Dict[str, PlasmidFeatures],
                   bins_table: Optional[pd.DataFrame],
                   coverage_df: Optional[pd.DataFrame],
                   scorer: PlasRiskScorer,
                   min_plasmid_score: float = 0.0,
                   complete_length: int = 10_000,
                   header_circular_map: Optional[Dict[str, bool]] = None
                   ) -> pd.DataFrame:
    """Build Layer 2 bin rows from Layer 1 contigs.

    Contigs without a bin assignment become singleton "unbinned:<id>" bins.
    Bins below min_plasmid_score are dropped.
    """
    header_circular_map = header_circular_map or {}
    assignments: Dict[str, List[str]] = {}
    bin_meta: Dict[str, Dict] = {}
    if bins_table is not None:
        for _, r in bins_table.iterrows():
            sid, bid = str(r["seq_id"]), str(r["bin_id"])
            if sid not in features:
                continue
            assignments.setdefault(bid, []).append(sid)
            bin_meta[bid] = {
                "circular": bool(r.get("circular", False)),
                "plasmid_score": r.get("plasmid_score", np.nan),
            }

    for sid in features:
        if not any(sid in members for members in assignments.values()):
            assignments["unbinned:%s" % sid] = [sid]

    coverages = None
    if coverage_df is not None and not coverage_df.empty:
        coverages = dict(zip(coverage_df["seq_id"].astype(str),
                             coverage_df["coverage"]))

    rows = []
    for bid, members in assignments.items():
        feats = [features[m] for m in members if m in features]
        if not feats:
            continue
        meta = bin_meta.get(bid, {"circular": False,
                                  "plasmid_score": np.nan})
        score = meta.get("plasmid_score", np.nan)
        if (min_plasmid_score and not pd.isna(score)
                and float(score) < min_plasmid_score):
            continue
        circ_contigs = {s for s in members
                        if bool(header_circular_map.get(s, False))}
        rows.append(aggregate_bin(
            bid, feats, scorer,
            circular=meta.get("circular", False),
            plasmid_score=score,
            coverages=coverages,
            complete_length=complete_length,
            circ_contigs=circ_contigs))

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(
            ["S_lite", "S_frag4", "S_cargo3"],
            ascending=False, na_position="last").reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)
    return df


# ---------------------------------------------------------------------------
# Layer 3: sample / community level
# ---------------------------------------------------------------------------

def aggregate_sample(sample_id: str, bin_df: pd.DataFrame,
                     contig_df: pd.DataFrame,
                     total_assembly_bp: Optional[int] = None) -> Dict:
    """Descriptive community statistics for one sample.

    Uses the best available per-bin score (S_lite, then S_frag4, then
    S_cargo3) so fragmented bins with detectable cargo still contribute to
    the community index; grades are reported only for fully assessable bins.
    """
    b = bin_df.copy()
    s_use = b["S_lite"] if "S_lite" in b.columns else pd.Series(
        np.nan, index=b.index)
    for col in ("S_frag4", "S_cargo3"):
        if col in b.columns:
            s_use = s_use.combine_first(b[col])
    b["S_use"] = s_use
    scored = b.dropna(subset=["S_use"]).copy()
    graded = scored[scored["grade"].fillna("") != ""].copy()
    n_scored = len(scored)
    n_graded = len(graded)

    if "coverage" in b.columns and b["coverage"].notna().any():
        w = b["coverage"].fillna(0.0).clip(lower=0).to_numpy(dtype=float)
        if w.sum() <= 0:
            w = b["total_length_bp"].to_numpy(dtype=float)
        weighting = "coverage_x_length"
    else:
        w = b["total_length_bp"].to_numpy(dtype=float)
        weighting = "length_only"
    if w.sum() <= 0:
        w = np.ones(len(b))
        weighting = "equal_weight"

    def wquantile(values, weights, q):
        values = np.asarray(values)
        order = np.argsort(values)
        values, weights = values[order], weights[order]
        cum = np.cumsum(weights) - 0.5 * weights
        return float(np.interp(q * weights.sum(), cum, values))

    vals = scored["S_use"].to_numpy(dtype=float)
    ww = w[scored.index.to_numpy()]
    if n_scored:
        cpri_mean = float(np.sum(vals * ww) / np.sum(ww))
        cpri_p90 = wquantile(vals, ww, 0.90)
        cpri_p99 = wquantile(vals, ww, 0.99)
        s_max = float(vals.max())
    else:
        cpri_mean = cpri_p90 = cpri_p99 = s_max = np.nan

    grade_counts = {g: int((graded["grade"] == g).sum())
                    for g in ("A", "B", "C", "D", "E")}
    ab = graded[graded["grade"].isin(("A", "B"))]
    # critical ARGs remain detectable on fragments; count them across all
    # scored bins, not only bins with a complete five-dimension grade
    hr = scored[scored["high_risk_genes"].fillna("") != ""]

    plasmid_bp = float(np.nansum(b["total_length_bp"].to_numpy(
        dtype=float)))
    denom_bp = float(total_assembly_bp) if total_assembly_bp else plasmid_bp

    def _load(df, score_col):
        if df.empty or "coverage" not in df.columns:
            return float(df["total_length_bp"].sum()) if not df.empty else 0.0
        c = df["coverage"].fillna(0.0).clip(lower=0).to_numpy(dtype=float)
        if c.sum() <= 0:
            return float(df["total_length_bp"].sum())
        read_bp = c * df["total_length_bp"].to_numpy(dtype=float)
        scale = denom_bp / read_bp.sum() if read_bp.sum() else 1.0
        return float(read_bp.sum() * scale)

    ab_load = _load(ab, "S_lite")
    hr_load = _load(hr, "S_use")

    n_complete = int((b["confidence"] == "complete").sum())
    n_probable = int((b["confidence"] == "probable").sum())
    n_fragmented = int((b["confidence"] == "fragmented").sum())

    row = {
        "sample_id": sample_id,
        "n_contig_elements": int(len(contig_df)),
        "n_plasmid_bins": int(len(b)),
        "n_scored_bins": int(n_scored),
        "n_graded_bins": int(n_graded),
        "n_complete": n_complete,
        "n_probable": n_probable,
        "n_fragmented": n_fragmented,
        "pct_ARG_positive": round(
            100.0 * (b["n_ARG"].fillna(0) > 0).mean(), 2),
        "pct_BMG_positive": round(
            100.0 * (b["n_BM"].fillna(0) > 0).mean(), 2),
        "pct_high_risk_ARG": round(
            100.0 * (b["high_risk_genes"].fillna("") != "").mean(), 2),
        "CPRI_mean": round(cpri_mean, 4),
        "CPRI_p90": round(cpri_p90, 4),
        "CPRI_p99": round(cpri_p99, 4),
        "S_max": round(s_max, 4),
        "pct_grade_A": round(100.0 * grade_counts["A"] / n_graded, 3)
        if n_graded else np.nan,
        "pct_grade_AB": round(
            100.0 * (grade_counts["A"] + grade_counts["B"]) / n_graded, 3)
        if n_graded else np.nan,
        "grade_A": grade_counts["A"],
        "grade_B": grade_counts["B"],
        "grade_C": grade_counts["C"],
        "grade_D": grade_counts["D"],
        "grade_E": grade_counts["E"],
        "risk_load_AB": round(ab_load, 1),
        "critical_ARG_load": round(hr_load, 1),
        "risk_load_AB_per_Mb": round(1e6 * ab_load / denom_bp, 3)
        if denom_bp else np.nan,
        "critical_ARG_load_per_Mb": round(1e6 * hr_load / denom_bp, 3)
        if denom_bp else np.nan,
        "weighting": weighting,
        "note": "Descriptive extension; CPRI was not trained or validated "
                "on metagenomic data",
    }
    return row


# ---------------------------------------------------------------------------
# CLI:  plasrisk meta ...
# ---------------------------------------------------------------------------

def _collect_fasta(paths: List[str]) -> List[str]:
    fasta = []
    exts = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
    for p in paths:
        if os.path.isdir(p):
            for root, _d, files in os.walk(p):
                for f in sorted(files):
                    if f.lower().endswith(exts):
                        fasta.append(os.path.join(root, f))
        elif any(ch in p for ch in "*?"):
            import glob
            fasta.extend(sorted(glob.glob(p)))
        else:
            fasta.append(p)
    return fasta


def _sample_id_from_path(path: str) -> str:
    base = os.path.basename(path)
    for suf in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"):
        if base.lower().endswith(suf):
            return base[: -len(suf)]
    return os.path.splitext(base)[0]


def build_meta_parser() -> "argparse.ArgumentParser":
    import argparse
    p = argparse.ArgumentParser(
        prog="plasrisk meta",
        description="Metagenomic three-layer PlasRisk scoring "
                    "(contig / plasmid-bin / community sample).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+",
                   help="assembled contig FASTA file(s) or a directory")
    p.add_argument("-o", "--output", default="plasrisk_meta",
                   help="output directory")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--min-id", type=float, default=75.0)
    p.add_argument("--min-cov", type=float, default=50.0)
    p.add_argument("--db", default=None,
                   help="abricate database root (default: bundled lookup)")
    p.add_argument("--no-abricate", action="store_true",
                   help="skip annotation (sequence / header evidence only)")
    p.add_argument("--aware-version", choices=["legacy", "2025"],
                   default="legacy",
                   help="WHO AWaRe version used by S_ARG")
    p.add_argument("--min-contig", type=int, default=1000,
                   help="ignore contigs shorter than this (bp)")
    p.add_argument("--complete-length", type=int, default=10000,
                   help="single contig length at/above which S_SIZE is "
                        "assessed without a circular flag (bp)")
    p.add_argument("--bins", default=None,
                   help="generic contig-to-bin TSV/CSV "
                        "(seq_id,bin_id[,circular,plasmid_score])")
    p.add_argument("--genomad-bin", default=None, metavar="DIR",
                   help="geNomad aggregate output directory "
                        "(*_plasmid_components.tsv / _plasmid_summary.tsv)")
    p.add_argument("--min-plasmid-score", type=float, default=0.0,
                   help="drop bins below this geNomad plasmid_score "
                        "(0 keeps all; geNomad recommends 0.7)")
    p.add_argument("--coverage", default=None,
                   help="per-contig coverage TSV/CSV (seq_id,coverage); "
                        "jgi_summarize_bam_contig_depths output accepted")
    p.add_argument("--coverage-column", default=None,
                   help="coverage column name for multi-sample coverage files")
    p.add_argument("--sample-id", default=None,
                   help="sample identifier (single FASTA input only)")
    p.add_argument("--total-assembly-bp", type=int, default=None,
                   help="total assembly length for per-Mb normalisation "
                        "(default: sum of plasmid bin lengths)")
    p.add_argument("--json", action="store_true",
                   help="also write meta_sample_summary.json")
    p.add_argument("-q", "--quiet", action="store_true")
    return p


def meta_main(argv: Optional[List[str]] = None) -> int:
    import json
    args = build_meta_parser().parse_args(argv)

    from .lookup import load_replicon_lookup
    lookup = None
    if not args.no_abricate:
        try:
            lookup = load_replicon_lookup()
        except Exception:
            lookup = None
    scorer = PlasRiskScorer(replicon_lookup=lookup, mode="lite",
                            aware_version=args.aware_version)

    fastas = _collect_fasta(args.inputs)
    if not fastas:
        print("[plasrisk meta] no FASTA inputs found", file=sys.stderr)
        return 2
    os.makedirs(args.output, exist_ok=True)

    bins_df = None
    if args.genomad_bin:
        bins_df = load_genomad(args.genomad_bin)
        if not args.quiet:
            print("[plasrisk meta] geNomad bins: %d contig-bin assignments"
                  % len(bins_df))
    if args.bins:
        generic = load_bins_table(args.bins)
        bins_df = generic if bins_df is None else pd.concat(
            [bins_df, generic], ignore_index=True).drop_duplicates(
            ["seq_id", "bin_id"])
    coverage_df = (load_coverage_table(args.coverage,
                                       args.coverage_column)
                   if args.coverage else None)

    all_contig, all_bin, sample_rows = [], [], []
    for fi, fasta in enumerate(fastas):
        sid = (args.sample_id if args.sample_id and len(fastas) == 1
               else _sample_id_from_path(fasta))
        if not args.quiet:
            print("[plasrisk meta] (%d/%d) %s" % (fi + 1, len(fastas), sid))

        if args.no_abricate:
            records = parse_fasta(fasta)
            from .annotate import PlasmidFeatures as _PF
            feats = []
            for seq_id, seq in records.items():
                feats.append(_PF(seq_id=seq_id, length_bp=len(seq)))
            annotation = AnnotationResult(file=fasta, features=feats)
        else:
            annotation = annotate_fasta(
                fasta, db_path=args.db, min_id=args.min_id,
                min_cov=args.min_cov, threads=args.threads,
                quiet=True)

        circ_header = circular_map_from_fasta(fasta)
        feats_kept = {f.seq_id: f for f in annotation.features
                      if f.length_bp >= args.min_contig}
        contig_df = score_contigs(annotation, scorer,
                                  complete_length=args.complete_length,
                                  min_length=args.min_contig,
                                  circular_map=circ_header)
        contig_df.insert(0, "sample_id", sid)

        bin_df = aggregate_bins(
            contig_df, feats_kept, bins_df, coverage_df, scorer,
            min_plasmid_score=args.min_plasmid_score,
            complete_length=args.complete_length,
            header_circular_map=circ_header)
        bin_df.insert(0, "sample_id", sid)

        total_bp = args.total_assembly_bp
        if total_bp is None:
            total_bp = int(sum(f.length_bp for f in feats_kept.values()))
        sample_row = aggregate_sample(sid, bin_df, contig_df, total_bp)
        all_contig.append(contig_df)
        all_bin.append(bin_df)
        sample_rows.append(sample_row)

    contig_all = pd.concat(all_contig, ignore_index=True)
    bin_all = pd.concat(all_bin, ignore_index=True)
    sample_all = pd.DataFrame(sample_rows)

    cpath = os.path.join(args.output, "meta_contig_scores.tsv")
    bpath = os.path.join(args.output, "meta_bin_scores.tsv")
    spath = os.path.join(args.output, "meta_sample_summary.tsv")
    contig_all.to_csv(cpath, sep="\t", index=False)
    bin_all.to_csv(bpath, sep="\t", index=False)
    sample_all.to_csv(spath, sep="\t", index=False)
    if args.json:
        sample_all.to_json(os.path.join(
            args.output, "meta_sample_summary.json"),
            orient="records", indent=2)

    if not args.quiet:
        print("[plasrisk meta] samples: %d | contig elements: %d | "
              "plasmid bins: %d" % (len(sample_all), len(contig_all),
                                    len(bin_all)))
        print("[plasrisk meta] wrote %s" % cpath)
        print("[plasrisk meta] wrote %s" % bpath)
        print("[plasrisk meta] wrote %s" % spath)
    return 0
