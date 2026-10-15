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
                  coverage: Optional[float],
                  complete_length: int = 10_000,
                  circ_contigs: bool = False) -> Dict:
    """Aggregate contig-level features of one predicted plasmid bin."""
    feats = sorted(feats, key=lambda f: f.seq_id)
    n_contigs = len(feats)
    total_len = int(sum(f.length_bp for f in feats))

    arg_names = _union([f.arg_names for f in feats])
    vf_names = _union([f.vf_names for f in feats])
    vf_cats = _union([f.vf_categories for f in feats])
    bm_names = _union([f.bm_gene_names for f in feats])
    replicons = _union([[r] for f in feats for r in
                        re.split(r"[;,/]", f.replicon) if r.strip()])

    # Merge per-gene product / resistance metadata
    arg_products, arg_resistance = {}, {}
    for f in feats:
        arg_products.update(f.metadata.get("arg_products", {}))
        arg_resistance.update(f.metadata.get("arg_resistance", {}))

    has_t4cp = any(bool(f.has_t4cp) for f in feats)
    has_rel = any(bool(f.has_relaxase) for f in feats)
    has_oriT = any(bool(f.has_oriT) for f in feats)
    has_aux = any(bool(f.has_auxiliary) for f in feats)
    has_int = any(bool(f.has_integron) for f in feats)
    # NOTE: IS copies split across contig boundaries are counted once per
    # contig conservatively; summing can slightly over-count fragmented ISs.
    n_is = int(sum(f.n_is for f in feats))
    any_mob_evidence = (has_t4cp or has_rel or has_oriT or has_aux or has_int
                        or n_is > 0)
    any_mob_annotated = any(f.mobility_annotated for f in feats)

    circ = (bool(circular) or bool(circ_contigs)
            or any(is_circular_header(f.seq_id) for f in feats))
    # S_SIZE is assessed only when the bin plausibly spans a whole plasmid:
    # circularised, or an aggregate long enough to be a near-complete
    # assembly. Fragment lengths never feed S_SIZE.
    size_assessable = circ or total_len >= complete_length
    complete = circ or (n_contigs == 1 and total_len >= complete_length)

    agg = PlasmidFeatures(
        seq_id=bin_id,
        length_bp=total_len,
        arg_names=arg_names,
        vf_names=vf_names,
        vf_categories=vf_cats,
        bm_gene_names=bm_names,
        replicon=";".join(replicons),
        has_t4cp=has_t4cp,
        has_relaxase=has_rel,
        has_oriT=has_oriT,
        has_auxiliary=has_aux,
        has_integron=has_int,
        n_is=n_is,
        mobility_annotated=any_mob_annotated or any_mob_evidence,
        metadata={"arg_products": arg_products,
                  "arg_resistance": arg_resistance},
    )

    s_arg = scorer.score_arg(agg)
    s_bmg = scorer.score_bmg(agg)
    s_vf = scorer.score_vf(agg)
    s_size = scorer.score_size(total_len) if size_assessable else None
    s_mob = mob_evidence_score(
        agg, complete=size_assessable or any_mob_annotated)
    s_rep = scorer.score_rep(agg) if agg.replicon else None

    cargo_vals = {"S_ARG": s_arg, "S_BMG": s_bmg, "S_VF": s_vf}
    s_cargo3, _dims3 = _renorm(cargo_vals)
    s_frag4, frag4_dims = None, []
    if s_mob is not None:
        f4 = dict(cargo_vals)
        f4["S_MOB"] = s_mob
        s_frag4, frag4_dims = _renorm(f4)
    s_lite, dims = None, []
    if s_size is not None and s_mob is not None:
        lv = dict(cargo_vals)
        lv["S_MOB"] = s_mob
        lv["S_SIZE"] = s_size
        s_lite, dims = _renorm(lv)

    if s_lite is not None and set(dims) == set(LITE_DIMENSIONS):
        grade, grade_label = _grade(s_lite)
    else:
        grade, grade_label = "", ""

    if circ:
        confidence = "complete"
    elif total_len >= 10_000:
        confidence = "probable"
    else:
        confidence = "fragmented"

    # Mobility class label (best evidence across contigs)
    mob_for_label = s_mob if s_mob is not None else 0.0
    if mob_for_label >= 0.85:
        mob_class = "conjugative_complete"
    elif mob_for_label >= 0.60:
        mob_class = "conjugative_likely"
    elif mob_for_label >= 0.35:
        mob_class = "mobilizable"
    else:
        mob_class = "non-mobilizable" if s_mob is not None else "unassessed"

    return {
        "bin_id": bin_id,
        "n_contigs": n_contigs,
        "total_length_bp": total_len,
        "total_length_kb": round(total_len / 1000.0, 2),
        "circular": circ,
        "confidence": confidence,
        "plasmid_score": plasmid_score,
        "coverage": coverage,
        "replicon": agg.replicon or "Unknown",
        "n_ARG": len(arg_names),
        "n_VF": len(vf_names),
        "n_BM": len(bm_names),
        "n_IS": n_is,
        "high_risk_genes": ";".join(agg.high_risk_args)
        if agg.high_risk_args else "",
        "mobility_class": mob_class,
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
    """Parse geNomad end-to-end outputs into the generic bins table.

    geNomad classifies each sequence independently, so every plasmid
    sequence is its own bin. Reads, searched anywhere under
    ``genomad_dir``:

      * ``*_plasmid_summary.tsv`` (preferred): one row per plasmid, with
        seq_name, topology and plasmid_score;
      * otherwise ``*_aggregated_classification.tsv``: per-sequence calls,
        filtered to the plasmid class;
      * legacy ``*_plasmid_components.tsv`` from older releases that
        grouped several contigs into one plasmid component.
    """
    comp_files, summ_files, agg_files = [], [], []
    for root, _dirs, files in os.walk(genomad_dir):
        for f in files:
            if f.endswith("_plasmid_components.tsv"):
                comp_files.append(os.path.join(root, f))
            elif f.endswith("_plasmid_summary.tsv"):
                summ_files.append(os.path.join(root, f))
            elif f.endswith("_aggregated_classification.tsv"):
                agg_files.append(os.path.join(root, f))

    def _summary_rows(sdf: "pd.DataFrame") -> List[dict]:
        name_col = _find_column(list(sdf.columns),
                                ("seq_name", "seq_id", "name", "bin_id"))
        if name_col is None:
            return []
        topo_col = _find_column(list(sdf.columns),
                                ("topology", "circular", "circ"))
        score_col = _find_column(list(sdf.columns),
                                 ("plasmid_score", "score"))
        out = []
        for _, r in sdf.iterrows():
            name = str(r[name_col])
            if name in ("", "nan"):
                continue
            circ = (topo_col is not None
                    and str(r[topo_col]).strip().lower()
                    in ("circular", "1", "true", "yes"))
            score = np.nan
            if score_col is not None and str(r[score_col]) not in ("", "nan"):
                try:
                    score = float(r[score_col])
                except (ValueError, TypeError):
                    score = np.nan
            out.append({"seq_id": name, "bin_id": name,
                        "circular": bool(circ), "plasmid_score": score})
        return out

    rows: List[dict] = []

    if summ_files:
        sdf = pd.read_csv(summ_files[0], sep="\t", dtype=str)
        rows = _summary_rows(sdf)

    if not rows and agg_files:
        adf = pd.read_csv(agg_files[0], sep="\t", dtype=str)
        lower_cols = {str(c).lower().strip(): c for c in adf.columns}
        name_col = lower_cols.get("seq_name") or lower_cols.get("seq_id") \
            or lower_cols.get("name")
        # Prefer an exact boolean "plasmid" column so it is never confused
        # with "plasmid_score"; then a categorical class column.
        bool_col = lower_cols.get("plasmid")
        class_col = None
        if bool_col is None:
            class_col = lower_cols.get("class") \
                or lower_cols.get("classification") \
                or lower_cols.get("prediction") or lower_cols.get("label")
        score_col = lower_cols.get("plasmid_score") \
            or lower_cols.get("score")
        chrom_col = lower_cols.get("chromosome_score")
        virus_col = lower_cols.get("virus_score")
        true_set = {"1", "true", "yes", "y"}
        if name_col is not None:
            for _, r in adf.iterrows():
                name = str(r[name_col])
                if name in ("", "nan"):
                    continue
                keep = False
                if bool_col is not None:
                    keep = str(r[bool_col]).strip().lower() in true_set
                elif class_col is not None:
                    keep = "plasmid" in str(r[class_col]).lower()
                elif score_col is not None:
                    try:
                        ps = float(r[score_col])
                        cs = float(r[chrom_col]) if chrom_col else -1.0
                        vs = float(r[virus_col]) if virus_col else -1.0
                        keep = ps >= max(cs, vs) and ps >= 0.5
                    except (ValueError, TypeError):
                        keep = False
                score = np.nan
                if score_col is not None \
                        and str(r[score_col]) not in ("", "nan"):
                    try:
                        score = float(r[score_col])
                    except (ValueError, TypeError):
                        score = np.nan
                if keep:
                    rows.append({"seq_id": name, "bin_id": name,
                                 "circular": False,
                                 "plasmid_score": score})

    # Legacy multi-contig components (older geNomad releases).
    if not rows and comp_files:
        cdf = pd.read_csv(comp_files[0], sep="\t", dtype=str)
        seq_col = _find_column(list(cdf.columns), ("seq_name", "seq_id"))
        comp_col = _find_column(list(cdf.columns),
                                ("components", "component"))
        if seq_col is not None and comp_col is not None:
            for _, r in cdf.iterrows():
                cell = str(r[comp_col])
                if cell in ("", "nan"):
                    continue
                for i, group in enumerate(cell.split(";")):
                    contigs = [x.strip() for x in group.split(",")
                               if x.strip()]
                    bin_id = "%s_component_%d" % (str(r[seq_col]), i + 1)
                    for contig in contigs:
                        rows.append({"seq_id": contig, "bin_id": bin_id,
                                     "circular": False,
                                     "plasmid_score": np.nan})

    if not rows:
        raise FileNotFoundError(
            "no geNomad plasmid output found under %s (expected "
            "*_plasmid_summary.tsv, *_aggregated_classification.tsv or "
            "*_plasmid_components.tsv)" % genomad_dir)

    bins = pd.DataFrame(rows)
    return bins.drop_duplicates(["seq_id", "bin_id"])


def _jgi_sample_columns(columns) -> List[str]:
    """Per-sample depth columns in a jgi_summarize table: contigName, length,
    totalAvgDepth and each ``<sample>-var`` variance column are excluded."""
    reserved = {"contigname", "length", "totalavgdepth"}
    out = []
    for c in columns:
        cs = str(c).strip().lower()
        if cs in reserved or cs.endswith("-var"):
            continue
        out.append(c)
    return out


def _cov_column_basename(name: str) -> str:
    """Normalize a jgi column label: drop directory and a trailing .bam."""
    base = os.path.basename(str(name)).strip()
    if base.endswith(".bam"):
        base = base[:-4]
    return base


def _resolve_coverage_column(columns, requested: Optional[str]) -> str:
    """Resolve a coverage column robustly.

    jgi_summarize_bam_contig_depths names the depth column after the BAM
    argument verbatim, which may be a full path or carry a ``.bam`` suffix.
    Accept an exact match, then a basename match, then the first jgi sample
    column, so callers can pass the bare accession without KeyError."""
    cols = [c for c in columns]
    if requested is not None:
        if requested in cols:
            return requested
        target = _cov_column_basename(requested)
        cand = [c for c in cols
                if not str(c).endswith("-var")
                and _cov_column_basename(c) == target]
        if cand:
            return cand[0]
    sample = _jgi_sample_columns(cols)
    if sample:
        return sample[0]
    raise KeyError(
        "could not find coverage column %r; available columns: %s"
        % (requested, list(cols)))


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
        cov_col = _resolve_coverage_column(df.columns, coverage_column)
    else:
        cov_col = _find_column(list(df.columns),
                               ("coverage", "cov", "mean_cov",
                                "meandepth", "totalavgdepth", "depth"))
        if cov_col is None:
            numeric = df.select_dtypes(include=[np.number]).columns
            if len(numeric) == 0:
                raise ValueError("no numeric coverage column found in %s"
                                 % path)
            cov_col = numeric[0]
    out = pd.DataFrame({"seq_id": df[seq_col].astype(str),
                        "coverage": pd.to_numeric(df[cov_col],
                                                   errors="coerce")})
    return out.dropna(subset=["coverage"])


def aggregate_bins(contig_df: pd.DataFrame,
                   features_by_id: Dict[str, PlasmidFeatures],
                   bins_df: Optional[pd.DataFrame],
                   coverage_df: Optional[pd.DataFrame],
                   scorer: PlasRiskScorer,
                   min_plasmid_score: float = 0.0,
                   complete_length: int = 10_000,
                   header_circular_map: Optional[Dict[str, bool]] = None
                   ) -> pd.DataFrame:
    """Layer 2: build one score row per predicted plasmid bin.

    Contigs absent from ``bins_df`` become singleton bins (flagged
    ``unbinned`` via the bin_id prefix)."""
    cov_map = (dict(zip(coverage_df["seq_id"], coverage_df["coverage"]))
               if coverage_df is not None else {})

    if bins_df is not None:
        bins_df = bins_df.copy()
        if min_plasmid_score > 0:
            bins_df = bins_df[
                bins_df["plasmid_score"].isna()
                | (bins_df["plasmid_score"] >= min_plasmid_score)]
        bin_map = dict(zip(bins_df["seq_id"], bins_df["bin_id"]))
        circ_map = bins_df.groupby("bin_id")["circular"].any().to_dict()
        score_map = (bins_df.dropna(subset=["plasmid_score"])
                     .groupby("bin_id")["plasmid_score"].max().to_dict())
    else:
        bin_map, circ_map, score_map = {}, {}, {}

    groups: Dict[str, List[PlasmidFeatures]] = {}
    for cid, feat in features_by_id.items():
        if cid in bin_map:
            groups.setdefault(bin_map[cid], []).append(feat)
        else:
            groups.setdefault("unbinned:%s" % cid, []).append(feat)

    rows = []
    header_map = header_circular_map or {}
    for bin_id, feats in groups.items():
        circ = circ_map.get(bin_id)
        circ_contigs = any(header_map.get(f.seq_id, False) for f in feats)
        pscore = score_map.get(bin_id)
        covs = [cov_map[f.seq_id] for f in feats if f.seq_id in cov_map]
        if covs:
            lens = np.array([f.length_bp for f in feats], dtype=float)
            cs = np.array(covs, dtype=float)
            coverage = float(np.sum(cs * lens) / np.sum(lens))
        else:
            coverage = None
        rows.append(aggregate_bin(bin_id, feats, scorer, circ, pscore,
                                  coverage,
                                  complete_length=complete_length,
                                  circ_contigs=circ_contigs))
    df = pd.DataFrame(rows)
    if not df.empty:
        # Bins from the external caller (geNomad/components) are predicted
        # plasmids; singleton bins built for contigs absent from the bin map
        # are an unbinned sensitivity layer (dominated by chromosome
        # fragments in a whole-community assembly).
        df["bin_origin"] = np.where(
            df["bin_id"].astype(str).str.startswith("unbinned:"),
            "unbinned", "predicted")
        sort_col = df["S_lite"].combine_first(df["S_frag4"]) \
            .combine_first(df["S_cargo3"]).fillna(-1)
        df = df.assign(_s=sort_col).sort_values("_s", ascending=False)
        df = df.drop(columns="_s").reset_index(drop=True)
        df.insert(0, "rank", range(1, len(df) + 1))
    return df


# ---------------------------------------------------------------------------
# Layer 3: sample / community summary
# ---------------------------------------------------------------------------

def aggregate_sample(sample_id: str, bin_df: pd.DataFrame,
                     contig_df: pd.DataFrame,
                     total_assembly_bp: Optional[int] = None) -> Dict:
    """Descriptive community-level statistics (CPRI = community plasmid
    risk index; not trained/validated in the isolate study).

    Unprefixed columns are computed over all scored elements (predicted
    bins plus unbinned singletons). Columns with the ``pred_`` prefix are
    computed over caller-predicted plasmid bins only (bin_origin ==
    'predicted'); these are the primary plasmid-level statistics and are
    not diluted by the typically much larger set of chromosome
    fragments."""
    b = bin_df.copy()
    s_use = b["S_lite"] if "S_lite" in b.columns else pd.Series(
        np.nan, index=b.index)
    for col in ("S_frag4", "S_cargo3"):
        if col in b.columns:
            s_use = s_use.combine_first(b[col])
    b["S_use"] = s_use

    n_bins = len(b)
    n_predicted = (
        int((b["bin_origin"] == "predicted").sum())
        if "bin_origin" in b.columns else n_bins)
    n_unbinned = n_bins - n_predicted
    n_complete = int((b["confidence"] == "complete").sum())
    n_probable = int((b["confidence"] == "probable").sum())
    n_fragmented = int((b["confidence"] == "fragmented").sum())

    denom_bp = (float(total_assembly_bp) if total_assembly_bp
                else float(np.nansum(b["total_length_bp"].to_numpy(
                    dtype=float))))

    def _block(d, prefix):
        """Statistics for one bin subset (all bins, or predicted-plasmid
        bins only). Key names carry ``prefix`` when it is non-empty."""
        sd = d.dropna(subset=["S_use"]).copy()
        gd = d[d["grade"] != ""]
        n_s = len(sd)
        n_g = len(gd)

        weights = pd.Series(1.0, index=sd.index)
        weighting = "unweighted"
        if "coverage" in sd.columns and sd["coverage"].notna().any():
            cov = sd["coverage"].fillna(0).to_numpy(dtype=float)
            length = sd["total_length_bp"].to_numpy(dtype=float)
            wv = cov * length
            if wv.sum() > 0:
                weights = pd.Series(wv, index=sd.index)
                weighting = "coverage_x_length"
        else:
            length = sd["total_length_bp"].to_numpy(dtype=float)
            if length.sum() > 0:
                weights = pd.Series(length, index=sd.index)
                weighting = "length"

        sv = sd["S_use"].to_numpy(dtype=float)
        wv = weights.to_numpy(dtype=float)

        def _wq(q):
            if len(sv) == 0:
                return np.nan
            order = np.argsort(sv)
            sw, ss = wv[order], sv[order]
            cw = np.cumsum(sw)
            if cw[-1] <= 0:
                return np.nan
            idx = np.searchsorted(cw, q * cw[-1])
            return ss[min(idx, len(ss) - 1)]

        cpri_mean = float(np.sum(wv * sv) / np.sum(wv)) if wv.sum() > 0 \
            else np.nan

        grade_counts = {g: int((gd["grade"] == g).sum())
                        for g in ("A", "B", "C", "D", "E")}
        ab = gd[gd["grade"].isin(("A", "B"))]
        # critical ARGs remain detectable on fragments; count them across
        # all scored bins, not only bins with a complete five-dimension grade
        hr = sd[sd["high_risk_genes"].fillna("") != ""]

        def _cov(dd):
            return np.nan_to_num(dd.get("coverage", pd.Series(
                [1.0] * len(dd))).to_numpy(dtype=float))

        ab_load = float(np.sum(
            ab["total_length_bp"].to_numpy(dtype=float) * _cov(ab)))
        hr_load = float(np.sum(
            hr["total_length_bp"].to_numpy(dtype=float) * _cov(hr)))

        p = prefix
        out = {
            "%sn_scored_bins" % p: int(n_s),
            "%sn_graded_bins" % p: int(n_g),
            "%spct_ARG_positive" % p: round(
                100.0 * (sd["n_ARG"] > 0).mean(), 2) if n_s else np.nan,
            "%spct_BMG_positive" % p: round(
                100.0 * (sd["n_BM"] > 0).mean(), 2) if n_s else np.nan,
            "%spct_high_risk_ARG" % p: round(
                100.0 * (sd["high_risk_genes"].fillna("") != "").mean(), 2)
            if n_s else np.nan,
            "%sCPRI_mean" % p: round(cpri_mean, 4)
            if not np.isnan(cpri_mean) else np.nan,
            "%sCPRI_p90" % p: round(float(_wq(0.90)), 4) if n_s else np.nan,
            "%sCPRI_p99" % p: round(float(_wq(0.99)), 4) if n_s else np.nan,
            "%sS_max" % p: round(float(np.nanmax(sv)), 4) if n_s else np.nan,
            "%spct_grade_A" % p: round(
                100.0 * grade_counts["A"] / n_g, 3) if n_g else np.nan,
            "%spct_grade_AB" % p: round(
                100.0 * (grade_counts["A"] + grade_counts["B"]) / n_g, 3)
            if n_g else np.nan,
            "%sgrade_A" % p: grade_counts["A"],
            "%sgrade_B" % p: grade_counts["B"],
            "%sgrade_C" % p: grade_counts["C"],
            "%sgrade_D" % p: grade_counts["D"],
            "%sgrade_E" % p: grade_counts["E"],
            "%srisk_load_AB" % p: round(ab_load, 1),
            "%scritical_ARG_load" % p: round(hr_load, 1),
            "%srisk_load_AB_per_Mb" % p: round(
                1e6 * ab_load / denom_bp, 3) if denom_bp else np.nan,
            "%scritical_ARG_load_per_Mb" % p: round(
                1e6 * hr_load / denom_bp, 3) if denom_bp else np.nan,
        }
        if not p:
            out["weighting"] = weighting
        return out

    row = {
        "sample_id": sample_id,
        "n_contig_elements": int(len(contig_df)),
        "n_plasmid_bins": int(n_bins),
        "n_predicted_bins": int(n_predicted),
        "n_unbinned_bins": int(n_unbinned),
        "n_complete": n_complete,
        "n_probable": n_probable,
        "n_fragmented": n_fragmented,
    }
    row.update(_block(b, ""))
    # Predicted-plasmid subset (caller-supplied bin map, e.g. geNomad):
    # primary statistics for plasmid risk, free of whole-community
    # chromosome-fragment dilution.
    if "bin_origin" in b.columns:
        bp = b[b["bin_origin"] == "predicted"]
        row.update(_block(bp, "pred_"))
    row["note"] = (
        "Descriptive extension; CPRI was not trained or validated on "
        "metagenomic data. Unprefixed metrics use all elements; pred_* "
        "metrics use only caller-predicted plasmid bins")
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
    p.add_argument("--strict-annotation", action="store_true",
                   help="exit 3 when abricate was requested but no "
                        "annotation database was available (prevents "
                        "silently producing all-zero ARG/VF results)")
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
                fasta, use_abricate=True, min_id=args.min_id,
                min_cov=args.min_cov, threads=args.threads,
                lookup=lookup)
            for w in annotation.warnings:
                print("[plasrisk meta][annotation] %s" % w,
                      file=sys.stderr)
            if not annotation.databases_used:
                print("[plasrisk meta] WARNING: no abricate annotation "
                      "database available; ARG/VF/BMG/replicon columns "
                      "will be empty for %s (check `which abricate`, "
                      "`abricate --list` and the BLAST dependency)"
                      % sid, file=sys.stderr)
                if args.strict_annotation:
                    return 3
            else:
                n_annotated = sum(
                    1 for f in annotation.features
                    if f.arg_names or f.vf_names or f.bm_gene_names
                    or f.replicon or f.n_is > 0 or f.has_integron
                    or f.has_relaxase or f.has_t4cp)
                if n_annotated == 0 and len(annotation.features) > 100:
                    print("[plasrisk meta] WARNING: databases %s were "
                          "queried but produced zero hits across %d "
                          "sequences; abricate likely failed at runtime "
                          "(check blastn and database installation)"
                          % (",".join(annotation.databases_used),
                             len(annotation.features)), file=sys.stderr)
                elif not args.quiet:
                    print("[plasrisk meta] annotation dbs: %s | annotated "
                          "elements: %d/%d"
                          % (",".join(annotation.databases_used),
                             n_annotated, len(annotation.features)))

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

