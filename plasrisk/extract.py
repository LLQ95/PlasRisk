"""
extract.py - Raw-assembly front end for PlasRisk (v1.6.0).

Most real-world FASTA inputs are not curated plasmid collections. They are
short-read whole-genome assemblies, long-read assemblies, or multi-strain
metagenomic assemblies in which the plasmid sequences have to be predicted
first. This module:

  1. classifies every contig as plasmid / chromosome / unclassified using
     PlasFlow and/or geNomad (external tools invoked through subprocess when
     available; precomputed output tables can also be supplied);
  2. extracts the predicted plasmid contigs to a FASTA file;
  3. routes long/circular plasmid contigs to the isolate scorer and short
     fragmented contigs / bins to the metagenomic three-layer;
  4. reports assemblies with no detectable plasmid as UNSCORABLE, with a
     reason, instead of forcing a chromosome through a plasmid score.

The frozen weights are never refit here.
"""

from __future__ import annotations

import glob as _glob
import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from .annotate import annotate_fasta, parse_fasta
from .lookup import load_replicon_lookup
from .scoring import PlasRiskScorer, PlasmidFeatures

# A contig is treated as a chromosome (rather than a small unclassifiable
# fragment) when the assembly or its longest contig reaches this length.
MIN_CHROMOSOME_TOTAL_BP = 500_000
LONG_CONTIG_BP = 500_000
DEFAULT_MIN_LENGTH_BP = 1_000
DEFAULT_COMPLETE_LENGTH_BP = 10_000


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ClassifiedContig:
    seq_id: str
    length_bp: int
    label: str  # 'plasmid', 'chromosome', 'unclassified'
    plasmid_prob: Optional[float] = None
    circular: bool = False
    source: str = ""  # classifier that produced the label


@dataclass
class ExtractionResult:
    status: str  # 'scored' or 'unscorable'
    input_fasta: str
    sample_id: str
    n_contigs: int = 0
    total_length_bp: int = 0
    longest_contig_bp: int = 0
    n_plasmid: int = 0
    n_chromosome: int = 0
    n_unclassified: int = 0
    classifier: str = ""
    threshold: float = 0.7
    route: str = ""  # 'isolate', 'meta', or '' when unscorable
    reason: str = ""
    plasmid_fasta: str = ""
    output_files: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Table parsing
# ---------------------------------------------------------------------------

def _read_table(path: str) -> pd.DataFrame:
    """Read a classifier table, tolerating BOMs, comments and any separator."""
    df = pd.read_csv(
        path, sep=None, engine="python", comment="#",
        dtype=str, encoding="utf-8-sig")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _normalize_label(value) -> str:
    if value is None:
        return "unclassified"
    v = str(value).strip().lower()
    if "plasmid" in v:
        return "plasmid"
    if "chromosom" in v:
        return "chromosome"
    return "unclassified"


def _bp(value) -> Optional[int]:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def _prob(value) -> Optional[float]:
    try:
        p = float(str(value).strip())
        if 0.0 <= p <= 1.0:
            return p
        if 0.0 < p <= 100.0:
            return p / 100.0
    except (TypeError, ValueError):
        pass
    return None


def _col(df: pd.DataFrame, candidates: Sequence[str],
         contains: Optional[Sequence[str]] = None) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    lower = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if contains:
        for c in df.columns:
            cl = str(c).lower()
            if any(tok in cl for tok in contains):
                return c
    return None


def _col_contains(df: pd.DataFrame, token: str,
                  exclude: Sequence[str] = ()) -> Optional[str]:
    for c in df.columns:
        cl = str(c).lower()
        if token in cl and not any(x in cl for x in exclude):
            return c
    return None


def parse_plasflow_table(path: str) -> List[ClassifiedContig]:
    """Parse a PlasFlow predictions TSV.

    Column names (not positions) are used. The contig id falls back across
    common names; the plasmid probability column is any column containing
    'plasmid' that is not the id/label/length column.
    """
    df = _read_table(path)
    id_c = _col(df, ("contig_name", "contig_id", "seq_name", "name", "id"))
    len_c = _col(df, ("contig_length", "length", "length_bp", "len"))
    label_c = _col(df, ("label", "prediction", "class"))
    prob_c = _col_contains(df, "plasmid",
                           exclude=("name", "label", "length", "id", "len"))
    if id_c is None:
        raise ValueError("PlasFlow table has no contig id column: %s"
                         % ", ".join(df.columns))
    out: List[ClassifiedContig] = []
    for _, row in df.iterrows():
        sid = str(row[id_c]).strip()
        length = _bp(row[len_c]) if len_c else None
        label = _normalize_label(row[label_c]) if label_c else "unclassified"
        prob = _prob(row[prob_c]) if prob_c else None
        if label == "unclassified" and prob is not None and prob >= 0.5:
            label = "plasmid"
        out.append(ClassifiedContig(
            seq_id=sid,
            length_bp=length if length is not None else 0,
            label=label, plasmid_prob=prob, source="plasflow"))
    return out


def parse_genomad_summary(path: str) -> List[ClassifiedContig]:
    """Parse a geNomad plasmid summary TSV, or a directory containing one."""
    if os.path.isdir(path):
        hits = sorted(_glob.glob(os.path.join(path, "**", "*_plasmid_summary.tsv"),
                                 recursive=True))
        if not hits:
            hits = sorted(_glob.glob(os.path.join(path, "**", "*_summary.tsv"),
                                     recursive=True))
        if not hits:
            raise FileNotFoundError(
                "no geNomad *_plasmid_summary.tsv under %s" % path)
        path = hits[0]
    df = _read_table(path)
    id_c = _col(df, ("seq_name", "seq_id", "name"))
    len_c = _col(df, ("length", "length_bp", "seq_length"))
    pscore_c = _col_contains(df, "plasmid_score")
    cscore_c = _col_contains(df, "chromosome_score")
    topo_c = _col(df, ("topology", "circular"))
    if id_c is None:
        raise ValueError("geNomad summary has no seq_name column")
    out: List[ClassifiedContig] = []
    for _, row in df.iterrows():
        sid = str(row[id_c]).strip()
        length = _bp(row[len_c]) if len_c else 0
        pscore = _prob(row[pscore_c]) if pscore_c else None
        cscore = _prob(row[cscore_c]) if cscore_c else None
        circular = False
        if topo_c is not None and str(row[topo_c]).strip().lower() in (
                "circular", "true", "yes", "y"):
            circular = True
        if pscore is not None and cscore is not None:
            label = "plasmid" if pscore >= cscore else "chromosome"
        elif pscore is not None:
            label = "plasmid" if pscore >= 0.5 else "unclassified"
        else:
            label = "unclassified"
        out.append(ClassifiedContig(
            seq_id=sid, length_bp=length or 0, label=label,
            plasmid_prob=pscore, circular=circular, source="genomad"))
    return out


def parse_generic_table(path: str) -> List[ClassifiedContig]:
    """Parse any id + label table, optionally with plasmid_score/length."""
    df = _read_table(path)
    id_c = _col(df, ("seq_id", "seq_name", "contig_id", "contig_name",
                      "id", "name"))
    label_c = _col(df, ("label", "prediction", "class"))
    len_c = _col(df, ("length", "length_bp", "contig_length", "len"))
    pscore_c = _col(df, ("plasmid_score", "plasmid_prob", "prob"))
    if id_c is None or label_c is None:
        raise ValueError(
            "generic table needs id and label columns: %s"
            % ", ".join(df.columns))
    out: List[ClassifiedContig] = []
    for _, row in df.iterrows():
        length = _bp(row[len_c]) if len_c else 0
        out.append(ClassifiedContig(
            seq_id=str(row[id_c]).strip(),
            length_bp=length or 0,
            label=_normalize_label(row[label_c]),
            plasmid_prob=_prob(row[pscore_c]) if pscore_c else None,
            source="generic"))
    return out


# ---------------------------------------------------------------------------
# External classifier runners (optional; tools are not bundled)
# ---------------------------------------------------------------------------

def plasflow_available() -> bool:
    return shutil.which("PlasFlow.py") is not None or shutil.which(
        "plasflow") is not None


def genomad_available() -> bool:
    return shutil.which("genomad") is not None


def _write_fasta(records: Dict[str, str], path: str) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for sid, seq in records.items():
            out.write(">%s\n" % sid)
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def _filter_fasta(input_fasta: str, output_fasta: str,
                  min_length: int = DEFAULT_MIN_LENGTH_BP) -> int:
    records = parse_fasta(input_fasta)
    kept = {sid: seq for sid, seq in records.items() if len(seq) >= min_length}
    _write_fasta(kept, output_fasta)
    return len(kept)


def run_plasflow(input_fasta: str, workdir: str,
                 threshold: float = 0.7,
                 min_length: int = DEFAULT_MIN_LENGTH_BP,
                 threads: int = 4) -> str:
    """Run PlasFlow on a length-filtered copy and return its TSV path."""
    os.makedirs(workdir, exist_ok=True)
    filtered = os.path.join(workdir, "plasflow_input.fasta")
    _filter_fasta(input_fasta, filtered, min_length)
    out_tsv = os.path.join(workdir, "plasflow_predictions.tsv")
    base = ["PlasFlow.py", "--input", filtered, "--output", out_tsv,
            "--threshold", str(threshold)]
    try:
        subprocess.run(base + ["--batch_size", str(max(1, threads * 8))],
                       check=True, capture_output=True, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Older PlasFlow releases do not support --batch_size
        subprocess.run(base, check=True, capture_output=True, text=True)
    if not os.path.exists(out_tsv):
        # Some versions write into the output directory rather than a file
        alt = os.path.join(workdir, "plasflow_predictions.tsv",
                          "predictions.tsv")
        if os.path.exists(alt):
            shutil.copy(alt, out_tsv)
    return out_tsv


def run_genomad(input_fasta: str, workdir: str, db_path: str,
                threshold: float = 0.7, threads: int = 4) -> str:
    """Run geNomad classify and return the output directory."""
    if not db_path or not os.path.isdir(db_path):
        raise ValueError(
            "geNomad requires a database path (--genomad-db)")
    os.makedirs(workdir, exist_ok=True)
    cmd = ["genomad", "classify", "--threads", str(threads),
           "--cleanup", input_fasta, workdir, db_path]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return workdir


# ---------------------------------------------------------------------------
# Combining classifiers
# ---------------------------------------------------------------------------

def combine_classifications(sets: Sequence[List[ClassifiedContig]],
                            vote: str = "majority") -> List[ClassifiedContig]:
    """Merge per-classifier contig labels by id.

    vote='majority' requires more than half of the classifiers seeing a
    contig to call it a plasmid; vote='any' calls it a plasmid when any
    classifier does (high recall).
    """
    by_id: Dict[str, List[ClassifiedContig]] = {}
    for cls in sets:
        for c in cls:
            by_id.setdefault(c.seq_id, []).append(c)
    out: List[ClassifiedContig] = []
    for sid, votes_ in by_id.items():
        n = len(votes_)
        n_plasmid = sum(1 for v in votes_ if v.label == "plasmid")
        n_chr = sum(1 for v in votes_ if v.label == "chromosome")
        probs = [v.plasmid_prob for v in votes_
                 if v.plasmid_prob is not None]
        if vote == "any":
            label = "plasmid" if n_plasmid >= 1 else (
                "chromosome" if n_chr == n else "unclassified")
        else:
            label = "plasmid" if n_plasmid > n / 2 else (
                "chromosome" if n_chr > n / 2 else "unclassified")
        out.append(ClassifiedContig(
            seq_id=sid,
            length_bp=max(v.length_bp for v in votes_),
            label=label,
            plasmid_prob=(sum(probs) / len(probs) if probs else None),
            circular=any(v.circular for v in votes_),
            source="+".join(sorted({v.source for v in votes_ if v.source}))))
    return out


# ---------------------------------------------------------------------------
# Unscorable logic and routing
# ---------------------------------------------------------------------------

def _unscorable_reason(contigs: Sequence[ClassifiedContig],
                       n_plasmid: int) -> str:
    n_chr = sum(1 for c in contigs if c.label == "chromosome")
    total = sum(c.length_bp for c in contigs)
    longest = max((c.length_bp for c in contigs), default=0)
    has_chr_label = n_chr > 0
    long_enough = (total >= MIN_CHROMOSOME_TOTAL_BP
                   or longest >= LONG_CONTIG_BP)
    if has_chr_label and long_enough:
        return ("No plasmid sequence detected. The input appears to contain "
                "chromosomal sequence only (total length %d bp, longest "
                "contig %d bp). PlasRisk scores plasmids rather than "
                "chromosomes." % (total, longest))
    if long_enough:
        return ("No plasmid sequence detected at the current threshold, and "
                "no contig was confidently classified. Lower the classifier "
                "threshold or supply a different classifier output.")
    return ("Input is too small or fragmented to classify as plasmid or "
            "chromosome (total length %d bp, longest contig %d bp)." % (
                total, longest))


def decide_route(plasmid_contigs: Sequence[ClassifiedContig],
                 complete_length: int = DEFAULT_COMPLETE_LENGTH_BP) -> str:
    """Decide whether predicted plasmids go to the isolate or meta scorer."""
    n = len(plasmid_contigs)
    if n == 0:
        return ""
    n_complete = sum(1 for c in plasmid_contigs
                     if c.circular or c.length_bp >= complete_length)
    n_circular = sum(1 for c in plasmid_contigs if c.circular)
    if (n_circular > 0 and n <= 25) or (n_complete / n) >= 0.5:
        return "isolate"
    return "meta"


# ---------------------------------------------------------------------------
# FASTA handling with full (untruncated) headers
# ---------------------------------------------------------------------------

def _read_full_fasta(path: str) -> Dict[str, str]:
    """Read FASTA keyed by the full, untruncated header line."""
    records: Dict[str, str] = {}
    current = None
    buf: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(buf)
                current = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if current is not None:
        records[current] = "".join(buf)
    return records


def _build_id_resolver(full_headers: Sequence[str]) -> Dict[str, str]:
    """Map classifier ids (full header or first token) onto full headers."""
    resolver: Dict[str, str] = {}
    for header in full_headers:
        resolver[header] = header
        first = header.split()[0] if header else header
        resolver.setdefault(first, header)
    return resolver


# ---------------------------------------------------------------------------
# Scoring routes
# ---------------------------------------------------------------------------

def _score_isolate(plasmid_fasta: str, outdir: str, mode: str,
                   aware_version: str, lookup,
                   min_id: float, min_cov: float, threads: int,
                   use_abricate: bool) -> List[str]:
    ann = annotate_fasta(
        plasmid_fasta, use_abricate=use_abricate, min_id=min_id,
        min_cov=min_cov, threads=threads, lookup=lookup)
    if not ann.features:
        raise ValueError("no annotatable plasmid features in extracted set")
    scorer = PlasRiskScorer(replicon_lookup=lookup, mode=mode,
                            aware_version=aware_version)
    df = scorer.score_dataframe(ann.features)
    path = os.path.join(outdir, "plasrisk_isolate_scores.tsv")
    df.to_csv(path, sep="\t", index=False, float_format="%.4f")
    return [path]


def _score_meta(plasmid_fasta: str, outdir: str, bins_table: Optional[str],
                coverage_table: Optional[str], aware_version: str, lookup,
                min_id: float, min_cov: float, threads: int,
                use_abricate: bool,
                complete_length: int) -> List[str]:
    from .metagenome import (
        score_contigs, aggregate_bins, aggregate_sample,
        load_bins_table, load_coverage_table, circular_map_from_fasta)
    ann = annotate_fasta(
        plasmid_fasta, use_abricate=use_abricate, min_id=min_id,
        min_cov=min_cov, threads=threads, lookup=lookup)
    scorer = PlasRiskScorer(replicon_lookup=lookup, mode="lite",
                            aware_version=aware_version)
    circ = circular_map_from_fasta(plasmid_fasta)
    kept = {f.seq_id: f for f in ann.features if f.length_bp >= 1000}
    contig_df = score_contigs(ann, scorer, complete_length=complete_length,
                              min_length=1000, circular_map=circ)
    bins_df = load_bins_table(bins_table) if bins_table else None
    cov_df = (load_coverage_table(coverage_table)
              if coverage_table else None)
    bin_df = aggregate_bins(contig_df, kept, bins_df, cov_df, scorer,
                            complete_length=complete_length,
                            header_circular_map=circ)
    total_bp = int(sum(f.length_bp for f in kept.values()))
    sample = aggregate_sample("plasmids", bin_df, contig_df, total_bp)
    cpath = os.path.join(outdir, "meta_contig_scores.tsv")
    bpath = os.path.join(outdir, "meta_bin_scores.tsv")
    spath = os.path.join(outdir, "meta_sample_summary.tsv")
    contig_df.to_csv(cpath, sep="\t", index=False)
    bin_df.to_csv(bpath, sep="\t", index=False)
    pd.DataFrame([sample]).to_csv(spath, sep="\t", index=False)
    return [cpath, bpath, spath]


# ---------------------------------------------------------------------------
# Orchestration per sample
# ---------------------------------------------------------------------------

def process_one(input_fasta: str, outdir: str,
                classifier: str = "plasflow",
                threshold: float = 0.7,
                plasflow_table: Optional[str] = None,
                genomad_table: Optional[str] = None,
                generic_table: Optional[str] = None,
                genomad_db: Optional[str] = None,
                bins_table: Optional[str] = None,
                coverage_table: Optional[str] = None,
                mode: str = "lite",
                aware_version: str = "legacy",
                route: str = "auto",
                vote: str = "majority",
                min_id: float = 75.0,
                min_cov: float = 50.0,
                threads: int = 4,
                use_abricate: bool = True,
                complete_length: int = DEFAULT_COMPLETE_LENGTH_BP,
                sample_id: Optional[str] = None,
                write_json: bool = False) -> ExtractionResult:
    sample_id = sample_id or os.path.splitext(
        os.path.basename(input_fasta.rstrip("/\\")))[0]
    sample_out = os.path.join(outdir, sample_id)
    os.makedirs(sample_out, exist_ok=True)

    full_records = _read_full_fasta(input_fasta)
    full_headers = list(full_records.keys())
    total_length = sum(len(s) for s in full_records.values())
    longest = max((len(s) for s in full_records.values()), default=0)
    resolver = _build_id_resolver(full_headers)

    classifications: List[List[ClassifiedContig]] = []
    classifier_names: List[str] = []

    if plasflow_table:
        classifications.append(parse_plasflow_table(plasflow_table))
        classifier_names.append("plasflow")
    if genomad_table:
        classifications.append(parse_genomad_summary(genomad_table))
        classifier_names.append("genomad")
    if generic_table:
        classifications.append(parse_generic_table(generic_table))
        classifier_names.append("generic")

    if not classifications:
        workdir = os.path.join(sample_out, "classifier_work")
        if classifier in ("plasflow", "both") and plasflow_available():
            tsv = run_plasflow(input_fasta, workdir, threshold=threshold,
                               threads=threads)
            classifications.append(parse_plasflow_table(tsv))
            classifier_names.append("plasflow")
        if classifier in ("genomad", "both") and genomad_available():
            gdir = run_genomad(input_fasta, workdir, genomad_db or "",
                               threshold=threshold, threads=threads)
            classifications.append(parse_genomad_summary(gdir))
            classifier_names.append("genomad")
        if not classifications:
            raise RuntimeError(
                "No classifier output available. Install PlasFlow or geNomad, "
                "or supply a precomputed table.")

    merged = combine_classifications(classifications, vote=vote)

    # Align classifier ids back to the full FASTA headers and fill lengths.
    aligned: List[ClassifiedContig] = []
    for c in merged:
        full_header = resolver.get(c.seq_id, c.seq_id)
        seq = full_records.get(full_header)
        length = len(seq) if seq is not None else c.length_bp
        aligned.append(ClassifiedContig(
            seq_id=full_header, length_bp=length or c.length_bp,
            label=c.label, plasmid_prob=c.plasmid_prob,
            circular=c.circular or _circular_header(full_header),
            source=c.source))

    manifest = pd.DataFrame({
        "id": [c.seq_id for c in aligned],
        "length_bp": [c.length_bp for c in aligned],
        "label": [c.label for c in aligned],
        "plasmid_prob": [c.plasmid_prob if c.plasmid_prob is not None
                         else "" for c in aligned],
        "circular": [c.circular for c in aligned],
    })
    manifest_path = os.path.join(sample_out, "extraction_manifest.tsv")
    manifest.to_csv(manifest_path, sep="\t", index=False)

    plasmid_contigs = [c for c in aligned if c.label == "plasmid"]
    n_plasmid = len(plasmid_contigs)
    n_chr = sum(1 for c in aligned if c.label == "chromosome")
    n_unc = sum(1 for c in aligned if c.label == "unclassified")

    result = ExtractionResult(
        status="scored", input_fasta=os.path.abspath(input_fasta),
        sample_id=sample_id, n_contigs=len(aligned),
        total_length_bp=total_length, longest_contig_bp=longest,
        n_plasmid=n_plasmid, n_chromosome=n_chr,
        n_unclassified=n_unc, classifier="+".join(classifier_names),
        threshold=threshold)
    result.output_files.append(manifest_path)

    if n_plasmid == 0:
        result.status = "unscorable"
        result.route = ""
        result.reason = _unscorable_reason(aligned, n_plasmid)
        upath = os.path.join(sample_out, "plasrisk_unscorable.tsv")
        with open(upath, "w", encoding="utf-8") as fh:
            fh.write("field\tvalue\n")
            fh.write("status\tunscorable\n")
            fh.write("reason\t%s\n" % result.reason)
            fh.write("total_length_bp\t%d\n" % total_length)
            fh.write("longest_contig_bp\t%d\n" % longest)
        result.output_files.append(upath)
        if write_json:
            jpath = os.path.join(sample_out, "extraction_result.json")
            with open(jpath, "w", encoding="utf-8") as fh:
                json.dump(asdict(result), fh, indent=2)
            result.output_files.append(jpath)
        return result

    plasmid_fasta = os.path.join(sample_out, "predicted_plasmids.fasta")
    selected = {c.seq_id: full_records.get(c.seq_id, "")
                for c in plasmid_contigs}
    _write_fasta(selected, plasmid_fasta)
    result.plasmid_fasta = plasmid_fasta
    result.output_files.append(plasmid_fasta)

    chosen_route = (decide_route(plasmid_contigs, complete_length)
                    if route == "auto" else route)
    result.route = chosen_route

    lookup = None
    if use_abricate:
        try:
            lookup = load_replicon_lookup()
        except Exception:
            lookup = None

    if chosen_route == "isolate":
        result.output_files.extend(_score_isolate(
            plasmid_fasta, sample_out, mode, aware_version, lookup,
            min_id, min_cov, threads, use_abricate))
    else:
        result.output_files.extend(_score_meta(
            plasmid_fasta, sample_out, bins_table, coverage_table,
            aware_version, lookup, min_id, min_cov, threads,
            use_abricate, complete_length))

    if write_json:
        jpath = os.path.join(sample_out, "extraction_result.json")
        with open(jpath, "w", encoding="utf-8") as fh:
            json.dump(asdict(result), fh, indent=2)
        result.output_files.append(jpath)
    return result


def _circular_header(header: str) -> bool:
    from .metagenome import is_circular_header
    try:
        return is_circular_header(header)
    except Exception:
        return False


# ---------------------------------------------------------------------------
# CLI:  plasrisk run ...
# ---------------------------------------------------------------------------

def _collect_fastas(inputs: Sequence[str]) -> List[str]:
    fastas: List[str] = []
    exts = (".fasta", ".fa", ".fna", ".fsa", ".fas")
    for item in inputs:
        if os.path.isdir(item):
            for root, _dirs, files in os.walk(item):
                for f in sorted(files):
                    if f.lower().endswith(exts):
                        fastas.append(os.path.join(root, f))
        elif any(ch in item for ch in "*?"):
            fastas.extend(sorted(_glob.glob(item)))
        elif os.path.isfile(item):
            fastas.append(item)
    seen, unique = set(), []
    for f in fastas:
        ap = os.path.abspath(f)
        if ap not in seen:
            seen.add(ap)
            unique.append(f)
    return unique


def build_run_parser():
    import argparse
    p = argparse.ArgumentParser(
        prog="plasrisk run",
        description="Predict plasmids in raw WGS/long-read/metagenomic "
                    "assemblies, extract and score them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+", help="assembly FASTA(s) or directory")
    p.add_argument("-o", "--output", default="plasrisk_run")
    p.add_argument("--classifier", choices=["plasflow", "genomad", "both"],
                   default="plasflow")
    p.add_argument("--threshold", type=float, default=0.7)
    p.add_argument("--plasflow-table", default=None)
    p.add_argument("--genomad-table", default=None)
    p.add_argument("--generic-table", default=None)
    p.add_argument("--genomad-db", default=None)
    p.add_argument("--bins", default=None)
    p.add_argument("--coverage", default=None)
    p.add_argument("--mode", choices=["lite", "full", "ordinal"],
                   default="lite")
    p.add_argument("--aware-version", choices=["legacy", "2025"],
                   default="legacy")
    p.add_argument("--route", choices=["auto", "isolate", "meta"],
                   default="auto")
    p.add_argument("--vote", choices=["majority", "any"], default="majority")
    p.add_argument("--min-id", type=float, default=75.0)
    p.add_argument("--min-cov", type=float, default=50.0)
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--no-abricate", action="store_true")
    p.add_argument("--complete-length", type=int,
                   default=DEFAULT_COMPLETE_LENGTH_BP)
    p.add_argument("--sample-id", default=None)
    p.add_argument("--json", action="store_true")
    return p


def run_main(argv: Optional[List[str]] = None) -> int:
    args = build_run_parser().parse_args(argv)
    fastas = _collect_fastas(args.inputs)
    if not fastas:
        print("[plasrisk run] no FASTA inputs found")
        return 2
    os.makedirs(args.output, exist_ok=True)
    status_rows = []
    for fasta in fastas:
        try:
            res = process_one(
                fasta, args.output, classifier=args.classifier,
                threshold=args.threshold,
                plasflow_table=args.plasflow_table,
                genomad_table=args.genomad_table,
                generic_table=args.generic_table,
                genomad_db=args.genomad_db, bins_table=args.bins,
                coverage_table=args.coverage, mode=args.mode,
                aware_version=args.aware_version, route=args.route,
                vote=args.vote, min_id=args.min_id, min_cov=args.min_cov,
                threads=args.threads, use_abricate=not args.no_abricate,
                complete_length=args.complete_length,
                sample_id=(args.sample_id if len(fastas) == 1 else None),
                write_json=args.json)
            status_rows.append({
                "sample_id": res.sample_id, "status": res.status,
                "route": res.route, "n_plasmid": res.n_plasmid,
                "n_chromosome": res.n_chromosome,
                "total_length_bp": res.total_length_bp,
                "reason": res.reason})
        except Exception as exc:  # keep processing the remaining samples
            status_rows.append({
                "sample_id": os.path.basename(fasta), "status": "error",
                "route": "", "n_plasmid": "", "n_chromosome": "",
                "total_length_bp": "", "reason": str(exc)})
    status_df = pd.DataFrame(status_rows)
    status_path = os.path.join(args.output, "run_status.tsv")
    status_df.to_csv(status_path, sep="\t", index=False)
    print("[plasrisk run] wrote %s" % status_path)
    print(status_df.to_string(index=False))
    return 0
