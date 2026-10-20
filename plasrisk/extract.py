"""
extract.py - Predict plasmid sequences from raw genomic data before scoring.

Most real-world FASTA files are not ready-made plasmid collections. They are
short-read whole-genome assemblies, long-read assemblies, or multi-strain
metagenomic assemblies in which the plasmid sequences have to be predicted
first. This module adds that front end:

    raw assembly FASTA
        -> length filter (>=1 kb)
        -> plasmid/chromosome classifier (PlasFlow or geNomad; tables can
           also be supplied after running the classifier separately)
        -> extract predicted plasmid contigs
        -> route to the appropriate scorer:
              long/circular plasmid contigs  -> isolate scorer
              short fragmented contigs/bins  -> metagenomic three-layer
        -> if no plasmid is found, the sample is NOT forced into a score.
           It is reported as unscorable with a reason, including a
           chromosome-only determination based on assembly length and the
           classifier labels.

The external classifiers are invoked through subprocess when available; they
are not bundled. PlasFlow depends on an older TensorFlow/Python stack, so it
is commonly run in its own conda environment. Its tabular output can be
passed with --plasflow-table. The frozen PlasRisk weights are never refit.
"""

from __future__ import annotations

import glob as _glob
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from shutil import which
from typing import Dict, List, Optional

import pandas as pd

from .annotate import annotate_fasta, parse_fasta, run_abricate
from .metagenome import (
    aggregate_bins,
    aggregate_sample,
    circular_map_from_fasta,
    is_circular_header,
    load_bins_table,
    load_coverage_table,
    load_genomad,
    score_contigs,
)
from .scoring import LITE_DIMENSIONS, PlasRiskScorer

# A bacterial chromosome is typically 0.5-10 Mb. An assembly whose total
# length (or longest contig) is at least this value and which contains no
# predicted plasmid is treated as chromosome-like.
MIN_CHROMOSOME_TOTAL_BP = 500_000
LONG_CONTIG_BP = 500_000
# A circular molecule no longer than this with no chromosome-scale length is
# treated as a plasmid in the replicon-based long-read classifier.
PLASMID_CIRC_MAX_BP = 500_000
# Median contig length above which a chromosome-containing assembly is treated
# as long-read (few, very long contigs) rather than short-read (many fragments).
LONGREAD_MEDIAN_BP = 20_000
LONGREAD_MAX_CONTIGS = 300

_FASTA_EXTS = (".fasta", ".fa", ".fna", ".fsa", ".fas")
_READ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ClassifiedContig:
    seq_id: str
    length_bp: int
    label: str                      # plasmid / chromosome / unclassified
    plasmid_prob: Optional[float]
    circular: bool = False
    source: str = ""


@dataclass
class ExtractionResult:
    status: str                     # scored / unscorable
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
    route: str = ""
    reason: str = ""
    plasmid_fasta: Optional[str] = None
    output_files: Dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _bp(value) -> int:
    m = re.search(r"\d+", str(value))
    return int(m.group()) if m else 0


def _prob(value) -> Optional[float]:
    try:
        f = float(str(value).strip())
    except (ValueError, TypeError):
        return None
    if f != f:  # NaN
        return None
    return f


def _col(columns: List[str], *wanted: str) -> Optional[str]:
    lower = {str(c).lower().strip(): c for c in columns}
    for w in wanted:
        if w in lower:
            return lower[w]
    return None


def _col_contains(columns: List[str], sub: str,
                  exclude: Optional[List[str]] = None) -> Optional[str]:
    exclude = [str(c).lower() for c in (exclude or [])]
    for c in columns:
        cl = str(c).lower()
        if sub in cl and cl not in exclude:
            return c
    return None


def _read_table(path: str) -> pd.DataFrame:
    # utf-8-sig strips a leading BOM (common in files exported on Windows)
    # and otherwise behaves as plain utf-8.
    df = pd.read_csv(path, sep=None, engine="python", dtype=str,
                     encoding="utf-8-sig", comment="#",
                     keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _normalize_label(text: str,
                     plasmid_prob: Optional[float] = None,
                     threshold: float = 0.5) -> str:
    t = str(text).lower().strip()
    if "plasmid" in t:
        return "plasmid"
    if "chromos" in t:
        return "chromosome"
    if "unclass" in t or t in ("", "nan", "none", "-"):
        if plasmid_prob is not None and plasmid_prob >= threshold:
            return "plasmid"
        return "unclassified"
    if t in ("0", "1") and plasmid_prob is not None:
        return "plasmid" if plasmid_prob >= threshold else "chromosome"
    return "unclassified"


def _wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def _write_fasta(path: str, records: Dict[str, str]) -> None:
    with open(path, "w", newline="\n") as fh:
        for sid, seq in records.items():
            fh.write(">%s\n%s\n" % (sid, _wrap(seq)))


def _filter_fasta(inp: str, outp: str, min_length: int) -> int:
    recs = parse_fasta(inp)
    kept = {sid: seq for sid, seq in recs.items() if len(seq) >= min_length}
    _write_fasta(outp, kept)
    return len(kept)


def _read_full_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA keyed by the FULL raw header (not truncated at spaces).

    Classifier tables (PlasFlow contig_name, SPAdes headers) carry the full
    header, while parse_fasta keeps only the first whitespace token, so the
    full header is needed to match them."""
    full: Dict[str, str] = {}
    sid: Optional[str] = None
    buf: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if sid is not None:
                    full[sid] = "".join(buf)
                sid = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if sid is not None:
        full[sid] = "".join(buf)
    return full


def _build_id_resolver(full: Dict[str, str]):
    """Map a classifier seq id to the stored full header, tolerating the
    first-token form when headers contain spaces."""
    by_token: Dict[str, str] = {}
    for header in full:
        tokens = header.split()
        token = tokens[0] if tokens else header
        by_token.setdefault(token, header)

    def resolve(candidate: str) -> Optional[str]:
        if candidate in full:
            return candidate
        tokens = candidate.split()
        token = tokens[0] if tokens else candidate
        return by_token.get(token)

    return resolve


# ---------------------------------------------------------------------------
# Classifier output parsers
# ---------------------------------------------------------------------------

def parse_plasflow_table(path: str, threshold: float = 0.7
                         ) -> List[ClassifiedContig]:
    """Parse a PlasFlow prediction TSV.

    Header layout: contig_id, contig_name, contig_length, id, label, then
    one probability column per class (unclassified, chromosome, plasmid)."""
    df = _read_table(path)
    cols = list(df.columns)
    name_c = _col(cols, "contig_name", "contig_id", "seq_name", "name")
    len_c = _col(cols, "contig_length", "length")
    label_c = _col(cols, "label", "prediction", "class")
    prob_c = _col_contains(cols, "plasmid",
                           exclude=[str(name_c), str(label_c), str(len_c)])

    out: List[ClassifiedContig] = []
    for _, row in df.iterrows():
        sid = str(row[name_c]).strip() if name_c else ""
        length = _bp(row[len_c]) if len_c else 0
        prob = _prob(row[prob_c]) if prob_c else None
        label = _normalize_label(row[label_c] if label_c else "", prob,
                                 threshold)
        # When the label is missing but probabilities are available, derive it.
        if label == "unclassified" and prob is not None and prob >= threshold:
            label = "plasmid"
        out.append(ClassifiedContig(sid, length, label, prob,
                                    is_circular_header(sid), "plasflow"))
    return out


def parse_genomad_summary(path: str, threshold: float = 0.7
                          ) -> List[ClassifiedContig]:
    """Parse a geNomad plasmid/summary TSV (a directory is also accepted)."""
    target = path
    if os.path.isdir(path):
        # geNomad >=1.7 end-to-end writes nested <stem>_summary/
        # <stem>_plasmid_summary.tsv, so search recursively first.
        hits = sorted(_glob.glob(os.path.join(path, "**",
                                              "*_plasmid_summary.tsv"),
                                 recursive=True))
        if not hits:
            hits = sorted(_glob.glob(os.path.join(path,
                                                  "*_plasmid_summary.tsv")))
        if not hits:
            hits = sorted(_glob.glob(os.path.join(path, "**", "*_summary.tsv"),
                                     recursive=True))
        if not hits:
            hits = sorted(_glob.glob(os.path.join(path, "*_summary.tsv")))
        if not hits:
            raise FileNotFoundError(
                "no geNomad *_plasmid_summary.tsv found under %s" % path)
        target = hits[0]

    df = _read_table(target)
    cols = list(df.columns)
    name_c = _col(cols, "seq_name", "contig_name", "seq_id", "name")
    len_c = _col(cols, "length", "contig_length")
    p_c = _col(cols, "plasmid_score", "prob_plasmid")
    c_c = _col(cols, "chromosome_score")
    topo_c = _col(cols, "topology")

    out: List[ClassifiedContig] = []
    for _, row in df.iterrows():
        sid = str(row[name_c]).strip() if name_c else ""
        length = _bp(row[len_c]) if len_c else 0
        pscore = _prob(row[p_c]) if p_c else None
        cscore = _prob(row[c_c]) if c_c else None
        circular = False
        if topo_c and str(row[topo_c]).lower().strip() in ("circular", "yes"):
            circular = True
        circular = circular or is_circular_header(sid)

        if pscore is not None:
            label = "plasmid" if pscore >= threshold else "chromosome"
        else:
            label = "unclassified"
        if label == "chromosome" and cscore is not None and pscore is not None \
                and pscore > cscore:
            label = "plasmid"
        out.append(ClassifiedContig(sid, length, label, pscore, circular,
                                    "genomad"))
    return out


def parse_generic_table(path: str, threshold: float = 0.5
                        ) -> List[ClassifiedContig]:
    """Parse any classifier table with a sequence id and a label column."""
    df = _read_table(path)
    cols = list(df.columns)
    name_c = _col(cols, "seq_id", "seq_name", "contig_name", "contig_id",
                  "name", "id")
    label_c = _col(cols, "label", "prediction", "predicted", "class")
    prob_c = _col(cols, "plasmid_score", "prob_plasmid", "plasmid_prob",
                  "plasmid")
    len_c = _col(cols, "length", "contig_length", "seq_length")
    out: List[ClassifiedContig] = []
    for _, row in df.iterrows():
        sid = str(row[name_c]).strip() if name_c else ""
        length = _bp(row[len_c]) if len_c else 0
        prob = _prob(row[prob_c]) if prob_c else None
        label = _normalize_label(row[label_c] if label_c else "", prob,
                                 threshold)
        out.append(ClassifiedContig(sid, length, label, prob,
                                    is_circular_header(sid), "generic"))
    return out


# ---------------------------------------------------------------------------
# Replicon-based classifier (long-read assemblies)
# ---------------------------------------------------------------------------

def parse_abricate_replicon_table(path: str) -> set:
    """Return the set of contig ids carrying a PlasmidFinder replicon hit.

    Accepts a raw abricate TSV (columns include SEQUENCE and DATABASE) or any
    table whose header contains a SEQUENCE/seq_id column and a GENE/name
    column mentioning a replicon (Inc*/Col/rep)."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        df = _read_table(path)
    if df is None or df.empty:
        return set()
    cols = list(df.columns)
    seq_c = _col(cols, "sequence", "seq_id", "seq_name", "contig_name",
                 "contig_id", "name")
    gene_c = _col(cols, "gene", "replicon", "name", "accession")
    db_c = _col(cols, "database")
    hits: set = set()
    if seq_c is None:
        return hits
    for _, row in df.iterrows():
        sid = str(row[seq_c]).strip()
        text = " ".join(str(row.get(c, "")) for c in
                        ([gene_c] if gene_c else []) +
                        ([db_c] if db_c else []))
        # A raw PlasmidFinder abricate table already reports only replicons;
        # for a generic table require an explicit replicon-like token.
        is_replicon = (db_c is not None and
                       "plasmidfinder" in str(row.get(db_c, "")).lower())
        if not is_replicon and gene_c is not None:
            g = str(row[gene_c])
            is_replicon = bool(re.search(r"(inc[a-z]?[a-z0-9_-]+|col\w*|"
                                         r"rep[a-z0-9_]*)", g, re.I))
        if is_replicon:
            hits.add(sid)
    return hits


def classify_by_replicon(fasta: str,
                         replicon_table: Optional[str] = None,
                         min_length: int = 1000,
                         threads: int = 4,
                         min_id: float = 75.0,
                         min_cov: float = 50.0
                         ) -> List[ClassifiedContig]:
    """Classify long-read assembly contigs by plasmid replicon evidence.

    A contig is called a plasmid when (i) it carries a PlasmidFinder replicon
    / replication gene (the defining, sequence-based plasmid marker that
    survives in complete long-read assemblies), or (ii) it is explicitly
    circularised and within the plasmid size range. Contigs at chromosome
    scale without either signal are called chromosome; the rest are left
    unclassified.
    """
    records = _read_full_fasta(fasta)
    resolve_id = _build_id_resolver(records)

    if replicon_table:
        raw_hits = parse_abricate_replicon_table(replicon_table)
    else:
        if not which("abricate"):
            raise RuntimeError(
                "abricate (with the plasmidfinder database) is required for "
                "the replicon classifier. Install it or run abricate "
                "separately and pass --replicon-table; alternatively use "
                "--classifier genomad.")
        hits_df = run_abricate(fasta, "plasmidfinder", min_id=min_id,
                               min_cov=min_cov, threads=threads)
        raw_hits = set()
        if hits_df is not None and not hits_df.empty and "SEQUENCE" in \
                hits_df.columns:
            raw_hits = {str(s).strip() for s in hits_df["SEQUENCE"]}

    rep_hits: set = set()
    for raw in raw_hits:
        resolved = resolve_id(raw)
        if resolved is not None:
            rep_hits.add(resolved)
        else:
            rep_hits.add(raw.split()[0] if raw.split() else raw)

    out: List[ClassifiedContig] = []
    for header, seq in records.items():
        length = len(seq)
        circular = is_circular_header(header)
        has_replicon = header in rep_hits or \
            (header.split()[0] if header.split() else header) in rep_hits
        is_circ_plasmid = (circular and min_length <= length <=
                           PLASMID_CIRC_MAX_BP)
        if has_replicon or is_circ_plasmid:
            label = "plasmid"
        elif length >= LONG_CONTIG_BP:
            label = "chromosome"
        else:
            label = "unclassified"
        out.append(ClassifiedContig(header, length, label, None,
                                    circular, "replicon"))
    return out


# ---------------------------------------------------------------------------
# Classifier runners (optional, via subprocess)
# ---------------------------------------------------------------------------

def plasflow_available() -> bool:
    return bool(which("PlasFlow.py") or which("plasflow"))


def genomad_available() -> bool:
    return bool(which("genomad"))


def _run(cmd: List[str]) -> None:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip()[-2000:]
        raise RuntimeError("command failed (%s): %s" % (" ".join(cmd), tail))


def run_plasflow(fasta: str, workdir: str, threshold: float = 0.7,
                 min_length: int = 1000, batch_size: int = 25000) -> str:
    os.makedirs(workdir, exist_ok=True)
    filtered = os.path.join(workdir, "plasflow_input.fasta")
    _filter_fasta(fasta, filtered, min_length)
    out_tsv = os.path.join(workdir, "plasflow_predictions.tsv")
    exe = which("PlasFlow.py") or which("plasflow")
    if not exe:
        raise RuntimeError(
            "PlasFlow not found in PATH (run it separately and pass "
            "--plasflow-table)")
    base = [exe, "--input", filtered, "--output", out_tsv,
            "--threshold", str(threshold)]
    try:
        _run(base + ["--batch_size", str(batch_size)])
    except RuntimeError:
        _run(base)  # older releases without --batch_size
    if not os.path.exists(out_tsv):
        raise RuntimeError("PlasFlow finished but no output TSV found")
    return out_tsv


def run_genomad(fasta: str, workdir: str, db_path: str,
                threads: int = 4, min_length: int = 1000,
                threshold: float = 0.7) -> str:
    """Run geNomad end-to-end (geNomad >=1.7; verified against v1.12).

    The legacy `genomad classify` subcommand was removed; v1.7+ uses
    `genomad end-to-end INPUT OUTPUT DATABASE` and writes the per-sequence
    plasmid table to <OUTPUT>/<stem>_summary/<stem>_plasmid_summary.tsv.
    """
    db_path = db_path or os.environ.get("GENOMAD_DB", "")
    if not db_path:
        raise RuntimeError(
            "geNomad requires a database path (--genomad-db or the "
            "GENOMAD_DB environment variable); download it with "
            "'genomad download-database <dir>'")
    os.makedirs(workdir, exist_ok=True)
    filtered = os.path.join(workdir, "genomad_input.fasta")
    _filter_fasta(fasta, filtered, min_length)
    cmd = ["genomad", "end-to-end", "--threads", str(threads),
           "--cleanup", "--min-score", str(threshold),
           filtered, workdir, db_path]
    _run(cmd)
    return workdir


# ---------------------------------------------------------------------------
# Combine classifiers
# ---------------------------------------------------------------------------

def combine_classifications(groups: List[List[ClassifiedContig]],
                            vote: str = "majority",
                            threshold: float = 0.5
                            ) -> List[ClassifiedContig]:
    """Combine several classifier outputs by sequence id.

    vote='majority': plasmid/chromosome only when more than half agree.
    vote='any':      plasmid if any classifier calls it plasmid (high recall).
    """
    if len(groups) == 1:
        return groups[0]

    votes: Dict[str, List[ClassifiedContig]] = {}
    order: List[str] = []
    for grp in groups:
        for c in grp:
            if c.seq_id not in votes:
                votes[c.seq_id] = []
                order.append(c.seq_id)
            votes[c.seq_id].append(c)

    out: List[ClassifiedContig] = []
    for sid in order:
        cs = votes[sid]
        n = len(cs)
        n_pl = sum(1 for c in cs if c.label == "plasmid")
        n_ch = sum(1 for c in cs if c.label == "chromosome")
        probs = [c.plasmid_prob for c in cs if c.plasmid_prob is not None]
        prob = sum(probs) / len(probs) if probs else None
        length = max(c.length_bp for c in cs)
        circular = any(c.circular for c in cs)
        sources = "+".join(sorted({c.source for c in cs}))

        if vote == "any":
            label = "plasmid" if n_pl >= 1 else (
                "chromosome" if n_ch == n else "unclassified")
        else:
            if n_pl > n / 2:
                label = "plasmid"
            elif n_ch > n / 2:
                label = "chromosome"
            else:
                label = "unclassified"
        if label == "unclassified" and prob is not None \
                and prob >= threshold:
            label = "plasmid"
        out.append(ClassifiedContig(sid, length, label, prob, circular,
                                    sources))
    return out


# ---------------------------------------------------------------------------
# Routing and unscorable determination
# ---------------------------------------------------------------------------

def decide_route(plasmid_contigs: List[ClassifiedContig],
                 complete_length: int = 10_000) -> str:
    n = len(plasmid_contigs)
    if n == 0:
        return "meta"
    n_complete = sum(1 for c in plasmid_contigs
                     if c.circular or c.length_bp >= complete_length)
    any_circular = any(c.circular for c in plasmid_contigs)
    if (any_circular and n <= 25) or n_complete / n >= 0.5:
        return "isolate"
    return "meta"


def _unscorable_reason(res: ExtractionResult) -> str:
    chromosome_like = (res.total_length_bp >= MIN_CHROMOSOME_TOTAL_BP or
                       res.longest_contig_bp >= LONG_CONTIG_BP)
    if res.n_chromosome > 0 and chromosome_like:
        return ("No plasmid sequence detected. The input appears to contain "
                "chromosomal sequence only (total length {tot:,} bp across "
                "{n} contig(s); longest {long:,} bp; {nc} contig(s) labeled "
                "chromosome). PlasRisk scores plasmids rather than "
                "chromosomes, so no risk score is returned.").format(
                    tot=res.total_length_bp, n=res.n_contigs,
                    long=res.longest_contig_bp, nc=res.n_chromosome)
    if chromosome_like:
        return ("No plasmid sequence detected. The assembly is long enough "
                "to be chromosome-like (total {tot:,} bp; longest {long:,} "
                "bp) but the classifier did not assign confident labels. "
                "Consider lowering --threshold, using another classifier, "
                "or checking that the assembly is not a single "
                "chromosome.").format(tot=res.total_length_bp,
                                      long=res.longest_contig_bp)
    return ("No plasmid sequence detected, and the assembly is too small or "
            "fragmented to classify (total {tot:,} bp across {n} contig(s); "
            "longest {long:,} bp). PlasRisk cannot produce a reliable score "
            "from this input.").format(tot=res.total_length_bp,
                                       n=res.n_contigs,
                                       long=res.longest_contig_bp)


# ---------------------------------------------------------------------------
# Scoring routes
# ---------------------------------------------------------------------------

def _make_scorer(mode: str, aware_version: str,
                 lookup) -> PlasRiskScorer:
    return PlasRiskScorer(replicon_lookup=lookup, mode=mode,
                          aware_version=aware_version)


def _score_isolate(plasmid_fasta: str, outdir: str, scorer: PlasRiskScorer,
                   use_abricate: bool, min_id: float, min_cov: float,
                   threads: int, lookup) -> Dict[str, str]:
    ann = annotate_fasta(plasmid_fasta, use_abricate=use_abricate,
                         min_id=min_id, min_cov=min_cov, threads=threads,
                         lookup=lookup)
    if not ann.features:
        raise RuntimeError("no features in extracted plasmid FASTA")
    df = scorer.score_dataframe(ann.features)
    tsv = os.path.join(outdir, "plasrisk_isolate_scores.tsv")
    df.to_csv(tsv, sep="\t", index=False, float_format="%.4f")
    return {"isolate_scores": tsv}


def _score_meta(plasmid_fasta: str, outdir: str, scorer: PlasRiskScorer,
                sample_id: str, bins_df, coverage_df,
                use_abricate: bool, min_id: float, min_cov: float,
                threads: int, min_contig: int, complete_length: int,
                lookup) -> Dict[str, str]:
    ann = annotate_fasta(plasmid_fasta, use_abricate=use_abricate,
                         min_id=min_id, min_cov=min_cov, threads=threads,
                         lookup=lookup)
    circ_header = circular_map_from_fasta(plasmid_fasta)
    feats_kept = {f.seq_id: f for f in ann.features
                  if f.length_bp >= min_contig}
    contig_df = score_contigs(ann, scorer,
                              complete_length=complete_length,
                              min_length=min_contig,
                              circular_map=circ_header)
    contig_df.insert(0, "sample_id", sample_id)
    bin_df = aggregate_bins(contig_df, feats_kept, bins_df, coverage_df,
                            scorer, complete_length=complete_length,
                            header_circular_map=circ_header)
    bin_df.insert(0, "sample_id", sample_id)
    total_bp = int(sum(f.length_bp for f in feats_kept.values()))
    sample_row = aggregate_sample(sample_id, bin_df, contig_df, total_bp)

    cpath = os.path.join(outdir, "meta_contig_scores.tsv")
    bpath = os.path.join(outdir, "meta_bin_scores.tsv")
    spath = os.path.join(outdir, "meta_sample_summary.tsv")
    contig_df.to_csv(cpath, sep="\t", index=False)
    bin_df.to_csv(bpath, sep="\t", index=False)
    pd.DataFrame([sample_row]).to_csv(spath, sep="\t", index=False)
    return {"contig_scores": cpath, "bin_scores": bpath,
            "sample_summary": spath}


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def process_one(fasta: str, outdir: str, *,
                classifier: str = "plasflow",
                plasflow_table: Optional[str] = None,
                genomad_table: Optional[str] = None,
                replicon_table: Optional[str] = None,
                generic_table: Optional[str] = None,
                bins_table: Optional[str] = None,
                coverage_table: Optional[str] = None,
                coverage_column: Optional[str] = None,
                threshold: float = 0.7,
                genomad_db: Optional[str] = None,
                route: str = "auto",
                mode: str = "lite",
                aware_version: str = "legacy",
                use_abricate: bool = True,
                threads: int = 4,
                min_id: float = 75.0,
                min_cov: float = 50.0,
                min_contig: int = 1000,
                complete_length: int = 10_000,
                vote: str = "majority",
                write_json: bool = False,
                quiet: bool = False) -> ExtractionResult:
    sample_id = os.path.splitext(os.path.basename(fasta))[0]
    os.makedirs(outdir, exist_ok=True)

    records = _read_full_fasta(fasta)
    lengths = {sid: len(seq) for sid, seq in records.items()}
    resolve_id = _build_id_resolver(records)
    res = ExtractionResult(
        status="unscorable", input_fasta=fasta, sample_id=sample_id,
        n_contigs=len(records),
        total_length_bp=int(sum(lengths.values())),
        longest_contig_bp=int(max(lengths.values()) if lengths else 0),
        classifier=classifier, threshold=threshold)

    # Gather classifier outputs (pre-existing tables take precedence).
    groups: List[List[ClassifiedContig]] = []
    classifiers_used: List[str] = []

    if plasflow_table:
        groups.append(parse_plasflow_table(plasflow_table, threshold))
        classifiers_used.append("plasflow")
    if genomad_table:
        groups.append(parse_genomad_summary(genomad_table, threshold))
        classifiers_used.append("genomad")
    if replicon_table:
        groups.append(classify_by_replicon(
            fasta, replicon_table=replicon_table, min_length=min_contig,
            threads=threads, min_id=min_id, min_cov=min_cov))
        classifiers_used.append("replicon")
    if generic_table:
        groups.append(parse_generic_table(generic_table, threshold))
        classifiers_used.append("generic")

    if not groups:
        if classifier == "plasflow":
            if not plasflow_available():
                raise RuntimeError(
                    "PlasFlow not found in PATH. Run PlasFlow separately "
                    "and pass its TSV with --plasflow-table.")
            tsv = run_plasflow(fasta, outdir, threshold, min_contig)
            groups.append(parse_plasflow_table(tsv, threshold))
            classifiers_used.append("plasflow")
        elif classifier == "genomad":
            if not genomad_available():
                raise RuntimeError(
                    "geNomad not found in PATH. Run it separately and pass "
                    "the summary/directory with --genomad-table.")
            gdir = run_genomad(fasta, outdir, genomad_db, threads,
                               min_contig, threshold)
            groups.append(parse_genomad_summary(gdir, threshold))
            classifiers_used.append("genomad")
        elif classifier == "replicon":
            groups.append(classify_by_replicon(
                fasta, replicon_table=None, min_length=min_contig,
                threads=threads, min_id=min_id, min_cov=min_cov))
            classifiers_used.append("replicon")
        else:
            raise RuntimeError("unknown classifier %r" % classifier)

    classified = combine_classifications(groups, vote=vote,
                                         threshold=threshold)
    # Align classifier ids to the input full headers; fill lengths and
    # circularity from the sequences when the classifier omitted them.
    for c in classified:
        resolved = resolve_id(c.seq_id)
        if resolved is not None:
            c.seq_id = resolved
        if not c.length_bp and c.seq_id in lengths:
            c.length_bp = lengths[c.seq_id]
        if c.seq_id in lengths and not c.circular:
            c.circular = is_circular_header(c.seq_id)

    plasmid = [c for c in classified if c.label == "plasmid"]
    res.n_plasmid = len(plasmid)
    res.n_chromosome = sum(1 for c in classified
                           if c.label == "chromosome")
    res.n_unclassified = sum(1 for c in classified
                             if c.label == "unclassified")
    res.classifier = "+".join(classifiers_used)

    # Manifest of every classified contig.
    manifest = os.path.join(outdir, "extraction_manifest.tsv")
    with open(manifest, "w", newline="\n") as fh:
        fh.write("seq_id\tlength_bp\tlabel\tplasmid_prob\tcircular\n")
        for c in classified:
            prob = "" if c.plasmid_prob is None else "%.4f" % c.plasmid_prob
            fh.write("%s\t%d\t%s\t%s\t%s\n" % (
                c.seq_id, c.length_bp, c.label, prob,
                "true" if c.circular else "false"))
    res.output_files["manifest"] = manifest

    if not plasmid:
        res.status = "unscorable"
        res.route = ""
        res.reason = _unscorable_reason(res)
        unsc = os.path.join(outdir, "plasrisk_unscorable.tsv")
        with open(unsc, "w", newline="\n") as fh:
            fh.write("field\tvalue\n")
            fh.write("status\tunscorable\n")
            fh.write("sample_id\t%s\n" % sample_id)
            fh.write("n_contigs\t%d\n" % res.n_contigs)
            fh.write("total_length_bp\t%d\n" % res.total_length_bp)
            fh.write("longest_contig_bp\t%d\n" % res.longest_contig_bp)
            fh.write("n_chromosome\t%d\n" % res.n_chromosome)
            fh.write("n_unclassified\t%d\n" % res.n_unclassified)
            fh.write("reason\t%s\n" % res.reason)
        res.output_files["unscorable"] = unsc
        if write_json:
            jpath = os.path.join(outdir, "plasrisk_unscorable.json")
            with open(jpath, "w") as jf:
                json.dump({k: v for k, v in res.__dict__.items()
                           if k != "output_files"}, jf, indent=2,
                          default=str)
            res.output_files["json"] = jpath
        return res

    # Extract predicted plasmid sequences.
    plasmid_fasta = os.path.join(outdir, "predicted_plasmids.fasta")
    kept = {c.seq_id: records[c.seq_id] for c in plasmid
            if c.seq_id in records}
    _write_fasta(plasmid_fasta, kept)
    res.plasmid_fasta = plasmid_fasta
    res.output_files["plasmid_fasta"] = plasmid_fasta

    chosen_route = route if route in ("isolate", "meta") else \
        decide_route(plasmid, complete_length)
    res.route = chosen_route

    from .lookup import load_replicon_lookup
    lookup = None
    if use_abricate:
        try:
            lookup = load_replicon_lookup()
        except Exception:
            lookup = None
    scorer = _make_scorer(mode if mode in ("lite", "full") else "lite",
                          aware_version, lookup)

    bins_df = load_bins_table(bins_table) if bins_table else None
    coverage_df = load_coverage_table(coverage_table, coverage_column) \
        if coverage_table else None

    if chosen_route == "isolate":
        files = _score_isolate(plasmid_fasta, outdir, scorer, use_abricate,
                               min_id, min_cov, threads, lookup)
    else:
        files = _score_meta(plasmid_fasta, outdir, scorer, sample_id,
                            bins_df, coverage_df, use_abricate, min_id,
                            min_cov, threads, min_contig, complete_length,
                            lookup)
    res.output_files.update(files)
    res.status = "scored"
    res.reason = ""

    if write_json:
        jpath = os.path.join(outdir, "extraction_result.json")
        with open(jpath, "w") as jf:
            json.dump({k: v for k, v in res.__dict__.items()
                       if k != "output_files"}, jf, indent=2, default=str)
        res.output_files["json"] = jpath
    return res


def _collect_fastas(inputs: List[str]) -> List[str]:
    found: List[str] = []
    for inp in inputs:
        matches = _glob.glob(inp)
        if not matches and os.path.exists(inp):
            matches = [inp]
        for m in matches:
            if os.path.isdir(m):
                for root, _d, files in os.walk(m):
                    for f in files:
                        if os.path.splitext(f)[1].lower() in _FASTA_EXTS:
                            found.append(os.path.join(root, f))
            elif os.path.isfile(m):
                found.append(m)
    seen, uniq = set(), []
    for f in found:
        ap = os.path.abspath(f)
        if ap not in seen:
            seen.add(ap)
            uniq.append(f)
    return uniq


def build_run_parser():
    import argparse
    p = argparse.ArgumentParser(
        prog="plasrisk run",
        description="Predict plasmids in a raw WGS/long-read/metagenomic "
                    "assembly and then score them.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+", help="assembly FASTA file(s) or a "
                                             "directory")
    p.add_argument("-o", "--output", default="plasrisk_run")
    p.add_argument("--classifier",
                   choices=["plasflow", "genomad", "replicon"],
                   default="plasflow",
                   help="plasflow/genomad for short-read assemblies; "
                        "'replicon' uses PlasmidFinder replicon evidence and "
                        "is designed for complete long-read assemblies")
    p.add_argument("--plasflow-table", default=None,
                   help="precomputed PlasFlow prediction TSV")
    p.add_argument("--genomad-table", default=None,
                   help="precomputed geNomad summary TSV or directory")
    p.add_argument("--replicon-table", default=None,
                   help="precomputed abricate plasmidfinder TSV for the "
                        "replicon classifier")
    p.add_argument("--generic-table", default=None,
                   help="generic classifier table (id + label columns)")
    p.add_argument("--bins", default=None,
                   help="contig-to-bin TSV (seq_id,bin_id[,circular,"
                        "plasmid_score]) for the meta route")
    p.add_argument("--coverage", default=None,
                   help="per-contig coverage TSV (seq_id,coverage)")
    p.add_argument("--coverage-column", default=None)
    p.add_argument("--threshold", type=float, default=0.7,
                   help="classifier probability threshold")
    p.add_argument("--genomad-db", default=None,
                   help="geNomad database path")
    p.add_argument("--route", choices=["auto", "isolate", "meta"],
                   default="auto")
    p.add_argument("--mode", choices=["lite", "full"], default="lite")
    p.add_argument("--aware-version", choices=["legacy", "2025"],
                   default="legacy")
    p.add_argument("--vote", choices=["majority", "any"], default="majority",
                   help="how to combine multiple classifiers")
    p.add_argument("--no-abricate", action="store_true",
                   help="skip abricate; sequence/header evidence only")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--min-id", type=float, default=75.0)
    p.add_argument("--min-cov", type=float, default=50.0)
    p.add_argument("--min-contig", type=int, default=1000)
    p.add_argument("--complete-length", type=int, default=10000)
    p.add_argument("--json", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    return p


def run_main(argv: Optional[List[str]] = None) -> int:
    args = build_run_parser().parse_args(argv)
    fastas = _collect_fastas(args.inputs)
    if not fastas:
        print("[plasrisk run] no FASTA inputs found", file=sys.stderr)
        return 2
    os.makedirs(args.output, exist_ok=True)

    status_rows = []
    for i, fasta in enumerate(fastas, 1):
        sample_id = os.path.splitext(os.path.basename(fasta))[0]
        sub = os.path.join(args.output, sample_id)
        if not args.quiet:
            print("[plasrisk run] (%d/%d) %s" % (i, len(fastas), sample_id))
        try:
            res = process_one(
                fasta, sub, classifier=args.classifier,
                plasflow_table=args.plasflow_table,
                genomad_table=args.genomad_table,
                replicon_table=args.replicon_table,
                generic_table=args.generic_table,
                bins_table=args.bins, coverage_table=args.coverage,
                coverage_column=args.coverage_column,
                threshold=args.threshold, genomad_db=args.genomad_db,
                route=args.route, mode=args.mode,
                aware_version=args.aware_version,
                use_abricate=not args.no_abricate, threads=args.threads,
                min_id=args.min_id, min_cov=args.min_cov,
                min_contig=args.min_contig,
                complete_length=args.complete_length, vote=args.vote,
                write_json=args.json, quiet=args.quiet)
        except Exception as exc:
            print("    ERROR: %s" % exc, file=sys.stderr)
            status_rows.append({"sample_id": sample_id, "status": "error",
                                "route": "", "reason": str(exc)})
            continue

        if not args.quiet:
            if res.status == "scored":
                print("    %s | route=%s | plasmids=%d" % (
                    res.status.upper(), res.route, res.n_plasmid))
            else:
                print("    UNSCORABLE | %s" % res.reason)
        status_rows.append({
            "sample_id": sample_id, "status": res.status,
            "route": res.route, "n_plasmid": res.n_plasmid,
            "n_chromosome": res.n_chromosome,
            "total_length_bp": res.total_length_bp,
            "reason": res.reason})

    status_df = pd.DataFrame(status_rows)
    spath = os.path.join(args.output, "run_status.tsv")
    status_df.to_csv(spath, sep="\t", index=False)
    if not args.quiet:
        print("[plasrisk run] status written to %s" % spath)
    return 0


# ---------------------------------------------------------------------------
# Unified auto entry point
# ---------------------------------------------------------------------------

def characterize_assembly(fasta: str) -> Dict[str, object]:
    """Summarise an assembly so the auto router can decide how to handle it.

    Pure plasmid collections contain no chromosome-scale molecule. Raw WGS /
    long-read / metagenomic assemblies do; they are then split by contig size
    and number into long-read (few, very long contigs) vs short-read (many
    fragmented contigs)."""
    records = _read_full_fasta(fasta)
    lengths = sorted((len(s) for s in records.values()), reverse=True)
    n = len(lengths)
    total = int(sum(lengths))
    longest = lengths[0] if lengths else 0
    median = lengths[n // 2] if n else 0
    n_circular = sum(1 for h in records if is_circular_header(h))
    chromosome_like = (total >= MIN_CHROMOSOME_TOTAL_BP or
                       longest >= LONG_CONTIG_BP)
    looks_long = ((longest >= LONG_CONTIG_BP and n <= LONGREAD_MAX_CONTIGS)
                  or median >= LONGREAD_MEDIAN_BP)
    return {
        "n_seqs": n,
        "total_bp": total,
        "longest_bp": longest,
        "median_bp": median,
        "n_circular": n_circular,
        "chromosome_like": chromosome_like,
        "read_type": "long" if looks_long else "short",
    }


def choose_auto_strategy(char: Dict[str, object],
                         input_type: str = "auto",
                         read_type: str = "auto") -> Dict[str, str]:
    """Map an assembly characterisation to an action/classifier/route.

    action='direct' scores the FASTA as a ready plasmid collection;
    action='extract' predicts and extracts plasmids first."""
    if input_type == "plasmid":
        action = "direct"
    elif input_type in ("genome", "metagenome"):
        action = "extract"
    else:  # auto: a chromosome-scale molecule means non-plasmid input
        action = "extract" if char["chromosome_like"] else "direct"

    rt = read_type if read_type in ("short", "long") else \
        str(char["read_type"])

    if action == "direct":
        return {"action": "direct", "classifier": "", "route": "direct",
                "read_type": rt,
                "rationale": "no chromosome-scale molecule; treated as a "
                             "plasmid collection and scored directly"}

    if input_type == "metagenome":
        # Both short- and long-read metagenomic assemblies are fragmented
        # multi-strain mixtures; geNomad + three-layer aggregation applies.
        return {"action": "extract", "classifier": "genomad",
                "route": "meta", "read_type": rt,
                "rationale": "metagenomic assembly; geNomad plasmid "
                             "prediction followed by three-layer scoring"}

    if rt == "long":
        return {"action": "extract", "classifier": "replicon",
                "route": "isolate", "read_type": rt,
                "rationale": "long-read assembly with chromosome-scale "
                             "contig(s); plasmids identified by PlasmidFinder "
                             "replicon/circularity and extracted"}

    return {"action": "extract", "classifier": "genomad", "route": "auto",
            "read_type": "short",
            "rationale": "short-read assembly; geNomad predicts plasmid "
                         "contigs, then isolate/meta routing is decided from "
                         "contig completeness"}


def _score_direct(fasta: str, outdir: str, *, mode: str,
                  aware_version: str, use_abricate: bool, min_id: float,
                  min_cov: float, threads: int, lookup) -> Dict[str, str]:
    """Score a ready-made plasmid FASTA exactly like the default command."""
    os.makedirs(outdir, exist_ok=True)
    ann = annotate_fasta(fasta, use_abricate=use_abricate, min_id=min_id,
                         min_cov=min_cov, threads=threads, lookup=lookup)
    if not ann.features:
        raise RuntimeError(
            "no sequences/features found; check that the FASTA is valid and "
            "contigs meet the length filter")
    scorer = _make_scorer(mode if mode in ("lite", "full") else "lite",
                          aware_version, lookup)
    df = scorer.score_dataframe(ann.features)
    tsv = os.path.join(outdir, "plasrisk_results.tsv")
    df.to_csv(tsv, sep="\t", index=False, float_format="%.4f")
    return {"scores": tsv}


def build_auto_parser():
    import argparse
    p = argparse.ArgumentParser(
        prog="plasrisk auto",
        description="Unified entry point: auto-detect plasmid collection vs "
                    "raw WGS/long-read/metagenomic assembly, short vs long "
                    "reads, extract plasmids when needed, and score them. "
                    "Reads (FASTQ) must be assembled first.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+",
                   help="FASTA file(s) or a directory of assemblies")
    p.add_argument("-o", "--output", default="plasrisk_auto")
    p.add_argument("--input-type", choices=["auto", "plasmid", "genome",
                                            "metagenome"], default="auto",
                   help="override the automatic input detection")
    p.add_argument("--read-type", choices=["auto", "short", "long"],
                   default="auto",
                   help="override the automatic short/long read detection")
    p.add_argument("--classifier",
                   choices=["auto", "plasflow", "genomad", "replicon"],
                   default="auto",
                   help="force a plasmid classifier instead of the auto "
                        "choice")
    p.add_argument("--genomad-db", default=None,
                   help="geNomad database path (or set GENOMAD_DB)")
    p.add_argument("--genomad-table", default=None,
                   help="precomputed geNomad summary TSV or directory")
    p.add_argument("--replicon-table", default=None,
                   help="precomputed abricate plasmidfinder TSV")
    p.add_argument("--plasflow-table", default=None,
                   help="precomputed PlasFlow prediction TSV")
    p.add_argument("--threshold", type=float, default=0.7)
    p.add_argument("--mode", choices=["lite", "full"], default="lite")
    p.add_argument("--aware-version", choices=["legacy", "2025"],
                   default="legacy")
    p.add_argument("--no-abricate", action="store_true")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--min-id", type=float, default=75.0)
    p.add_argument("--min-cov", type=float, default=50.0)
    p.add_argument("--min-contig", type=int, default=1000)
    p.add_argument("--complete-length", type=int, default=10000)
    p.add_argument("--json", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    return p


def _is_reads(path: str) -> bool:
    low = path.lower()
    return low.endswith(_READ_EXTS)


def auto_main(argv: Optional[List[str]] = None) -> int:
    args = build_auto_parser().parse_args(argv)
    fastas = _collect_fastas(args.inputs)
    if not fastas:
        print("[plasrisk auto] no FASTA inputs found", file=sys.stderr)
        return 2

    reads = [f for f in fastas if _is_reads(f)]
    if reads:
        print("[plasrisk auto] FASTQ reads are not assembled in-place. "
              "Assemble first, then re-run on the contigs FASTA:\n"
              "  short reads : megahit -1 R1.fq.gz -2 R2.fq.gz -o asm\n"
              "                (or metaSPAdes --meta)\n"
              "  long reads  : flye --meta --nano-raw reads.fq.gz -o asm\n"
              "                (or metaFlye)\n"
              "offending file(s): %s" % "; ".join(reads[:5]),
              file=sys.stderr)
        return 2

    os.makedirs(args.output, exist_ok=True)
    from .lookup import load_replicon_lookup
    lookup = None
    if not args.no_abricate:
        try:
            lookup = load_replicon_lookup()
        except Exception:
            lookup = None

    rows = []
    for i, fasta in enumerate(fastas, 1):
        sample_id = os.path.splitext(os.path.basename(fasta))[0]
        sub = os.path.join(args.output, sample_id)
        char = characterize_assembly(fasta)
        strat = choose_auto_strategy(char, args.input_type, args.read_type)
        if args.classifier != "auto":
            strat["classifier"] = args.classifier
            if args.classifier in ("genomad", "plasflow"):
                strat["route"] = "auto"
        if not args.quiet:
            print("[plasrisk auto] (%d/%d) %s | %s" % (
                i, len(fastas), sample_id, strat["rationale"]))

        row = {"sample_id": sample_id, "n_seqs": char["n_seqs"],
               "total_bp": char["total_bp"], "longest_bp": char["longest_bp"],
               "read_type": strat["read_type"], "action": strat["action"],
               "classifier": strat["classifier"], "route": strat["route"],
               "status": "", "n_plasmid": 0, "reason": ""}
        try:
            if strat["action"] == "direct":
                files = _score_direct(
                    fasta, sub, mode=args.mode,
                    aware_version=args.aware_version,
                    use_abricate=not args.no_abricate, min_id=args.min_id,
                    min_cov=args.min_cov, threads=args.threads, lookup=lookup)
                row.update(status="scored", route="direct")
                if not args.quiet:
                    print("    SCORED directly | %s" % files["scores"])
            else:
                res = process_one(
                    fasta, sub, classifier=strat["classifier"],
                    plasflow_table=args.plasflow_table,
                    genomad_table=args.genomad_table,
                    replicon_table=args.replicon_table,
                    threshold=args.threshold, genomad_db=args.genomad_db,
                    route=strat["route"], mode=args.mode,
                    aware_version=args.aware_version,
                    use_abricate=not args.no_abricate,
                    threads=args.threads, min_id=args.min_id,
                    min_cov=args.min_cov, min_contig=args.min_contig,
                    complete_length=args.complete_length,
                    write_json=args.json, quiet=args.quiet)
                # Long-read replicon scan found nothing: fall back to geNomad.
                if (res.status != "scored" and strat["classifier"] ==
                        "replicon" and genomad_available() and
                        (args.genomad_db or os.environ.get("GENOMAD_DB"))):
                    if not args.quiet:
                        print("    replicon scan found no plasmid; retrying "
                              "with geNomad")
                    res = process_one(
                        fasta, sub, classifier="genomad",
                        genomad_table=args.genomad_table,
                        threshold=args.threshold, genomad_db=args.genomad_db,
                        route="auto", mode=args.mode,
                        aware_version=args.aware_version,
                        use_abricate=not args.no_abricate,
                        threads=args.threads, min_id=args.min_id,
                        min_cov=args.min_cov, min_contig=args.min_contig,
                        complete_length=args.complete_length,
                        write_json=args.json, quiet=args.quiet)
                    row["classifier"] = res.classifier or "genomad"
                row.update(status=res.status, route=res.route,
                           n_plasmid=res.n_plasmid, reason=res.reason)
                if not args.quiet:
                    if res.status == "scored":
                        print("    SCORED | route=%s | plasmids=%d" % (
                            res.route, res.n_plasmid))
                    else:
                        print("    UNSCORABLE | %s" % res.reason)
        except Exception as exc:
            row.update(status="error", reason=str(exc))
            print("    ERROR: %s" % exc, file=sys.stderr)
        rows.append(row)

    status_df = pd.DataFrame(rows)
    spath = os.path.join(args.output, "auto_status.tsv")
    status_df.to_csv(spath, sep="\t", index=False)
    if not args.quiet:
        print("[plasrisk auto] status written to %s" % spath)
    return 0
