"""
conjscan.py - CONJScan/MacSyFinder system-level conjugative annotation.

Gene-based mobility annotation (oriTfinder/oriTDB/MOB-typer or abricate gene
names) detects T4CP, relaxase and oriT homologues one gene at a time. CONJScan
(Cury et al., 2020; in MacSyFinder v2, Neron et al., 2023) instead uses
curated HMM profiles together with genetic-organisation rules to classify
complete or likely conjugative systems, putative conjugative systems and
putative mobilizable systems.

This module provides the system-level, independent evidence stream used by
``annotate.py``. For complete plasmids the search is run in ordered-replicon
mode with the appropriate topology (linear for linear molecules and
short-read fragments, circular for closed circular plasmids); for fragmented
or mixed collections the content-only unordered mode can be selected.

Reference command (single replicon):

    macsyfinder --db-type ordered_replicon \\
                --sequence-db <proteins.faa> \\
                --models CONJScan/Plasmids all \\
                --replicon-topology {circular,linear}
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd

# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

CLASS_NONE = "none"
CLASS_MOBILIZABLE = "mobilizable"
CLASS_PCONJ = "putative_conjugative"
CLASS_LIKELY = "conjugative_likely"
CLASS_CONJ = "conjugative_complete"

CLASS_RANK = {
    CLASS_NONE: 0,
    CLASS_MOBILIZABLE: 1,
    CLASS_PCONJ: 2,
    CLASS_LIKELY: 3,
    CLASS_CONJ: 4,
}

# PlasRisk mobility vocabulary used by scoring.py / the locked catalogue.
PLASRISK_CLASSES = {
    CLASS_CONJ: "conjugative_complete",
    CLASS_LIKELY: "conjugative_likely",
    CLASS_PCONJ: "conjugative_putative",
    CLASS_MOBILIZABLE: "mobilizable",
    CLASS_NONE: "non-mobilizable",
}


def classify_model(model_fqn: str) -> str:
    """Map a CONJScan model fully-qualified name to a CONJScan class."""
    name = str(model_fqn).split("/")[-1].strip()
    low = name.lower()
    if low.startswith("likely_conj"):
        return CLASS_LIKELY
    if low.startswith("pconj"):
        return CLASS_PCONJ
    if low.startswith("pmob") or low.startswith("mob"):
        return CLASS_MOBILIZABLE
    if low.startswith("conj"):
        return CLASS_CONJ
    return CLASS_NONE


def _mpf_of_model(model_fqn: str) -> str:
    """Extract the MPF/type suffix, e.g. Conj_T4SS_typeG -> T4SS_typeG."""
    name = str(model_fqn).split("/")[-1].strip()
    low = name.lower()
    for prefix in ("likely_conj_", "pconj_", "conj_", "pmob_"):
        if low.startswith(prefix):
            return name[len(prefix):]
    return name


def union_mobility_class(gene_class: str, conj_class: str) -> str:
    """Return the higher-confidence class of the two evidence streams.

    Parameters
    ----------
    gene_class : str
        mobility class inferred from gene-level evidence
        (conjugative_complete / conjugative_likely / conjugative_putative /
        mobilizable / non-mobilizable).
    conj_class : str
        CONJScan class (see PLASRISK_CLASSES).
    """
    g = _normalise(gene_class)
    c = _normalise(conj_class)
    best = g if CLASS_RANK[g] >= CLASS_RANK[c] else c
    return PLASRISK_CLASSES[best]


def rank_of(cls: str) -> int:
    """Confidence rank of any PlasRisk/CONJScan class label."""
    return CLASS_RANK[_normalise(cls)]


def _normalise(cls: str) -> str:
    low = str(cls or "").lower()
    if "complete" in low:
        return CLASS_CONJ
    if "non-mobil" in low or "nonmobil" in low:
        return CLASS_NONE
    if "likely" in low:
        return CLASS_LIKELY
    if "putative" in low or "pconj" in low:
        return CLASS_PCONJ
    if "mobiliz" in low or "pmob" in low:
        return CLASS_MOBILIZABLE
    return CLASS_NONE


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class ConjScanResult:
    """System-level conjugative annotation for one replicon."""

    seq_id: str
    conjscan_class: str = CLASS_NONE
    mpf_type: str = ""
    wholeness: float = 0.0
    n_loci: Optional[int] = None
    has_relaxase: bool = False
    has_t4cp: bool = False
    has_mpf: bool = False
    unordered_class: str = ""
    raw_models: List[str] = field(default_factory=list)

    @property
    def is_conjugative(self) -> bool:
        """Complete or likely conjugative system."""
        return self.conjscan_class in (CLASS_CONJ, CLASS_LIKELY)

    def plasrisk_class(self) -> str:
        return PLASRISK_CLASSES.get(self.conjscan_class, "non-mobilizable")


# ---------------------------------------------------------------------------
# Availability
# ---------------------------------------------------------------------------

def conjscan_available() -> bool:
    """True when both macsyfinder and prodigal are on PATH."""
    return shutil.which("macsyfinder") is not None and shutil.which("prodigal") is not None


def missing_tools() -> List[str]:
    return [t for t in ("macsyfinder", "prodigal") if shutil.which(t) is None]


# ---------------------------------------------------------------------------
# Protein prediction
# ---------------------------------------------------------------------------

def _prodigal_proteins(seq_id: str, seq: str, workdir: str) -> Optional[str]:
    """Predict proteins for one replicon and return gembase-formatted FAA text.

    Prodigal runs in single mode per replicon, which is the correct mode for
    a complete plasmid; very short fragments fall back to meta mode.
    """
    nucl = os.path.join(workdir, f"{_safe(seq_id)}.fn")
    faa = os.path.join(workdir, f"{_safe(seq_id)}.faa")
    with open(nucl, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n{seq}\n")

    def run(mode: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["prodigal", "-i", nucl, "-a", faa,
             "-o", os.devnull, "-p", mode, "-q", "-f", "gff"],
            capture_output=True, text=True,
        )

    run("single" if len(seq) >= 20000 else "meta")
    if not os.path.exists(faa) or os.path.getsize(faa) == 0:
        run("meta")
    if not os.path.exists(faa) or os.path.getsize(faa) == 0:
        return None

    lines: List[str] = []
    n = 0
    with open(faa, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
                desc = line.rstrip("\n").split(" ", 1)
                extra = desc[1] if len(desc) > 1 else ""
                lines.append(f">{seq_id}_{n:05d} {extra}\n")
            else:
                lines.append(line)
    return "".join(lines)


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


# ---------------------------------------------------------------------------
# FASTA streaming (no biopython dependency)
# ---------------------------------------------------------------------------

def _iter_fasta(path: str):
    sid, chunks = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(chunks)
                sid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if sid is not None:
            yield sid, "".join(chunks)


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def _recover_replicon(value: str, known: set) -> str:
    v = str(value).strip()
    if v in known:
        return v
    if v.rsplit("_", 1)[0] in known:
        return v.rsplit("_", 1)[0]
    return v


def parse_systems(tsv_path: str, known_ids: set) -> Dict[str, ConjScanResult]:
    """Parse a best_solution/all_systems TSV into per-replicon results."""
    out: Dict[str, ConjScanResult] = {}
    if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
        return out
    try:
        # MacSyFinder writes several '#' comment lines before the header;
        # strip them manually so a '#'-prefixed header (older releases) is
        # still parsed as the column names.
        import io
        lines = []
        for line in open(tsv_path, encoding="utf-8", errors="replace"):
            if line.startswith("#") and not line[1:].lstrip().startswith(
                    "model_fqn"):
                continue
            lines.append(line.lstrip("# ") if line.startswith("#") else line)
        df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    except Exception:
        return out
    if df.empty or "model_fqn" not in df.columns:
        return out
    id_col = "replicon" if "replicon" in df.columns else df.columns[0]
    df["_rid"] = df[id_col].apply(lambda x: _recover_replicon(x, known_ids))

    for rid, grp in df.groupby("_rid", sort=False):
        best, best_rank = CLASS_NONE, 0
        mpf = ""
        whole, loci = 0.0, None
        genes = " ".join(grp.get("gene_name", pd.Series(dtype=str)).astype(str)).lower()
        for _, r in grp.iterrows():
            cls = classify_model(r.get("model_fqn", ""))
            if CLASS_RANK[cls] > best_rank:
                best_rank = CLASS_RANK[cls]
                best = cls
                mpf = _mpf_of_model(r.get("model_fqn", ""))
            try:
                whole = max(whole, float(r.get("sys_wholeness", 0)))
            except (TypeError, ValueError):
                pass
            try:
                loci = int(r.get("sys_loci")) if loci is None else min(loci, int(r.get("sys_loci")))
            except (TypeError, ValueError):
                pass
        out[rid] = ConjScanResult(
            seq_id=rid,
            conjscan_class=best,
            mpf_type=mpf,
            wholeness=round(whole, 3),
            n_loci=loci,
            has_relaxase=("relaxase" in genes or "mob" in genes),
            has_t4cp=("t4cp" in genes),
            has_mpf=("mpf" in genes or "t4ss" in genes),
            raw_models=sorted(set(grp.get("model_fqn", pd.Series(dtype=str)).astype(str))),
        )
    return out


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def annotate_conjugation(
    fasta_path: str,
    topology: str = "linear",
    threads: int = 4,
    db_type: str = "auto",
    keep_workdir: bool = False,
) -> Dict[str, ConjScanResult]:
    """Run CONJScan on every replicon in a nucleotide FASTA.

    Parameters
    ----------
    fasta_path : str
        Nucleotide FASTA; one or more plasmid/contig sequences.
    topology : str
        "linear" (default, conservative for short-read fragments) or
        "circular" (closed circular molecules).
    threads : int
        CPU threads for MacSyFinder.
    db_type : str
        "auto" uses ordered_replicon (gembase batch with a topology file)
        for ordered searches, which applies the gene-content and
        gene-organisation rules; "unordered" evaluates gene content only
        and is appropriate for highly fragmented/mixed metagenomic
        contigs.
    keep_workdir : bool
        Keep the temporary Prodigal/MacSyFinder working directory.

    Returns
    -------
    dict mapping sequence id to ConjScanResult. Replicons without any
    detected system are present with conjscan_class == "none".
    """
    if not conjscan_available():
        raise RuntimeError(
            "macsyfinder/prodigal not found. Install them with "
            "conda install -c conda-forge -c bioconda macsyfinder prodigal "
            "and the CONJScan models (pip install conjscan-models)."
        )
    if topology not in ("linear", "circular"):
        raise ValueError("topology must be 'linear' or 'circular'")
    if db_type == "auto":
        db_type = "ordered"

    records = list(_iter_fasta(fasta_path))
    if not records:
        return {}

    workdir = tempfile.mkdtemp(prefix="conjscan_")
    try:
        faa = os.path.join(workdir, "replicons.faa")
        known = set()
        with open(faa, "w", encoding="utf-8") as out:
            for sid, seq in records:
                known.add(sid)
                text = _prodigal_proteins(sid, seq, workdir)
                if text:
                    out.write(text)

        topo = os.path.join(workdir, "topology.txt")
        with open(topo, "w", encoding="utf-8") as fh:
            fh.write("# replicon: topology\n")
            for sid, _ in records:
                fh.write(f"{sid}: {topology}\n")

        if db_type == "unordered":
            results = _run_macsyfinder(
                workdir, faa, "unordered", threads, topology, topo, known)
        else:
            results = _run_macsyfinder(
                workdir, faa, "gembase", threads, topology, topo, known)

        # ensure every input replicon is represented
        for sid, _ in records:
            results.setdefault(sid, ConjScanResult(seq_id=sid))
        return results
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)


def _run_macsyfinder(workdir: str, faa: str, mode: str, threads: int,
                     topology: str, topo: str, known: set
                     ) -> Dict[str, ConjScanResult]:
    if mode == "unordered":
        outdir = os.path.join(workdir, "unordered")
        cmd = [
            "macsyfinder", "--db-type", "unordered",
            "--sequence-db", faa,
            "--models", "CONJScan/Plasmids", "all",
            "--out-dir", outdir, "-w", str(threads),
        ]
    else:
        outdir = os.path.join(workdir, "ordered")
        # gembase is the batch form of ordered_replicon; the single-replicon
        # equivalent uses --db-type ordered_replicon --replicon-topology.
        cmd = [
            "macsyfinder", "--db-type", "gembase",
            "--sequence-db", faa,
            "--models", "CONJScan/Plasmids", "all",
            "--topology-file", topo,
            "--out-dir", outdir, "-w", str(threads),
        ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log_path = os.path.join(outdir, "macsyfinder.stderr.log")
    try:
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(proc.stdout or "")
            fh.write("\n")
            fh.write(proc.stderr or "")
    except OSError:
        pass
    if proc.returncode != 0:
        raise RuntimeError(
            f"macsyfinder failed (rc={proc.returncode}): "
            f"{(proc.stderr or '')[-500:]}"
        )

    candidates = [
        "gembase_best_solution.tsv", "best_solution.tsv",
        "gembase_all_systems.tsv", "all_systems.tsv",
    ]
    for name in candidates:
        path = os.path.join(outdir, name)
        if os.path.exists(path):
            return parse_systems(path, known)
    return {}
