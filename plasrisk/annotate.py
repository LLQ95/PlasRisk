"""
annotate.py - FASTA annotation using abricate and sequence analysis.

Wraps abricate for ARG/VF/replicon/BMRG annotation when available,
with graceful fallback to sequence-only scoring.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .scoring import (
    AUX_PATTERN,
    ORIT_PATTERN,
    RELAXASE_PATTERN,
    T4CP_PATTERN,
    PlasmidFeatures,
)


# ---------------------------------------------------------------------------
# FASTA parsing (lightweight, no biopython hard dependency)
# ---------------------------------------------------------------------------

def parse_fasta(filepath: str) -> Dict[str, str]:
    """Parse a FASTA file, return dict of {id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                # Take first word after '>' as ID
                current_id = line[1:].split()[0] if len(line) > 1 else "seq"
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)
    return sequences


def fasta_stats(seq: str) -> Dict:
    """Compute basic sequence statistics."""
    length = len(seq)
    gc = (seq.count("G") + seq.count("C")) / length * 100 if length > 0 else 0
    return {"length": length, "gc_percent": round(gc, 2)}


# ---------------------------------------------------------------------------
# abricate wrapper
# ---------------------------------------------------------------------------

ABRICATE_DATABASES = {
    "card": "ARG (CARD)",
    "resfinder": "ARG (ResFinder)",
    "vfdb": "Virulence factors (VFDB)",
    "plasmidfinder": "Replicon typing (PlasmidFinder)",
    "bacmet": "Biocide/metal resistance (BacMet)",
    "ecoli_vf": "E. coli virulence factors",
    "ncbi": "ARG (NCBI AMRFinder)",
    "megares": "ARG (MEGARes)",
    "argannot": "ARG (ARG-ANNOT)",
}


def abricate_available() -> bool:
    """Check if abricate is in PATH."""
    return shutil.which("abricate") is not None


def abricate_databases() -> List[str]:
    """List abricate databases currently installed."""
    try:
        result = subprocess.run(
            ["abricate", "--list"],
            capture_output=True, text=True, timeout=30,
        )
        dbs = []
        for line in result.stdout.strip().split("\n")[1:]:  # skip header
            parts = line.split()
            if parts:
                dbs.append(parts[0])
        return dbs
    except Exception:
        return []


def run_abricate(fasta_path: str, db: str, min_id: float = 75.0,
                min_cov: float = 50.0, threads: int = 4) -> pd.DataFrame:
    """
    Run abricate on a FASTA file with a specified database.

    Returns
    -------
    pd.DataFrame with abricate output columns, or empty DataFrame on failure.
    """
    if not abricate_available():
        return pd.DataFrame()
    try:
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False, mode="w") as tmp:
            tmp_path = tmp.name
        cmd = [
            "abricate",
            "--db", db,
            "--minid", str(min_id),
            "--mincov", str(min_cv),
            "--threads", str(threads),
            "--quiet",
            fasta_path,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0 or not result.stdout.strip():
            return pd.DataFrame()
        from io import StringIO
        df = pd.read_csv(StringIO(result.stdout), sep="\t")
        return df
    except Exception:
        return pd.DataFrame()
    finally:
        if "tmp_path" in locals():
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


# ---------------------------------------------------------------------------
# Annotation result container
# ---------------------------------------------------------------------------

@dataclass
class AnnotationResult:
    """Aggregated annotation for one FASTA file (may contain multiple plasmids)."""
    file: str
    features: List[PlasmidFeatures] = field(default_factory=list)
    databases_used: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Main annotation function
# ---------------------------------------------------------------------------

def annotate_fasta(
    fasta_path: str,
    use_abricate: bool = True,
    abricate_dbs: Optional[List[str]] = None,
    min_id: float = 75.0,
    min_cov: float = 50.0,
    threads: int = 4,
    lookup: Optional[pd.DataFrame] = None,
) -> AnnotationResult:
    """
    Annotate a plasmid FASTA file and return PlasmidFeatures for each sequence.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file (one or more plasmid sequences).
    use_abricate : bool
        Whether to use abricate for annotation (default True).
    abricate_dbs : list of str, optional
        abricate databases to query. If None, auto-detects available databases.
    min_id, min_cov : float
        Minimum identity and coverage thresholds for abricate.
    threads : int
        Number of threads for abricate (default 4).
    lookup : pd.DataFrame, optional
        Replicon lookup table for prior scores.

    Returns
    -------
    AnnotationResult
    """
    result = AnnotationResult(file=os.path.basename(fasta_path))

    # Parse sequences
    try:
        seqs = parse_fasta(fasta_path)
    except Exception as exc:
        result.warnings.append(f"Failed to parse FASTA: {exc}")
        return result

    if not seqs:
        result.warnings.append("No sequences found in FASTA file.")
        return result

    # Determine which databases to use
    available_dbs = []
    abricate_hits = {}
    if use_abricate and abricate_available():
        installed = abricate_databases()
        if abricate_dbs is None:
            # Priority order: card + vfdb + plasmidfinder + bacmet
            preferred = ["card", "vfdb", "plasmidfinder", "bacmet",
                         "resfinder", "ecoli_vf", "ncbi"]
            available_dbs = [d for d in preferred if d in installed]
        else:
            available_dbs = [d for d in abricate_dbs if d in installed]

        if not available_dbs:
            result.warnings.append(
                "abricate found but no standard databases installed. "
                "Run 'abricate-get_db --db card --force' etc."
            )
        else:
            result.databases_used = available_dbs
            for db in available_dbs:
                abricate_hits[db] = run_abricate(fasta_path, db, min_id, min_cov, threads)
    elif use_abricate and not abricate_available():
        result.warnings.append(
            "abricate not found in PATH. Install with: "
            "conda install -c bioconda abricate. "
            "Proceeding with sequence-only scoring."
        )

    # Build lookup dict for priors
    lookup_dict = {}
    if lookup is not None:
        for _, row in lookup.iterrows():
            lookup_dict[row["replicon_primary"]] = row.to_dict()

    # Process each sequence
    for seq_id, seq in seqs.items():
        feat = PlasmidFeatures(seq_id=seq_id, length_bp=len(seq))

        # Extract hits for this sequence
        arg_names = []
        vf_names = []
        vf_categories = []
        bm_names = []
        replicons = []
        all_gene_names = []

        for db, df in abricate_hits.items():
            if df.empty:
                continue
            # abricate output: #FILE SEQUENCE START END GENE COVERAGE ...
            seq_hits = df[df["SEQUENCE"].astype(str) == seq_id] if "SEQUENCE" in df.columns else df
            if seq_hits.empty:
                continue

            for _, hit in seq_hits.iterrows():
                gene = str(hit.get("GENE", "")).strip()
                product = str(hit.get("PRODUCT", "")).strip()
                accession = str(hit.get("ACCESSION", "")).strip()
                name = gene if gene and gene != "nan" else product
                if not name or name == "nan":
                    continue
                all_gene_names.append(name)

                if db in ("card", "resfinder", "ncbi", "argannot", "megares"):
                    # Clean gene name: take first allele before parenthesis
                    clean = re.split(r"[(_;]", name)[0].strip()
                    if clean:
                        arg_names.append(name)
                elif db in ("vfdb", "ecoli_vf"):
                    vf_names.append(name)
                    # VFDB PRODUCT often contains category info
                    if product and product != "nan":
                        vf_categories.append(product)
                elif db == "plasmidfinder":
                    replicons.append(name)
                elif db == "bacmet":
                    bm_names.append(name)

        feat.arg_names = arg_names
        feat.vf_names = vf_names
        feat.vf_categories = vf_categories
        feat.bm_gene_names = bm_names

        # Replicon: take first PlasmidFinder hit, normalize
        if replicons:
            feat.replicon = _normalize_replicon(replicons[0])
        else:
            feat.replicon = ""

        # Mobility: infer from gene names in all abricate hits
        all_text = ",".join(all_gene_names).upper()
        if all_text:
            feat.has_t4cp = bool(T4CP_PATTERN.search(all_text))
            feat.has_relaxase = bool(RELAXASE_PATTERN.search(all_text))
            feat.has_oriT = bool(ORIT_PATTERN.search(all_text))
            feat.has_auxiliary = bool(AUX_PATTERN.search(all_text))
        else:
            feat.has_t4cp = False
            feat.has_relaxase = False
            feat.has_oriT = False
            feat.has_auxiliary = False

        # Apply replicon lookup priors
        if feat.replicon and feat.replicon in lookup_dict:
            priors = lookup_dict[feat.replicon]
            feat.s_rep_prior = priors.get("S_REP", 0.3)
            feat.s_geo_prior = priors.get("S_GEO", 0.3)
            feat.s_hab_prior = priors.get("S_HAB", 0.3)
            feat.s_grow_prior = priors.get("S_GROW", 0.3)
            feat.s_host_prior = priors.get("S_HOST", 0.5)

        result.features.append(feat)

    return result


def _normalize_replicon(name: str) -> str:
    """
    Normalize PlasmidFinder replicon names to match PIPdb naming.
    E.g., 'IncFII_1__AP001918' -> 'IncFII', 'Col440I_1__CP014494_1' -> 'Col440I'
    """
    # Remove allele numbers and accession suffixes
    clean = re.sub(r"_\d+__.*$", "", name)
    clean = re.sub(r"_\d+$", "", clean)
    clean = re.sub(r"\(.*?\)", lambda m: m.group(0), clean)  # keep parentheses
    return clean.strip()
