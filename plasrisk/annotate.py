"""
annotate.py - FASTA annotation using abricate and sequence analysis.

Wraps abricate for ARG/VF/replicon/BMRG/IS/transposon annotation when available,
with graceful fallback to replicon-empirical priors when databases are missing.
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
    CARD_BM_PATTERN,
    GLOBAL_MEDIANS,
    HOUSEKEEPING_GENES,
    ORIT_PATTERN,
    RELAXASE_PATTERN,
    T4CP_PATTERN,
    WHO_PATHOGEN_PATTERN,
    PlasmidFeatures,
)

# IS detection: gene names starting with "IS" followed by digits/letters
IS_PATTERN = re.compile(r"\bIS[A-Z]?[\d_]+", re.I)
# Integron detection: integrase genes
INTEGRON_PATTERN = re.compile(
    r"\bintI[1-4]\b|integrase|integron", re.I)


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
    "bacmet2": "Biocide/metal resistance (BacMet2, protein)",
    "ecoli_vf": "E. coli virulence factors",
    "ncbi": "ARG (NCBI AMRFinder)",
    "megares": "ARG (MEGARes)",
    "argannot": "ARG (ARG-ANNOT)",
    "ISfinder": "Insertion sequences (ISfinder)",
    "Tn": "Transposons (Tn)",
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
        cmd = [
            "abricate",
            "--db", db,
            "--minid", str(min_id),
            "--mincov", str(min_cov),
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
        Replicon lookup table for prior scores (all 10 dimensions).

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
    mobility_dbs_run = False
    if use_abricate and abricate_available():
        installed = abricate_databases()
        if abricate_dbs is None:
            # Priority order: card + resfinder + vfdb + plasmidfinder + bacmet2
            # + ISfinder + Tn for mobility/IS detection
            preferred = ["card", "resfinder", "vfdb", "plasmidfinder", "bacmet2",
                         "ISfinder", "Tn", "ecoli_vf", "ncbi"]
            available_dbs = [d for d in preferred if d in installed]
        else:
            available_dbs = [d for d in abricate_dbs if d in installed]

        if not available_dbs:
            result.warnings.append(
                "abricate found but no standard databases installed. "
                "Run 'abricate-get_db --db card --force' etc. "
                "Proceeding with replicon-empirical priors."
            )
        else:
            result.databases_used = available_dbs
            for db in available_dbs:
                abricate_hits[db] = run_abricate(fasta_path, db, min_id, min_cov, threads)
            # Check if mobility-specific databases were run.
            # CARD/ResFinder/NCBI do not reliably detect T4CP/relaxase/oriT,
            # so only ISfinder and Tn count as mobility-specific.
            mobility_dbs_run = bool(set(available_dbs) & {"ISfinder", "Tn"})
    elif use_abricate and not abricate_available():
        result.warnings.append(
            "abricate not found in PATH. Install with: "
            "conda install -c bioconda abricate. "
            "Proceeding with replicon-empirical priors."
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
        arg_products = {}
        arg_resistance = {}
        vf_names = []
        vf_categories = []
        bm_names = []
        replicons = []
        all_gene_names = []
        is_from_isfinder = set()

        for db, df in abricate_hits.items():
            if df.empty:
                continue
            # abricate SEQUENCE column may have version suffix (e.g., NC_019050.1)
            if "SEQUENCE" in df.columns:
                seq_col = df["SEQUENCE"].astype(str).apply(
                    lambda x: re.sub(r"\.\d+$", "", x))
                seq_hits = df[seq_col == seq_id]
            else:
                seq_hits = df
            if seq_hits.empty:
                continue

            for _, hit in seq_hits.iterrows():
                gene = str(hit.get("GENE", "")).strip()
                product = str(hit.get("PRODUCT", "")).strip()
                accession = str(hit.get("ACCESSION", "")).strip()
                resistance = str(hit.get("RESISTANCE", "")).strip()
                name = gene if gene and gene != "nan" else product
                if not name or name == "nan":
                    continue
                all_gene_names.append(name)

                if db in ("card", "resfinder", "ncbi", "argannot", "megares"):
                    # Clean gene name: take first allele before parenthesis
                    clean = re.split(r"[(_;]", name)[0].strip()
                    if clean:
                        # Exclude chromosomal housekeeping genes
                        if not HOUSEKEEPING_GENES.match(clean):
                            arg_names.append(name)
                            arg_products[name] = product
                            arg_resistance[name] = resistance
                elif db in ("vfdb", "ecoli_vf"):
                    vf_names.append(name)
                    if product and product != "nan":
                        vf_categories.append(product)
                elif db == "plasmidfinder":
                    replicons.append(name)
                elif db == "bacmet2":
                    bm_names.append(name)
                elif db == "ISfinder":
                    is_from_isfinder.add(name)
                    all_gene_names.append(name)

        # CARD-based fallback for biocide/metal genes when BacMet is
        # unavailable or returns no hits.
        if not bm_names:
            for name in arg_names:
                if CARD_BM_PATTERN.search(name):
                    bm_names.append(name)

        feat.arg_names = arg_names
        feat.metadata["arg_products"] = arg_products
        feat.metadata["arg_resistance"] = arg_resistance
        feat.vf_names = vf_names
        feat.vf_categories = vf_categories
        feat.bm_gene_names = bm_names

        # Replicon: collect ALL PlasmidFinder hits, normalize each,
        # join with ';' for multi-replicon plasmids (e.g., IncFII;IncFIA;IncR)
        if replicons:
            normalized = [_normalize_replicon(r) for r in replicons]
            # Deduplicate while preserving order
            seen = set()
            unique_reps = []
            for r in normalized:
                if r and r not in seen:
                    seen.add(r)
                    unique_reps.append(r)
            feat.replicon = ";".join(unique_reps)
        else:
            feat.replicon = ""

        # Mobility: infer from gene names in all abricate hits
        all_text = ",".join(all_gene_names).upper()
        if all_text:
            feat.has_t4cp = bool(T4CP_PATTERN.search(all_text))
            feat.has_relaxase = bool(RELAXASE_PATTERN.search(all_text))
            feat.has_oriT = bool(ORIT_PATTERN.search(all_text))
            feat.has_auxiliary = bool(AUX_PATTERN.search(all_text))
            feat.has_integron = bool(INTEGRON_PATTERN.search(all_text))
            # IS count: combine ISfinder hits with pattern detection
            pattern_is = set(IS_PATTERN.findall(all_text))
            feat.n_is = len(is_from_isfinder | pattern_is)
            # mobility_annotated = True only if we have actual mobility
            # evidence OR mobility-specific databases were run.
            # abricate CARD/VFDB alone do not reliably detect T4CP/relaxase,
            # so we should not claim mobility was assessed if no mobility
            # genes were found.
            has_mob_evidence = (feat.has_t4cp or feat.has_relaxase or
                                feat.has_oriT or feat.has_auxiliary or
                                feat.has_integron or feat.n_is > 0)
            feat.mobility_annotated = has_mob_evidence or mobility_dbs_run
        else:
            feat.has_t4cp = False
            feat.has_relaxase = False
            feat.has_oriT = False
            feat.has_auxiliary = False
            feat.has_integron = False
            feat.n_is = 0
            # No genes found at all; only claim mobility was annotated if
            # mobility-specific databases (ISfinder/Tn) were actually run.
            feat.mobility_annotated = mobility_dbs_run

        # Detect WHO priority pathogen from sequence header or metadata
        header_text = seq_id.upper()
        if WHO_PATHOGEN_PATTERN.search(header_text):
            feat.host_is_who_priority = True

        # Apply replicon lookup priors (all 10 dimensions)
        # For multi-replicon plasmids, take the max prior across replicons
        if feat.replicon and lookup_dict:
            rep_list = [r.strip() for r in re.split(r"[;,/]", feat.replicon) if r.strip()]
            best_priors = {}
            for dim in ("S_ARG", "S_VF", "S_MOB", "S_HOST", "S_REP",
                        "S_SIZE", "S_BM", "S_GEO", "S_HAB", "S_GROW"):
                best_val = None
                for rep in rep_list:
                    # Exact match
                    if rep in lookup_dict:
                        val = lookup_dict[rep].get(dim)
                        if val is not None and val == val:  # not NaN
                            if best_val is None or val > best_val:
                                best_val = val
                    else:
                        # Prefix match
                        for rep_key in lookup_dict:
                            if rep.startswith(rep_key) or rep_key.startswith(rep):
                                val = lookup_dict[rep_key].get(dim)
                                if val is not None and val == val:
                                    if best_val is None or val > best_val:
                                        best_val = val
                                break
                if best_val is not None:
                    best_priors[dim] = best_val

            feat.s_rep_prior = best_priors.get("S_REP", GLOBAL_MEDIANS["S_REP"])
            feat.s_geo_prior = best_priors.get("S_GEO", GLOBAL_MEDIANS["S_GEO"])
            feat.s_hab_prior = best_priors.get("S_HAB", GLOBAL_MEDIANS["S_HAB"])
            feat.s_grow_prior = best_priors.get("S_GROW", GLOBAL_MEDIANS["S_GROW"])
            feat.s_host_prior = best_priors.get("S_HOST", GLOBAL_MEDIANS["S_HOST"])
            feat.s_mob_prior = best_priors.get("S_MOB", GLOBAL_MEDIANS["S_MOB"])

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
