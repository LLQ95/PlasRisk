"""
PlasRisk - Plasmid Risk Assessment Tool
========================================

Data-driven weighted risk scoring framework for bacterial plasmids,
derived on 792,964 PIPdb PSCs (Random Forest MDG, LASSO, and grid-search
consensus weights).

Default output is the 5-dimension FASTA-only lite core:
S = 0.253*S_ARG + 0.237*S_BM + 0.201*S_MOB + 0.175*S_SIZE + 0.135*S_VF

The full 10-dimension model (--mode full) adds S_HOST, S_REP, S_GEO,
S_HAB and S_GROW, which require PIPdb-style metadata and are otherwise
imputed from replicon priors.

Reference: [to be updated upon publication]
"""

__version__ = "1.1.1"
__author__ = "PlasRisk Team"

from .scoring import (PlasRiskScorer, PIPdbScorer, PlasmidFeatures,
                      RISK_WEIGHTS, RISK_WEIGHTS_LITE, RISK_GRADES,
                      WEIGHT_SUM, get_scorer, size_rank_concordance)
from .annotate import annotate_fasta, AnnotationResult
from .lookup import load_replicon_lookup

__all__ = [
    "PlasRiskScorer",
    "PIPdbScorer",
    "PlasmidFeatures",
    "RISK_WEIGHTS",
    "RISK_WEIGHTS_LITE",
    "RISK_GRADES",
    "WEIGHT_SUM",
    "get_scorer",
    "size_rank_concordance",
    "annotate_fasta",
    "AnnotationResult",
    "load_replicon_lookup",
]
