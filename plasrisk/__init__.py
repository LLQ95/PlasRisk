"""
PlasRisk - Plasmid Risk Assessment Tool
========================================

A 10-dimension weighted risk scoring framework for bacterial plasmids.
Weights derived by data-driven consensus (Random Forest MDG, LASSO, and
grid-search optimization) on 792,964 PIPdb PSCs.

S = 0.245*S_ARG + 0.110*S_VF + 0.204*S_MOB + 0.028*S_HOST
  + 0.003*S_REP + 0.181*S_SIZE + 0.211*S_BM
  + 0.002*S_GEO + 0.002*S_HAB + 0.015*S_GROW

Two scoring models are available:
  - PlasRiskScorer: 10-dimension weighted continuous model (default)
  - PIPdbScorer: original PIPdb 8-item ordinal model (use --model pipdb)

Reference: [to be updated upon publication]
"""

__version__ = "1.1.0"
__author__ = "PlasRisk Team"

from .scoring import (PlasRiskScorer, PIPdbScorer, PlasmidFeatures,
                      RISK_WEIGHTS, RISK_WEIGHTS_LITE, RISK_GRADES,
                      WEIGHT_SUM, get_scorer)
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
    "annotate_fasta",
    "AnnotationResult",
    "load_replicon_lookup",
]
