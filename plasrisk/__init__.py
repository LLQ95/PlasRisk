"""
PlasRisk - Plasmid Risk Assessment Tool
========================================

A 10-dimension weighted risk scoring framework for bacterial plasmids.
Weights derived by data-driven consensus (Random Forest MDG, LASSO, and
grid-search optimization) on 792,964 PIPdb PSCs.

S = 0.245*S_ARG + 0.110*S_VF + 0.204*S_MOB + 0.028*S_HOST
  + 0.003*S_REP + 0.181*S_SIZE + 0.211*S_BM
  + 0.002*S_GEO + 0.002*S_HAB + 0.015*S_GROW

Reference: [to be updated upon publication]
"""

__version__ = "1.0.0"
__author__ = "PlasRisk Team"

from .scoring import PlasRiskScorer, PlasmidFeatures, RISK_WEIGHTS, RISK_GRADES, WEIGHT_SUM
from .annotate import annotate_fasta, AnnotationResult
from .lookup import load_replicon_lookup

__all__ = [
    "PlasRiskScorer",
    "PlasmidFeatures",
    "RISK_WEIGHTS",
    "RISK_GRADES",
    "WEIGHT_SUM",
    "annotate_fasta",
    "AnnotationResult",
    "load_replicon_lookup",
]
