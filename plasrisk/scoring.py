"""
scoring.py - Ten-dimension plasmid risk scoring model.

Each component S_i is normalized to [0, 1].
Composite score:  S = sum(w_i * S_i), normalized by sum(w_i).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Weight scheme (data-driven consensus weights from RF-MDG / LASSO / grid-search
# optimization on 792,964 PIPdb PSCs; normalized to sum = 1.0)
# ---------------------------------------------------------------------------
RISK_WEIGHTS: Dict[str, float] = {
    "S_ARG":  0.2448,  # Antimicrobial resistance gene burden
    "S_VF":   0.1096,  # Virulence factor burden
    "S_MOB":  0.2041,  # Mobility / conjugation potential
    "S_HOST": 0.0282,  # Host range breadth
    "S_REP":  0.0030,  # Replicon-specific risk prior
    "S_SIZE": 0.1808,  # Plasmid size (cargo capacity)
    "S_BM":   0.2112,  # Biocide / metal resistance (co-selection)
    "S_GEO":  0.0015,  # Geographic spread
    "S_HAB":  0.0022,  # Habitat breadth (One Health)
    "S_GROW": 0.0147,  # Temporal growth rate
}

WEIGHT_SUM = sum(RISK_WEIGHTS.values())  # 1.0001 ~ 1.0

# ---------------------------------------------------------------------------
# Risk grade thresholds (on normalized S_norm in [0, 1])
# ---------------------------------------------------------------------------
RISK_GRADES = [
    (0.60, "A", "Very High"),
    (0.45, "B", "High"),
    (0.30, "C", "Moderate"),
    (0.15, "D", "Low"),
    (0.00, "E", "Minimal"),
]

# ---------------------------------------------------------------------------
# High-risk ARG patterns (clinically critical gene families)
# ---------------------------------------------------------------------------
HIGH_RISK_ARG_PATTERNS = {
    "mcr":       re.compile(r"mcr", re.I),
    "blaNDM":    re.compile(r"NDM", re.I),
    "blaKPC":    re.compile(r"KPC", re.I),
    "blaOXA48":  re.compile(r"OXA-?(48|181|232|244)", re.I),
    "blaCTX-M":  re.compile(r"CTX-M", re.I),
    "tetX":      re.compile(r"tet\(?X\)?", re.I),
    "van":       re.compile(r"van[A-Z]", re.I),
    "cfr":       re.compile(r"\bcfr\b", re.I),
    "optrA":     re.compile(r"optrA", re.I),
    "poxtA":     re.compile(r"poxtA", re.I),
    "blaIMP":    re.compile(r"\bIMP", re.I),
    "blaVIM":    re.compile(r"\bVIM", re.I),
    "rmt":       re.compile(r"rmt[A-H]", re.I),
}

# WHO-priority ARG indicator (broad set of clinically important families)
WHO_ARG_PATTERNS = re.compile(
    r"mcr|NDM|KPC|OXA-?48|CTX-M|CMY|DHA|VIM|IMP|tet\(?X\)?|van|"
    r"cfr|optrA|poxtA|qnr|aac.6.-Ib-cr|rmt|fos[A-Z]|SHV|TEM|PER|VEB|GES",
    re.I,
)

# ---------------------------------------------------------------------------
# Mobility gene markers (keyword-based detection in annotation)
# ---------------------------------------------------------------------------
T4CP_PATTERN = re.compile(
    r"\b(traD|traG|virD4|trwB|T4CP|coupling)\b", re.I)
RELAXASE_PATTERN = re.compile(
    r"\b(traI|mobA|nikB|traA|relaxase|traH|traJ_rel|mobC)\b", re.I)
ORIT_PATTERN = re.compile(
    r"\b(oriT|nic|traI_rel|mob)\b", re.I)
AUX_PATTERN = re.compile(
    r"\b(traJ|traK|traM|traY|traN|traO|traP|trb[A-Z]|Pil|taxC)\b", re.I)

# ---------------------------------------------------------------------------
# Biocide / metal resistance patterns
# ---------------------------------------------------------------------------
MER_PATTERN = re.compile(r"\bmer[A-Z]?\b|mercury", re.I)
QAC_PATTERN = re.compile(r"qac|quaternary|disinfectant", re.I)
ARS_COP_PATTERN = re.compile(
    r"\bars[A-Z]?\b|\bcop[A-Z]?\b|\bsil[A-Z]?\b|czc[A-Z]?|cad[A-Z]?|pco[A-Z]",
    re.I,
)

# VF category bonuses
EXOTOXIN_PATTERN = re.compile(r"exotoxin|toxin|enterotoxin|hemolysin|cytolysin", re.I)
EFFECTOR_PATTERN = re.compile(r"effector delivery|type III|type IV|T3SS|T4SS|secretion", re.I)


# ---------------------------------------------------------------------------
# Data class for plasmid features
# ---------------------------------------------------------------------------
@dataclass
class PlasmidFeatures:
    """Container for all features needed to score a single plasmid."""

    # Sequence
    seq_id: str = ""
    length_bp: int = 0

    # Annotation counts
    arg_names: List[str] = field(default_factory=list)
    vf_names: List[str] = field(default_factory=list)
    vf_categories: List[str] = field(default_factory=list)
    bm_gene_names: List[str] = field(default_factory=list)

    # Replicon
    replicon: str = ""

    # Mobility elements (booleans)
    has_t4cp: Optional[bool] = None
    has_relaxase: Optional[bool] = None
    has_oriT: Optional[bool] = None
    has_auxiliary: Optional[bool] = None
    mobility_class: str = ""  # if pre-classified

    # Host range (if known; otherwise None -> use replicon lookup)
    n_host_genera: Optional[int] = None
    n_countries: Optional[int] = None
    n_habitats: Optional[int] = None
    annual_growth_rate: Optional[float] = None

    # Replicon lookup priors (filled by scorer if not provided)
    s_rep_prior: Optional[float] = None
    s_geo_prior: Optional[float] = None
    s_hab_prior: Optional[float] = None
    s_grow_prior: Optional[float] = None
    s_host_prior: Optional[float] = None

    # Convenience: any extra metadata
    metadata: Dict = field(default_factory=dict)

    @property
    def n_arg(self) -> int:
        return len(self.arg_names)

    @property
    def n_vf(self) -> int:
        return len(self.vf_names)

    @property
    def n_bm(self) -> int:
        return len(self.bm_gene_names)

    @property
    def high_risk_args(self) -> List[str]:
        hits = []
        for name in self.arg_names:
            for label, pat in HIGH_RISK_ARG_PATTERNS.items():
                if pat.search(name):
                    hits.append(label)
                    break
        return hits

    @property
    def has_who_arg(self) -> bool:
        return any(WHO_ARG_PATTERNS.search(n) for n in self.arg_names)

    @property
    def has_high_risk_arg(self) -> bool:
        return len(self.high_risk_args) > 0


# ---------------------------------------------------------------------------
# Scorer
# ---------------------------------------------------------------------------
class PlasRiskScorer:
    """Compute 10-dimension PlasRisk scores for plasmids."""

    def __init__(self, replicon_lookup: Optional[pd.DataFrame] = None):
        """
        Parameters
        ----------
        replicon_lookup : pd.DataFrame, optional
            Lookup table with columns:
            replicon_primary, S_REP, S_GEO, S_HAB, S_GROW, S_HOST
            If None, built-in defaults are used.
        """
        if replicon_lookup is not None:
            self.lookup = replicon_lookup.set_index("replicon_primary")
        else:
            self.lookup = None

        # Default priors for unknown replicons (moderate / neutral)
        self.defaults = {
            "S_REP":  0.30,
            "S_GEO":  0.30,
            "S_HAB":  0.30,
            "S_GROW": 0.30,
            "S_HOST": 0.50,
        }

    # ------------------------------------------------------------------
    # Individual component scorers
    # ------------------------------------------------------------------

    @staticmethod
    def score_arg(feat: PlasmidFeatures) -> float:
        """
        S_ARG: ARG burden.
        0 if no ARGs.
        0.25 base + min(n_arg*0.05, 0.35) + 0.20*WHO + 0.20*high-risk, cap 1.0.
        """
        n = feat.n_arg
        if n == 0:
            return 0.0
        base = 0.25
        count_bonus = min(n * 0.05, 0.35)
        who_bonus = 0.20 if feat.has_who_arg else 0.0
        hr_bonus = 0.20 if feat.has_high_risk_arg else 0.0
        return min(base + count_bonus + who_bonus + hr_bonus, 1.0)

    @staticmethod
    def score_vf(feat: PlasmidFeatures) -> float:
        """
        S_VF: virulence factor burden.
        0 if no VFs.
        0.30 base + min(n_vf*0.03, 0.40) + 0.15*exotoxin + 0.15*effector, cap 1.0.
        """
        n = feat.n_vf
        if n == 0:
            return 0.0
        base = 0.30
        count_bonus = min(n * 0.03, 0.40)
        cat_text = ",".join(feat.vf_categories) if feat.vf_categories else ""
        name_text = ",".join(feat.vf_names)
        all_text = cat_text + "," + name_text
        exo_bonus = 0.15 if EXOTOXIN_PATTERN.search(all_text) else 0.0
        eff_bonus = 0.15 if EFFECTOR_PATTERN.search(all_text) else 0.0
        return min(base + count_bonus + exo_bonus + eff_bonus, 1.0)

    @staticmethod
    def score_mob(feat: PlasmidFeatures) -> float:
        """
        S_MOB: mobility potential based on conjugation elements.
        0.10 base if any element + 0.35*T4CP + 0.25*relaxase
        + 0.15*oriT + 0.15*auxiliary, cap 1.0.
        """
        # If mobility_class is pre-classified, use it
        if feat.mobility_class:
            mc = feat.mobility_class.lower()
            if "complete" in mc:
                return 1.00
            if "likely" in mc or "conjugative" in mc:
                return 0.75
            if "mobilizable" in mc:
                return 0.45
            return 0.10

        # Otherwise infer from element booleans
        elements = [feat.has_t4cp, feat.has_relaxase,
                    feat.has_oriT, feat.has_auxiliary]
        if all(e is None for e in elements):
            return 0.10  # unknown -> low default

        has_t4cp = bool(feat.has_t4cp)
        has_rel = bool(feat.has_relaxase)
        has_oriT = bool(feat.has_oriT)
        has_aux = bool(feat.has_auxiliary)

        if not any([has_t4cp, has_rel, has_oriT, has_aux]):
            return 0.05

        score = 0.10
        if has_t4cp:
            score += 0.35
        if has_rel:
            score += 0.25
        if has_oriT:
            score += 0.15
        if has_aux:
            score += 0.15
        return min(score, 1.0)

    @staticmethod
    def score_size(length_bp: int) -> float:
        """
        S_SIZE: sigmoidal function of plasmid length.
        Midpoint at 30 kb, steepness parameter 15 kb.
        Small plasmids (<10 kb) score low; large (>80 kb) score high.
        """
        if length_bp <= 0:
            return 0.0
        length_kb = length_bp / 1000.0
        return 1.0 / (1.0 + math.exp(-(length_kb - 30.0) / 15.0))

    @staticmethod
    def score_bm(feat: PlasmidFeatures) -> float:
        """
        S_BM: biocide/metal resistance (co-selection potential).
        0 if no BMRGs.
        0.25 base + min(n_bm*0.04, 0.35) + 0.15*mer + 0.15*qac + 0.10*ars/cop/sil.
        """
        n = feat.n_bm
        if n == 0:
            return 0.0
        base = 0.25
        count_bonus = min(n * 0.04, 0.35)
        bm_text = ",".join(feat.bm_gene_names)
        mer_bonus = 0.15 if MER_PATTERN.search(bm_text) else 0.0
        qac_bonus = 0.15 if QAC_PATTERN.search(bm_text) else 0.0
        ars_bonus = 0.10 if ARS_COP_PATTERN.search(bm_text) else 0.0
        return min(base + count_bonus + mer_bonus + qac_bonus + ars_bonus, 1.0)

    def score_host(self, feat: PlasmidFeatures) -> float:
        """
        S_HOST: host range breadth.
        If n_host_genera known, map to [0,1].
        Otherwise use replicon lookup prior.
        """
        if feat.n_host_genera is not None:
            n = feat.n_host_genera
            if n <= 1:
                return 0.20
            elif n <= 3:
                return 0.50
            elif n <= 10:
                return 0.75
            else:
                return 1.00
        if feat.s_host_prior is not None:
            return float(feat.s_host_prior)
        return self._lookup_prior(feat.replicon, "S_HOST")

    def score_rep(self, feat: PlasmidFeatures) -> float:
        if feat.s_rep_prior is not None:
            return float(feat.s_rep_prior)
        return self._lookup_prior(feat.replicon, "S_REP")

    def score_geo(self, feat: PlasmidFeatures) -> float:
        if feat.n_countries is not None:
            n = feat.n_countries
            return min(n / 45.0, 1.0)  # 45+ countries -> 1.0
        if feat.s_geo_prior is not None:
            return float(feat.s_geo_prior)
        return self._lookup_prior(feat.replicon, "S_GEO")

    def score_hab(self, feat: PlasmidFeatures) -> float:
        if feat.n_habitats is not None:
            n = feat.n_habitats
            return min(n / 8.0, 1.0)  # 8+ habitats -> 1.0
        if feat.s_hab_prior is not None:
            return float(feat.s_hab_prior)
        return self._lookup_prior(feat.replicon, "S_HAB")

    def score_grow(self, feat: PlasmidFeatures) -> float:
        if feat.annual_growth_rate is not None:
            g = feat.annual_growth_rate
            # Map growth rate to [0,1]: negative -> 0.1, 0 -> 0.3,
            # 0.1 -> 0.7, 0.3+ -> 1.0
            return float(np.clip(0.3 + g * 2.0, 0.0, 1.0))
        if feat.s_grow_prior is not None:
            return float(feat.s_grow_prior)
        return self._lookup_prior(feat.replicon, "S_GROW")

    # ------------------------------------------------------------------
    # Lookup helpers
    # ------------------------------------------------------------------

    def _lookup_prior(self, replicon: str, column: str) -> float:
        """Retrieve a prior from the replicon lookup table."""
        if not replicon or self.lookup is None:
            return self.defaults.get(column, 0.3)
        # Try exact match first, then prefix match (e.g., "IncFII" for "IncFII(K)")
        if replicon in self.lookup.index:
            val = self.lookup.loc[replicon, column]
            return float(val) if pd.notna(val) else self.defaults[column]
        for idx in self.lookup.index:
            if replicon.startswith(idx) or idx.startswith(replicon):
                val = self.lookup.loc[idx, column]
                return float(val) if pd.notna(val) else self.defaults[column]
        return self.defaults[column]

    # ------------------------------------------------------------------
    # Main scoring
    # ------------------------------------------------------------------

    def score(self, feat: PlasmidFeatures) -> Dict:
        """
        Compute all 10 component scores and the composite.

        Returns
        -------
        dict with keys: seq_id, length_bp, replicon, n_ARG, n_VF, n_BM,
                        S_ARG ... S_GROW, S_total, S_norm, grade, grade_label,
                        high_risk_genes, mobility_class
        """
        components = {
            "S_ARG":  self.score_arg(feat),
            "S_VF":   self.score_vf(feat),
            "S_MOB":  self.score_mob(feat),
            "S_HOST": self.score_host(feat),
            "S_REP":  self.score_rep(feat),
            "S_SIZE": self.score_size(feat.length_bp),
            "S_BM":   self.score_bm(feat),
            "S_GEO":  self.score_geo(feat),
            "S_HAB":  self.score_hab(feat),
            "S_GROW": self.score_grow(feat),
        }

        s_total = sum(RISK_WEIGHTS[k] * components[k] for k in components)
        s_norm = s_total / WEIGHT_SUM

        grade, grade_label = self._grade(s_norm)

        # Determine mobility class if not pre-set
        mob_class = feat.mobility_class or self._infer_mob_class(components["S_MOB"])

        result = {
            "seq_id": feat.seq_id,
            "length_bp": feat.length_bp,
            "length_kb": round(feat.length_bp / 1000.0, 2),
            "replicon": feat.replicon or "Unknown",
            "n_ARG": feat.n_arg,
            "n_VF": feat.n_vf,
            "n_BM": feat.n_bm,
            "high_risk_genes": ";".join(feat.high_risk_args) if feat.high_risk_args else "",
            "mobility_class": mob_class,
        }
        result.update({k: round(v, 4) for k, v in components.items()})
        result["S_total"] = round(s_total, 4)
        result["S_norm"] = round(s_norm, 4)
        result["grade"] = grade
        result["grade_label"] = grade_label
        return result

    def score_dataframe(self, features: List[PlasmidFeatures]) -> pd.DataFrame:
        """Score multiple plasmids and return a sorted DataFrame."""
        rows = [self.score(f) for f in features]
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values("S_norm", ascending=False).reset_index(drop=True)
            df.insert(0, "rank", range(1, len(df) + 1))
        return df

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _grade(s_norm: float) -> Tuple[str, str]:
        for threshold, grade, label in RISK_GRADES:
            if s_norm >= threshold:
                return grade, label
        return "E", "Minimal"

    @staticmethod
    def _infer_mob_class(s_mob: float) -> str:
        if s_mob >= 0.85:
            return "conjugative_complete"
        elif s_mob >= 0.55:
            return "conjugative_likely"
        elif s_mob >= 0.25:
            return "mobilizable"
        else:
            return "non-mobilizable"
