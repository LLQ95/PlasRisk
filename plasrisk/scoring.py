#!/usr/bin/env python3
"""
scoring.py - Risk scoring functions for PlasRisk.

Implements 10 dimensions:
  S_ARG  - Antimicrobial resistance gene burden (WHO AWaRe weighted)
  S_VF   - Virulence factor burden
  S_MOB  - Mobility / conjugative potential
  S_HOST - Host range
  S_REP  - Replicon type
  S_SIZE - Plasmid size
  S_BMG  - Biocide/metal resistance
  S_GEO  - Geographic spread
  S_HAB  - Habitat breadth
  S_GROW - Temporal growth rate

Default scorer is the 5-dimension FASTA-only lite core; the full 10-dimension
model is selected with mode="full". The original PIPdb 8-item ordinal index
is provided as PIPdbScorer / mode="ordinal".
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# =====================================================================
# Final consensus weights (locked; manuscript Table 1)
# Order: descending weight
# =====================================================================
RISK_WEIGHTS: Dict[str, float] = {
    "S_ARG":  0.237,
    "S_BMG":  0.222,
    "S_MOB":  0.189,
    "S_SIZE": 0.164,
    "S_VF":   0.126,
    "S_HOST": 0.030,
    "S_HAB":  0.019,
    "S_GEO":  0.005,
    "S_REP":  0.004,
    "S_GROW": 0.002,
}

# Lite 5-dimension FASTA-only weights (renormalized from full weights over
# the five sequence-derived dimensions; manuscript Table 6)
RISK_WEIGHTS_LITE: Dict[str, float] = {
    "S_ARG":  0.253,
    "S_BMG":  0.237,
    "S_MOB":  0.201,
    "S_SIZE": 0.175,
    "S_VF":   0.135,
}

WEIGHT_SUM = sum(RISK_WEIGHTS.values())

# Risk grade thresholds (score, grade, label); manuscript Table 5
RISK_GRADES = [
    (0.60, "A", "Very High"),
    (0.45, "B", "High"),
    (0.30, "C", "Moderate"),
    (0.15, "D", "Low"),
    (0.00, "E", "Minimal"),
]

# =====================================================================
# WHO AWaRe classification patterns (legacy mapping used for the frozen
# component matrix; gene-name based)
# =====================================================================
AWARE_ACCESS = re.compile(
    r"^(bla)?(TEM|SHV|OXA-1|blaZ|LEN|LAP)"
    r"|cat[ABC]?|cmlA|floR"
    r"|tet[ABCDEFGHIJKLMYZ]?|tet[S]?"
    r"|dfr[A-Z]?|sul[123]?"
    r"|aac(3|6\x27)?-?(I[ab]?|II|III|VI|APH|ANT)"
    r"|str[AB]|ant(3|6)?|aph(3|6)?"
    r"|erm[A-F]?|mef[AB]?|msr[ABCD]?|mph[ABCDE]?"
    r"|lnu[A-F]?|vga[ABCDE]?|vgb[ABC]?"
    r"|fus[ABC]?|mec[ABC]?"
    r"|qnr[A-DS]?|qepA|oqxA?B?"
    r"|aac\(6\x27\)-Ib"
    r"|bla[A-Z]{3,}?",
    re.IGNORECASE,
)

AWARE_RESERVE = re.compile(
    r"mcr-?[0-9]|tetX[0-9]?|cfr|optrA|poxtA|van[A-Z]?"
    r"|rmt[A-H]?|armA|npmA"
    r"|fos[A-Z]?|blaKPC|blaNDM|blaVIM|blaIMP|blaGES"
    r"|blaOXA-?(23|24|48|58)|blaSME|blaIMI|blaSPM|blaSIM|blaGIM"
    r"|vat[A-E]?|vgb[ABC]?",
    re.IGNORECASE,
)

# =====================================================================
# WHO AWaRe 2025 category table (B09489; sensitivity option).
# Returns 0 = not an antibiotic (scored under S_BMG), 1 = Access,
# 2 = Watch, 3 = Reserve.
# =====================================================================
A25_ACCESS = re.compile(
    r"^(bla)?(TEM|SHV|OXA-1|blaZ|LEN|LAP)"
    r"|^cat[ABC]?[0-9]*$|^cmlA|^floR"
    r"|^tet[A-M][0-9]*$|^dfr[A-Z]?[0-9]*$|^sul[123]?[0-9]*$"
    r"|^aac3|^aac6\x27?(?!.*cr)|^ant(3|6)|^aph(3|6)"
    r"|^str[AB][0-9]*$|^blaCM[YT]|^blaROB|^blaBRO|^blaOXA[2-9]?(?!48|23|24|58)",
    re.IGNORECASE,
)
A25_RESERVE = re.compile(
    r"mcr-?[0-9]|tetX[0-9]?|^cfr|optrA|poxtA"
    r"|rmt[A-H]?|armA|npmA|^fos[A-Z]?[0-9]*$"
    r"|^vat[A-E]?|^vgb[ABC]?|linezolid|^optr|^cfr[ABC]?",
    re.IGNORECASE,
)
A25_WATCH = re.compile(
    r"blaKPC|blaNDM|blaVIM|blaIMP|blaGES|blaSME|blaIMI|blaSPM|blaSIM|blaGIM"
    r"|blaOXA-?(23|24|48|58)|^van[A-Z]|^erm[A-F]?[0-9]*$|^mef[AB]?"
    r"|^msr[ABCD]?|^mph[A-E]?|^lnu[A-F]?|^vga[ABCDE]?|^fus[ABC]?"
    r"|^mec[ABC]?|qnr[A-DS]?|qepA|^oqxA?B?|aac\(6\x27\)-Ib-cr"
    r"|^bla(CTX|CMY|DHA|FOX|MOX|ACC|ACT|MIR|CFM|CRO|CAZ|CTX|FEP|CPO)"
    r"|^aac(?!3)|^aad|^aph(?!3)|^ant(?!3)|^sxt|^bla[A-Z]{3,}",
    re.IGNORECASE,
)
A25_NON_ANTIBIOTIC = re.compile(
    r"qac|^mer[A-Z]?[0-9]*$|^ars[A-Z]?[0-9]*$|^cop[A-Z]?[0-9]*$"
    r"|^sil[A-Z]?[0-9]*$|^czc[A-Z]?[0-9]*$|^cad[A-Z]?[0-9]*$"
    r"|^pco[A-Z]?[0-9]*$|^ter[A-Z]?[0-9]*$|^zin[T]?"
    r"|^bfp|^iut|^iro|^fim|^pil|^pap|^fep|^ent|^ybt|^irp"
    r"|^clp|^hly|^stx|^ipa|^vir|^esp|^eae|^tir|^lee|^chu|^sit",
    re.IGNORECASE,
)


def _aware_category_2025(gene_name: str) -> int:
    """Classify one gene name under the WHO AWaRe 2025 table.

    Returns 0 (not an antibiotic), 1 (Access), 2 (Watch), 3 (Reserve).
    Non-antibiotic resistance genes (biocide/metal) are scored under
    S_BMG and contribute 0 to S_ARG under the 2025 mapping.
    """
    g = gene_name.strip()
    if not g:
        return 0
    if A25_NON_ANTIBIOTIC.search(g):
        return 0
    if A25_RESERVE.search(g):
        return 3
    if A25_ACCESS.search(g):
        return 1
    if A25_WATCH.search(g):
        return 2
    # Unrecognized but plausible clinical resistance gene: default Watch
    if re.search(r"bla|aac|aph|ant|aad|erm|tet|dfr|sul|cat|cml|flo|mcr|van"
                 r"|rRNA|methyl|oxa|carbapenem", g, re.IGNORECASE):
        return 2
    return 0


# =====================================================================
# High-risk gene patterns (manuscript Method S2; hazard multipliers)
# =====================================================================
HIGH_RISK_PATTERNS = re.compile(
    r"mcr-?[0-9]|blaKPC|blaNDM|blaVIM|blaIMP|blaGES"
    r"|blaOXA-?(23|24|48|58)|blaSME|blaIMI|blaSPM|blaSIM|blaGIM"
    r"|tetX[0-9]?|^cfr|optrA|poxtA|van[A-Z]?"
    r"|qnr[A-DS]?|aac\(6\x27\)-Ib-cr|qepA|oqxA?B?"
    r"|rmt[A-H]?|armA|npmA|^fos[A-Z]?|^CTX-M|mec[ABC]?",
    re.IGNORECASE,
)

# Housekeeping / chromosomal genes excluded from FASTA-only ARG hazard
HOUSEKEEPING_GENES = re.compile(
    r"^(omp|por|fhu|sdh|mdh|gyr|par|rpo|rpl|rps|inf|dna|rec|uvr)"
    r"|AcrAB-TolC|MexAB-OprM|AdeFGH|EmhABC|MdtEF|RosAB|Fsr|KdpE|TelAB"
    r"|MFS|SMC|EmrE|QacE$|sugE|NorM|DinF|H-NS|arcA|marA|soxS|ramA",
    re.IGNORECASE,
)

# =====================================================================
# Mobility element patterns
# =====================================================================
T4CP_PATTERNS = re.compile(
    r"T4CP|T4CP\.|coupling protein|VirD4|TraD|Tcp"
    r"|relaxosome coupling|Tic?P|Tra[GS]|Cpl?",
    re.IGNORECASE,
)
RELAXASE_PATTERNS = re.compile(
    r"relaxase|helicase|TraI|VirD2|Mob[ABCDEFGH]?|NikB"
    r"|Tral|TraI|relaxase|MobF|MobC|TraW|TraH|TraJ|Sog",
    re.IGNORECASE,
)
ORIT_PATTERNS = re.compile(
    r"oriT|origin of transfer|nic site|bom|rlx",
    re.IGNORECASE,
)
MPF_PATTERNS = re.compile(
    r"T4SS|type IV secretion|VirB[1-9]?|Tra[ABCDEFGHKLMNUVWXY]?"
    r"|Trb[A-Z]?|T4[PB]|Pil[LTWXYZ]?|Cpa|Cag|Dot/Icm|VirD4"
    r"|mpf|T2SS|conjugative|mating pair",
    re.IGNORECASE,
)

# =====================================================================
# VF category weights
# =====================================================================
VF_CATEGORY_WEIGHTS = {
    "Exotoxin": 0.45,
    "Type III secretion system": 0.40,
    "T3SS effector": 0.35,
    "Type IV secretion system": 0.30,
    "Adherence": 0.25,
    "Invasion": 0.30,
    "Serum resistance": 0.30,
    "Nutritional/Metabolic factor": 0.20,
    "Immune evasion": 0.30,
    "Biofilm": 0.20,
    "Capsule": 0.20,
    "Protease": 0.20,
    "Toxin": 0.40,
    "Hemolysin": 0.35,
}

# =====================================================================
# Biocide/metal family bonuses
# =====================================================================
BMG_FAMILY_BONUS = {
    "mer": 0.35,   # mercury operon
    "qac": 0.25,   # quaternary ammonium compounds
    "ars": 0.25,   # arsenic
    "cop": 0.20,   # copper
    "sil": 0.20,   # silver
    "pco": 0.20,   # copper plasmid-borne
    "czc": 0.20,   # cadmium/zinc/cobalt
    "cad": 0.15,
    "ter": 0.15,   # tellurium
    "zin": 0.10,
}

ARG_HAZARD_CAP = 30.0
BMG_HAZARD_CAP = 20.0


@dataclass
class PlasmidFeatures:
    """Container for plasmid annotation features."""
    seq_id: str
    length_bp: int
    arg_names: List[str] = field(default_factory=list)
    vf_names: List[str] = field(default_factory=list)
    vf_categories: List[str] = field(default_factory=list)
    bm_gene_names: List[str] = field(default_factory=list)
    replicon: str = ""
    host_genera: Optional[List[str]] = None
    countries: Optional[List[str]] = None
    habitats: Optional[List[str]] = None
    collection_year: Optional[int] = None
    has_t4cp: bool = False
    has_relaxase: bool = False
    has_oriT: bool = False
    has_auxiliary: bool = False


class PlasRiskScorer:
    """Compute PlasRisk component scores and the composite risk score."""

    def __init__(
        self,
        mode: str = "lite",
        replicon_lookup: Optional[Dict] = None,
        aware_version: str = "legacy",
    ):
        if mode not in ("lite", "full"):
            raise ValueError("mode must be 'lite' or 'full' for PlasRiskScorer")
        if aware_version not in ("legacy", "2025"):
            raise ValueError("aware_version must be 'legacy' or '2025'")
        self.mode = mode
        self.aware_version = aware_version
        self.weights = RISK_WEIGHTS if mode == "full" else RISK_WEIGHTS_LITE
        self.replicon_lookup = replicon_lookup or {}

    # -----------------------------------------------------------------
    # S_ARG
    # -----------------------------------------------------------------
    def classify_aware(self, gene_name: str) -> int:
        """Legacy AWaRe classification: 1 Access, 2 Watch, 3 Reserve."""
        g = gene_name.strip()
        if AWARE_RESERVE.search(g):
            return 3
        if AWARE_ACCESS.search(g):
            return 1
        return 2

    def score_arg(self, feat: PlasmidFeatures) -> float:
        """WHO AWaRe-weighted ARG burden with high-risk multipliers."""
        genes = [g for g in feat.arg_names if not HOUSEKEEPING_GENES.search(g)]
        if not genes:
            return 0.0
        hazard = 0.0
        n_highrisk = 0
        for g in genes:
            if self.aware_version == "2025":
                cat = _aware_category_2025(g)
            else:
                cat = self.classify_aware(g)
            hazard += float(cat)
            if HIGH_RISK_PATTERNS.search(g):
                n_highrisk += 1
        if n_highrisk:
            hazard *= 1.3 ** n_highrisk
        if self.aware_version == "legacy":
            hazard = min(hazard, ARG_HAZARD_CAP)
            s = math.log10(1 + hazard) / math.log10(1 + ARG_HAZARD_CAP)
        else:
            # 2025: non-antibiotic genes already contribute 0; Reserve is
            # the top category so no separate 1.5x last-resort multiplier
            hazard = min(hazard, ARG_HAZARD_CAP)
            s = math.log10(1 + hazard) / math.log10(1 + ARG_HAZARD_CAP)
        return float(min(1.0, s))

    # -----------------------------------------------------------------
    # S_VF
    # -----------------------------------------------------------------
    def score_vf(self, feat: PlasmidFeatures) -> float:
        if not feat.vf_names:
            return 0.0
        hazard = 0.0
        for name in feat.vf_names:
            hazard += 0.3
        for cat in set(feat.vf_categories):
            for key, w in VF_CATEGORY_WEIGHTS.items():
                if key.lower() in cat.lower():
                    hazard += w
                    break
        return float(min(1.0, hazard / 3.0))

    # -----------------------------------------------------------------
    # S_MOB
    # -----------------------------------------------------------------
    def score_mob(self, feat: PlasmidFeatures) -> float:
        score = 0.0
        if feat.has_relaxase:
            score += 0.35
        if feat.has_t4cp:
            score += 0.30
        if feat.has_oriT:
            score += 0.15
        if feat.has_auxiliary:
            score += 0.20
        return float(min(1.0, score))

    def score_mob_from_names(self, gene_names: List[str]) -> float:
        """Infer mobility score from annotated gene names."""
        text = " ".join(gene_names)
        score = 0.0
        if RELAXASE_PATTERNS.search(text):
            score += 0.35
        if T4CP_PATTERNS.search(text):
            score += 0.30
        if ORIT_PATTERNS.search(text):
            score += 0.15
        if MPF_PATTERNS.search(text):
            score += 0.20
        return float(min(1.0, score))

    # -----------------------------------------------------------------
    # S_HOST
    # -----------------------------------------------------------------
    def score_host(self, feat: PlasmidFeatures) -> float:
        if feat.host_genera:
            n = len(feat.host_genera)
            return float(min(1.0, n / 10.0))
        if feat.replicon and feat.replicon in self.replicon_lookup:
            info = self.replicon_lookup[feat.replicon]
            if isinstance(info, dict):
                return float(info.get("host_range_score", 0.5))
        return 0.5

    # -----------------------------------------------------------------
    # S_REP
    # -----------------------------------------------------------------
    def score_rep(self, feat: PlasmidFeatures) -> float:
        if feat.replicon and feat.replicon in self.replicon_lookup:
            info = self.replicon_lookup[feat.replicon]
            if isinstance(info, dict):
                return float(info.get("backbone_risk", 0.5))
        return 0.5

    # -----------------------------------------------------------------
    # S_SIZE
    # -----------------------------------------------------------------
    def score_size(self, length_bp: int) -> float:
        """Monotonic sigmoid on log10 length, midpoint 30 kb."""
        x = 4.0 * (math.log10(max(length_bp, 1)) - math.log10(30000))
        return 1.0 / (1.0 + math.exp(-x))

    @staticmethod
    def score_size_logratio(length_bp: int) -> float:
        """Discovery-pipeline log-ratio transform (weight-fitting analyses)."""
        s = math.log10(max(length_bp, 1)) / math.log10(1_000_000)
        return float(min(1.0, max(0.0, s)))

    # -----------------------------------------------------------------
    # S_BMG
    # -----------------------------------------------------------------
    def score_bmg(self, feat: PlasmidFeatures) -> float:
        if not feat.bm_gene_names:
            return 0.0
        hazard = 0.0
        families_hit = set()
        for g in feat.bm_gene_names:
            gl = g.lower()
            hazard += 0.15
            for fam, bonus in BMG_FAMILY_BONUS.items():
                if gl.startswith(fam) or fam in gl:
                    families_hit.add(fam)
        for fam in families_hit:
            hazard += BMG_FAMILY_BONUS[fam]
        return float(min(1.0, hazard / 2.0))

    # -----------------------------------------------------------------
    # S_GEO
    # -----------------------------------------------------------------
    def score_geo(self, feat: PlasmidFeatures) -> float:
        if feat.countries:
            return float(min(1.0, len(feat.countries) / 50.0))
        if feat.replicon and feat.replicon in self.replicon_lookup:
            info = self.replicon_lookup[feat.replicon]
            if isinstance(info, dict):
                return float(info.get("geo_spread_score", 0.1))
        return 0.1

    # -----------------------------------------------------------------
    # S_HAB
    # -----------------------------------------------------------------
    def score_hab(self, feat: PlasmidFeatures) -> float:
        if feat.habitats:
            return float(min(1.0, len(feat.habitats) / 3.0))
        if feat.replicon and feat.replicon in self.replicon_lookup:
            info = self.replicon_lookup[feat.replicon]
            if isinstance(info, dict):
                return float(info.get("habitat_breadth_score", 0.1))
        return 0.1

    # -----------------------------------------------------------------
    # S_GROW
    # -----------------------------------------------------------------
    def score_grow(self, feat: PlasmidFeatures) -> float:
        if feat.replicon and feat.replicon in self.replicon_lookup:
            info = self.replicon_lookup[feat.replicon]
            if isinstance(info, dict):
                return float(min(1.0, max(0.0,
                                          info.get("growth_rate_score", 0.1))))
        return 0.1

    # -----------------------------------------------------------------
    # Grade
    # -----------------------------------------------------------------
    def _grade(self, score: float):
        for threshold, grade, label in RISK_GRADES:
            if score >= threshold:
                return grade, label
        return "E", "Minimal"

    # -----------------------------------------------------------------
    # Full scoring
    # -----------------------------------------------------------------
    def score(self, feat: PlasmidFeatures) -> Dict:
        s_arg = self.score_arg(feat)
        s_vf = self.score_vf(feat)
        s_mob = self.score_mob(feat)
        s_size = self.score_size(feat.length_bp)
        s_bmg = self.score_bmg(feat)

        components = {
            "S_ARG": s_arg,
            "S_VF": s_vf,
            "S_MOB": s_mob,
            "S_SIZE": s_size,
            "S_BMG": s_bmg,
        }

        if self.mode == "full":
            components.update({
                "S_HOST": self.score_host(feat),
                "S_REP": self.score_rep(feat),
                "S_GEO": self.score_geo(feat),
                "S_HAB": self.score_hab(feat),
                "S_GROW": self.score_grow(feat),
            })
        else:
            components.update({
                "S_HOST": None, "S_REP": None, "S_GEO": None,
                "S_HAB": None, "S_GROW": None,
            })

        s_total = sum(self.weights[k] * components[k] for k in self.weights)
        s_norm = s_total / WEIGHT_SUM
        grade, grade_label = self._grade(s_norm)

        hr_genes = sorted({
            g for g in feat.arg_names if HIGH_RISK_PATTERNS.search(g)
        })

        return {
            "seq_id": feat.seq_id,
            "length_bp": feat.length_bp,
            "replicon": feat.replicon,
            **components,
            "S_total": round(s_total, 6),
            "S_norm": round(s_norm, 6),
            "grade": grade,
            "grade_label": grade_label,
            "high_risk_genes": ";".join(hr_genes),
            "n_ARG": len(feat.arg_names),
            "n_VF": len(feat.vf_names),
            "n_BM": len(feat.bm_gene_names),
            "model_mode": self.mode,
            "aware_version": self.aware_version,
        }

    def score_dataframe(self, features: List[PlasmidFeatures]) -> pd.DataFrame:
        rows = [self.score(f) for f in features]
        df = pd.DataFrame(rows)
        df = df.sort_values("S_norm", ascending=False).reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)
        return df


# =====================================================================
# Original PIPdb 8-item ordinal index (Zhu et al. 2025, Table 1)
# =====================================================================
PIPdb_BINS = {
    "pathogenic_phylum": [(1, 1), (5, 5)],
    "pathogenic_species": [(1, 1), (5, 5)],
    "habitats": [(0, 1), (1, 2), (2, 3), (3, 5)],
    "args": [(0, 1), (2, 2), (5, 3), (10, 4), (20, 5)],
    "vfgs": [(0, 1), (2, 2), (5, 3), (10, 4), (20, 5)],
    "who_args": [(0, 1), (1, 3), (3, 4), (5, 5)],
    "iss": [(0, 1), (1, 2), (3, 3), (5, 4), (10, 5)],
    "growth_rate": [(-0.1, 1), (0, 2), (0.1, 3), (0.3, 4), (0.5, 5)],
}


def _bin_value(value: float, bins: List) -> int:
    result = bins[0][1]
    for threshold, score in bins:
        if value >= threshold:
            result = score
    return result


class PIPdbScorer:
    """Reproduce the published PIPdb combined risk index (ordinal 1-5)."""

    def __init__(self, **kwargs):
        self.mode = "ordinal"

    def score(self, feat: PlasmidFeatures) -> Dict:
        pathogenic_phylum = 5 if feat.host_genera else 1
        pathogenic_species = 5 if feat.vf_names or feat.arg_names else 1
        habitats = _bin_value(len(feat.habitats or []), PIPdb_BINS["habitats"])
        args = _bin_value(len(feat.arg_names), PIPdb_BINS["args"])
        vfgs = _bin_value(len(feat.vf_names), PIPdb_BINS["vfgs"])
        who_args = _bin_value(
            len([g for g in feat.arg_names
                 if HIGH_RISK_PATTERNS.search(g)]),
            PIPdb_BINS["who_args"],
        )
        iss = _bin_value(0, PIPdb_BINS["iss"])
        growth = 1

        raw = (pathogenic_phylum + pathogenic_species + habitats
               + args + vfgs + 2 * who_args + iss + growth)
        index = int(min(5, max(1, round(raw / 8.0 + 0.6))))
        normalized = (index - 1) / 4.0

        return {
            "seq_id": feat.seq_id,
            "length_bp": feat.length_bp,
            "replicon": feat.replicon,
            "combined_risk_index": index,
            "risk_index_normalized": round(normalized, 6),
            "grade": str(index),
            "grade_label": "",
            "n_ARG": len(feat.arg_names),
            "n_VF": len(feat.vf_names),
            "n_BM": len(feat.bm_gene_names),
            "model_mode": "ordinal",
        }

    def score_dataframe(self, features: List[PlasmidFeatures]) -> pd.DataFrame:
        rows = [self.score(f) for f in features]
        df = pd.DataFrame(rows)
        df = df.sort_values("risk_index_normalized",
                            ascending=False).reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)
        return df


# =====================================================================
# S_SIZE rank concordance
# =====================================================================
def size_rank_concordance(lengths: List[int]) -> Dict:
    """Verify sigmoid and log-ratio S_SIZE rank plasmids identically."""
    s_sig = [PlasRiskScorer().score_size(int(L)) for L in lengths]
    s_log = [PlasRiskScorer.score_size_logratio(int(L)) for L in lengths]
    from scipy.stats import spearmanr, kendalltau
    rho = float(spearmanr(s_sig, s_log).statistic)
    tau = float(kendalltau(s_sig, s_log).statistic)
    return {
        "n": len(lengths),
        "spearman_rho": rho,
        "kendall_tau": tau,
        "min_length_bp": int(min(lengths)),
        "max_length_bp": int(max(lengths)),
    }


# =====================================================================
# Factory
# =====================================================================
def get_scorer(
    model: str = "plasrisk",
    mode: Optional[str] = None,
    replicon_lookup: Optional[Dict] = None,
    aware_version: str = "legacy",
):
    if model == "pipdb" or mode == "ordinal":
        return PIPdbScorer()
    effective_mode = mode or "lite"
    return PlasRiskScorer(
        mode=effective_mode,
        replicon_lookup=replicon_lookup,
        aware_version=aware_version,
    )
