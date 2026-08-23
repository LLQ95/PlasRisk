"""
scoring.py - Ten-dimension plasmid risk scoring model.

Each component S_i is normalized to [0, 1].
Composite score:  S = sum(w_i * S_i), normalized by sum(w_i).

When sequence-derived dimensions cannot be computed (e.g., abricate not
available or no mobility genes detected), replicon-specific empirical
medians from 792,964 PIPdb PSCs are used as fallback priors.
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
# Lite mode: 5-dimension core (S_ARG + S_VF + S_MOB + S_SIZE + S_BM)
# Renormalized from the full consensus weights; captures >=99.9% of mean AUC
# (0.920 vs 0.920 for 10-dim at 3 decimal places).
# Derived from all-subsets dimensionality analysis (pipdb_20, 5-fold CV):
#   k=5 mean CV AUC = 0.9201 vs k=10 = 0.9203 (difference 0.0002, NS).
# Includes S_VF because it raises MDR-VF fusion AUC from 0.943 to 0.963.
# Overfitting analysis (pipdb_21) confirmed no train-test gap for either model.
# ---------------------------------------------------------------------------
LITE_DIMENSIONS = ("S_ARG", "S_VF", "S_MOB", "S_SIZE", "S_BM")
LITE_WEIGHTS_RAW = {k: RISK_WEIGHTS[k] for k in LITE_DIMENSIONS}
LITE_WEIGHT_SUM = sum(LITE_WEIGHTS_RAW.values())  # 0.9505
RISK_WEIGHTS_LITE: Dict[str, float] = {
    k: v / LITE_WEIGHT_SUM for k, v in LITE_WEIGHTS_RAW.items()
}
# S_ARG=0.2576, S_VF=0.1153, S_MOB=0.2147, S_SIZE=0.1902, S_BM=0.2222

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
# AWaRe classification mapping (WHO Access, Watch, Reserve groups)
# Used for S_ARG hazard weighting, matching pipdb_02_risk_score.py logic.
# Keywords in gene/product/resistance names map to AWaRe categories.
# ---------------------------------------------------------------------------
AWARE_ACCESS = re.compile(
    r"blaTEM|blaSHV|blaOXA-?1\b|ampC|cat|cml|floR|tet[A-M]|dfrA|sul[123]|"
    r"erm[A-T]|mef[A-E]|msr[A-E]|mph[A-E]|lnu[A-C]|vat|vga|vgb|"
    r"aac\(3\)|aac\(6'\)-I[abd]|ant\(3\)|aph\(3'\)-III|strA|strB|"
    r"blaCMY-?2\b|blaDHA-?1\b|fosA[0-9]?|qnrD",
    re.I,
)
AWARE_RESERVE = re.compile(
    r"mcr|NDM|KPC|OXA-?48|OXA-?181|OXA-?232|OXA-?244|VIM|IMP|GES|SPM|"
    r"GIM|SIM|DIM|AIM|SMB|KHM|FRI|TMB|FRI-|tet\(?X\)?|van[ABC]|cfr|optrA|poxtA|"
    r"rmt[A-H]|armA|npmA",
    re.I,
)
# Everything else clinically relevant defaults to Watch
AWARE_WATCH_KEYWORDS = re.compile(
    r"CTX-M|CMY|DHA|FOX|MOX|LAT|ACC|MIR|ACT|MOR|CFE|ESBL|carbapenem|"
    r"qnr[A-CS]|aac\(6'\)-Ib-cr|qepA|oqxAB|gyrA|parC|"
    r"blaPER|blaVEB|blaGES|blaBEL|blaTLA|blaSCO|blaBIC|blaIMI|blaSME|blaNMC|"
    r"blaZ|mecA|mecC|pbp2a|pbp2|mupA|fusB|ileS|"
    r"blaOXA-?(2|9|10|23|24|40|58|72|113|143|181|235|239|244|245|247|253|255|258|259|276|347|437|488|514|515|517|546|547|552|553|554|555|556|557|558|559|560|561|562|563|564|565|566|567|568|569|570|571|572|573|574|575|576|577|578|579|580|581|582|583|584|585|586|587|588|589|590|591|592|593|594|595|596|597|598|599|600|601|602|603|604|605|606|607|608|609|610|611|612|613|614|615|616|617|618|619|620|621|622|623|624|625|626|627|628|629|630|631|632|633|634|635|636|637|638|639|640|641|642|643|644|645|646|647|648|649|650|651|652|653|654|655|656|657|658|659|660|661|662|663|664|665|666|667|668|669|670|671|672|673|674|675|676|677|678|679|680|681|682|683|684|685|686|687|688|689|690|691|692|693|694|695|696|697|698|699|700|701|702|703|704|705|706|707|708|709|710|711|712|713|714|715|716|717|718|719|720|721|722|723|724|725|726|727|728|729|730|731|732|733|734|735|736|737|738|739|740|741|742|743|744|745|746|747|748|749|750|751|752|753|754|755|756|757|758|759|760|761|762|763|764|765|766|767|768|769|770|771|772|773|774|775|776|777|778|779|780|781|782|783|784|785|786|787|788|789|790|791|792|793|794|795|796|797|798|799|800|801|802|803|804|805|806|807|808|809|810|811|812|813|814|815|816|817|818|819|820|821|822|823|824|825|826|827|828|829|830|831|832|833|834|835|836|837|838|839|840|841|842|843|844|845|846|847|848|849|850|851|852|853|854|855|856|857|858|859|860|861|862|863|864|865|866|867|868|869|870|871|872|873|874|875|876|877|878|879|880|881|882|883|884|885|886|887|888|889|890|891|892|893|894|895|896|897|898|899|900)",
    re.I,
)

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

# Last-resort antibiotic classes (1.5x multiplier)
LAST_RESORT_PATTERN = re.compile(
    r"carbapenem|colistin|polymyxin|glycopeptide|vancomycin|tigecycline|"
    r"linezolid|ceftazidime-avibactam|meropenem-vaborbactam|NDM|KPC|"
    r"mcr|van[ABC]|tet\(?X\)?|cfr|optrA|poxtA|oxazolidinone|phenicol.*last",
    re.I,
)

# WHO critical-priority pathogen association (1.2x multiplier)
WHO_PATHOGEN_PATTERN = re.compile(
    r"escherichia|klebsiella|enterobacter|salmonella|acinetobacter|"
    r"pseudomonas|staphylococcus|enterococcus|neisseria|haemophilus",
    re.I,
)

# Chromosomal housekeeping genes that should NOT count as plasmid-borne ARGs
# These are multi-drug efflux pumps / porins that abricate CARD may annotate
# but are not horizontally transferred resistance determinants.
HOUSEKEEPING_GENES = re.compile(
    r"^(acr[ABEF]|tolC|mdfA|mdt[A-Z0-9]*|emr[A-Z0-9]*|msbA|acrD|acrEF|"
    r"ampC|ampH|ampG|nfxB|mex[A-Z0-9]*|oprM|oprJ|mexR|nalC|nalD|"
    r"ade[A-Z0-9]*|adeR|adeS|cmeABC|cmeR|ceoAB|opcR|amrAB|amrR|"
    r"farAB|rosAB|mtrRCDE|norM|norA|bmr|bmrA|bmr3|blt|bmrR|"
    r"patAB|pmrA|carO|oprD|omp[A-Z0-9]*|ompF|ompC|ompK3[56]|"
    r"gyrA|parC|parE|folA|folP|rpoB|rpsL|kasA|inhA|embB|pncA|"
    r"ald|aldH|gidB|tlyA|rrs|rrl|rpl[A-Z]|rps[A-Z]|rpm[A-Z])$",
    re.I,
)

# ---------------------------------------------------------------------------
# Mobility gene markers (keyword-based detection in annotation)
# ---------------------------------------------------------------------------
T4CP_PATTERN = re.compile(
    r"\b(traD|traG|virD4|trwB|T4CP|coupling)\b", re.I)
RELAXASE_PATTERN = re.compile(
    r"\b(traI|mobA|nikB|traA|relaxase|traH|traJ_rel|mobC|virD2)\b", re.I)
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
# CARD-based biocide/metal gene fallback (when BacMet database unavailable)
CARD_BM_PATTERN = re.compile(
    r"qac[A-Z]?|mer[A-Z]?|ars[A-Z]?|cop[A-Z]?|sil[A-Z]?|czc[A-Z]?|"
    r"cad[A-Z]?|pco[A-Z]|ter[A-Z]?|zin[TAB]|nim[A-Z]?|ncc[A-Z]?|"
    r"chr[A-Z]?|act[A-Z]?|mdt[A-Z]?|emr[A-Z]?|mex[A-Z]?|ade[A-Z]?|"
    r"acr[ABEF]|tolC|mdfA|ceo[A-Z]|amr[A-Z]|far[A-Z]|ros[AB]|"
    r"bmr|blt|nor[AM]|patAB|mep[A-Z]|sugE|betA|opr|mtrR|"
    r"biocide|disinfectant|metal|copper|silver|arsenic|cadmium|zinc",
    re.I,
)

# VF category bonuses
EXOTOXIN_PATTERN = re.compile(r"exotoxin|toxin|enterotoxin|hemolysin|cytolysin", re.I)
EFFECTOR_PATTERN = re.compile(r"effector delivery|type III|type IV|T3SS|T4SS|secretion", re.I)

# PIPdb 99th percentile ARG hazard cap (from pipdb_02 analysis)
ARG_HAZARD_CAP = 30.0
# Maximum plasmid length for S_SIZE normalization (1 Mb)
MAX_PLASMID_LENGTH = 1_000_000

# PIPdb global medians (used as defaults when replicon is unknown)
# Computed from 792,964 PSCs in tab_psc_final_scores.csv
GLOBAL_MEDIANS = {
    "S_ARG":  0.000,
    "S_VF":   0.000,
    "S_MOB":  0.085,
    "S_HOST": 0.760,
    "S_REP":  0.300,
    "S_SIZE": 0.609,
    "S_BM":   0.000,
    "S_GEO":  0.000,
    "S_HAB":  0.000,
    "S_GROW": 0.000,
}


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
    has_integron: Optional[bool] = None
    n_is: int = 0  # insertion sequence count
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
    s_mob_prior: Optional[float] = None

    # Whether mobility was directly assessed from sequence
    # If False and no mobility genes found, S_MOB will use replicon prior
    mobility_annotated: bool = False

    # WHO pathogen association (for S_ARG 1.2x bonus)
    host_is_who_priority: bool = False

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
    """Compute PlasRisk scores for plasmids (10-dim full or 5-dim lite)."""

    def __init__(self, replicon_lookup: Optional[pd.DataFrame] = None,
                 mode: str = "full"):
        """
        Parameters
        ----------
        replicon_lookup : pd.DataFrame, optional
            Lookup table with columns:
            replicon_primary, S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE,
            S_BM, S_GEO, S_HAB, S_GROW, n_PSC
            Empirical medians from 792,964 PIPdb PSCs.
            If None, built-in global medians are used.
        mode : str
            "full" (default): 10-dimension model.
            "lite": 5-dimension core (S_ARG, S_VF, S_MOB, S_SIZE, S_BM),
                    achieving equivalent mean AUC (0.920) to the full model.
                    Renormalized weights: 0.258/0.115/0.215/0.190/0.222.
        """
        if mode not in ("full", "lite"):
            raise ValueError("mode must be 'full' or 'lite', got %r" % mode)
        self.mode = mode
        self.weights = RISK_WEIGHTS_LITE if mode == "lite" else RISK_WEIGHTS
        self.weight_sum = sum(self.weights.values())
        self.active_dims = tuple(self.weights.keys())

        if replicon_lookup is not None:
            self.lookup = replicon_lookup.set_index("replicon_primary")
        else:
            self.lookup = None

        # Default priors for unknown replicons: PIPdb global medians
        self.defaults = dict(GLOBAL_MEDIANS)

    # ------------------------------------------------------------------
    # AWaRe classification helper
    # ------------------------------------------------------------------

    @staticmethod
    def classify_aware(gene_name: str, product: str = "",
                       resistance: str = "") -> int:
        """
        Classify an ARG into WHO AWaRe category.
        Returns 1 (Access), 2 (Watch), or 3 (Reserve).
        """
        text = f"{gene_name},{product},{resistance}"
        if AWARE_RESERVE.search(text):
            return 3
        if AWARE_ACCESS.search(text):
            return 1
        # Default to Watch for unrecognized but clinically relevant genes
        return 2

    # ------------------------------------------------------------------
    # Individual component scorers
    # ------------------------------------------------------------------

    def score_arg(self, feat: PlasmidFeatures) -> float:
        """
        S_ARG: AWaRe-weighted ARG hazard with log10 compression.

        Each ARG contributes:
          base weight by AWaRe class (Access=1, Watch=2, Reserve=3)
          x1.5 if last-resort antibiotic class
          x1.3 per high-risk ARG (carbapenemase, mcr, tetX, van, cfr, optrA)
          x1.2 if WHO priority pathogen host

        Raw hazard w is log10-compressed and normalized to 99th percentile:
          S_ARG = log10(1 + min(w, 30)) / log10(1 + 30)

        Housekeeping chromosomal genes (acrAB, tolC, etc.) are excluded.
        """
        # Filter out housekeeping genes
        arg_names = [g for g in feat.arg_names
                     if not HOUSEKEEPING_GENES.match(g.split("_")[0].split("-")[0])]
        n = len(arg_names)
        if n == 0:
            return 0.0

        hazard = 0.0
        n_highrisk = len(feat.high_risk_args)
        for name in arg_names:
            # Get product/resistance info from metadata if available
            product = feat.metadata.get("arg_products", {}).get(name, "")
            resistance = feat.metadata.get("arg_resistance", {}).get(name, "")
            aware = self.classify_aware(name, product, resistance)
            weight = float(aware)

            # Last-resort multiplier
            if LAST_RESORT_PATTERN.search(f"{name},{product},{resistance}"):
                weight *= 1.5

            hazard += weight

        # High-risk ARG multiplier (1.3x per high-risk gene)
        if n_highrisk > 0:
            hazard *= (1.3 ** n_highrisk)

        # WHO priority pathogen bonus
        if feat.host_is_who_priority:
            hazard *= 1.2

        # Log10 compression to 99th percentile cap
        hazard = min(hazard, ARG_HAZARD_CAP)
        return float(math.log10(1 + hazard) / math.log10(1 + ARG_HAZARD_CAP))

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

    def score_mob(self, feat: PlasmidFeatures) -> float:
        """
        S_MOB: mobility and MGE plasticity score.
        Base conjugative apparatus class + integron bonus,
        plus IS density bonus (up to +0.30).

        Class bases: conjugative_complete=1.0, conjugative_likely=0.85,
                     mobilizable=0.6, non-mobilizable=0.1.
        Integron bonus: +0.10 if integron present.
        IS bonus: 0.30 * min(n_is / (length_kb), 2.0) / 2.0.

        When mobility cannot be assessed from sequence (no mobility genes
        detected and no IS elements, or mobility databases not run), the
        replicon-specific empirical median S_MOB from PIPdb is used as a
        prior. This is critical for external validation where MOB-suite is
        not available and abricate CARD does not detect T4CP/relaxase.
        """
        # If mobility was not annotated AND no mobility evidence exists,
        # use replicon prior. But if T4CP/relaxase/oriT/IS are explicitly
        # set (even without abricate), use that evidence.
        has_any_evidence = (bool(feat.has_t4cp) or bool(feat.has_relaxase) or
                            bool(feat.has_oriT) or bool(feat.has_auxiliary) or
                            bool(feat.has_integron) or feat.n_is > 0 or
                            bool(feat.mobility_class))
        if not feat.mobility_annotated and not has_any_evidence:
            if feat.s_mob_prior is not None:
                return float(feat.s_mob_prior)
            return self._lookup_prior(feat.replicon, "S_MOB")

        # Determine mobility class base score
        if feat.mobility_class:
            mc = feat.mobility_class.lower()
            if "complete" in mc:
                mob_base = 1.0
            elif "likely" in mc:
                mob_base = 0.85
            elif "mobilizable" in mc:
                mob_base = 0.6
            else:
                mob_base = 0.1
        else:
            # Infer from element booleans
            has_t4cp = bool(feat.has_t4cp)
            has_rel = bool(feat.has_relaxase)
            has_oriT = bool(feat.has_oriT)

            if has_t4cp and has_rel and has_oriT:
                mob_base = 1.0  # conjugative_complete
            elif has_t4cp and has_rel:
                mob_base = 0.85  # conjugative_likely
            elif has_rel or has_oriT:
                mob_base = 0.6  # mobilizable
            else:
                # No mobility genes detected from sequence.
                # If ISfinder found IS elements, give partial credit for
                # MGE plasticity; otherwise fall back to replicon prior
                # rather than assuming non-mobilizable (abricate CARD/VFDB
                # do not detect T4CP/relaxase/oriT reliably).
                if feat.n_is > 0:
                    mob_base = 0.1  # non-mobilizable but has IS plasticity
                elif feat.s_mob_prior is not None:
                    return float(feat.s_mob_prior)
                else:
                    return self._lookup_prior(feat.replicon, "S_MOB")

        # Integron bonus
        integron_bonus = 0.10 if feat.has_integron else 0.0
        mob_component = min(mob_base + integron_bonus, 1.0)

        # IS density bonus (up to +0.30)
        length_kb = feat.length_bp / 1000.0 if feat.length_bp > 0 else 1.0
        is_density = feat.n_is / length_kb if length_kb > 0 else 0.0
        is_bonus = 0.30 * min(is_density, 2.0) / 2.0

        return float(np.clip(mob_component + is_bonus, 0.0, 1.0))

    @staticmethod
    def _has_mobility_evidence(feat: PlasmidFeatures) -> bool:
        """Check if any mobility-related annotation was attempted."""
        return any(v is not None for v in [
            feat.has_t4cp, feat.has_relaxase, feat.has_oriT,
            feat.has_auxiliary, feat.has_integron
        ]) or feat.n_is > 0

    @staticmethod
    def score_size(length_bp: int) -> float:
        """
        S_SIZE: sigmoidal plasmid length score.
        Uses a logistic function on log10(length), centered at 30 kb with
        steepness k=4. This gives:
          ~1 kb  -> ~0.05
          ~30 kb -> ~0.50
          ~120 kb-> ~0.92
          ~200 kb-> ~0.97
        """
        if length_bp <= 0:
            return 0.0
        log_len = math.log10(length_bp)
        log_center = math.log10(30000)  # 30 kb midpoint
        return float(1.0 / (1.0 + math.exp(-4.0 * (log_len - log_center))))

    @staticmethod
    def score_bm(feat: PlasmidFeatures) -> float:
        """
        S_BM: biocide/metal resistance (co-selection potential).
        0 if no BMRGs.
        0.25 base + min(n_bm*0.04, 0.35) + 0.15*mer + 0.15*qac + 0.10*ars/cop/sil.

        When BacMet annotation is unavailable (bm_gene_names empty), falls back
        to detecting biocide/metal genes from ARG names (CARD/ResFinder hits).
        """
        # Collect BM gene names; if empty, try CARD-based fallback
        bm_names = list(feat.bm_gene_names)
        if not bm_names and feat.arg_names:
            # Fallback: scan ARG names for biocide/metal resistance genes
            for name in feat.arg_names:
                if CARD_BM_PATTERN.search(name):
                    bm_names.append(name)

        n = len(bm_names)
        if n == 0:
            return 0.0
        base = 0.25
        count_bonus = min(n * 0.04, 0.35)
        bm_text = ",".join(bm_names)
        mer_bonus = 0.15 if MER_PATTERN.search(bm_text) else 0.0
        qac_bonus = 0.15 if QAC_PATTERN.search(bm_text) else 0.0
        ars_bonus = 0.10 if ARS_COP_PATTERN.search(bm_text) else 0.0
        return min(base + count_bonus + mer_bonus + qac_bonus + ars_bonus, 1.0)

    def score_host(self, feat: PlasmidFeatures) -> float:
        """
        S_HOST: host range breadth.
        If n_host_genera known, map to [0,1].
        Otherwise use replicon lookup prior (PIPdb empirical median).
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
            return float(np.clip(0.3 + g * 2.0, 0.0, 1.0))
        if feat.s_grow_prior is not None:
            return float(feat.s_grow_prior)
        return self._lookup_prior(feat.replicon, "S_GROW")

    # ------------------------------------------------------------------
    # Lookup helpers
    # ------------------------------------------------------------------

    def _lookup_prior(self, replicon: str, column: str) -> float:
        """Retrieve a prior from the replicon lookup table.

        Handles multi-replicon strings (e.g., "IncFII;IncFIA;IncR") by
        splitting on ';' and taking the maximum prior across all replicons.
        """
        default = self.defaults.get(column, 0.3)
        if not replicon or self.lookup is None:
            return default

        # Split multi-replicon strings and evaluate each
        rep_list = [r.strip() for r in re.split(r"[;,/]", str(replicon)) if r.strip()]
        if not rep_list:
            return default

        best_val = None
        for rep in rep_list:
            val = self._lookup_single(rep, column)
            if val is not None:
                if best_val is None or val > best_val:
                    best_val = val
        return float(best_val) if best_val is not None else default

    def _lookup_single(self, replicon: str, column: str) -> Optional[float]:
        """Lookup prior for a single replicon name (exact + prefix match)."""
        if not replicon or self.lookup is None:
            return None
        # Exact match
        if replicon in self.lookup.index:
            val = self.lookup.loc[replicon, column]
            return float(val) if pd.notna(val) else None
        # Prefix match (e.g., "IncFII" for "IncFII(K)")
        for idx in self.lookup.index:
            if replicon.startswith(idx) or idx.startswith(replicon):
                val = self.lookup.loc[idx, column]
                return float(val) if pd.notna(val) else None
        return None

    # ------------------------------------------------------------------
    # Main scoring
    # ------------------------------------------------------------------

    def score(self, feat: PlasmidFeatures) -> Dict:
        """
        Compute component scores and the composite.

        In full mode, all 10 components are computed.
        In lite mode, only the 5 core components (S_ARG, S_VF, S_MOB, S_SIZE, S_BM)
        are computed; other S_* fields are reported as None.

        When sequence-derived dimensions cannot be determined (abricate not
        available, databases missing, or no genes detected), replicon-specific
        empirical medians from PIPdb are used as fallback priors.

        Returns
        -------
        dict with keys: seq_id, length_bp, replicon, n_ARG, n_VF, n_BM,
                        S_ARG ... S_GROW, S_total, S_norm, grade, grade_label,
                        high_risk_genes, mobility_class, model_mode,
                        imputed_dimensions
        """
        all_components = {
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

        # Track which dimensions were imputed from replicon priors
        imputed = []
        if not feat.mobility_annotated or not self._has_mobility_evidence(feat):
            imputed.append("S_MOB")
        if feat.n_host_genera is None:
            imputed.append("S_HOST")
        if feat.n_countries is None:
            imputed.append("S_GEO")
        if feat.n_habitats is None:
            imputed.append("S_HAB")
        if feat.annual_growth_rate is None:
            imputed.append("S_GROW")
        if not feat.replicon:
            imputed.append("S_REP")

        if self.mode == "lite":
            # Report only active dimensions; others as None
            components = {k: (all_components[k] if k in self.weights else None)
                          for k in RISK_WEIGHTS}
        else:
            components = all_components

        s_total = sum(self.weights[k] * all_components[k] for k in self.active_dims)
        s_norm = s_total / self.weight_sum

        grade, grade_label = self._grade(s_norm)

        # Determine mobility class if not pre-set
        mob_class = feat.mobility_class or self._infer_mob_class(components["S_MOB"])

        # Count BM genes (including CARD fallback for reporting)
        bm_names = list(feat.bm_gene_names)
        bm_fallback_names = []
        if not bm_names and feat.arg_names:
            for name in feat.arg_names:
                if CARD_BM_PATTERN.search(name):
                    bm_fallback_names.append(name)
        n_bm_total = len(bm_names) + len(bm_fallback_names)

        result = {
            "seq_id": feat.seq_id,
            "length_bp": feat.length_bp,
            "length_kb": round(feat.length_bp / 1000.0, 2),
            "replicon": feat.replicon or "Unknown",
            "n_ARG": feat.n_arg,
            "n_VF": feat.n_vf,
            "n_BM": n_bm_total,
            "high_risk_genes": ";".join(feat.high_risk_args) if feat.high_risk_args else "",
            "mobility_class": mob_class,
            "model_mode": self.mode,
            "imputed_dimensions": ";".join(imputed) if imputed else "",
        }
        for k, v in components.items():
            result[k] = round(v, 4) if v is not None else None
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
        # Thresholds for additive formula (mob_base + IS bonus)
        if s_mob >= 0.85:
            return "conjugative_complete"
        elif s_mob >= 0.60:
            return "conjugative_likely"
        elif s_mob >= 0.35:
            return "mobilizable"
        else:
            return "non-mobilizable"


# Backward-compatible alias (tests reference score_mob)
PlasRiskScorer.score_mob = PlasRiskScorer.score_mob
