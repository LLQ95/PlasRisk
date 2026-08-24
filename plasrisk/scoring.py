"""
scoring.py - Ten-dimension plasmid risk scoring model.

Each component S_i is normalized to [0, 1].
Composite score:  S = sum(w_i * S_i), normalized by sum(w_i).

When sequence-derived dimensions cannot be computed (e.g., abricate not
available or no mobility genes detected), replicon-specific empirical
medians from 792,964 PIPdb PSCs are used as fallback priors.

Two scoring models are provided:
  - PlasRiskScorer: 10-dimension weighted continuous model (default)
  - PIPdbScorer: original PIPdb 8-item ordinal model (use get_scorer('pipdb'))
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

# WHO-priority ARG indicator
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
# Mobility gene markers
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
CARD_BM_PATTERN = re.compile(
    r"qac[A-Z]?|mer[A-Z]?|ars[A-Z]?|cop[A-Z]?|sil[A-Z]?|czc[A-Z]?|"
    r"cad[A-Z]?|pco[A-Z]|ter[A-Z]?|zin[TAB]|nim[A-Z]?|ncc[A-Z]?|"
    r"chr[A-Z]?|act[A-Z]?|mdt[A-Z]?|emr[A-Z]?|mex[A-Z]?|ade[A-Z]?|"
    r"acr[ABEF]|tolC|mdfA|ceo[A-Z]|amr[A-Z]|far[A-Z]|ros[AB]|"
    r"bmr|blt|nor[AM]|patAB|mep[A-Z]|sugE|betA|opr|mtrR|"
    r"biocide|disinfectant|metal|copper|silver|arsenic|cadmium|zinc",
    re.I,
)

EXOTOXIN_PATTERN = re.compile(r"exotoxin|toxin|enterotoxin|hemolysin|cytolysin", re.I)
EFFECTOR_PATTERN = re.compile(r"effector delivery|type III|type IV|T3SS|T4SS|secretion", re.I)

ARG_HAZARD_CAP = 30.0
MAX_PLASMID_LENGTH = 1_000_000

# PIPdb global medians (from 792,964 PSCs)
GLOBAL_MEDIANS = {
    "S_ARG":  0.000, "S_VF":   0.000, "S_MOB":  0.085, "S_HOST": 0.760,
    "S_REP":  0.300, "S_SIZE": 0.609, "S_BM":   0.000, "S_GEO":  0.000,
    "S_HAB":  0.000, "S_GROW": 0.000,
}


# ---------------------------------------------------------------------------
# Data class for plasmid features
# ---------------------------------------------------------------------------
@dataclass
class PlasmidFeatures:
    """Container for all features needed to score a single plasmid."""

    seq_id: str = ""
    length_bp: int = 0

    arg_names: List[str] = field(default_factory=list)
    vf_names: List[str] = field(default_factory=list)
    vf_categories: List[str] = field(default_factory=list)
    bm_gene_names: List[str] = field(default_factory=list)

    replicon: str = ""

    has_t4cp: Optional[bool] = None
    has_relaxase: Optional[bool] = None
    has_oriT: Optional[bool] = None
    has_auxiliary: Optional[bool] = None
    has_integron: Optional[bool] = None
    n_is: int = 0
    mobility_class: str = ""

    n_host_genera: Optional[int] = None
    n_countries: Optional[int] = None
    n_habitats: Optional[int] = None
    annual_growth_rate: Optional[float] = None

    # PIPdb original model fields
    n_pathogenic_phylum: Optional[int] = None
    n_pathogenic_species: Optional[int] = None

    s_rep_prior: Optional[float] = None
    s_geo_prior: Optional[float] = None
    s_hab_prior: Optional[float] = None
    s_grow_prior: Optional[float] = None
    s_host_prior: Optional[float] = None
    s_mob_prior: Optional[float] = None

    mobility_annotated: bool = False
    host_is_who_priority: bool = False
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
    def n_who_arg(self) -> int:
        """Count of ARGs on the WHO priority pathogen list."""
        return sum(1 for n in self.arg_names if WHO_ARG_PATTERNS.search(n))

    @property
    def has_high_risk_arg(self) -> bool:
        return len(self.high_risk_args) > 0


# ---------------------------------------------------------------------------
# PlasRisk 10-dimension weighted scorer
# ---------------------------------------------------------------------------
class PlasRiskScorer:
    """Compute PlasRisk scores for plasmids (10-dim full or 5-dim lite)."""

    def __init__(self, replicon_lookup: Optional[pd.DataFrame] = None,
                 mode: str = "full"):
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

        self.defaults = dict(GLOBAL_MEDIANS)

    @staticmethod
    def classify_aware(gene_name: str, product: str = "",
                       resistance: str = "") -> int:
        text = f"{gene_name},{product},{resistance}"
        if AWARE_RESERVE.search(text):
            return 3
        if AWARE_ACCESS.search(text):
            return 1
        return 2

    def score_arg(self, feat: PlasmidFeatures) -> float:
        arg_names = [g for g in feat.arg_names
                     if not HOUSEKEEPING_GENES.match(g.split("_")[0].split("-")[0])]
        n = len(arg_names)
        if n == 0:
            return 0.0

        hazard = 0.0
        n_highrisk = len(feat.high_risk_args)
        for name in arg_names:
            product = feat.metadata.get("arg_products", {}).get(name, "")
            resistance = feat.metadata.get("arg_resistance", {}).get(name, "")
            aware = self.classify_aware(name, product, resistance)
            weight = float(aware)
            if LAST_RESORT_PATTERN.search(f"{name},{product},{resistance}"):
                weight *= 1.5
            hazard += weight

        if n_highrisk > 0:
            hazard *= (1.3 ** n_highrisk)
        if feat.host_is_who_priority:
            hazard *= 1.2

        hazard = min(hazard, ARG_HAZARD_CAP)
        return float(math.log10(1 + hazard) / math.log10(1 + ARG_HAZARD_CAP))

    @staticmethod
    def score_vf(feat: PlasmidFeatures) -> float:
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
        has_any_evidence = (bool(feat.has_t4cp) or bool(feat.has_relaxase) or
                            bool(feat.has_oriT) or bool(feat.has_auxiliary) or
                            bool(feat.has_integron) or feat.n_is > 0 or
                            bool(feat.mobility_class))
        if not feat.mobility_annotated and not has_any_evidence:
            if feat.s_mob_prior is not None:
                return float(feat.s_mob_prior)
            return self._lookup_prior(feat.replicon, "S_MOB")

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
            has_t4cp = bool(feat.has_t4cp)
            has_rel = bool(feat.has_relaxase)
            has_oriT = bool(feat.has_oriT)

            if has_t4cp and has_rel and has_oriT:
                mob_base = 1.0
            elif has_t4cp and has_rel:
                mob_base = 0.85
            elif has_rel or has_oriT:
                mob_base = 0.6
            else:
                if feat.n_is > 0:
                    mob_base = 0.1
                elif feat.s_mob_prior is not None:
                    return float(feat.s_mob_prior)
                else:
                    return self._lookup_prior(feat.replicon, "S_MOB")

        integron_bonus = 0.10 if feat.has_integron else 0.0
        mob_component = min(mob_base + integron_bonus, 1.0)

        length_kb = feat.length_bp / 1000.0 if feat.length_bp > 0 else 1.0
        is_density = feat.n_is / length_kb if length_kb > 0 else 0.0
        is_bonus = 0.30 * min(is_density, 2.0) / 2.0

        return float(np.clip(mob_component + is_bonus, 0.0, 1.0))

    @staticmethod
    def _has_mobility_evidence(feat: PlasmidFeatures) -> bool:
        return any(v is not None for v in [
            feat.has_t4cp, feat.has_relaxase, feat.has_oriT,
            feat.has_auxiliary, feat.has_integron
        ]) or feat.n_is > 0

    @staticmethod
    def score_size(length_bp: int) -> float:
        if length_bp <= 0:
            return 0.0
        log_len = math.log10(length_bp)
        log_center = math.log10(30000)
        return float(1.0 / (1.0 + math.exp(-4.0 * (log_len - log_center))))

    @staticmethod
    def score_bm(feat: PlasmidFeatures) -> float:
        bm_names = list(feat.bm_gene_names)
        if not bm_names and feat.arg_names:
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
            return min(n / 45.0, 1.0)
        if feat.s_geo_prior is not None:
            return float(feat.s_geo_prior)
        return self._lookup_prior(feat.replicon, "S_GEO")

    def score_hab(self, feat: PlasmidFeatures) -> float:
        if feat.n_habitats is not None:
            n = feat.n_habitats
            return min(n / 8.0, 1.0)
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

    def _lookup_prior(self, replicon: str, column: str) -> float:
        default = self.defaults.get(column, 0.3)
        if not replicon or self.lookup is None:
            return default
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
        if not replicon or self.lookup is None:
            return None
        if replicon in self.lookup.index:
            val = self.lookup.loc[replicon, column]
            return float(val) if pd.notna(val) else None
        for idx in self.lookup.index:
            if replicon.startswith(idx) or idx.startswith(replicon):
                val = self.lookup.loc[idx, column]
                return float(val) if pd.notna(val) else None
        return None

    def score(self, feat: PlasmidFeatures) -> Dict:
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
            components = {k: (all_components[k] if k in self.weights else None)
                          for k in RISK_WEIGHTS}
        else:
            components = all_components

        s_total = sum(self.weights[k] * all_components[k] for k in self.active_dims)
        s_norm = s_total / self.weight_sum

        grade, grade_label = self._grade(s_norm)
        mob_class = feat.mobility_class or self._infer_mob_class(components["S_MOB"])

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
        rows = [self.score(f) for f in features]
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values("S_norm", ascending=False).reset_index(drop=True)
            df.insert(0, "rank", range(1, len(df) + 1))
        return df

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
        elif s_mob >= 0.60:
            return "conjugative_likely"
        elif s_mob >= 0.35:
            return "mobilizable"
        else:
            return "non-mobilizable"


# Backward-compatible alias
PlasRiskScorer.score_mob = PlasRiskScorer.score_mob


# ---------------------------------------------------------------------------
# PIPdb original ordinal scoring model
# ---------------------------------------------------------------------------
# Reference: PIPdb: a comprehensive plasmid sequence resource for tracking
# the horizontal transfer of pathogenic factors and antimicrobial resistance
# genes. Table 1: Risk scoring system for PSCs.
#
# Combined risk index:
#   Round((Pathogenic_phylum_score + Pathogenic_species_score + Habitats_score
#          + ARGs_score + VFGs_score + 2*WHO_ARGs_score + ISs_score
#          + Annual_average_growth_rate) / 8 + 0.6)
# ---------------------------------------------------------------------------

def _ordinal_score(value: float, bins: List[Tuple[int, float, float]]) -> int:
    """Assign an ordinal score based on bin thresholds (high to low)."""
    for score, lower, upper in bins:
        if lower <= value < upper:
            return score
    return bins[-1][0]


PIPDB_BINS = {
    "pathogenic_phylum": [
        (5, 5, float("inf")), (4, 4, 5), (3, 3, 4), (2, 2, 3), (1, 1, 2),
    ],
    "pathogenic_species": [
        (5, 8, float("inf")), (4, 6, 8), (3, 4, 6), (2, 2, 4), (1, 1, 2),
    ],
    "habitats": [
        (5, 8, float("inf")), (4, 6, 8), (3, 4, 6), (2, 2, 4), (1, 1, 2),
    ],
    "args": [
        (5, 10, float("inf")), (4, 5, 10), (3, 2, 5), (2, 1, 2), (1, 0, 1),
    ],
    "vfgs": [
        (5, 6, float("inf")), (4, 4, 6), (3, 2, 4), (2, 1, 2), (1, 0, 1),
    ],
    "who_args": [
        (5, 3, float("inf")), (4, 2, 3), (3, 1, 2), (1, 0, 1),
    ],
    "iss": [
        (5, 30, float("inf")), (4, 15, 30), (3, 5, 15), (2, 2, 5), (1, 1, 2),
    ],
    "growth_rate": [
        (5, 0.2, float("inf")), (4, 0.1, 0.2), (3, 0.05, 0.1),
        (2, 0.01, 0.05), (1, 0.0, 0.01),
    ],
}

PIPDB_GRADES = {
    5: ("5", "Very High"), 4: ("4", "High"), 3: ("3", "Moderate"),
    2: ("2", "Low"), 1: ("1", "Minimal"),
}


class PIPdbScorer:
    """
    Original PIPdb ordinal risk scoring model.

    Implements the 8-item additive scoring system from the PIPdb paper
    (Table 1): Pathogenic phylum, Pathogenic species, Habitats, ARGs,
    VFGs, WHO ARGs (2x weight), Insertion sequences, and Annual average
    growth rate. Each item is binned into an ordinal 1-5 score, then:

        combined_risk_index = Round((sum + 2*WHO_ARGs) / 8 + 0.6)

    Use as an alternative to PlasRiskScorer for comparison or when
    reproducing the original PIPdb results.

    Notes
    -----
    - Pathogenic phylum/species counts require host metadata not available
      from sequence alone; if not provided, they default to 1 (minimum).
    - Habitats and growth rate similarly default to 1 if unknown.
    - ARG/VFG/WHO-ARG/IS counts can be derived from abricate annotations.
    """

    def __init__(self):
        self.model_name = "pipdb"

    @staticmethod
    def score_pathogenic_phylum(n: Optional[int]) -> int:
        if n is None:
            return 1
        return _ordinal_score(float(n), PIPDB_BINS["pathogenic_phylum"])

    @staticmethod
    def score_pathogenic_species(n: Optional[int]) -> int:
        if n is None:
            return 1
        return _ordinal_score(float(n), PIPDB_BINS["pathogenic_species"])

    @staticmethod
    def score_habitats(n: Optional[int]) -> int:
        if n is None:
            return 1
        return _ordinal_score(float(n), PIPDB_BINS["habitats"])

    @staticmethod
    def score_args(n: int) -> int:
        return _ordinal_score(float(n), PIPDB_BINS["args"])

    @staticmethod
    def score_vfgs(n: int) -> int:
        return _ordinal_score(float(n), PIPDB_BINS["vfgs"])

    @staticmethod
    def score_who_args(n: int) -> int:
        return _ordinal_score(float(n), PIPDB_BINS["who_args"])

    @staticmethod
    def score_iss(n: int) -> int:
        return _ordinal_score(float(n), PIPDB_BINS["iss"])

    @staticmethod
    def score_growth_rate(rate: Optional[float]) -> int:
        if rate is None:
            return 1
        return _ordinal_score(float(max(rate, 0.0)), PIPDB_BINS["growth_rate"])

    def score(self, feat: PlasmidFeatures) -> Dict:
        """Compute the PIPdb combined risk index for one plasmid."""
        s_phylum = self.score_pathogenic_phylum(feat.n_pathogenic_phylum)
        s_species = self.score_pathogenic_species(feat.n_pathogenic_species)
        s_habitats = self.score_habitats(feat.n_habitats)
        s_args = self.score_args(feat.n_arg)
        s_vfgs = self.score_vfgs(feat.n_vf)
        s_who = self.score_who_args(feat.n_who_arg)
        s_iss = self.score_iss(feat.n_is)
        s_growth = self.score_growth_rate(feat.annual_growth_rate)

        raw_sum = (s_phylum + s_species + s_habitats + s_args +
                   s_vfgs + 2 * s_who + s_iss + s_growth)
        combined_index = int(round(raw_sum / 8.0 + 0.6))
        combined_index = max(1, min(5, combined_index))

        normalized = (combined_index - 1) / 4.0
        grade, grade_label = PIPDB_GRADES.get(combined_index, ("1", "Minimal"))

        return {
            "seq_id": feat.seq_id,
            "length_bp": feat.length_bp,
            "replicon": feat.replicon or "Unknown",
            "n_ARG": feat.n_arg,
            "n_VF": feat.n_vf,
            "n_IS": feat.n_is,
            "n_WHO_ARG": feat.n_who_arg,
            "score_phylum": s_phylum,
            "score_species": s_species,
            "score_habitats": s_habitats,
            "score_args": s_args,
            "score_vfgs": s_vfgs,
            "score_who_args": s_who,
            "score_iss": s_iss,
            "score_growth": s_growth,
            "combined_risk_index": combined_index,
            "risk_index_normalized": round(normalized, 4),
            "grade": grade,
            "grade_label": grade_label,
            "model_mode": "pipdb",
        }

    def score_dataframe(self, features: List[PlasmidFeatures]) -> pd.DataFrame:
        rows = [self.score(f) for f in features]
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(
                "combined_risk_index", ascending=False
            ).reset_index(drop=True)
            df.insert(0, "rank", range(1, len(df) + 1))
        return df


def get_scorer(model: str = "plasrisk", mode: str = "full",
               replicon_lookup: Optional[pd.DataFrame] = None):
    """
    Factory function to get a scorer by model name.

    Parameters
    ----------
    model : str
        ``"plasrisk"`` (default) - 10-dimension weighted continuous model.
        ``"pipdb"`` - original PIPdb 8-item ordinal model.
    mode : str
        For PlasRisk: ``"full"`` (10-dim) or ``"lite"`` (5-dim core).
        Ignored for PIPdb.
    replicon_lookup : pd.DataFrame, optional
        Replicon empirical prior table (PlasRisk only).

    Returns
    -------
    PlasRiskScorer or PIPdbScorer
    """
    model = model.lower().strip()
    if model in ("plasrisk", "plasrisk10", "10dim", "10-dim"):
        return PlasRiskScorer(replicon_lookup=replicon_lookup, mode=mode)
    elif model in ("pipdb", "ordinal", "pipdb-original"):
        return PIPdbScorer()
    else:
        raise ValueError(
            f"Unknown model '{model}'. Choose 'plasrisk' or 'pipdb'.")
