"""Tests for the CONJScan/MacSyFinder integration.

The actual MacSyFinder/Prodigal execution is skipped when those binaries
are not installed; parsing, classification and union logic are tested with
static fixtures.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plasrisk.conjscan import (  # noqa: E402
    CLASS_CONJ,
    CLASS_LIKELY,
    CLASS_MOBILIZABLE,
    CLASS_PCONJ,
    annotate_conjugation,
    classify_model,
    conjscan_available,
    parse_systems,
    rank_of,
    union_mobility_class,
)
from plasrisk.scoring import PlasmidFeatures, PlasRiskScorer  # noqa: E402

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Minimal MacSyFinder best_solution.tsv shape (two replicons, one empty).
SYSTEMS_TSV = """# This file is part of MacSyFinder
# model_fqn	sys_id	replicon	hit_id	gene_name	hit_pos	sys_loci	sys_wholeness	sys_score
CONJScan/Plasmids/Conj_T4SS	Conj1	psc1	psc1_00001	relaxase	1	6	0.95	12.5
CONJScan/Plasmids/Conj_T4SS	Conj1	psc1	psc1_00002	T4CP	2	6	0.95	12.5
CONJScan/Plasmids/Conj_T4SS	Conj1	psc1	psc1_00003	MPF_T4SS	3	6	0.95	12.5
CONJScan/Plasmids/pConj_T4SS_typeG	PC1	psc2	psc2_00001	T4CP	1	4	0.80	8.1
CONJScan/Plasmids/pConj_T4SS_typeG	PC1	psc2	psc2_00002	MPF_T4SS	2	4	0.80	8.1
CONJScan/Plasmids/pMob_MobP	Mob1	psc3	psc3_00001	relaxase	1	2	0.75	5.0
"""


@pytest.fixture
def systems_tsv(tmp_path):
    path = tmp_path / "best_solution.tsv"
    path.write_text(SYSTEMS_TSV, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Model classification
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model,expected", [
    ("CONJScan/Plasmids/Conj_T4SS", CLASS_CONJ),
    ("CONJScan/Plasmids/Likely_Conj_T4SS", CLASS_LIKELY),
    ("CONJScan/Plasmids/pConj_T4SS_typeG", CLASS_PCONJ),
    ("CONJScan/Plasmids/pMob_MobP", CLASS_MOBILIZABLE),
    ("CONJScan/Plasmids/Mob_P", CLASS_MOBILIZABLE),
])
def test_classify_model(model, expected):
    assert classify_model(model) == expected


def test_union_mobility_class():
    # gene-only positive dominates
    assert union_mobility_class("conjugative_complete", "non-mobilizable") \
        == "conjugative_complete"
    # system-level pConj upgrades gene-only mobilizable
    assert union_mobility_class("mobilizable", "conjugative_putative") \
        == "conjugative_putative"
    # no upgrade when gene evidence is stronger
    assert union_mobility_class("conjugative_likely", "mobilizable") \
        == "conjugative_likely"
    # both negative
    assert union_mobility_class("non-mobilizable", "none") \
        == "non-mobilizable"


def test_rank_of():
    assert rank_of("conjugative_putative") > rank_of("mobilizable")
    assert rank_of("conjugative_likely") > rank_of("conjugative_putative")
    assert rank_of("conjugative_complete") > rank_of("conjugative_likely")


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------

def test_parse_systems(systems_tsv):
    results = parse_systems(str(systems_tsv), {"psc1", "psc2", "psc3"})
    assert set(results) == {"psc1", "psc2", "psc3"}
    assert results["psc1"].conjscan_class == CLASS_CONJ
    assert results["psc1"].is_conjugative
    assert results["psc1"].has_t4cp and results["psc1"].has_relaxase
    assert results["psc2"].conjscan_class == CLASS_PCONJ
    assert not results["psc2"].is_conjugative
    assert results["psc2"].has_t4cp
    assert results["psc3"].conjscan_class == CLASS_MOBILIZABLE
    assert results["psc3"].has_relaxase
    assert results["psc1"].plasrisk_class() == "conjugative_complete"
    assert results["psc2"].plasrisk_class() == "conjugative_putative"


def test_parse_missing_file(tmp_path):
    assert parse_systems(str(tmp_path / "nope.tsv"), {"x"}) == {}


# ---------------------------------------------------------------------------
# Scoring of putative conjugative class
# ---------------------------------------------------------------------------

def _mobility_score(mobility_class, integron=False, n_is=0, length=30000):
    scorer = PlasRiskScorer(mode="lite")
    feat = PlasmidFeatures(
        seq_id="x", length_bp=length,
        mobility_class=mobility_class,
        mobility_annotated=True,
        has_integron=integron, n_is=n_is,
    )
    return scorer.score_mob(feat)


def test_putative_base_score():
    complete = _mobility_score("conjugative_complete")
    likely = _mobility_score("conjugative_likely")
    putative = _mobility_score("conjugative_putative")
    mobilizable = _mobility_score("mobilizable")
    assert complete > likely > putative > mobilizable
    assert abs(putative - 0.7) < 1e-9


# ---------------------------------------------------------------------------
# End-to-end (skipped without macsyfinder/prodigal)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not conjscan_available(),
                    reason="macsyfinder/prodigal not installed")
def test_annotate_conjugation(tmp_path):
    fasta = tmp_path / "small.fasta"
    # 6 kb synthetic fragment is enough for prodigal to run
    fasta.write_text(">c1\n" + "ATG" * 2000 + "\n", encoding="utf-8")
    results = annotate_conjugation(str(fasta), topology="linear",
                                   threads=1)
    assert "c1" in results
    assert results["c1"].conjscan_class in (
        "none", "mobilizable", "putative_conjugative",
        "conjugative_likely", "conjugative_complete")
