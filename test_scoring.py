"""
tests/test_scoring.py - Unit tests for PlasRisk scoring functions.
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from plasrisk.scoring import (
    PlasRiskScorer,
    PlasmidFeatures,
    RISK_WEIGHTS,
    WEIGHT_SUM,
)


def test_weight_sum():
    """Weights should sum to ~1.0."""
    assert abs(sum(RISK_WEIGHTS.values()) - 1.0) < 0.01
    assert abs(WEIGHT_SUM - 1.0) < 0.01


def test_empty_plasmid():
    """Plasmid with no genes should score very low."""
    scorer = PlasRiskScorer()
    feat = PlasmidFeatures(seq_id="test", length_bp=5000)
    result = scorer.score(feat)
    assert result["S_ARG"] == 0.0
    assert result["S_VF"] == 0.0
    assert result["S_BM"] == 0.0
    assert result["grade"] in ("D", "E")
    assert result["S_norm"] < 0.3


def test_high_risk_plasmid():
    """Plasmid with NDM, mcr, VFs, and conjugation should score high."""
    scorer = PlasRiskScorer()
    feat = PlasmidFeatures(
        seq_id="pHighRisk",
        length_bp=120000,
        arg_names=["NDM-1", "CTX-M-15", "mcr-1", "TEM-1", "AAC(6')-Ib-cr"],
        vf_names=["aerobactin", "iroN", "iutA"],
        vf_categories=["Exotoxin", "Nutritional/Metabolic factor"],
        bm_gene_names=["merA", "merR", "qacEdelta1", "arsC"],
        replicon="IncX3",
        has_t4cp=True,
        has_relaxase=True,
        has_oriT=True,
        has_auxiliary=True,
    )
    result = scorer.score(feat)
    assert result["S_ARG"] > 0.7
    assert result["S_VF"] > 0.5
    assert result["S_BM"] > 0.5
    assert result["S_MOB"] >= 0.9
    assert result["S_SIZE"] > 0.9
    assert result["grade"] in ("A", "B")
    assert result["S_norm"] > 0.4
    assert "blaNDM" in result["high_risk_genes"]
    assert "mcr" in result["high_risk_genes"]


def test_size_sigmoid():
    """S_SIZE should be monotonic and sigmoidal."""
    scorer = PlasRiskScorer()
    s_small = scorer.score_size(1000)
    s_med = scorer.score_size(30000)
    s_large = scorer.score_size(200000)
    assert s_small < s_med < s_large
    assert 0.0 < s_small < 0.2
    assert 0.4 < s_med < 0.6
    assert s_large > 0.95


def test_mobility_classes():
    """S_MOB should distinguish complete conjugative from non-mobilizable."""
    scorer = PlasRiskScorer()
    complete = PlasmidFeatures(
        length_bp=50000,
        has_t4cp=True, has_relaxase=True, has_oriT=True, has_auxiliary=True,
    )
    none = PlasmidFeatures(
        length_bp=5000,
        has_t4cp=False, has_relaxase=False, has_oriT=False, has_auxiliary=False,
    )
    assert scorer.score_mob(complete) > 0.9
    assert scorer.score_mob(none) < 0.15


def test_bm_scoring():
    """S_BM should detect mer and qac bonuses."""
    scorer = PlasRiskScorer()
    no_bm = PlasmidFeatures(length_bp=5000, bm_gene_names=[])
    mer_only = PlasmidFeatures(length_bp=5000, bm_gene_names=["merA", "merR", "merD"])
    assert scorer.score_bm(no_bm) == 0.0
    assert scorer.score_bm(mer_only) >= 0.5


def test_grade_thresholds():
    """Grade assignment should follow thresholds."""
    scorer = PlasRiskScorer()
    assert scorer._grade(0.70) == ("A", "Very High")
    assert scorer._grade(0.50) == ("B", "High")
    assert scorer._grade(0.35) == ("C", "Moderate")
    assert scorer._grade(0.20) == ("D", "Low")
    assert scorer._grade(0.05) == ("E", "Minimal")


def test_dataframe_output():
    """score_dataframe should return sorted DataFrame with rank."""
    import pandas as pd
    scorer = PlasRiskScorer()
    feats = [
        PlasmidFeatures(seq_id="low", length_bp=2000),
        PlasmidFeatures(seq_id="high", length_bp=100000,
                        arg_names=["NDM-1"], has_t4cp=True, has_relaxase=True),
    ]
    df = scorer.score_dataframe(feats)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert df.iloc[0]["seq_id"] == "high"  # highest risk first
    assert "rank" in df.columns


if __name__ == "__main__":
    test_weight_sum()
    test_empty_plasmid()
    test_high_risk_plasmid()
    test_size_sigmoid()
    test_mobility_classes()
    test_bm_scoring()
    test_grade_thresholds()
    test_dataframe_output()
    print("All tests passed.")
