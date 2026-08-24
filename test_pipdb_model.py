#!/usr/bin/env python3
"""Tests for PIPdbScorer (original PIPdb ordinal model)."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plasrisk import PIPdbScorer, PlasmidFeatures, get_scorer, PlasRiskScorer


def test_basic_pipdb():
    """Low-risk plasmid should score index=2."""
    scorer = PIPdbScorer()
    f = PlasmidFeatures(seq_id="pBR322", length_bp=4361)
    r = scorer.score(f)
    assert r["combined_risk_index"] == 2
    assert r["grade"] == "2"
    assert r["risk_index_normalized"] == 0.25
    print(f"  pBR322: index={r['combined_risk_index']}, grade={r['grade']} OK")


def test_high_risk_pipdb():
    """High-risk plasmid with many ARGs/VFs/IS should score 5."""
    scorer = PIPdbScorer()
    f = PlasmidFeatures(
        seq_id="pNDM-KPC", length_bp=120000,
        arg_names=["blaNDM-1", "blaKPC-2", "blaCTX-M-15", "qnrS1",
                   "aac6-Ib", "tetA", "sul1", "catA1", "mphA", "ermB", "dfrA12"],
        vf_names=["iroN", "iutA", "hlyF", "sitA", "sitB", "ompA"],
        n_is=35,
        n_pathogenic_phylum=5, n_pathogenic_species=8, n_habitats=8,
        annual_growth_rate=0.25,
    )
    r = scorer.score(f)
    assert r["combined_risk_index"] == 5
    assert r["score_args"] == 5
    assert r["score_vfgs"] == 5
    assert r["score_who_args"] == 5
    assert r["score_iss"] == 5
    assert r["n_WHO_ARG"] >= 3
    print(f"  pNDM-KPC: index={r['combined_risk_index']}, grade={r['grade']} OK")


def test_who_double_weight():
    """WHO ARGs are double-weighted in the formula."""
    scorer = PIPdbScorer()
    f = PlasmidFeatures(seq_id="pMCR", length_bp=30000,
                        arg_names=["mcr-1", "blaNDM-1", "blaKPC-2"], n_is=5)
    r = scorer.score(f)
    assert r["score_who_args"] == 5
    assert r["n_WHO_ARG"] == 3
    print(f"  pMCR+NDM+KPC: WHO score={r['score_who_args']}, n_WHO={r['n_WHO_ARG']} OK")


def test_formula_min_max():
    """Verify formula at boundaries."""
    scorer = PIPdbScorer()
    f_min = PlasmidFeatures(seq_id="min")
    r_min = scorer.score(f_min)
    assert r_min["combined_risk_index"] == 2
    print(f"  Minimum: 9/8+0.6=1.725 -> round={r_min['combined_risk_index']} OK")

    f_max = PlasmidFeatures(
        seq_id="max",
        arg_names=[f"arg{i}" for i in range(15)],
        vf_names=[f"vf{i}" for i in range(8)],
        n_is=35, n_pathogenic_phylum=6, n_pathogenic_species=10,
        n_habitats=10, annual_growth_rate=0.3)
    r_max = scorer.score(f_max)
    assert r_max["combined_risk_index"] == 5
    print(f"  Maximum: 45/8+0.6=6.225 -> clamp={r_max['combined_risk_index']} OK")


def test_ordinal_bins():
    """Test ordinal bin thresholds."""
    scorer = PIPdbScorer()
    assert scorer.score_args(0) == 1
    assert scorer.score_args(1) == 2
    assert scorer.score_args(2) == 3
    assert scorer.score_args(4) == 3
    assert scorer.score_args(5) == 4
    assert scorer.score_args(9) == 4
    assert scorer.score_args(10) == 5
    assert scorer.score_vfgs(0) == 1
    assert scorer.score_vfgs(1) == 2
    assert scorer.score_vfgs(3) == 3
    assert scorer.score_vfgs(5) == 4
    assert scorer.score_vfgs(6) == 5
    assert scorer.score_iss(0) == 1
    assert scorer.score_iss(1) == 1
    assert scorer.score_iss(2) == 2
    assert scorer.score_iss(14) == 3
    assert scorer.score_iss(15) == 4
    assert scorer.score_iss(29) == 4
    assert scorer.score_iss(30) == 5
    assert scorer.score_growth_rate(0.0) == 1
    assert scorer.score_growth_rate(0.005) == 1
    assert scorer.score_growth_rate(0.01) == 2
    assert scorer.score_growth_rate(0.05) == 3
    assert scorer.score_growth_rate(0.1) == 4
    assert scorer.score_growth_rate(0.2) == 5
    assert scorer.score_pathogenic_phylum(None) == 1
    assert scorer.score_pathogenic_species(None) == 1
    assert scorer.score_habitats(None) == 1
    assert scorer.score_growth_rate(None) == 1
    print("  All ordinal bin thresholds OK")


def test_factory():
    """Test get_scorer factory function."""
    s1 = get_scorer("plasrisk")
    s2 = get_scorer("pipdb")
    s3 = get_scorer("plasrisk", mode="lite")
    assert isinstance(s1, PlasRiskScorer)
    assert isinstance(s2, PIPdbScorer)
    assert s1.mode == "full"
    assert s3.mode == "lite"
    print("  get_scorer factory OK")


def test_plasrisk_still_works():
    """Verify PlasRiskScorer is not broken by changes."""
    pr = PlasRiskScorer()
    f = PlasmidFeatures(seq_id="test", length_bp=30000,
                        arg_names=["blaNDM-1", "mcr-1"], n_is=10,
                        has_t4cp=True, has_relaxase=True)
    r = pr.score(f)
    assert r["model_mode"] == "full"
    assert 0 <= r["S_norm"] <= 1
    assert r["grade"] in ("A", "B", "C", "D", "E")
    print(f"  PlasRisk: S_norm={r['S_norm']:.3f}, grade={r['grade']} OK")


def test_score_dataframe():
    """Test batch scoring with PIPdbScorer."""
    scorer = PIPdbScorer()
    features = [
        PlasmidFeatures(seq_id="p1", length_bp=5000),
        PlasmidFeatures(seq_id="p2", length_bp=100000,
                        arg_names=["blaNDM-1", "mcr-1", "blaKPC-2"],
                        n_is=20, n_habitats=6, annual_growth_rate=0.15),
    ]
    df = scorer.score_dataframe(features)
    assert len(df) == 2
    assert "combined_risk_index" in df.columns
    assert "risk_index_normalized" in df.columns
    print(f"  Batch scoring: {len(df)} plasmids OK")


if __name__ == "__main__":
    print("=== PIPdbScorer Tests ===\n")
    tests = [
        test_basic_pipdb,
        test_high_risk_pipdb,
        test_who_double_weight,
        test_formula_min_max,
        test_ordinal_bins,
        test_factory,
        test_plasrisk_still_works,
        test_score_dataframe,
    ]
    passed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            print(f"  FAILED: {t.__name__}: {e}")
            import traceback
            traceback.print_exc()
    print(f"\n{passed}/{len(tests)} tests passed")
    sys.exit(0 if passed == len(tests) else 1)
