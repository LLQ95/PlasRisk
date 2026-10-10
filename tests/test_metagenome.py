"""Tests for the metagenomic three-layer extension (v1.5.0)."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from plasrisk.metagenome import (
    aggregate_bin,
    aggregate_sample,
    aggregate_bins,
    is_circular_header,
    load_bins_table,
    load_coverage_table,
    mob_evidence_score,
    score_contig,
)
from plasrisk.fragment_sim import run_simulation, tile_fragments
from plasrisk.scoring import PlasRiskScorer, PlasmidFeatures


def make_feat(seq_id, length, args=None, bms=None, vfs=None,
              t4cp=False, relaxase=False, oriT=False, integron=False,
              n_is=0, replicon="", mobility_class=None):
    return PlasmidFeatures(
        seq_id=seq_id, length_bp=length,
        arg_names=list(args or []), bm_gene_names=list(bms or []),
        vf_names=list(vfs or []), replicon=replicon,
        has_t4cp=t4cp, has_relaxase=relaxase, has_oriT=oriT,
        has_integron=integron, n_is=n_is,
        mobility_class=mobility_class,
        mobility_annotated=bool(t4cp or relaxase or oriT or integron
                                or n_is or mobility_class))


@pytest.fixture(scope="module")
def scorer():
    return PlasRiskScorer(replicon_lookup=None, mode="lite")


# --------------------------------------------------------------------------
# Layer 1
# --------------------------------------------------------------------------

def test_small_fragment_with_critical_arg_gets_cargo_score(scorer):
    feat = make_feat("node_1", 1500, args=["mcr-9.1"])
    row = score_contig(feat, scorer, complete_length=10000)
    assert row.S_SIZE is None            # fragment length never used
    assert row.S_lite5 is None           # no 5-dim grade on fragments
    assert row.S_ARG > 0.4               # one critical ARG gives strong signal
    assert row.S_cargo3 is not None and row.S_cargo3 > 0.15
    assert row.grade == ""
    assert "S_ARG" in row.assessable_dims


def test_circular_contig_gets_full_lite_score(scorer):
    feat = make_feat("contig_circ circular=true", 45000,
                     args=["blaNDM-1"], bms=["sul1"],
                     t4cp=True, relaxase=True, oriT=True, n_is=3)
    row = score_contig(feat, scorer, complete_length=10000)
    assert row.circular and row.complete
    assert row.S_SIZE is not None and row.S_MOB is not None
    assert row.S_lite5 is not None
    assert row.grade in {"A", "B", "C", "D", "E"}


def test_mob_evidence_none_on_bare_fragment(scorer):
    feat = make_feat("node_2", 2000)
    assert mob_evidence_score(feat, complete=False) is None
    # complete assembly-wide absence is informative
    assert mob_evidence_score(feat, complete=True) == pytest.approx(0.1)


# --------------------------------------------------------------------------
# Layer 2
# --------------------------------------------------------------------------

def test_bin_aggregation_unions_genes_and_lengths(scorer):
    f1 = make_feat("c1", 6000, args=["mcr-9.1", "sul1"], bms=["arsR"],
                   t4cp=True, relaxase=True)
    f2 = make_feat("c2", 5000, args=["mcr-9.1", "tet(A)"], bms=["copA"],
                   oriT=True, n_is=2)
    row = aggregate_bin("bin_1", [f1, f2], scorer, circular=False,
                        plasmid_score=None, coverage=None)
    assert row["total_length_bp"] == 11000          # lengths summed
    assert row["n_ARG"] == 3                        # mcr-9.1 deduped
    assert row["n_BM"] == 2
    assert row["S_lite"] is not None
    assert row["confidence"] == "probable"          # >= 10 kb aggregate
    assert row["grade"] in {"A", "B", "C", "D", "E"}


def test_circular_bin_is_complete(scorer):
    f = make_feat("c1", 3000)
    row = aggregate_bin("b", [f], scorer, circular=True,
                        plasmid_score=0.99, coverage=12.0)
    assert row["confidence"] == "complete"
    assert row["coverage"] == pytest.approx(12.0)


def test_aggregate_bins_singletons_and_coverage(scorer):
    feats = {
        "c1": make_feat("c1", 40000, args=["blaNDM-1"], t4cp=True,
                        relaxase=True, oriT=True),
        "c2": make_feat("c2", 3000),
    }
    bins = pd.DataFrame({"seq_id": ["c1"], "bin_id": ["binA"],
                         "circular": [True], "plasmid_score": [0.95]})
    cov = pd.DataFrame({"seq_id": ["c1", "c2"], "coverage": [10.0, 2.0]})
    contig_df = pd.DataFrame({"seq_id": ["c1", "c2"]})
    out = aggregate_bins(contig_df, feats, bins, cov, scorer)
    ids = set(out["bin_id"])
    assert "binA" in ids and "unbinned:c2" in ids
    bina = out[out["bin_id"] == "binA"].iloc[0]
    assert bina["coverage"] == pytest.approx(10.0)


def test_bins_and_coverage_parsers(tmp_path):
    bp = tmp_path / "bins.tsv"
    bp.write_text("contig\tcomponent\ttopology\n"
                  "c1\tbinX\tcircular\nc2\tbinX\tlinear\n",
                  encoding="utf-8")
    bins = load_bins_table(str(bp))
    assert list(bins["bin_id"]) == ["binX", "binX"]
    assert bins.loc[0, "circular"] and not bins.loc[1, "circular"]

    cp = tmp_path / "cov.tsv"
    cp.write_text("contigName\tsampleA\tsampleB\nc1\t5.0\t9.0\nc2\t1.0\t2.0\n",
                  encoding="utf-8")
    cov = load_coverage_table(str(cp))
    assert set(cov["seq_id"]) == {"c1", "c2"}
    cov2 = load_coverage_table(str(cp), coverage_column="sampleB")
    assert cov2.loc[cov2.seq_id == "c1", "coverage"].iloc[0] == 9.0


# --------------------------------------------------------------------------
# Layer 3
# --------------------------------------------------------------------------

def test_sample_summary_weighted_math(scorer):
    bin_df = pd.DataFrame([
        {"bin_id": "b1", "total_length_bp": 10000, "coverage": 10.0,
         "S_lite": 0.80, "grade": "A", "confidence": "complete",
         "n_ARG": 1, "n_BM": 1, "high_risk_genes": "mcr-9.1"},
        {"bin_id": "b2", "total_length_bp": 30000, "coverage": 1.0,
         "S_lite": 0.20, "grade": "D", "confidence": "fragmented",
         "n_ARG": 0, "n_BM": 0, "high_risk_genes": ""},
    ])
    contig_df = pd.DataFrame({"x": [1, 2]})
    row = aggregate_sample("s1", bin_df, contig_df)
    # weights: 10000*10 = 100000 vs 30000*1 = 30000
    expect = (100000 * 0.80 + 30000 * 0.20) / 130000
    assert row["CPRI_mean"] == round(expect, 4)
    assert row["grade_A"] == 1 and row["pct_grade_AB"] == 50.0
    assert row["critical_ARG_load"] == pytest.approx(100000.0)
    assert "not trained or validated" in row["note"]


# --------------------------------------------------------------------------
# Fragmentation simulation
# --------------------------------------------------------------------------

def test_tile_fragments_cover_sequence():
    rng = np.random.default_rng(0)
    seq = "ACGT" * 5000  # 20 kb
    pieces = tile_fragments(seq, 5000, rng)
    joined = "".join(p[2] for pieces in [pieces] for p in pieces)
    assert joined == seq.upper()
    assert all(p[1] - p[0] == len(p[2]) for p in pieces)


def test_simulation_plumbing_no_abricate(tmp_path):
    fp = tmp_path / "plasmids.fasta"
    # one 60 kb + one 8 kb synthetic plasmid
    fp.write_text(">p60kb\n%s\n>p8kb\n%s\n"
                  % ("A" * 60000, "C" * 8000), encoding="utf-8")
    outdir = str(tmp_path / "sim")
    elem, summary = run_simulation(
        str(fp), frag_lengths=(2000, 5000), n_replicates=2,
        min_fragment=1000, use_abricate=False, quiet=True)
    assert not elem.empty
    assert set(summary["frag_len_bp"]) == {2000, 5000}
    # oracle bin of size-only fragments must recover the intact S_SIZE
    p60 = elem[elem["plasmid_id"] == "p60kb"]
    assert (p60["delta_bin"].abs().max()) < 0.02
