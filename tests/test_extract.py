"""Tests for the raw-assembly prediction/extraction front end."""

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from plasrisk.extract import (
    ClassifiedContig,
    combine_classifications,
    decide_route,
    parse_generic_table,
    parse_genomad_summary,
    parse_plasflow_table,
    process_one,
)


def _write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


PLASFLOW_TSV = (
    "contig_name\tcontig_length\tlabel\tplasmid\tchromosome\tunclassified\n"
    "NODE_1\t700000\tchromosome\t0.05\t0.94\t0.01\n"
    "NODE_2\t55000\tplasmid\t0.91\t0.08\t0.01\n"
)

GENOMAD_TSV = (
    "seq_name\tlength\tplasmid_score\tchromosome_score\ttopology\n"
    "c1\t700000\t0.04\t0.95\tLinear\n"
    "c2\t55000\t0.92\t0.07\tCircular\n"
)


def test_parse_plasflow_table():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "pf.tsv")
        _write(p, PLASFLOW_TSV)
        rows = parse_plasflow_table(p)
        assert len(rows) == 2
        by_id = {r.seq_id: r for r in rows}
        assert by_id["NODE_1"].label == "chromosome"
        assert by_id["NODE_2"].label == "plasmid"
        assert by_id["NODE_2"].plasmid_prob == pytest.approx(0.91)


def test_parse_genomad_summary_file():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "asm_plasmid_summary.tsv")
        _write(p, GENOMAD_TSV)
        rows = parse_genomad_summary(p)
        by_id = {r.seq_id: r for r in rows}
        assert by_id["c2"].label == "plasmid"
        assert by_id["c2"].circular is True


def test_parse_genomad_summary_directory():
    with tempfile.TemporaryDirectory() as d:
        sub = os.path.join(d, "asm")
        os.makedirs(sub)
        _write(os.path.join(sub, "x_plasmid_summary.tsv"), GENOMAD_TSV)
        rows = parse_genomad_summary(d)
        assert len(rows) == 2


def test_parse_generic_table():
    text = ("seq_id\tlength\tlabel\tplasmid_score\n"
            "a\t100\tunclassified\t0.1\n"
            "b\t40000\tplasmid\t0.88\n")
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "g.tsv")
        _write(p, text)
        rows = parse_generic_table(p)
        by_id = {r.seq_id: r for r in rows}
        assert by_id["b"].label == "plasmid"
        assert by_id["b"].length_bp == 40000


def test_combine_majority_agreement():
    a = [ClassifiedContig("x", 50000, "plasmid", 0.9),
         ClassifiedContig("y", 700000, "chromosome", 0.05)]
    b = [ClassifiedContig("x", 50000, "plasmid", 0.8),
         ClassifiedContig("y", 700000, "chromosome", 0.06)]
    out = {c.seq_id: c for c in combine_classifications([a, b], "majority")}
    assert out["x"].label == "plasmid"
    assert out["y"].label == "chromosome"


def test_combine_majority_disagreement_is_unclassified():
    a = [ClassifiedContig("x", 50000, "plasmid", 0.9)]
    b = [ClassifiedContig("x", 50000, "chromosome", 0.2)]
    out = combine_classifications([a, b], "majority")
    assert out[0].label == "unclassified"


def test_combine_any_high_recall():
    a = [ClassifiedContig("x", 50000, "unclassified", 0.4)]
    b = [ClassifiedContig("x", 50000, "plasmid", 0.6)]
    out = combine_classifications([a, b], "any")
    assert out[0].label == "plasmid"


def test_decide_route():
    circular = [ClassifiedContig("p", 55000, "plasmid", circular=True)]
    assert decide_route(circular) == "isolate"
    frags = [ClassifiedContig("f%d" % i, 1500, "plasmid")
             for i in range(5)]
    assert decide_route(frags) == "meta"


def _fasta(chrom_len, plasmid_specs):
    parts = [">chrom\n" + "A" * chrom_len]
    for name, length, circular in plasmid_specs:
        tag = " circular=true" if circular else ""
        parts.append(">%s%s\n" % (name, tag) + "C" * length)
    return "\n".join(parts)


def _generic_table_for(specs):
    lines = ["seq_id\tlength\tlabel\tplasmid_score"]
    lines.append("chrom\t%d\tchromosome\t0.02" % specs["chrom"])
    for name, length, label, score in specs["rows"]:
        lines.append("%s\t%d\t%s\t%s" % (name, length, label, score))
    return "\n".join("\n".join(lines).splitlines()) + "\n"


def test_process_chromosome_only_is_unscorable():
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "asm.fasta")
        table = os.path.join(d, "cls.tsv")
        _write(fasta, _fasta(700000, []))
        _write(table, _generic_table_for({
            "chrom": 700000,
            "rows": [("chrom", 700000, "chromosome", "0.02")]}))
        res = process_one(fasta, d, generic_table=table, use_abricate=False)
        assert res.status == "unscorable"
        assert "chromosomal sequence only" in res.reason


def test_process_small_ambiguous_is_unscorable():
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "asm.fasta")
        table = os.path.join(d, "cls.tsv")
        _write(fasta, ">x\n" + "A" * 500)
        _write(table, ("seq_id\tlength\tlabel\tplasmid_score\n"
                       "x\t500\tunclassified\t0.3\n"))
        res = process_one(fasta, d, generic_table=table, use_abricate=False)
        assert res.status == "unscorable"
        assert "too small or fragmented" in res.reason


def test_process_longread_circular_plasmid_scored_isolate():
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "asm.fasta")
        table = os.path.join(d, "cls.tsv")
        _write(fasta, _fasta(1_000_000, [("plas", 55000, True)]))
        _write(table, _generic_table_for({
            "chrom": 1_000_000,
            "rows": [("chrom", 1_000_000, "chromosome", "0.02"),
                     ("plas", 55000, "plasmid", "0.9")]}))
        res = process_one(fasta, d, generic_table=table, use_abricate=False)
        assert res.status == "scored"
        assert res.route == "isolate"
        extracted = open(res.plasmid_fasta).read()
        assert "chrom" not in extracted
        assert "plas" in extracted


def test_process_shortread_fragmented_plasmid_scored_meta():
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "asm.fasta")
        table = os.path.join(d, "cls.tsv")
        specs = [("f%d" % i, 1500, "plasmid", "0.8") for i in range(5)]
        _write(fasta, _fasta(700000,
                             [("f%d" % i, 1500, False) for i in range(5)]))
        _write(table, _generic_table_for({
            "chrom": 700000,
            "rows": [("chrom", 700000, "chromosome", "0.02")] + specs}))
        res = process_one(fasta, d, generic_table=table, use_abricate=False)
        assert res.status == "scored"
        assert res.route == "meta"
