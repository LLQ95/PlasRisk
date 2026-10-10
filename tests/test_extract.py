"""Tests for the raw-assembly prediction front end (plasrisk.extract)."""

import os

import pytest

from plasrisk.extract import (
    ClassifiedContig,
    combine_classifications,
    decide_route,
    parse_genomad_summary,
    parse_generic_table,
    parse_plasflow_table,
    process_one,
)


def _write(path, text):
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


def _fasta(path, records):
    lines = []
    for sid, seq in records.items():
        lines.append(">%s" % sid)
        lines.append(seq)
    _write(path, "\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Parser tests
# ---------------------------------------------------------------------------

def test_parse_plasflow_table(tmp_path):
    tsv = tmp_path / "pf.tsv"
    _write(tsv, "contig_id\tcontig_name\tcontig_length\tid\tlabel\t"
                "unclassified\tchromosome\tplasmid\n"
                "1\tchrA\t5000000\t0\tchromosome\t0.01\t0.98\t0.01\n"
                "2\tplasB\t45000\t1\tplasmid\t0.02\t0.05\t0.93\n")
    cc = parse_plasflow_table(str(tsv))
    assert len(cc) == 2
    assert cc[0].seq_id == "chrA"
    assert cc[0].label == "chromosome"
    assert cc[1].label == "plasmid"
    assert cc[1].plasmid_prob == pytest.approx(0.93)
    assert cc[1].length_bp == 45000
    assert cc[1].source == "plasflow"


def test_parse_genomad_summary(tmp_path):
    tsv = tmp_path / "s_plasmid_summary.tsv"
    _write(tsv, "seq_name\tlength\ttopology\tn_genes\tplasmid_score\t"
                "chromosome_score\n"
                "ctg1\t48000\tCircular\t40\t0.92\t0.05\n"
                "ctg2\t3000000\t-\t3000\t0.02\t0.97\n")
    cc = parse_genomad_summary(str(tmp_path))
    labels = {c.seq_id: c for c in cc}
    assert labels["ctg1"].label == "plasmid"
    assert labels["ctg1"].circular is True
    assert labels["ctg2"].label == "chromosome"


def test_parse_generic_table(tmp_path):
    tsv = tmp_path / "g.tsv"
    _write(tsv, "seq_id\tlength\tlabel\tplasmid_score\n"
                "a\t20000\tplasmid\t0.8\n"
                "b\t900000\tchromosome\t0.1\n")
    cc = parse_generic_table(str(tsv))
    assert cc[0].label == "plasmid"
    assert cc[1].label == "chromosome"
    assert cc[1].length_bp == 900000


# ---------------------------------------------------------------------------
# Combine and routing
# ---------------------------------------------------------------------------

def test_combine_majority_and_any():
    g1 = [ClassifiedContig("x", 1000, "plasmid", 0.9),
          ClassifiedContig("y", 1000, "plasmid", 0.9)]
    g2 = [ClassifiedContig("x", 1000, "chromosome", 0.1),
          ClassifiedContig("y", 1000, "plasmid", 0.8)]
    g3 = [ClassifiedContig("x", 1000, "chromosome", 0.1)]
    maj = {c.seq_id: c for c in combine_classifications([g1, g2, g3])}
    assert maj["x"].label == "chromosome"   # 2 of 3
    assert maj["y"].label == "plasmid"
    anyv = {c.seq_id: c for c in combine_classifications([g1, g2, g3],
                                                         vote="any")}
    assert anyv["x"].label == "plasmid"


def test_decide_route():
    circular = [ClassifiedContig("p", 60000, "plasmid", 0.9,
                                 circular=True)]
    fragments = [ClassifiedContig("f%d" % i, 1500, "plasmid", 0.6)
                 for i in range(6)]
    assert decide_route(circular) == "isolate"
    assert decide_route(fragments) == "meta"


# ---------------------------------------------------------------------------
# End-to-end process_one
# ---------------------------------------------------------------------------

def _generic_table(path, rows):
    out = ["seq_id\tlength\tlabel\tplasmid_score"]
    for sid, length, label, prob in rows:
        out.append("%s\t%d\t%s\t%s" % (sid, length, label, prob))
    _write(path, "\n".join(out) + "\n")


def test_chromosome_only_unscorable(tmp_path):
    fasta = tmp_path / "asm.fasta"
    _fasta(fasta, {"chr1": "C" * 600_000})
    table = tmp_path / "cls.tsv"
    _generic_table(str(table), [("chr1", 600_000, "chromosome", 0.05)])
    res = process_one(str(fasta), str(tmp_path / "out"),
                      generic_table=str(table), use_abricate=False)
    assert res.status == "unscorable"
    assert "chromosomal sequence only" in res.reason
    assert os.path.exists(res.output_files["unscorable"])
    assert os.path.exists(res.output_files["manifest"])


def test_small_ambiguous_unscorable(tmp_path):
    fasta = tmp_path / "tiny.fasta"
    _fasta(fasta, {"c1": "A" * 800, "c2": "G" * 1200})
    table = tmp_path / "cls.tsv"
    _generic_table(str(table), [
        ("c1", 800, "unclassified", 0.4),
        ("c2", 1200, "unclassified", 0.45)])
    res = process_one(str(fasta), str(tmp_path / "out"),
                      generic_table=str(table), use_abricate=False)
    assert res.status == "unscorable"
    assert "too small or fragmented" in res.reason


def test_longread_isolate_route(tmp_path):
    fasta = tmp_path / "lr.fasta"
    _fasta(fasta, {
        "chr1": "C" * 1_000_000,
        "plas1 circular=true": "A" * 60_000})
    table = tmp_path / "cls.tsv"
    _generic_table(str(table), [
        ("chr1", 1_000_000, "chromosome", 0.05),
        ("plas1 circular=true", 60_000, "plasmid", 0.95)])
    out = tmp_path / "out"
    res = process_one(str(fasta), str(out), generic_table=str(table),
                      use_abricate=False)
    assert res.status == "scored"
    assert res.route == "isolate"
    assert os.path.exists(res.plasmid_fasta)
    # Only the plasmid sequence is carried into the extracted FASTA.
    from plasrisk.annotate import parse_fasta
    extracted = parse_fasta(res.plasmid_fasta)
    assert list(extracted) == ["plas1"]
    assert os.path.exists(res.output_files["isolate_scores"])


def test_shortread_meta_route(tmp_path):
    fasta = tmp_path / "sr.fasta"
    recs = {"chr1": "C" * 800_000}
    rows = [("chr1", 800_000, "chromosome", 0.05)]
    for i in range(5):
        sid = "frag%d" % i
        recs[sid] = "A" * 1500
        rows.append((sid, 1500, "plasmid", 0.8))
    _fasta(fasta, recs)
    table = tmp_path / "cls.tsv"
    _generic_table(str(table), rows)
    out = tmp_path / "out"
    res = process_one(str(fasta), str(out), generic_table=str(table),
                      use_abricate=False)
    assert res.status == "scored"
    assert res.route == "meta"
    assert os.path.exists(res.output_files["bin_scores"])
    assert os.path.exists(res.output_files["sample_summary"])
