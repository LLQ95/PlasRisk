"""Tests for the unified auto entry point and replicon classifier."""

import os

from plasrisk.extract import (
    characterize_assembly,
    choose_auto_strategy,
    classify_by_replicon,
    auto_main,
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
# Assembly characterisation and routing
# ---------------------------------------------------------------------------

def test_plasmid_collection_routes_direct(tmp_path):
    fasta = tmp_path / "plasmids.fasta"
    _fasta(fasta, {"p1": "A" * 20_000, "p2": "C" * 90_000})
    char = characterize_assembly(str(fasta))
    assert char["chromosome_like"] is False
    strat = choose_auto_strategy(char)
    assert strat["action"] == "direct"


def test_longread_assembly_routes_replicon(tmp_path):
    fasta = tmp_path / "lr.fasta"
    _fasta(fasta, {
        "chr1": "C" * 2_000_000,
        "plas1 circular=true": "A" * 90_000})
    char = characterize_assembly(str(fasta))
    assert char["chromosome_like"] is True
    assert char["read_type"] == "long"
    strat = choose_auto_strategy(char)
    assert strat["action"] == "extract"
    assert strat["classifier"] == "replicon"
    assert strat["route"] == "isolate"


def test_shortread_fragmented_assembly_routes_genomad(tmp_path):
    fasta = tmp_path / "sr.fasta"
    recs = {"frag%d" % i: "A" * 1_500 for i in range(600)}
    _fasta(fasta, recs)
    char = characterize_assembly(str(fasta))
    assert char["chromosome_like"] is True
    assert char["read_type"] == "short"
    strat = choose_auto_strategy(char)
    assert strat["classifier"] == "genomad"
    assert strat["route"] == "auto"


def test_metagenome_forces_genomad_meta(tmp_path):
    fasta = tmp_path / "meta.fasta"
    _fasta(fasta, {"c1": "C" * 1_000_000, "p": "A" * 40_000})
    char = characterize_assembly(str(fasta))
    strat = choose_auto_strategy(char, input_type="metagenome")
    assert strat["classifier"] == "genomad"
    assert strat["route"] == "meta"


# ---------------------------------------------------------------------------
# Replicon classifier
# ---------------------------------------------------------------------------

def _abricate_table(path, seq, gene="IncFII(S)"):
    header = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\t"
              "COVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\t"
              "PRODUCT\tRESISTANCE")
    row = ("in.fasta\t%s\t1\t1000\t+\t%s\t1-1000\t=\t0\t90\t99\tplasmidfinder\t"
           "X\treplicon\t-" % (seq, gene))
    _write(path, header + "\n" + row + "\n")


def test_classify_by_replicon(tmp_path):
    fasta = tmp_path / "lr.fasta"
    _fasta(fasta, {
        "chr1": "C" * 2_000_000,
        "plas1": "A" * 60_000,
        "other": "G" * 80_000})
    table = tmp_path / "rep.tsv"
    _abricate_table(str(table), "plas1")
    cc = {c.seq_id: c for c in classify_by_replicon(
        str(fasta), replicon_table=str(table))}
    assert cc["plas1"].label == "plasmid"
    assert cc["plas1"].source == "replicon"
    assert cc["chr1"].label == "chromosome"
    assert cc["other"].label == "unclassified"


def test_circular_small_contig_is_plasmid_without_replicon(tmp_path):
    fasta = tmp_path / "lr.fasta"
    _fasta(fasta, {
        "chr1": "C" * 2_000_000,
        "plas_c circular=true": "A" * 45_000})
    table = tmp_path / "empty.tsv"
    _write(table, "#FILE\tSEQUENCE\tGENE\tDATABASE\n")
    cc = {c.seq_id: c for c in classify_by_replicon(
        str(fasta), replicon_table=str(table))}
    assert cc["plas_c circular=true"].label == "plasmid"
    assert cc["chr1"].label == "chromosome"


# ---------------------------------------------------------------------------
# auto_main end to end (no external binaries; tables drive classifiers)
# ---------------------------------------------------------------------------

def test_auto_longread_extracts_and_scores(tmp_path):
    fasta = tmp_path / "lr.fasta"
    _fasta(fasta, {
        "chr1": "C" * 2_000_000,
        "plas1": "A" * 60_000})
    table = tmp_path / "rep.tsv"
    _abricate_table(str(table), "plas1")
    out = tmp_path / "out"
    rc = auto_main([str(fasta), "-o", str(out),
                    "--replicon-table", str(table),
                    "--no-abricate", "-q"])
    assert rc == 0
    status = _read_status(out / "auto_status.tsv")
    row = status[0]
    assert row["status"] == "scored"
    assert row["classifier"] == "replicon"
    assert row["route"] == "isolate"
    assert int(row["n_plasmid"]) == 1


def test_auto_direct_scores_plasmid_collection(tmp_path):
    fasta = tmp_path / "plasmids.fasta"
    _fasta(fasta, {"p1": "A" * 30_000, "p2": "C" * 70_000})
    out = tmp_path / "out"
    rc = auto_main([str(fasta), "-o", str(out), "--no-abricate", "-q"])
    assert rc == 0
    status = _read_status(out / "auto_status.tsv")
    assert status[0]["status"] == "scored"
    assert status[0]["route"] == "direct"
    assert os.path.exists(out / "plasmids" / "plasrisk_results.tsv")


def test_auto_rejects_reads(tmp_path):
    reads = tmp_path / "sample_R1.fastq.gz"
    _write(reads, "@r1\nACGT\n+\nIIII\n")
    rc = auto_main([str(reads), "-o", str(tmp_path / "out")])
    assert rc == 2


def _read_status(path):
    import csv
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
