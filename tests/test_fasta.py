import pytest
import struct
from openpyxl import Workbook

from plasmid_verify.fasta import infer_target_id, parse_fasta, parse_sequence_file, reverse_complement
from mac_app.main import _table_insert_rows, _xlsx_insert_rows


def test_filename_prefix_inference():
    assert infer_target_id("T01_consensus.fasta") == "T01"
    assert infer_target_id("clone-7__insert_02.fa") == "clone-7"


def test_parse_fasta_and_reverse_complement():
    records = parse_fasta(">sample notes\nacgt\nNN\n")
    assert records[0].name == "sample"
    assert records[0].sequence == "ACGTNN"
    assert reverse_complement("ACGTRY") == "RYACGT"


def test_parse_fasta_rejects_non_dna():
    with pytest.raises(ValueError, match="Unsupported"):
        parse_fasta(">bad\nACGTZ")


def _snapgene_packet(packet_type: int, payload: bytes) -> bytes:
    return struct.pack(">BI", packet_type, len(payload)) + payload


def test_parse_snapgene_dna_sequence():
    cookie = _snapgene_packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 1, 1))
    dna = _snapgene_packet(0x00, b"\x01ACGTNN")
    records = parse_sequence_file(cookie + dna, "T01_insert_01.dna")
    assert len(records) == 1
    assert records[0].name == "T01_insert_01"
    assert records[0].sequence == "ACGTNN"


def test_reject_invalid_snapgene_file():
    with pytest.raises(ValueError, match="Invalid SnapGene"):
        parse_sequence_file(b">not a binary SnapGene file\nACGT", "T01.dna")


def test_parse_csv_and_tsv_insert_tables():
    csv_rows = _table_insert_rows(
        "target_id,name,dna_sequence,order\nT01,part_a,ACGT,2\nT01,part_b,TTAA,1\n",
        "inserts.csv",
    )
    tsv_rows = _table_insert_rows("clone_id\tsequence\nT02\tggcc\n", "inserts.tsv")
    assert [(row["targetId"], row["recordName"], row["sequence"], row["order"]) for row in csv_rows] == [
        ("T01", "part_a", "ACGT", 2), ("T01", "part_b", "TTAA", 1)
    ]
    assert tsv_rows[0]["targetId"] == "T02"
    assert tsv_rows[0]["sequence"] == "GGCC"


def test_parse_xlsx_insert_table(tmp_path):
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.append(["target_id", "name", "dna_sequence", "order"])
    worksheet.append(["T01", "part_a", "ACGT", 2])
    worksheet.append(["T02", "part_b", "TTAA", 1])
    path = tmp_path / "inserts.xlsx"
    workbook.save(path)

    rows = _xlsx_insert_rows(str(path))
    assert [(row["targetId"], row["recordName"], row["sequence"], row["order"]) for row in rows] == [
        ("T01", "part_a", "ACGT", 2), ("T02", "part_b", "TTAA", 1)
    ]
