import pytest

from plasmid_verify.assembly import AssemblyError, assemble_gibson, assemble_golden_gate
from plasmid_verify.fasta import reverse_complement


def test_gibson_replaces_bracketed_interval_and_marks_regions():
    left = "ACGTACGTACGTTAG"
    right = "TTCGGAACCTTAGGC"
    parent = "GGGCCC" + left + "GATTACAAA" + right + "AATTCC"
    insert = left + "ATGGCTTAA" + right

    assembly = assemble_gibson(parent, [("T01_insert.fasta", insert)], "T01", arm_length=15)

    assert assembly.sequence == "GGGCCC" + insert + "AATTCC"
    assert [a.region for a in assembly.annotations[6 : 6 + 15]] == ["junction"] * 15
    assert [a.region for a in assembly.annotations[21:30]] == ["coding"] * 9
    assert assembly.annotations[21].feature_start == 21
    assert assembly.annotations[29].feature_end == 30


def test_gibson_accepts_reverse_complement_insert():
    left = "ACGTACGTACGTTAG"
    right = "TTCGGAACCTTAGGC"
    parent = "GGGCCC" + left + "GATTACAAA" + right + "AATTCC"
    forward_insert = left + "ATGGCTTAA" + right

    assembly = assemble_gibson(
        parent, [("reverse.fasta", reverse_complement(forward_insert))], "T01", 15
    )

    assert assembly.sequence == "GGGCCC" + forward_insert + "AATTCC"
    assert "reverse complement" in assembly.notes[0]


def _gg_fragment(left: str, payload: str, right: str) -> str:
    return "GGTCTC" + "A" + left + payload + right + "T" + "GAGACC"


def test_golden_gate_single_insert_removes_sites_and_drop_out():
    parent = "GGTCTC" + "A" + "AATG" + "TTAACCAA" + "CGCT" + "T" + "GAGACC" + "TTCCAAGG"
    insert = _gg_fragment("AATG", "ATGGCTTAA", "CGCT")

    assembly = assemble_golden_gate(parent, [("T01_insert.fasta", insert)], "T01", "BsaI")

    assert assembly.sequence == "AATG" + "ATGGCTTAA" + "CGCT" + "TTCCAAGG"
    assert "GGTCTC" not in assembly.sequence
    assert "GAGACC" not in assembly.sequence


def test_golden_gate_multiple_inserts_deduplicates_fusion_overhang():
    parent = "GGTCTC" + "A" + "AATG" + "TTAACCAA" + "GGCC" + "T" + "GAGACC" + "TTCCAAGG"
    first = _gg_fragment("AATG", "ATGAAA", "CGCT")
    second = _gg_fragment("CGCT", "GGGTAA", "GGCC")

    assembly = assemble_golden_gate(
        parent,
        [("T01_insert_01.fasta", first), ("T01_insert_02.fasta", second)],
        "T01",
        "BsaI",
    )

    assert assembly.sequence == "AATGATGAAACGCTGGGTAAGGCCTTCCAAGG"


def test_golden_gate_rejects_bad_overhang_chain():
    parent = "GGTCTC" + "A" + "AATG" + "TTAACCAA" + "GGCC" + "T" + "GAGACC" + "TTCCAAGG"
    first = _gg_fragment("AATG", "ATGAAA", "CGCT")
    second = _gg_fragment("TTTT", "GGGTAA", "GGCC")

    with pytest.raises(AssemblyError, match="Incompatible fusion overhangs"):
        assemble_golden_gate(
            parent,
            [("first.fasta", first), ("second.fasta", second)],
            "T01",
            "BsaI",
        )
