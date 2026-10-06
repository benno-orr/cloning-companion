from plasmid_verify.assembly import assemble_gibson
from plasmid_verify.fasta import reverse_complement
from plasmid_verify.verify import verify_consensus


def _assembly():
    left = "ACGTACGTACGTTAG"
    right = "TTCGGAACCTTAGGC"
    parent = "GGGCCC" + left + "GATTACAAA" + right + "AATTCC"
    insert = left + "ATGGCTGACTAA" + right
    return assemble_gibson(parent, [("insert.fasta", insert)], "T01", arm_length=15)


def test_exact_rotated_reverse_complement_is_pass():
    assembly = _assembly()
    rotated = assembly.sequence[11:] + assembly.sequence[:11]
    result = verify_consensus(assembly, reverse_complement(rotated))

    assert result.verdict == "PASS"
    assert result.orientation == "reverse complement"
    assert result.mutations == []


def test_silent_substitution_is_review():
    assembly = _assembly()
    coding_start = next(i for i, a in enumerate(assembly.annotations) if a.region == "coding")
    observed = list(assembly.sequence)
    # GCT (Ala) -> GCC (Ala)
    observed[coding_start + 5] = "C"
    result = verify_consensus(assembly, "".join(observed))

    assert result.verdict == "REVIEW"
    assert len(result.mutations) == 1
    assert result.mutations[0].category == "silent"
    assert result.mutations[0].protein_change == "A2="


def test_nonsense_substitution_is_fail():
    assembly = _assembly()
    coding_start = next(i for i, a in enumerate(assembly.annotations) if a.region == "coding")
    observed = list(assembly.sequence)
    # GAC (Asp) -> TAA (stop), represented as three adjacent SNVs. At least one
    # call must resolve to a nonsense codon using the full observed codon.
    observed[coding_start + 6 : coding_start + 9] = "TAA"
    result = verify_consensus(assembly, "".join(observed))

    assert result.verdict == "FAIL"
    assert "nonsense" in {mutation.category for mutation in result.mutations}


def test_one_base_coding_deletion_is_frameshift():
    assembly = _assembly()
    coding_start = next(i for i, a in enumerate(assembly.annotations) if a.region == "coding")
    observed = assembly.sequence[: coding_start + 2] + assembly.sequence[coding_start + 3 :]
    result = verify_consensus(assembly, observed)

    assert result.verdict == "FAIL"
    assert result.mutations[0].category == "frameshift"
