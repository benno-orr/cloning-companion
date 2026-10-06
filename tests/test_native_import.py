import struct

import pytest

pytest.importorskip("webview")

from mac_app.main import NativeAPI, _serialize_result
from plasmid_verify.assembly import assemble_gibson
from plasmid_verify.verify import verify_consensus
from plasmid_verify.golden_gate_design import _write_snapgene_map


def test_design_drop_registration_reports_missing_target_and_is_idempotent():
    class Element:
        def __init__(self):
            self.listeners = []
        def on(self, event, callback):
            self.listeners.append((event, callback))
    class DOM:
        elements = []
        def get_elements(self, selector):
            assert selector == "#pick-design-map"
            return self.elements
    class Window:
        dom = DOM()
    api = NativeAPI()
    assert not api.register_design_drop_target()
    window = Window()
    api.bind(window)
    assert not api.register_design_drop_target()
    element = Element()
    window.dom.elements = [element]
    assert api.register_design_drop_target()
    assert api.register_design_drop_target()
    assert len(element.listeners) == 1
    assert element.listeners[0][0] == "drop"


def test_design_defaults_save_each_run_and_embed_all_graphics(tmp_path, monkeypatch):
    monkeypatch.setattr("mac_app.main.application_support", lambda: tmp_path / "support")
    source = tmp_path / "design.dna"
    _write_snapgene_map(source, "CACC" + "A" * 20 + "TGAAATGGCC", [
        {"name": "[Backbone]", "start": 0, "end": 28},
        {"name": "pUC ori", "start": 4, "end": 20},
        {"name": "{Insert}", "start": 28, "end": 34},
        {"name": "-", "start": 24, "end": 28},
        {"name": "-", "start": 0, "end": 4},
    ], "Test design")
    api = NativeAPI()
    first = api.run_annotated_golden_gate_design(str(source))
    second = api.run_annotated_golden_gate_design(str(source), "")
    assert first["ok"] and second["ok"], (first, second)
    assert first["outputDir"] != second["outputDir"]
    assert len(first["graphics"]) == 3  # linear map and two fragment-end views
    assert all(graphic["src"].startswith("data:image/") for graphic in first["graphics"])
    assert {graphic["group"] for graphic in first["graphics"]} == {"Linear map", "Fragment ends"}
    assert api.annotated_design_map_from_path(str(source))["fragmentCount"] == 2
    class Window:
        def evaluate_js(self, script):
            self.script = script
    window = Window()
    api.bind(window)
    api._receive_dropped_design({"dataTransfer": {"files": [{"pywebviewFullPath": str(source)}]}})
    assert window.script.startswith("window.nativeDesignDropped(")
    api._receive_dropped_design({"dataTransfer": {"files": [{"pywebviewFullPath": str(tmp_path / 'missing.dna')}]}})
    assert window.script.startswith("window.nativeDesignDropFailed(")


def _packet(packet_type: int, payload: bytes) -> bytes:
    return struct.pack(">BI", packet_type, len(payload)) + payload


def _snapgene(sequence: str) -> bytes:
    cookie = _packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 1, 1))
    dna = _packet(0x00, b"\x01" + sequence.encode("ascii"))
    return cookie + dna


def test_native_verification_accepts_snapgene_for_all_input_roles(tmp_path):
    left = "ACGTACGTACGTTAG"
    right = "TTCGGAACCTTAGGC"
    parent = "GGGCCC" + left + "GATTACAAA" + right + "AATTCC"
    insert = left + "ATGGCTTAA" + right
    target = "GGGCCC" + insert + "AATTCC"
    for filename, sequence in (
        ("parent.dna", parent),
        ("T01_insert_01.dna", insert),
        ("T01_consensus.dna", target),
    ):
        (tmp_path / filename).write_bytes(_snapgene(sequence))

    api = NativeAPI()
    parent_info = api._file_info(str(tmp_path / "parent.dna"))
    insert_info = api._file_info(str(tmp_path / "T01_insert_01.dna"), 1)
    insert_info["backbone"] = parent_info
    consensus_info = api._file_info(str(tmp_path / "T01_consensus.dna"))
    assert parent_info["length"] == len(parent)
    assert insert_info["targetId"] == "T01"
    assert consensus_info["targetId"] == "T01"

    result = api.run_verification(
        {
            "parentPath": "",
            "method": "Gibson",
            "armLength": 15,
            "inserts": [insert_info],
            "consensuses": [consensus_info],
        }
    )
    assert result["ok"] is True
    assert result["failures"] == []
    assert result["results"][0]["verdict"] == "PASS"


def test_serialized_coding_mutation_includes_codon_substitution():
    left = "ACGTACGTACGTTAG"
    right = "TTCGGAACCTTAGGC"
    parent = "GGGCCC" + left + "GATTACAAA" + right + "AATTCC"
    insert = left + "ATGGCTGACTAA" + right
    assembly = assemble_gibson(parent, [("T01_insert.fasta", insert)], "T01", arm_length=15)
    coding_start = next(index for index, annotation in enumerate(assembly.annotations) if annotation.region == "coding")
    observed = list(assembly.sequence)
    observed[coding_start + 5] = "C"  # GCT (Ala) -> GCC (Ala)
    result = verify_consensus(assembly, "".join(observed))

    mutation = _serialize_result(assembly, result)["mutations"][0]
    assert mutation["codon"] == {
        "expected": "GCT", "observed": "GCC", "expectedAa": "A", "observedAa": "A",
        "aaPosition": 2, "codonStart": coding_start + 4,
    }
