import struct
import json
import zipfile
from pathlib import Path

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


def test_design_library_persists_runs_without_export_folders(tmp_path, monkeypatch):
    monkeypatch.setattr("mac_app.main.application_support", lambda: tmp_path / "support")
    monkeypatch.setattr("tempfile.tempdir", str(tmp_path))
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
    old_id = first['plasmids'][0]['downloadId']
    first_bytes = api._design_files[old_id]['data']
    # Even a legacy project with an output folder must not silently save.
    second = api.run_annotated_golden_gate_design(str(source), str(tmp_path / 'old-save-location'))
    assert first["ok"] and second["ok"], (first, second)
    assert 'outputDir' not in first and 'outputDir' not in second
    assert set(tmp_path.iterdir()) == {source, tmp_path / 'support'}
    assert not first['temporary'] and not second['temporary']
    assert first['savedRun']['id'] != second['savedRun']['id']
    assert len(api.list_saved_designs()['runs']) == 2
    assert not api.download_design_file(old_id)['ok']
    assert len(first["graphics"]) == 3  # linear map and two fragment-end views
    assert all(graphic["src"].startswith("data:image/") for graphic in first["graphics"])
    assert {graphic["group"] for graphic in first["graphics"]} == {"Linear map", "Fragment ends"}
    assert all(g['downloadId'] in api._design_files for g in second['graphics'])
    assert {'synthesis_order.tsv', 'synthesis_order.fasta', 'pcr_primers.tsv', 'assembly_schematic.dna'} <= {d['filename'] for d in second['downloads']}
    # Cancelling a download creates nothing and keeps it available.
    monkeypatch.setattr(api, '_save_dialog', lambda *args: None)
    assert api.download_design_file(second['plasmids'][0]['downloadId'])['cancelled']
    assert api.download_design_bundle()['cancelled']
    assert set(tmp_path.iterdir()) == {source, tmp_path / 'support'}
    # A selected file is byte-identical; full and filtered archives work.
    saved = tmp_path / 'selected.dna'
    monkeypatch.setattr(api, '_save_dialog', lambda *args: str(saved))
    key = second['plasmids'][0]['downloadId']
    assert api.download_design_file(key)['ok']
    assert saved.read_bytes() == api._design_files[key]['data']
    for group in ('all', 'plasmids', 'orders', 'graphics'):
        bundle = tmp_path / f'{group}.zip'
        monkeypatch.setattr(api, '_save_dialog', lambda *args: str(bundle))
        assert api.download_design_bundle(group)['ok']
        with zipfile.ZipFile(bundle) as archive:
            expected = [f for f in api._design_files.values() if group == 'all' or f['group'] == group]
            assert set(archive.namelist()) == {f['relativePath'] for f in expected}
            for f in expected:
                assert archive.read(f['relativePath']) == f['data']
            if group == 'all':
                report = json.loads(archive.read('design.json'))
                assert not Path(report['assembledMap']).is_absolute()
                assert report['assembledMap'] in archive.namelist()
    assert not api.download_design_file('/etc/passwd')['ok']
    assert not api.download_design_bundle('../other')['ok']
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
    # Simulate relaunch, without the original input map or in-memory cache.
    source.unlink()
    restarted = NativeAPI()
    assert len(restarted.list_saved_designs()['runs']) == 2
    restored = restarted.open_saved_design(first['savedRun']['id'])
    assert restored['ok'] and restored['graphics'] == first['graphics']
    assert restored['plasmids'] == first['plasmids']
    restored_id = restored['plasmids'][0]['downloadId']
    monkeypatch.setattr(restarted, '_save_dialog', lambda *args: str(tmp_path / 'reopened.dna'))
    assert restarted.download_design_file(restored_id)['ok']
    assert (tmp_path / 'reopened.dna').read_bytes() == first_bytes
    assert not restarted.open_saved_design("' OR 1=1 --")['ok']
    assert restored_id in restarted._design_files  # Failed loads preserve active run.


def test_design_library_reports_corruption_without_overwriting_it(tmp_path, monkeypatch):
    monkeypatch.setattr('mac_app.main.application_support', lambda: tmp_path)
    folder = tmp_path / 'DesignLibrary'
    folder.mkdir()
    database = folder / 'designs.sqlite3'
    database.write_bytes(b'invalid database')
    api = NativeAPI()
    assert not api.list_saved_designs()['ok']
    assert not api.open_saved_design('missing')['ok']
    assert database.read_bytes() == b'invalid database'


def test_design_library_failed_transaction_keeps_previous_runs(tmp_path):
    import sqlite3
    from mac_app import design_library
    result = {'ok': True, 'project': 'saved', 'downloads': [], 'graphics': []}
    run = design_library.save(tmp_path, result, {})
    with pytest.raises(sqlite3.IntegrityError):
        design_library.save(tmp_path, {**result, 'project': None}, {})
    assert design_library.list_runs(tmp_path) == [run]
    assert design_library.load(tmp_path, run['id'])[0]['project'] == 'saved'


def test_failed_design_cleans_staging_and_preserves_existing_downloads(tmp_path, monkeypatch):
    monkeypatch.setattr('tempfile.tempdir', str(tmp_path))
    api = NativeAPI()
    api._design_files = {'existing': {'data': b'previous result'}}
    def fail(root):
        (root / 'partial.txt').write_text('partial render')
        raise ValueError('renderer failed')
    with pytest.raises(ValueError, match='renderer failed'):
        api._temporary_design(fail)
    assert list(tmp_path.iterdir()) == []
    assert api._design_files['existing']['data'] == b'previous result'


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
