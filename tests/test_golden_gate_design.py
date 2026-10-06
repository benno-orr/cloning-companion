from pathlib import Path
import csv
import json
import struct
from xml.etree.ElementTree import fromstring, Element, SubElement

import yaml
from Bio import SeqIO

import pytest

from plasmid_verify.golden_gate_design import Scaffold, plan_from_annotated_snapgene, run_plan, parse_variable_elements, run_annotated_snapgene_design, _write_snapgene_map


def _native_linker_map(tmp_path):
    sequence = "CACC" + "AAAA" * 5 + "TGAA" + "ATGGCC" + "GGATCTGGATCTGGA" + "GCGGCG"
    path = tmp_path / "native_linker.dna"
    _write_snapgene_map(path, sequence, [
        {"name": "[Backbone]", "start": 0, "end": 28},
        {"name": "ori pUC", "start": 4, "end": 20},
        {"name": "{scFv}", "start": 28, "end": 34},
        {"name": "{TM-tail}", "start": 49, "end": 55},
        {"name": "5aa linker", "start": 34, "end": 49, "type": "CDS"},
        {"name": "-", "start": 24, "end": 28},
        {"name": "-", "start": 34, "end": 49},
        {"name": "-", "start": 0, "end": 4},
    ], "Native sequence with internal linker")
    return path, sequence


def test_native_map_export_preserves_every_base_and_linker(tmp_path):
    path, source = _native_linker_map(tmp_path)
    result = run_annotated_snapgene_design(path, tmp_path / "exports")
    assembled = str(SeqIO.read(result["assembledMap"], "snapgene").seq)
    assert len(assembled) == len(source)
    assert assembled in source + source
    assert "ATGGCCGGATCTGGATCTGGAGCGGCG" in assembled + assembled
    report = json.loads((Path(result["outputDir"]) / "design.json").read_text())
    assert all(row["left_fusion_owner"] == row["right_fusion_owner"] == "core" for row in report["fragments"])
    schematic = str(SeqIO.read(result["schematicMap"], "snapgene").seq)
    for row in report["fragments"]:
        core = row["core_sequence"]
        assert core in schematic
        assert row["final_fragment_sequence"] == "GGCTACCGTCTCA" + core + "AGAGACGGGCTAC"
    record = SeqIO.read(result["schematicMap"], "snapgene")
    tail = next(f for f in record.features if f.qualifiers.get("label") == ["{TM-tail}"])
    assert str(tail.extract(record.seq)) == "GCGGCG"
    assert schematic[int(tail.location.start) - 15:int(tail.location.start)] == "GGATCTGGATCTGGA"
    from plasmid_verify.golden_gate_design import _snapgene_features_xml, _qualifier
    tail_xml = next(f for f in _snapgene_features_xml(Path(result["schematicMap"])) if f.get("name") == "{TM-tail}")
    extensions = json.loads(_qualifier(tail_xml, "assigned_extensions"))
    assert extensions == [[int(tail.location.start) - 15, int(tail.location.start)]]


@pytest.mark.parametrize("offset", range(12))
def test_hatched_linker_is_partitioned_once_between_physical_fragments(offset):
    from plasmid_verify.golden_gate_design import _fragment_extensions
    # The selected four bases are present on both ends, displayed downstream.
    common = {"left_fusion_owner": "core", "right_fusion_owner": "core", "left_fusion": "CACC", "right_fusion": "GGAT"}
    upstream = dict(common, variant_window={"prefix": "", "suffix": "A" * (offset + 4)})
    downstream = dict(common, variant_window={"prefix": "A" * (15 - offset), "suffix": ""})
    left = _fragment_extensions(upstream, 100 + offset + 4)
    right = _fragment_extensions(downstream, 15 - offset + 100)
    assert left == ([(100, 100 + offset)] if offset else [])
    assert right == [(0, 15 - offset)]
    assert sum(b - a for a, b in left + right) == 15


def test_hatching_excludes_fixed_bases_inside_core_and_terminal_overlap():
    from plasmid_verify.golden_gate_design import _fragment_extensions, _fragment_core_bounds
    row = {"left_fusion_owner": "extension", "right_fusion_owner": "extension",
           "left_fusion": "CACC", "right_fusion": "TGAA",
           "variant_window": {"prefix": "AAA", "suffix": "GG", "fixed_core_prefix": "AT", "fixed_core_suffix": "GC"}}
    assert _fragment_core_bounds(row, 25) == (7, 19)
    assert _fragment_extensions(row, 25) == [(0, 7), (19, 21)]


def test_diagonal_hatching_is_clipped_to_its_rectangle():
    from PIL import Image, ImageDraw
    from plasmid_verify.golden_gate_design import _draw_hatched_region
    image = Image.new("RGB", (80, 40), "white")
    _draw_hatched_region(ImageDraw.Draw(image), (20, 10, 60, 30), "blue")
    changed = [(x, y) for y in range(40) for x in range(80) if image.getpixel((x, y)) != (255, 255, 255)]
    assert changed
    assert all(20 <= x <= 60 and 10 <= y <= 30 for x, y in changed)


@pytest.mark.parametrize("linker_offset", range(12))
@pytest.mark.parametrize("use_variants", [False, True])
def test_any_cut_within_linker_preserves_complete_uploaded_sequence(tmp_path, linker_offset, use_variants):
    from plasmid_verify.golden_gate_design import Candidate, _coordinate_fragment_plan, _expanded_fragment_specs, _assembled_circular_map
    path, source = _native_linker_map(tmp_path)
    scaffold = Scaffold.from_path(path)
    plan = plan_from_annotated_snapgene(path)
    if use_variants:
        plan["fragments"][1]["variants"] = [{"name": "changed", "sequence": "ATGAGC"}]
        plan["fragments"][2]["variants"] = [{"name": "changed", "sequence": "GCTGCT"}]
    chosen = {}
    for spec in plan["junction_candidates"]:
        start, end = spec["scaffold_range"]
        offset = linker_offset if start == 34 else 0
        chosen[spec["id"]] = Candidate(spec["id"], source[start + offset:start + offset + 4], offset, 0, "test")
    resolved, expected = _coordinate_fragment_plan(plan, scaffold, chosen)
    rows = []
    for spec in _expanded_fragment_specs(resolved):
        if "sequence" in spec["core"]:
            core = spec["core"]["sequence"]
        else:
            a, b = spec["core"]["scaffold_range"]
            core = source[a:b] if b > a else source[a:] + source[:b]
        rows.append({"fragment_id": spec["id"], "core_sequence": core, "left_fusion": chosen[spec["left_junction"]].fusion,
                     "right_fusion": chosen[spec["right_junction"]].fusion, "left_fusion_owner": "core", "right_fusion_owner": "core"})
    actual, _, _ = _assembled_circular_map(rows)
    assert len(actual) == len(expected) == len(source)
    target = source[:28] + "ATGAGC" + source[34:49] + "GCTGCT" if use_variants else source
    assert actual in target + target


def test_pasted_variants_replace_only_core_and_keep_native_linker(tmp_path):
    path, source = _native_linker_map(tmp_path)
    result = run_annotated_snapgene_design(path, tmp_path / "exports", variable_texts={
        "fragment_2": "changed,ATGAGC\nlonger,ATGGCCGCC", "fragment_3": "tail,GCTGCT"})
    assembled = str(SeqIO.read(result["assembledMap"], "snapgene").seq)
    expected = source[:28] + "ATGAGC" + source[34:49] + "GCTGCT"
    assert len(assembled) == len(expected)
    assert assembled in expected + expected
    assert "ATGAGCGGATCTGGATCTGGAGCTGCT" in assembled + assembled
    record = SeqIO.read(result["schematicMap"], "snapgene")
    tail = next(f for f in record.features if f.qualifiers.get("label") == ["{TM-tail}"])
    assert str(tail.extract(record.seq)) == "GCTGCT"
    for junction in result["automaticJunctions"]:
        assert junction["left_end"][-4:] == junction["right_end"][:4] == junction["shared_fusion"]
        assert junction["left_end"][:-4] + junction["right_end"] == junction["scaffold_region"]


def test_invalid_reconstruction_writes_no_order_files(tmp_path, monkeypatch):
    from plasmid_verify import golden_gate_design as designer
    path, source = _native_linker_map(tmp_path)
    real_assembly = designer._assembled_circular_map
    def broken(rows):
        sequence, starts, retained = real_assembly(rows)
        return sequence + "GGAT", starts, retained
    monkeypatch.setattr(designer, "_assembled_circular_map", broken)
    with pytest.raises(ValueError, match="does not reproduce"):
        run_annotated_snapgene_design(path, tmp_path / "exports")
    assert not (tmp_path / "exports").exists()


def test_explicit_plan_rejects_added_fusion_already_inside_core(tmp_path):
    scaffold = tmp_path / "source.fasta"
    scaffold.write_text(">source\nAAAATGAA\n")
    plan = {"scaffold": str(scaffold), "junction_candidates": [{"id": "a", "sequence": "TGAA"}],
            "fragments": [{"id": "bad", "left_junction": "a", "right_junction": "a", "core": {"sequence": "TGAAGCGGCC"}}]}
    path = tmp_path / "bad_repeat.yaml"
    path.write_text(yaml.safe_dump(plan))
    with pytest.raises(ValueError, match="adding it again would duplicate 4 bp"):
        run_plan(path, tmp_path / "exports")
    assert not (tmp_path / "exports").exists()


@pytest.mark.parametrize("junction_labels", [("(junction after backbone)", "(junction after payload)"), ("-", "-"), ("-", "(junction after payload)")])
def test_annotated_snapgene_map_becomes_design_plan(tmp_path: Path, monkeypatch, junction_labels):
    """The desktop workflow uses semantic map labels, not prescribed IDs."""
    map_path = tmp_path / "annotated_design.dna"
    map_path.write_bytes(b"SnapGene placeholder")
    scaffold = Scaffold(
        map_path,
        "CACCGTACAAAGGCCT",
        {
            "[Constant backbone]": [(0, 4)],
            "pUC ori": [(1, 3)],
            "{Payload}": [(8, 12)],
            "ordinary plasmid annotation": [(5, 7)],
            "ordinary-hyphenated-name": [(0, 16)],
        },
    )
    for label, span in zip(junction_labels, [(4, 8), (12, 16)]):
        scaffold.features.setdefault(label, []).append(span)
    monkeypatch.setattr(Scaffold, "from_path", classmethod(lambda cls, path: scaffold))

    plan = plan_from_annotated_snapgene(map_path, "BsaI")

    assert plan["enzyme"] == "BsaI"
    assert [fragment["name"] for fragment in plan["fragments"]] == ["Constant backbone", "Payload"]
    assert [fragment["source"] for fragment in plan["fragments"]] == ["PCR", "SYN"]
    assert plan["fragments"][0]["left_junction"] == "junction_2"
    assert plan["fragments"][0]["right_junction"] == "junction_1"
    assert [item["id"] for item in plan["junction_candidates"]] == ["junction_1", "junction_2"]


def test_annotated_map_rejects_overlapping_fragment_cores(tmp_path: Path, monkeypatch):
    map_path = tmp_path / "invalid_design.dna"
    map_path.write_bytes(b"SnapGene placeholder")
    scaffold = Scaffold(
        map_path,
        "CACCGTACAAAGGCCT",
        {
            "[Backbone]": [(0, 8)], "pUC ori": [(1, 3)],
            "{Overlapping insert}": [(6, 12)],
            "(junction one)": [(8, 12)], "(junction two)": [(12, 16)],
        },
    )
    monkeypatch.setattr(Scaffold, "from_path", classmethod(lambda cls, path: scaffold))

    with pytest.raises(ValueError, match="Fragment core annotations overlap"):
        plan_from_annotated_snapgene(map_path)


def test_design_writes_order_package_and_respects_core_owned_fusions(tmp_path: Path):
    scaffold = tmp_path / "scaffold.fasta"
    scaffold.write_text(">backbone\nCACCAAAAGCCT\n", encoding="utf-8")
    plan = {
        "project": "test design",
        "scaffold": str(scaffold),
        "enzyme": "BsmBI",
        "defaults": {"pcr_anneal_bp": 6},
        "junction_candidates": [
            {"id": "backbone_to_insert", "sequence": "CACC", "preferred": "CACC"},
            {"id": "insert_to_backbone", "sequence": "GCCT", "preferred": "GCCT"},
        ],
        "fragments": [
            {"id": "backbone", "source": "PCR", "core": {"scaffold_range": [0, 12]},
             "left_junction": "backbone_to_insert", "right_junction": "insert_to_backbone",
             "left_fusion_owner": "core", "right_fusion_owner": "core"},
            {"id": "insert", "source": "SYN", "core": {"sequence": "ATGGCC"},
             "left_junction": "insert_to_backbone", "right_junction": "backbone_to_insert",
             "left_fusion_owner": "extension", "right_fusion_owner": "extension"},
        ],
    }
    plan_path = tmp_path / "plan.yaml"
    plan_path.write_text(yaml.safe_dump(plan), encoding="utf-8")

    result = run_plan(plan_path, tmp_path / "exports")

    assert result["ok"] is True
    assert [item["selected_fusion"] for item in result["junctions"]] == ["CACC", "GCCT"]
    output = Path(result["outputDir"])
    assert (output / "fragments.tsv").exists()
    assert (output / "pcr_primers.tsv").exists()
    assert (output / "synthesis_order.fasta").exists()
    assert Path(result["assembledMap"]).exists()
    schematic = SeqIO.read(result["schematicMap"], "snapgene")
    assert str(schematic.seq).count("N") == 12  # six visual spacer bases per gap
    assert "N" not in (output / "synthesis_order.fasta").read_text()
    assert len(result["junctionScreencaps"]) == 2
    assert all(Path(image).exists() for image in result["junctionScreencaps"])
    assert all(Path(image).with_name(Path(image).stem + "_assembled.png").exists() for image in result["junctionScreencaps"])
    assembled = SeqIO.read(result["assembledMap"], "snapgene")
    assert assembled.annotations["topology"] == "circular"
    assert len(assembled) == 18
    fragments = (output / "fragments.tsv").read_text(encoding="utf-8")
    # Core-owned CACC and GCCT are retained once; they are not repeated in backbone adapters.
    assert "GGCTACCGTCTCA\tCACCAAAAGCCT\tAGAGACGGGCTAC" in fragments


@pytest.mark.parametrize("text", [
    "name\tdna_sequence\nparental\tatggcc\nsweeping\tATGAGC",
    "parental,ATGGCC\nsweeping,ATGAGC",
    ">parental\nATG\nGCC\n>sweeping\nATGAGC",
    "name,aa,dna\nparental,MA,ATGGCC\nsweeping,MS,ATGAGC",
])
def test_pasted_variants(text):
    assert parse_variable_elements(text) == [
        {"name": "parental", "sequence": "ATGGCC"},
        {"name": "sweeping", "sequence": "ATGAGC"},
    ]


@pytest.mark.parametrize("text, message", [
    ("same,ATG\nSAME,GCC", "duplicate variant"),
    ("bad,ATN", "unambiguous"),
    ("bad,XYZ", "Unsupported"),
    (">empty\n>valid\nATG", "Sequence is empty"),
    ("name,dna_sequence", "no named DNA"),
    ("nameonly", "needs a variant name"),
])
def test_pasted_variants_reject_invalid(text, message):
    with pytest.raises(ValueError, match=message):
        parse_variable_elements(text)


def test_pasted_libraries_export_one_order_per_option_not_per_combination(tmp_path, monkeypatch):
    scaffold_path = tmp_path / "library.dna"
    scaffold_path.write_bytes(b"placeholder")
    plan = {
        "project": "library", "enzyme": "BsmBI", "defaults": {"pcr_anneal_bp": 4},
        "junction_candidates": [{"id": name, "sequence": seq} for name, seq in [("a", "CACC"), ("b", "GCCT"), ("c", "TGAA")]],
        "fragments": [
            {"id": "bb", "name": "backbone", "source": "PCR", "core": {"sequence": "AAAAAA"}, "left_junction": "c", "right_junction": "a", "annotation_kind": "fixed"},
            {"id": "scfv", "name": "scFv", "source": "SYN", "core": {"sequence": "NNN"}, "left_junction": "a", "right_junction": "b", "annotation_kind": "variable"},
            {"id": "tail", "name": "TM-tail", "source": "SYN", "core": {"sequence": "NNN"}, "left_junction": "b", "right_junction": "c", "annotation_kind": "variable"},
        ],
    }
    monkeypatch.setattr("plasmid_verify.golden_gate_design.plan_from_annotated_snapgene", lambda *args: plan)
    monkeypatch.setattr("plasmid_verify.golden_gate_design._snapgene_features_xml", lambda *args: [])
    monkeypatch.setattr(Scaffold, "from_path", classmethod(lambda cls, path: Scaffold(path, "AAAAAA", {})))
    result = run_annotated_snapgene_design(scaffold_path, tmp_path / "exports", variable_texts={
        "scfv": "parental,ATGGCC\nsweeping,ATGAGCGCC",
        "tail": ">WT\nGGT\n>dead\nGCGGCG",
    })
    assert len(result["fragments"]) == 5  # One backbone + two scFvs + two tails.
    assert [row["assembly_position"] for row in result["fragments"]] == [1, 2, 2, 3, 3]
    output = Path(result["outputDir"])
    with (output / "synthesis_order.tsv").open() as handle:
        orders = list(csv.DictReader(handle, delimiter="\t"))
    assert len(orders) == 4
    assert len(list(SeqIO.parse(output / "synthesis_order.fasta", "fasta"))) == 4
    for row in orders:
        seq = row["final_fragment_sequence"]
        assert seq.startswith("GGCTACCGTCTCA" + row["left_fusion"])
        assert seq.endswith(row["right_fusion"] + "AGAGACGGGCTAC")
    assert {row["left_fusion"] for row in orders if row["element_id"] == "scfv"} == {"CACC"}
    assembled = SeqIO.read(result["assembledMap"], "snapgene")
    assert str(assembled.seq) == "AAAAAACACCATGGCCGCCTGGTTGAA"
    assert [row["variant_name"] for row in result["mapVariants"]] == ["", "parental", "WT"]
    report = json.loads((output / "design.json").read_text())
    assert report["mapVariants"] == result["mapVariants"]
    assert result["plasmidCount"] == result["combinationCount"] == 4
    assert len(list(Path(result["plasmidFolder"]).glob('*.dna'))) == 4
    assert {str(SeqIO.read(p['path'], 'snapgene').seq) for p in result['plasmids']} == {
        'AAAAAACACC' + scfv + 'GCCT' + tail + 'TGAA'
        for scfv in ('ATGGCC', 'ATGAGCGCC') for tail in ('GGT', 'GCGGCG')}
    assert "variants" not in plan["fragments"][1]  # Input plan remains untouched.
    with pytest.raises(ValueError, match="Unknown variable fragment"):
        run_annotated_snapgene_design(scaffold_path, tmp_path, variable_texts={"bb": "x,ATG"})


def test_variants_enforce_owned_fusions(tmp_path):
    scaffold = tmp_path / "scaffold.fasta"
    scaffold.write_text(">scaffold\nAAAAAA\n")
    plan = {"scaffold": str(scaffold), "junction_candidates": [{"id": "a", "sequence": "CACC"}],
            "fragments": [{"id": "insert", "left_junction": "a", "right_junction": "a",
                           "left_fusion_owner": "core", "variants": [{"name": "bad", "sequence": "ATGGCC"}]}]}
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(plan))
    with pytest.raises(ValueError, match="owned left fusion CACC"):
        run_plan(path, tmp_path / "exports")


def test_schematic_preserves_native_annotations_and_staggered_strand_colors(tmp_path):
    source = tmp_path / "native.dna"
    coding = Element("Feature", {"name": "test protein", "directionality": "1", "type": "CDS"})
    SubElement(coding, "Segment", {"range": "5-7", "color": "#99cc00", "translated": "1"})
    _write_snapgene_map(source, "CACCAAAAGCCT", [
        {"name": "promoter", "type": "promoter", "start": 4, "end": 8, "color": "#99cc00"},
        {"xml": coding},
    ], "Native scaffold")
    plan = {"scaffold": str(source), "defaults": {"pcr_anneal_bp": 4}, "project": "colors & annotations",
            "junction_candidates": [{"id": "a", "sequence": "CACC"}, {"id": "b", "sequence": "GCCT"}],
            "fragments": [
                {"id": "backbone", "name": "pCAG", "annotation_kind": "fixed", "is_backbone": True,
                 "source": "PCR", "core": {"scaffold_range": [0, 12]}, "left_junction": "a", "right_junction": "b",
                 "left_fusion_owner": "core", "right_fusion_owner": "core"},
                {"id": "insert", "name": "scFv", "annotation_kind": "variable", "core": {"sequence": "ATGGCC"},
                 "left_junction": "b", "right_junction": "a"},
            ]}
    path = tmp_path / "colors.yaml"
    path.write_text(yaml.safe_dump(plan))
    result = run_plan(path, tmp_path / "exports")
    raw = Path(result["schematicMap"]).read_bytes()
    packets, offset = {}, 0
    while offset < len(raw):
        kind, size = struct.unpack(">BI", raw[offset:offset + 5])
        packets[kind] = raw[offset + 5:offset + 5 + size]
        offset += 5 + size
    features = fromstring(packets[10])
    native = next(f for f in features if f.get("name") == "promoter")
    assert native.get("type") == "promoter"
    assert native.find("Segment").get("range") == "5-8"
    assert native.find("Segment").get("color") == "#99cc00"
    assert "[pCAG]" in [f.get("name") for f in features]
    assert "{scFv}" in [f.get("name") for f in features]
    block = next(f for f in features if f.get("name") == "[pCAG]")
    assert block.get("directionality") == "0"
    assert block.find("Segment").get("range") == "1-12"
    gaps = [f for f in features if f.get("name") == "└┐"]
    assert len(gaps) == 2
    assert [s.get("name") for s in gaps[0].findall("Segment")] == ["←BsmBI", "//", "BsmBI→"]
    colors = fromstring(packets[20])
    top = colors.find("TopStrand").findall("ColorRange")
    bottom = colors.find("BottomStrand").findall("ColorRange")
    assert top[0].attrib == {"range": "0..3", "colors": "#ff8080"}
    assert bottom[0].attrib == {"range": "0..3", "colors": "gray - 50"}
    assert top[2].attrib == {"range": "8..11", "colors": "gray - 50"}
    assert bottom[2].attrib == {"range": "8..11", "colors": "#ff8080"}
    assert str(SeqIO.read(result["assembledMap"], "snapgene").seq) == "AAAAGCCTATGGCCCACC"


@pytest.mark.parametrize("with_variant", [False, True])
def test_overhang_crossing_core_boundary_retains_per_base_ownership(tmp_path, with_variant):
    source = tmp_path / "crossing.dna"
    sequence = "CACC" + "A" * 20 + "TGAAATGGCC"
    _write_snapgene_map(source, sequence, [
        {"name": "[Backbone]", "start": 0, "end": 28},
        {"name": "pUC ori", "start": 4, "end": 20},
        {"name": "{Insert}", "start": 28, "end": 34},
        {"name": "-", "start": 26, "end": 30},  # AA|AT crosses the core boundary
        {"name": "-", "start": 0, "end": 4},
    ], "Shared fusion spans two native cores")
    variable_id = next(f["id"] for f in plan_from_annotated_snapgene(source)["fragments"] if f.get("annotation_kind") == "variable")
    # The leading AT is covered by '-', so it is supplied by the scaffold.
    result = run_annotated_snapgene_design(source, tmp_path / "out", variable_texts={variable_id: ">longer\nGAAAGCC"} if with_variant else None)
    raw, offset, packets = Path(result["schematicMap"]).read_bytes(), 0, {}
    while offset < len(raw):
        kind, size = struct.unpack(">BI", raw[offset:offset + 5])
        packets[kind] = raw[offset + 5:offset + 5 + size]
        offset += size + 5
    features, colors = fromstring(packets[10]), fromstring(packets[20])
    gap = next(feature for feature in features if feature.get("name") == "└┐")
    segments = gap.findall("Segment")
    left_end = int(segments[0].get("range").split("-")[0]) - 1
    right_start = int(segments[-1].get("range").split("-")[-1])
    def color_at(strand, position):
        for item in colors.find(strand):
            a, b = map(int, item.get("range").split(".."))
            if a <= position <= b:
                return item.get("colors")
    expected = ["#ff8080", "#ff8080", "#ffbf80", "#ffbf80"]
    assert [color_at("BottomStrand", left_end - 4 + i) for i in range(4)] == expected
    assert [color_at("TopStrand", right_start + i) for i in range(4)] == expected
    target = sequence[:28] + ("ATGAAAGCC" if with_variant else "ATGGCC")
    actual = str(SeqIO.read(result["assembledMap"], "snapgene").seq)
    assert len(actual) == len(target) and actual in target + target


def test_fragment_palette_uses_backbone_relative_hsv_and_wraps_hue():
    from plasmid_verify.golden_gate_design import _fragment_color, _fragment_palette, _display_color
    assert [_fragment_color(n) for n in range(5)] == ["#ff8080", "#ffbf80", "#ffff80", "#bfff80", "#80ff80"]
    assert _fragment_color(12) == _fragment_color(0)
    assert _fragment_color(13) == _fragment_color(1)
    assert _fragment_palette([{}, {}, {"is_backbone": True}, {}]) == ["#ffff80", "#bfff80", "#ff8080", "#ffbf80"]
    assert _fragment_palette([]) == []
    assert _display_color("#ffbf80") == "#ffbf80"
    assert _display_color("gray - 50") == "#a0a0a0"
    assert _display_color("red") == "#ff0000"


def test_unique_plasmids_keep_complete_linker_all_combinations_and_aliases(tmp_path):
    path, source = _native_linker_map(tmp_path)
    variables = [f for f in plan_from_annotated_snapgene(path)['fragments'] if f.get('annotation_kind') == 'variable']
    texts = {variables[0]['id']: 'one,ATGGCC\nalias,ATGGCC\nlonger,ATGAGCGCC',
             variables[1]['id']: 'tail/1,GCGGCG\ntail?1,GCTGCT'}
    result = run_annotated_snapgene_design(path, tmp_path/'out', variable_texts=texts)
    assert result['plasmidCount'] == 4 and result['combinationCount'] == 6
    assert sorted(p['combinationCount'] for p in result['plasmids']) == [1, 1, 2, 2]
    assert len({p['filename'] for p in result['plasmids']}) == 4
    expected = {source[:28]+a+'GGATCTGGATCTGGA'+b for a in ('ATGGCC','ATGAGCGCC') for b in ('GCGGCG','GCTGCT')}
    from plasmid_verify.golden_gate_design import _circular_identity
    assert {_circular_identity(str(SeqIO.read(p['path'],'snapgene').seq)) for p in result['plasmids']} == {_circular_identity(s) for s in expected}
    with (Path(result['outputDir'])/'plasmid_combinations.tsv').open() as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    assert len(records) == 6
    assert all((Path(result['outputDir'])/r['file']).exists() for r in records)
    repeated = run_annotated_snapgene_design(path, tmp_path/'out', variable_texts=texts)
    assert result['plasmidFolder'] != repeated['plasmidFolder']
    assert all(Path(p['path']).exists() for p in result['plasmids'])


def test_circular_identity_handles_origins_and_reverse_complements():
    from plasmid_verify.golden_gate_design import _circular_identity, reverse_complement
    import itertools
    for chars in itertools.product('ACGT', repeat=5):
        sequence=''.join(chars)
        variants=[s[i:]+s[:i] for s in (sequence, reverse_complement(sequence)) for i in range(len(sequence))]
        assert _circular_identity(sequence) == min(variants)


def test_library_size_limit_is_explicit_not_silent():
    from plasmid_verify.golden_gate_design import _unique_products
    with pytest.raises(ValueError, match='10,201 combinations'):
        _unique_products([{'assembly_position': i} for i in range(2) for _ in range(101)], None)


def test_input_core_and_strand_colors_survive_each_plasmid(tmp_path):
    from plasmid_verify.golden_gate_design import _snapgene_features_xml, _source_strand_colors
    path, sequence = _native_linker_map(tmp_path)
    features = _snapgene_features_xml(path)
    custom = {'[Backbone]':'#112233','{scFv}':'#ab4321','{TM-tail}':'#7788aa'}
    for feature in features:
        if feature.get('name') in custom:
            for segment in feature.findall('Segment'):
                segment.set('color', custom[feature.get('name')])
    colors=Element('StrandColors')
    top=SubElement(colors,'TopStrand')
    SubElement(top,'ColorRange',{'range':'0..33','colors':'#123456'})
    SubElement(top,'ColorRange',{'range':'34..54','colors':'#00cc22'})
    SubElement(SubElement(colors,'BottomStrand'),'ColorRange',{'range':'0..54','colors':'#aa00cc'})
    _write_snapgene_map(path, sequence, [{'xml': f} for f in features], 'Colored input', colors)
    variables=[f for f in plan_from_annotated_snapgene(path)['fragments'] if f.get('annotation_kind')=='variable']
    result=run_annotated_snapgene_design(path,tmp_path/'out',variable_texts={variables[0]['id']:'native,ATGGCC\nlong,ATGAGCGCC'})
    for product in result['plasmids']:
        output=Path(product['path'])
        exported=_snapgene_features_xml(output)
        selected=product['choices'][1]['variant_name']
        feature=next(f for f in exported if f.get('name')=='{'+selected+'}')
        assert feature.find('Segment').get('color')=='#ab4321'
        dna=str(SeqIO.read(output,'snapgene').seq)
        palette=_source_strand_colors(output,len(dna))
        assert set(palette['BottomStrand'])=={'#aa00cc'}
        # Both different-length products retain exact strand colors for linker.
        start=(dna+dna).index('GGATCTGGATCTGGA')
        assert [palette['TopStrand'][(start+i)%len(dna)] for i in range(15)]==['#00cc22']*15
        assert set(palette['TopStrand'])=={'#123456','#00cc22'}
        native_tail=next(f for f in exported if f.get('name')=='{TM-tail}')
        assert native_tail.find('Segment').get('color')=='#7788aa'


def test_nonrepresentative_failure_prevents_all_outputs(tmp_path, monkeypatch):
    import plasmid_verify.golden_gate_design as design
    path, _ = _native_linker_map(tmp_path)
    variable=next(f for f in plan_from_annotated_snapgene(path)['fragments'] if f.get('annotation_kind')=='variable')
    real=design._assembled_circular_map
    def assembled(rows):
        dna, starts, retained=real(rows)
        return (dna+'A' if any(r.get('variant_name')=='second' for r in rows) else dna),starts,retained
    monkeypatch.setattr(design,'_assembled_circular_map',assembled)
    with pytest.raises(ValueError, match='no outputs written'):
        run_annotated_snapgene_design(path,tmp_path/'out',variable_texts={variable['id']:'first,ATGGCC\nsecond,ATGAGC'})
    assert not (tmp_path/'out').exists()


def test_eight_by_six_library_exports_forty_eight_plasmids(tmp_path):
    path, _ = _native_linker_map(tmp_path)
    variables = [f for f in plan_from_annotated_snapgene(path)['fragments'] if f.get('annotation_kind') == 'variable']
    result = run_annotated_snapgene_design(path, tmp_path/'out', variable_texts={
        variables[0]['id']: '\n'.join(f'scFv-{i},ATG'+ 'GCC'*i for i in range(1,9)),
        variables[1]['id']: '\n'.join(f'tail-{i},'+'GCG'*i for i in range(1,7)),
    })
    assert result['plasmidCount'] == result['combinationCount'] == 48
    assert len(list(Path(result['plasmidFolder']).glob('*.dna'))) == 48
    assert len({str(SeqIO.read(p['path'],'snapgene').seq) for p in result['plasmids']}) == 48
    assert len(result['fragments']) == 1+8+6  # Orders are still per fragment, not per plasmid.
