from xml.etree.ElementTree import Element, SubElement
from pathlib import Path
import pytest

from Bio import SeqIO

from plasmid_verify.golden_gate_design import (
    _feature_codons, _visible_codons, _project_annotations, _qualifier,
    _write_snapgene_map, _snapgene_features_xml, run_annotated_snapgene_design,
    _junction_display_window,
    _write_schematic_junction_image, _HighResolutionDraw,
)


def feature(name="ORF 1", spans=("1-9",), reverse=False, phase=1, kind="misc_feature"):
    f = Element("Feature", {"name": name, "type": kind, "directionality": "2" if reverse else "1"})
    for span in spans:
        SubElement(f, "Segment", {"range": span, "color": "#aaccdd"})
    q = SubElement(f, "Q", {"name": "codon_start"})
    SubElement(q, "V", {"int": str(phase)})
    return f


def test_orf_without_snapgene_translation_flag_uses_annotated_phase():
    f = feature(spans=("1-10",), phase=2)
    codons = _feature_codons(f, "AATGGCCTAA")
    assert "".join(c["aa"] for c in codons) == "MA*"
    assert codons[0]["positions"] == [1, 2, 3]
    assert [c["number"] for c in _visible_codons(f, "AATGGCCTAA", 3, 10)] == [2, 3]


def test_junction_context_is_twelve_core_bases_plus_entire_intervening_linker():
    # Four linker bases on the left; all fifteen on the right; shared four
    # count once in the joined view. The 32 bp adapter gap is not product DNA.
    cores = [feature('[upstream]', ('1-30',)), feature('{downstream}', ('82-120',))]
    first, last = _junction_display_window(cores, 200, 34, 32)
    assert first == 18  # last 12 of upstream core 0..30
    assert last == 93   # first 12 of downstream core 81..120
    assert last - first - 32 - 4 == 12 + 15 + 12


def test_junction_context_handles_origin_and_short_cores():
    cores = [feature('[upstream]', ('121-140',)), feature('{downstream}', ('5-10',))]
    first, last = _junction_display_window(cores, 176, 144, 32)
    assert first == 128
    assert last == 186  # full six-base short core after the circular origin


def test_reverse_orf_uses_reverse_complement_and_biological_numbering():
    codons = _feature_codons(feature(reverse=True), "TTAGGCCAT")
    assert "".join(c["aa"] for c in codons) == "MA*"
    assert codons[0]["positions"] == [8, 7, 6]
    assert codons[-1]["number"] == 3


def test_multisegment_cds_carries_frame_across_segment_boundary():
    f = feature(name="protein", spans=("1-4", "5-9"), kind="CDS")
    assert "".join(c["aa"] for c in _feature_codons(f, "ATGGCCTAA")) == "MA*"


def test_explicit_noncoding_segment_is_not_translated():
    f = feature(spans=("1-3", "4-6", "7-9"), kind="CDS")
    for s, value in zip(f.findall('Segment'), ('1', '0', '1')):
        s.set('translated', value)
    assert "".join(c["aa"] for c in _feature_codons(f, "ATGAAATAA")) == "M*"
    assert _feature_codons(feature(name="promoter"), "ATGGCCTAA") == []


def test_partial_projection_preserves_source_frame_and_residue_numbers():
    f = feature()
    projected = _project_annotations([f], "ATGGCCTAA", {p: p + 18 for p in range(2, 9)})[0]['xml']
    assert projected.find('Segment').get('range') == '21-27'
    assert _qualifier(projected, 'codon_start') == '2'
    assert [(c['aa'], c['number'], c['positions']) for c in _feature_codons(projected, '')] == [
        ('A', 2, [21, 22, 23]), ('*', 3, [24, 25, 26])]
    assert _project_annotations([f], "ATGGCCTAA", {p: p for p in range(2, 9)}, require_complete=True) == []


def test_reverse_partial_projection_preserves_phase():
    projected = _project_annotations([feature(reverse=True)], "TTAGGCCAT", {p: p for p in range(7)})[0]['xml']
    assert _qualifier(projected, 'codon_start') == '2'
    assert [(c['aa'], c['number']) for c in _feature_codons(projected, '')] == [('A', 2), ('*', 3)]


def test_origin_spanning_orf_and_translation_table():
    f = feature(spans=("7-3",))
    codons = _feature_codons(f, 'TAACCCATG')
    assert [c['aa'] for c in codons] == ['M', '*']
    assert [c['start'] for c in _visible_codons(f, 'TAACCCATG', 6, 12)] == [6, 9]
    f = feature(spans=("1-3",))
    q = SubElement(f, 'Q', {'name': 'transl_table'})
    SubElement(q, 'V', {'int': '2'})
    assert _feature_codons(f, 'TGA')[0]['aa'] == 'W'


@pytest.mark.parametrize('reverse,phase,replacement,expected', [
    (False, 1, 'ATGAGCGCC', 'MSA'),
    (True, 1, 'TTAGGCCAT', 'MA*'),
    (False, 2, 'AATGAGCGCC', 'MSA'),
])
def test_variant_translation_uses_new_dna_frame_name_and_native_color(tmp_path, reverse, phase, replacement, expected):
    source = tmp_path / 'source.dna'
    sequence = 'CACC' + 'A' * 20 + 'TGAA' + 'ATGGCC' + 'GGATCTGGATCTGGA' + 'GCGGCG'
    protein = feature(name='Original protein', spans=('29-31', '32-34'), reverse=reverse, phase=phase, kind='CDS')
    SubElement(SubElement(protein, 'Q', {'name': 'translation'}), 'V', {'text': 'STALE'})
    _write_snapgene_map(source, sequence, [
        {'name': '[Backbone]', 'start': 0, 'end': 28},
        {'name': 'ori pUC', 'start': 4, 'end': 20},
        {'name': '{Insert}', 'start': 28, 'end': 34},
        {'name': '{Tail}', 'start': 49, 'end': 55},
        {'xml': protein},
        {'xml': feature(name='Tail protein', spans=('50-55',), kind='CDS')},
        {'name': 'Old mutation', 'start': 30, 'end': 31},
        {'name': '-', 'start': 24, 'end': 28},
        {'name': '-', 'start': 34, 'end': 49},
        {'name': '-', 'start': 0, 'end': 4},
    ], 'Synthetic replacement test')
    result = run_annotated_snapgene_design(source, tmp_path / 'out', variable_texts={
        'fragment_2': f'new_insert,{replacement}', 'fragment_3': 'new_tail,GCTGCTGCT'})
    for path in (result['assembledMap'], result['schematicMap']):
        annotations = _snapgene_features_xml(Path(path))
        dna = str(SeqIO.read(path, 'snapgene').seq)
        for name, amino in [('new_insert', expected), ('new_tail', 'AAA')]:
            f = next(f for f in annotations if f.get('name') == name)
            assert ''.join(c['aa'] for c in _feature_codons(f, dna)) == amino
            assert f.find('Segment').get('color') == '#aaccdd'
            assert len(f.findall('Segment')) == 1  # No invented domain boundaries.
            assert 'recalculated' in _qualifier(f, 'note')
        assert not any(f.get('name') in {'Original protein', 'Old mutation', 'Tail protein'} for f in annotations)
        insert = next(f for f in annotations if f.get('name') == 'new_insert')
        assert len(_feature_codons(insert, dna)) == 3


def test_linear_map_wedges_follow_individual_images_and_have_no_callout_headings(tmp_path, monkeypatch):
    from PIL import Image, ImageDraw
    from plasmid_verify.golden_gate_design import _write_plasmid_junction_figure
    images = []
    for i, size in enumerate([(600, 300), (900, 440)]):
        path = tmp_path / f'{i}.png'
        Image.new('RGB', size, 'white').save(path)
        images.append(path)
    rows = [dict(element_name=name, is_backbone=i == 0, left_fusion='CACC', right_fusion='TGAA',
                 left_fusion_owner='core', right_fusion_owner='core', variant_window={})
            for i, name in enumerate(['backbone', 'insert'])]
    texts, polygons = [], []
    original_text, original_polygon = ImageDraw.ImageDraw.text, ImageDraw.ImageDraw.polygon
    def text(self, xy, value, *args, **kwargs):
        texts.append(value)
        return original_text(self, xy, value, *args, **kwargs)
    def polygon(self, points, *args, **kwargs):
        polygons.append(points)
        return original_polygon(self, points, *args, **kwargs)
    monkeypatch.setattr(ImageDraw.ImageDraw, 'text', text)
    monkeypatch.setattr(ImageDraw.ImageDraw, 'polygon', polygon)
    _write_plasmid_junction_figure(tmp_path / 'map.png', rows, [0, 100], ['A' * 104] * 2, 'A' * 200,
        [{'junction_id': 'a', 'selected_fusion': 'CACC'}, {'junction_id': 'b', 'selected_fusion': 'TGAA'}],
        {'a': (96, 0, 1), 'b': (196, 1, 0)}, images)
    assert [p[1][0] - p[0][0] for p in polygons] == [600 - 168, 900 - 168]
    assert not any('→' in t or 'CACC' in t or 'TGAA' in t for t in texts)
    assert {'1', '2'} <= set(texts)  # Central map badges remain.


def test_annotations_are_exported_to_assembled_map_and_orf_images(tmp_path, monkeypatch):
    sequence = 'CACC' + 'A' * 20 + 'TGAA' + 'ATGGCC' + 'GGATCTGGATCTGGA' + 'GCGGCG'
    source = tmp_path / 'source.dna'
    _write_snapgene_map(source, sequence, [
        {'name': '[Backbone]', 'start': 0, 'end': 28},
        {'name': 'pUC ori', 'start': 4, 'end': 20},
        {'name': '{Insert}', 'start': 28, 'end': 34},
        {'name': '{Tail}', 'start': 49, 'end': 55},
        {'xml': feature(name='ORF insert', spans=('29-34',))},
        {'xml': feature(name='ORF linker', spans=('35-49',))},
        {'name': '-', 'start': 24, 'end': 28},
        {'name': '-', 'start': 34, 'end': 49},
        {'name': '-', 'start': 0, 'end': 4},
    ], 'Annotated ORF input')
    result = run_annotated_snapgene_design(source, tmp_path / 'out')
    result_path = Path(result['assembledMap'])
    features = _snapgene_features_xml(result_path)
    assembled = str(SeqIO.read(result_path, 'snapgene').seq)
    orf = next(f for f in features if f.get('name') == 'ORF linker')
    assert ''.join(c['aa'] for c in _feature_codons(orf, assembled)) == 'GSGSG'
    assert all('N' not in ''.join(assembled[p] for p in c['positions']) for c in _feature_codons(orf, assembled))
    assert len(assembled) == len(sequence)
    # Fragment ends keep the actual adapters, omit only the fictitious N gap,
    # and draw two identically oriented, four-base staggered cuts.
    schematic = Path(result['schematicMap'])
    gap = next(f for f in _snapgene_features_xml(schematic) if f.get('name') == '└┐')
    gap_start = int(gap.find('Segment').get('range').split('-')[0]) - 1
    gap_end = int(gap.findall('Segment')[-1].get('range').split('-')[1])
    texts, cuts = [], []
    real_text, real_line = _HighResolutionDraw.text, _HighResolutionDraw.line
    def text(self, xy, value, **kwargs):
        texts.append(value)
        return real_text(self, xy, value, **kwargs)
    def line(self, xy, **kwargs):
        if len(xy) == 4 and isinstance(xy[0], tuple) and xy[0][1] == 7:
            cuts.append(xy)
        return real_line(self, xy, **kwargs)
    monkeypatch.setattr(_HighResolutionDraw, 'text', text)
    monkeypatch.setattr(_HighResolutionDraw, 'line', line)
    _write_schematic_junction_image(tmp_path/'ends.png', schematic, gap_start, gap_end-gap_start, fragment_ends=True)
    assert texts.count('BsmBI') == 2
    assert 'N' not in texts and 'n' not in texts and '//' not in texts
    assert len(cuts) == 2
    assert all(c[2][0] - c[1][0] == 4 * 16 for c in cuts)
    assert all([p[1] for p in c] == [7, 43, 43, 76] for c in cuts)
