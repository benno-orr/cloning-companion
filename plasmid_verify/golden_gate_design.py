"""Scaffold-aware Golden Gate fragment design used by the native app.

Plans can be YAML files or annotated SnapGene maps. Both declare core
fragments, allowed junction regions, and whether each fusion belongs to the
core or its Type IIS adapter. The engine selects non-palindromic,
non-conflicting 4 bp fusions and writes a flat order/assembly package.
"""

from __future__ import annotations

import csv
import copy
import colorsys
import io
import json
import re
import struct
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import Element, SubElement, tostring, fromstring

import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from PIL import Image, ImageDraw, ImageFont

from .fasta import clean_sequence, parse_sequence_file, reverse_complement, to_fasta


ENZYMES = {
    "BsmBI": {"forward_site": "CGTCTC", "reverse_site": "GAGACG"},
    "BsaI": {"forward_site": "GGTCTC", "reverse_site": "GAGACC"},
}

_FIXED_FRAGMENT_LABEL = re.compile(r"^\[([^\[\]]+)\]$")
_VARIABLE_FRAGMENT_LABEL = re.compile(r"^\{([^{}]+)\}$")
_PUC_ORIGIN_LABEL = re.compile(r"^(?:puc\s*ori|ori\s*puc)$", re.IGNORECASE)
_JUNCTION_PLACEMENT_LABEL = re.compile(r"^(?:\((.+)\)|-)$")
_JUNCTION_COLOR = "#e7a93b"


def _fragment_color(n: int) -> str:
    """HSV (30*n degrees, 50%, 100%), wrapping hue after twelve fragments."""
    rgb = colorsys.hsv_to_rgb((30 * n % 360) / 360, .5, 1)
    return "#" + "".join(f"{round(channel * 255):02x}" for channel in rgb)


def _fragment_palette(rows: list[dict[str, Any]]) -> list[str]:
    backbone = next((i for i, row in enumerate(rows) if row.get("is_backbone")), 0)
    return [_fragment_color((i - backbone) % len(rows)) for i in range(len(rows))]


def _display_color(value: str, default: str = "#555555") -> str:
    # Accept custom RGB colors and older named SnapGene strand colors.
    if re.fullmatch(r"#[0-9a-fA-F]{6}", value):
        return value
    return {"red": "#ff0000", "orange": "#ff8800", "green": "#009900", "blue": "#164cff",
            "violet": "#8844cc", "gray - 50": "#a0a0a0"}.get(value, default)


def _dna(value: str) -> str:
    sequence = clean_sequence(value)
    if set(sequence) - set("ACGT"):
        raise ValueError("Golden Gate design inputs must use unambiguous A/C/G/T DNA")
    return sequence


def _range(sequence: str, start: int, end: int) -> str:
    if start < 0 or end < 0 or start >= len(sequence) or end > len(sequence):
        raise ValueError(f"Range [{start}, {end}) is outside a {len(sequence)} bp scaffold")
    return sequence[start:end] if end >= start else sequence[start:] + sequence[:end]


def _feature_name(feature) -> str:
    for key in ("label", "gene", "note"):
        values = feature.qualifiers.get(key, [])
        if values and str(values[0]).strip():
            return str(values[0]).strip()
    return ""


@dataclass
class Scaffold:
    path: Path
    sequence: str
    features: dict[str, list[tuple[int, int]]]

    @classmethod
    def from_path(cls, path: Path) -> "Scaffold":
        if path.suffix.lower() != ".dna":
            record = parse_sequence_file(path.read_bytes(), path.name)[0]
            return cls(path, clean_sequence(record.sequence), {})
        record = SeqIO.read(path, "snapgene")
        features: dict[str, list[tuple[int, int]]] = {}
        for feature in record.features:
            name = _feature_name(feature)
            if not name:
                continue
            parts = getattr(feature.location, "parts", [feature.location])
            spans = [(int(part.start), int(part.end)) for part in parts]
            if spans:
                features.setdefault(name, []).extend(spans)
        return cls(path, clean_sequence(str(record.seq)), features)

    def feature(self, name: str) -> tuple[str, list[tuple[int, int]]]:
        spans = self.features.get(name)
        if not spans:
            available = ", ".join(sorted(self.features)[:30])
            raise ValueError(f"Scaffold feature {name!r} not found. Available features include: {available}")
        return "".join(_range(self.sequence, start, end) for start, end in spans), spans


@dataclass(frozen=True)
class Candidate:
    junction_id: str
    fusion: str
    offset: int
    score: float
    source: str


def _score(sequence: str, preferred: str | None) -> float:
    score = 1000.0 if preferred and sequence == _dna(preferred) else 0.0
    gc = (sequence.count("G") + sequence.count("C")) / len(sequence)
    score += 10 - abs(gc - 0.5) * 10 + len(set(sequence)) * 1.5
    if len(set(sequence)) == 1:
        score -= 25
    if sequence[:2] == sequence[2:]:
        score -= 4
    return score


def _region(spec: dict[str, Any], scaffold: Scaffold) -> tuple[str, str]:
    choices = [key for key in ("sequence", "scaffold_feature", "scaffold_range") if key in spec]
    if len(choices) != 1:
        raise ValueError(f"Junction {spec.get('id')!r} needs exactly one candidate region source")
    if "sequence" in spec:
        return _dna(str(spec["sequence"])), "inline sequence"
    if "scaffold_feature" in spec:
        sequence, spans = scaffold.feature(str(spec["scaffold_feature"]))
        if len(spans) != 1:
            raise ValueError(f"Junction {spec.get('id')!r} feature must be one contiguous highlighted region")
        return sequence, f"feature {spec['scaffold_feature']} ({spans[0][0]}..{spans[0][1]})"
    start, end = spec["scaffold_range"]
    label = str(spec.get("label", "")).strip()
    source = f"{label} (scaffold range {start}..{end})" if label else f"scaffold range {start}..{end}"
    return _range(scaffold.sequence, int(start), int(end)), source


def _select_junctions(plan: dict[str, Any], scaffold: Scaffold) -> tuple[dict[str, Candidate], list[dict[str, Any]]]:
    settings = plan.get("selection", {})
    if int(settings.get("fusion_length", 4)) != 4:
        raise ValueError("The desktop designer currently supports 4 bp fusion overhangs")
    avoid = {_dna(item) for item in settings.get("avoid", [])}
    reject_palindromes = bool(settings.get("reject_palindromes", True))
    reject_rc = bool(settings.get("reject_reverse_complement_pairs", True))
    options: dict[str, list[Candidate]] = {}
    rows: list[dict[str, Any]] = []
    for spec in plan.get("junction_candidates", []):
        junction_id = str(spec["id"])
        region, source = _region(spec, scaffold)
        if len(region) < 4:
            raise ValueError(f"Junction {junction_id!r} has fewer than 4 allowed bases")
        candidates = []
        for offset in range(len(region) - 3):
            fusion = region[offset : offset + 4]
            reason = "accepted"
            if set(fusion) - set("ACGT"):
                reason = "non-DNA"
            elif fusion in avoid:
                reason = "user avoid list"
            elif reject_palindromes and fusion == reverse_complement(fusion):
                reason = "palindromic"
            score = _score(fusion, spec.get("preferred")) if reason == "accepted" else None
            rows.append({"junction_id": junction_id, "fusion": fusion, "offset_in_allowed_region": offset,
                         "accepted": reason == "accepted", "reason": reason, "score": "" if score is None else f"{score:.2f}",
                         "candidate_region_source": source})
            if score is not None:
                candidates.append(Candidate(junction_id, fusion, offset, score, source))
        candidates.sort(key=lambda item: (-item.score, item.offset, item.fusion))
        if not candidates:
            raise ValueError(f"No usable 4 bp fusion remains for junction {junction_id!r}")
        options[junction_id] = candidates
    selected: dict[str, Candidate] = {}
    used: set[str] = set()
    for junction_id in sorted(options, key=lambda key: (len(options[key]), key)):
        pick = next((item for item in options[junction_id]
                     if item.fusion not in used and (not reject_rc or reverse_complement(item.fusion) not in used)), None)
        if pick is None:
            raise ValueError(f"No non-conflicting fusion remains for junction {junction_id!r}")
        selected[junction_id] = pick
        used.add(pick.fusion)
    return selected, rows


def _read_fasta(path: Path, record_name: str | None) -> str:
    records = parse_sequence_file(path.read_bytes(), path.name)
    if record_name is None:
        if len(records) != 1:
            raise ValueError(f"{path.name} contains multiple records; specify core.record")
        return _dna(records[0].sequence)
    for record in records:
        if record.name == record_name:
            return _dna(record.sequence)
    raise ValueError(f"FASTA record {record_name!r} was not found in {path.name}")


def _core(spec: dict[str, Any], scaffold: Scaffold, plan_dir: Path) -> str:
    choices = [key for key in ("sequence", "fasta", "scaffold_feature", "scaffold_range") if key in spec]
    if len(choices) != 1:
        raise ValueError("Each fragment core needs exactly one of sequence, fasta, scaffold_feature, scaffold_range")
    if "sequence" in spec:
        return _dna(str(spec["sequence"]))
    if "fasta" in spec:
        path = Path(str(spec["fasta"])).expanduser()
        return _read_fasta(path if path.is_absolute() else plan_dir / path, spec.get("record"))
    if "scaffold_feature" in spec:
        return _dna(scaffold.feature(str(spec["scaffold_feature"]))[0])
    start, end = spec["scaffold_range"]
    return _dna(_range(scaffold.sequence, int(start), int(end)))


def parse_variable_elements(text: str, element_name: str = "variable fragment") -> list[dict[str, str]]:
    """Read a named core DNA list pasted from a sheet or as FASTA."""
    text = text.strip().lstrip("\ufeff")
    if not text:
        return []
    if text.startswith(">"):
        rows = []
        for line in text.splitlines():
            if line.startswith(">"):
                rows.append({"name": line[1:].strip(), "sequence": ""})
            elif line.strip():
                rows[-1]["sequence"] += line.strip()
    else:
        delimiter = "\t" if "\t" in text.splitlines()[0] else ","
        table = list(csv.reader(io.StringIO(text), delimiter=delimiter))
        headers = [value.strip().lower().replace(" ", "_").replace("-", "_") for value in table[0]]
        name_keys = {"name", "variant", "variant_name", "element", "id"}
        dna_keys = {"sequence", "dna", "dna_sequence", "core_sequence", "dna_seq"}
        named_header = next((index for index, value in enumerate(headers) if value in name_keys), None)
        dna_header = next((index for index, value in enumerate(headers) if value in dna_keys), None)
        if named_header is not None and dna_header is not None:
            name_column, dna_column, table = named_header, dna_header, table[1:]
        else:
            name_column, dna_column = 0, 1
        rows = []
        for index, values in enumerate(table, 1):
            if not any(value.strip() for value in values):
                continue
            if len(values) <= max(name_column, dna_column):
                raise ValueError(f"{element_name}: row {index} needs a variant name and core DNA sequence")
            rows.append({"name": values[name_column].strip(), "sequence": values[dna_column]})
    seen = set()
    for row in rows:
        name = row["name"].strip()
        if not name:
            raise ValueError(f"{element_name}: every variant needs a name")
        if name.casefold() in seen:
            raise ValueError(f"{element_name}: duplicate variant name {name!r}")
        seen.add(name.casefold())
        try:
            row["sequence"] = _dna(row["sequence"])
        except ValueError as exc:
            raise ValueError(f"{element_name} / {name}: {exc}") from exc
        row["name"] = name
    if not rows:
        raise ValueError(f"{element_name}: no named DNA variants found")
    return rows


def _expanded_fragment_specs(plan: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand options within their assembly slot, without duplicating fixed parts."""
    expanded = []
    used_ids = set()
    for position, fragment in enumerate(plan.get("fragments", []), 1):
        element_id = str(fragment["id"])
        variants = fragment.get("variants")
        if variants is not None and (not isinstance(variants, list) or not variants):
            raise ValueError(f"{element_id}: variants must be a nonempty list")
        seen_names = set()
        for option_index, variant in enumerate(variants or [None], 1):
            spec = dict(fragment)
            spec.update(element_id=element_id, element_name=fragment.get("name", element_id),
                        assembly_position=position, option_index=option_index, variant_name="")
            if variant is not None:
                name = str(variant.get("name", "")).strip()
                if not name or name.casefold() in seen_names:
                    raise ValueError(f"{element_id}: variant names must be nonempty and unique")
                seen_names.add(name.casefold())
                slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "variant"
                variant_core = _dna(str(variant["sequence"]))
                window = fragment.get("variant_window")
                if window:
                    variant_core = window.get("fixed_core_prefix", "") + variant_core + window.get("fixed_core_suffix", "")
                    trim_left, trim_right = window["trim_left"], window["trim_right"]
                    if len(variant_core) <= trim_left + trim_right:
                        raise ValueError(f"{element_id} / {name}: variant is shorter than its junction trims")
                    if ((trim_left and variant_core[:trim_left] != window["discarded_prefix"]) or
                            (trim_right and variant_core[-trim_right:] != window["discarded_suffix"])):
                        raise ValueError(f"{element_id} / {name}: variant changes core bases outside its selected fragment window; choose a wider junction placement")
                    variant_core = window["prefix"] + variant_core[trim_left:len(variant_core) - trim_right if trim_right else None] + window["suffix"]
                spec.update(id=f"{element_id}__{option_index:02d}_{slug}", name=name,
                            variant_name=name, core={"sequence": variant_core})
            if str(spec["id"]) in used_ids:
                raise ValueError(f"Duplicate fragment ID: {spec['id']}")
            used_ids.add(str(spec["id"]))
            expanded.append(spec)
    return expanded


def _tsv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    fields = fields or (list(rows[0]) if rows else [])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def _packet(packet_type: int, payload: bytes) -> bytes:
    return struct.pack(">BI", packet_type, len(payload)) + payload


def _write_snapgene_map(path: Path, sequence: str, features: list[dict[str, Any]], title: str, strand_colors: Element | None = None) -> None:
    """Write a compact circular SnapGene map with colored, named features."""
    root = Element("Features", {"nextValidID": str(len(features) + 1)})
    for index, feature in enumerate(features, 1):
        if "xml" in feature:
            item = copy.deepcopy(feature["xml"])
            item.set("recentID", str(index))
            root.append(item)
            continue
        item = SubElement(root, "Feature", {
            "recentID": str(index), "name": str(feature["name"]), "directionality": str(feature.get("directionality", "0")),
            "type": str(feature.get("type", "misc_feature")), "swappedSegmentNumbering": "1",
            "allowSegmentOverlaps": "0", "readingFrame": "-1", "consecutiveTranslationNumbering": "1",
        })
        SubElement(item, "Segment", {
            "range": f"{int(feature['start']) + 1}-{int(feature['end'])}",
            "color": str(feature.get("color", "#808080")), "type": "standard",
        })
        qualifier = SubElement(item, "Q", {"name": "note"})
        SubElement(qualifier, "V", {"text": str(feature.get("note", ""))})
        if feature.get("owner_color"):
            qualifier = SubElement(item, "Q", {"name": "fragment_color"})
            SubElement(qualifier, "V", {"text": feature["owner_color"]})
        if feature.get("assigned_extensions"):
            qualifier = SubElement(item, "Q", {"name": "assigned_extensions"})
            SubElement(qualifier, "V", {"text": json.dumps(feature["assigned_extensions"])})
    notes = Element("Notes")
    for name, value in (("UUID", str(uuid.uuid4())), ("Type", "Synthetic"), ("ConfirmedExperimentally", "0"),
                        ("Comments", title), ("SequenceClass", "UNA"), ("TransformedInto", "unspecified")):
        SubElement(notes, name).text = value
    cookie = b"SnapGene" + struct.pack(">HHH", 1, 15, 20)
    path.write_bytes(b"".join((
        _packet(0x09, cookie), _packet(0x00, b"\x01" + sequence.encode("ascii")),
        _packet(0x0A, b'<?xml version="1.0"?>' + tostring(root, encoding="utf-8")), _packet(0x06, tostring(notes, encoding="utf-8")),
        _packet(0x14, tostring(strand_colors, encoding="utf-8")) if strand_colors is not None else b"",
    )))


def _snapgene_features_xml(path: Path) -> list[Element]:
    if path.suffix.lower() != ".dna":
        return []
    data, offset = path.read_bytes(), 0
    while offset < len(data):
        if offset + 5 > len(data):
            raise ValueError("Truncated SnapGene packet header")
        packet_type, size = struct.unpack(">BI", data[offset:offset + 5])
        offset += 5
        if offset + size > len(data):
            raise ValueError("Truncated SnapGene packet payload")
        payload = data[offset:offset + size]
        offset += size
        if packet_type == 0x0A:
            return list(fromstring(payload).findall("Feature"))
    return []


def _qualifier(feature: Element, name: str, default: str = "") -> str:
    value = feature.find(f"Q[@name='{name}']/V")
    return next((value.get(key) for key in ("text", "int", "predef") if value.get(key) is not None), default) if value is not None else default


def _coding_segment(feature: Element, segment: Element) -> bool:
    """Recognize explicit translation tracks and annotated CDS/ORF features.

    This does not find new ORFs or infer coding from arbitrary feature names.
    Explicit per-segment flags take precedence over the feature type.
    """
    if any(s.get("translated") is not None for s in feature.findall("Segment")):
        return segment.get("translated") == "1"
    return (feature.get("type", "").casefold() in {"cds", "orf", "open_reading_frame"} or
            bool(re.search(r"\borf(?:\b|\d)", feature.get("name", ""), re.I))) and segment.get("type") != "gap"


def _segment_positions(segment: Element, length: int) -> list[int]:
    first, last = map(int, segment.get("range").split("-"))
    return list(range(first - 1, last)) if last >= first else list(range(first - 1, length)) + list(range(last))


def _feature_codons(feature: Element, sequence: str) -> list[dict[str, Any]]:
    """Return annotated codons in biological order, with source residue numbers."""
    saved = _qualifier(feature, "cloning_companion_codons")
    if saved:
        return json.loads(saved)
    positions = [p for segment in feature.findall("Segment") if _coding_segment(feature, segment)
                 for p in _segment_positions(segment, len(sequence))]
    reverse = feature.get("directionality") == "2"
    if reverse:
        positions.reverse()
    phase = int(_qualifier(feature, "codon_start", "1")) - 1
    if phase not in (0, 1, 2):
        raise ValueError(f"{feature.get('name')}: codon_start must be 1, 2, or 3")
    table = int(_qualifier(feature, "transl_table", "1"))
    codons = []
    for offset in range(phase, len(positions) - 2, 3):
        bases = positions[offset:offset + 3]
        dna = "".join(sequence[p] for p in bases)
        if reverse:
            dna = str(Seq(dna).complement())
        codons.append({"positions": bases, "aa": str(Seq(dna).translate(table=table)), "number": (offset - phase) // 3 + 1})
    return codons


def _visible_codons(feature: Element, sequence: str, first: int, last: int) -> list[dict[str, Any]]:
    result = []
    for codon in _feature_codons(feature, sequence):
        positions = [(p - first) % len(sequence) + first for p in codon["positions"]]
        if max(positions) < last and max(positions) - min(positions) == 2:
            result.append({**codon, "start": min(positions), "end": max(positions) + 1})
    return result


def _native_row_mapping(row: dict[str, Any], scaffold: Scaffold, target_start: int) -> dict[int, int]:
    """Map unchanged native bases only; never reuse proteins from changed cores."""
    spec = row.get("core_source", {})
    if row.get("native_window") or "scaffold_range" in spec:
        start, end = row.get("native_window") or spec["scaffold_range"]
        spans = [(start, end)] if end >= start else [(start, len(scaffold.sequence)), (0, end)]
    elif "scaffold_feature" in spec:
        _, spans = scaffold.feature(spec["scaffold_feature"])
    else:
        return {}
    positions = [p for a, b in spans for p in range(a, b)]
    original = "".join(scaffold.sequence[p] for p in positions)
    if original != row["core_sequence"]:
        window = row.get("variant_window") or {}
        prefix, suffix = len(window.get("prefix", "")), len(window.get("suffix", ""))
        return {**{positions[i]: target_start + i for i in range(prefix)},
                **{positions[len(positions) - suffix + i]: target_start + len(row["core_sequence"]) - suffix + i for i in range(suffix)}}
    return {p: target_start + i for i, p in enumerate(positions)}


def _project_annotations(native_features: list[Element], sequence: str, mapping: dict[int, int], require_complete: bool = False) -> list[dict[str, Any]]:
    """Project native features and complete codons without translating adapters.

    Clipped schematic features keep original residue numbering. Assembled maps
    require complete annotations, so pasted replacement cores cannot inherit a
    stale protein annotation from the source plasmid.
    """
    mapped = []
    for feature in native_features:
        name = feature.get("name", "")
        if _FIXED_FRAGMENT_LABEL.fullmatch(name) or _VARIABLE_FRAGMENT_LABEL.fullmatch(name) or _JUNCTION_PLACEMENT_LABEL.fullmatch(name):
            continue
        segments = []
        complete = True
        for segment in feature.findall("Segment"):
            original_positions = _segment_positions(segment, len(sequence))
            complete = complete and all(p in mapping for p in original_positions)
            runs = []
            for p in original_positions:
                if p not in mapping:
                    continue
                target = mapping[p]
                if runs and target == runs[-1][-1] + 1:
                    runs[-1].append(target)
                else:
                    runs.append([target])
            for run in runs:
                target = copy.deepcopy(segment)
                target.set("range", f"{run[0] + 1}-{run[-1] + 1}")
                if _coding_segment(feature, segment):
                    target.set("translated", "1")
                segments.append(target)
        if not segments or (require_complete and not complete):
            continue
        clone = copy.deepcopy(feature)
        # Coordinate-dependent run-on/cleavage hints must not retain old positions.
        for attr in ("maxRunOn", "maxFusedRunOn", "cleavageArrows"):
            clone.attrib.pop(attr, None)
        for segment in clone.findall("Segment"):
            clone.remove(segment)
        for index, segment in enumerate(segments):
            clone.insert(index, segment)
        codons = [{**codon, "positions": [mapping[p] for p in codon["positions"]]}
                  for codon in _feature_codons(feature, sequence) if all(p in mapping for p in codon["positions"])]
        if any(_coding_segment(feature, s) for s in feature.findall("Segment")):
            # Store positions as well as AA so clipping never restarts a frame
            # or numbers residues from the edge of the displayed window.
            for q in list(clone.findall("Q")):
                if q.get("name") in {"translation", "cloning_companion_codons"}:
                    clone.remove(q)
            q = SubElement(clone, "Q", {"name": "cloning_companion_codons"})
            SubElement(q, "V", {"text": json.dumps(codons)})
            # SnapGene itself also needs the clipped 5′ phase, independently
            # of the original-frame codon coordinates used by our images.
            if codons:
                translated_positions = [p for s in clone.findall("Segment") if _coding_segment(clone, s)
                                        for p in _segment_positions(s, max(mapping.values()) + 1)]
                if clone.get("directionality") == "2":
                    translated_positions.reverse()
                phase = translated_positions.index(codons[0]["positions"][0]) % 3
                for q in list(clone.findall("Q")):
                    if q.get("name") == "codon_start":
                        clone.remove(q)
                q = SubElement(clone, "Q", {"name": "codon_start"})
                SubElement(q, "V", {"int": str(phase + 1)})
        mapped.append({"xml": clone})
    return mapped


def _project_core_annotations(row: dict[str, Any], scaffold: Scaffold, core_start: int, native_features: list[Element]) -> list[dict[str, Any]]:
    return _project_annotations(native_features, scaffold.sequence, _native_row_mapping(row, scaffold, core_start))


def _fragment_core_bounds(row: dict[str, Any], retained_length: int) -> tuple[int, int]:
    """Intrinsic annotation bounds within a retained physical fragment."""
    window = row.get("variant_window") or {}
    left = len(row["left_fusion"]) if row["left_fusion_owner"] == "extension" else 0
    right = len(row["right_fusion"]) if row["right_fusion_owner"] == "extension" else 0
    return left + len(window.get("prefix", "")), retained_length - right - len(window.get("suffix", ""))


def _fragment_extensions(row: dict[str, Any], retained_length: int) -> list[tuple[int, int]]:
    """Non-core display spans, assigning the shared four bases downstream once.

    Fixed junction bases already inside the intrinsic annotation are core, not
    extensions. This affects graphics only, never the assembly or order DNA.
    """
    first, last = _fragment_core_bounds(row, retained_length)
    end = retained_length - 4
    return [(a, b) for a, b in ((0, min(first, end)), (max(0, last), end)) if a < b]


def _draw_hatched_region(draw, bounds, color, background="white", spacing=7):
    """Draw clipped diagonal //// strokes without bleeding into adjacent core."""
    left, top, right, bottom = bounds
    if right <= left or bottom <= top:
        return
    draw.rectangle(bounds, fill=background)
    height = bottom - top
    origin = left - height
    while origin < right:
        a, b = max(left, origin), min(right, origin + height)
        if a < b:
            draw.line((a, bottom - (a - origin), b, bottom - (b - origin)), fill=color, width=1)
        origin += spacing
    draw.rectangle(bounds, fill=None, outline=color)


def _write_assembly_schematic(path: Path, rows: list[dict[str, Any]], all_rows: list[dict[str, Any]], scaffold: Scaffold, enzyme: dict[str, Any], padding: str, spacer: str) -> dict[str, int]:
    """Write a planning map: retained pieces separated by outward Type IIS sites.

    The N gaps are diagram spacers, never synthesis/order DNA. Both strands use
    independent colors so only the exposed staggered fusion is fragment-colored.
    """
    gap_left = spacer + enzyme["reverse_site"] + padding
    gap_right = padding + enzyme["forward_site"] + spacer
    gap = gap_left + "N" * 6 + gap_right
    sequence = ""
    features = []
    strand_colors = Element("StrandColors")
    top, bottom = SubElement(strand_colors, "TopStrand"), SubElement(strand_colors, "BottomStrand")
    native_features = _snapgene_features_xml(scaffold.path)
    palette = _fragment_palette(rows)
    junction_positions = {}

    # Native core membership is independent of which physical fragment carries
    # a copied overlap. Unlabelled connecting DNA belongs to the next core.
    owners = [None] * len(scaffold.sequence)
    for owner_index, row in enumerate(rows):
        if row.get("native_core_feature"):
            _, spans = scaffold.feature(row["native_core_feature"])
            for a, b in spans:
                owners[a:b] = [palette[owner_index]] * (b - a)
    if any(owners):
        next_owner = next(value for value in owners if value is not None)
        for position in range(len(owners) * 2 - 1, -1, -1):
            offset = position % len(owners)
            if owners[offset] is not None:
                next_owner = owners[offset]
            elif position < len(owners):
                owners[offset] = next_owner

    def retained_colors(row, index, left, right):
        current = palette[index]
        native_window = row.get("native_window")
        if native_window:
            a, b = native_window
            original = [owners[(a + n) % len(owners)] for n in range((b - a) % len(owners))]
            if row.get("variant_name"):
                window = row["variant_window"]
                prefix, suffix = len(window["prefix"]), len(window["suffix"])
                return original[:prefix] + [current] * (len(row["core_sequence"]) - prefix - suffix) + (original[-suffix:] if suffix else [])
            return original
        return ([palette[index - 1]] * len(left) + [current] * len(row["core_sequence"]) +
                [palette[(index + 1) % len(rows)]] * len(right))

    def color(parent, start, end, value):
        if end > start:
            SubElement(parent, "ColorRange", {"range": f"{start}..{end - 1}", "colors": value})

    def color_bases(parent, start, values):
        run = 0
        for end in range(1, len(values) + 1):
            if end == len(values) or values[end] != values[run]:
                color(parent, start + run, start + end, values[run])
                run = end

    for index, row in enumerate(rows):
        left = row["left_fusion"] if row["left_fusion_owner"] == "extension" else ""
        right = row["right_fusion"] if row["right_fusion_owner"] == "extension" else ""
        retained = left + row["core_sequence"] + right
        start, end = len(sequence), len(sequence) + len(retained)
        sequence += retained
        kind = row.get("annotation_kind")
        name = row["element_name"]
        label = f"{{{name}}}" if kind == "variable" or row["variant_name"] else f"[{name}]"
        options = [item["name"] for item in all_rows if item["element_id"] == row["element_id"]]
        # A physical insert may carry native linker/flanking DNA without that
        # DNA belonging to its labelled core. Keep feature bounds independent
        # of the strand provenance colors below (also for pasted variants).
        core_start, core_end = _fragment_core_bounds(row, len(retained))
        feature_start, feature_end = start + core_start, start + core_end
        features.append({"name": label, "start": feature_start, "end": feature_end,
                         "owner_color": palette[index],
                         "assigned_extensions": [(start + a, start + b) for a, b in _fragment_extensions(row, len(retained))],
                         "color": palette[index],
                         "note": f"Shown: {row['name']}. Options: {', '.join(options)}. Planning schematic, not order DNA."})
        features.extend(_project_core_annotations(row, scaffold, start + len(left), native_features))
        provenance = retained_colors(row, index, left, right)
        if len(provenance) != len(retained) or any(value is None for value in provenance):
            raise ValueError("Fragment base provenance does not match its retained sequence")
        color_bases(top, start, provenance[:4])
        color_bases(top, start + 4, provenance[4:-4])
        color(top, end - 4, end, "gray - 50")
        color(bottom, start, start + 4, "gray - 50")
        color_bases(bottom, start + 4, provenance[4:-4])
        color_bases(bottom, end - 4, provenance[-4:])
        for fusion, a, b in ((row["left_fusion"], start, start + 4), (row["right_fusion"], end - 4, end)):
            features.append({"name": fusion, "start": a, "end": b, "color": "#bfbfbf", "note": "Selected 4 bp Golden Gate fusion"})
        gap_start = len(sequence)
        sequence += gap
        color(top, gap_start, len(sequence), "gray - 50")
        color(bottom, gap_start, len(sequence), "gray - 50")
        junction_positions[str(row["right_junction"])] = gap_start
        composite = Element("Feature", {"name": "└┐", "type": "misc_feature", "allowSegmentOverlaps": "0", "consecutiveTranslationNumbering": "1"})
        for segment_name, a, b in ((f"←{enzyme.get('name', 'Type IIS')}", gap_start, gap_start + len(gap_left)),
                                   ("//", gap_start + len(gap_left), gap_start + len(gap_left) + 6),
                                   (f"{enzyme.get('name', 'Type IIS')}→", gap_start + len(gap_left) + 6, len(sequence))):
            SubElement(composite, "Segment", {"name": segment_name, "range": f"{a + 1}-{b}", "color": "noColor", "type": "standard"})
        features.append({"xml": composite})
    _write_snapgene_map(path, sequence, features, "Planning schematic only: duplicated fusions and N gaps are NOT an assembled plasmid or order sequence.", strand_colors)
    if str(SeqIO.read(path, "snapgene").seq).upper() != sequence:
        raise ValueError("Generated schematic did not round-trip correctly")
    return junction_positions


def _font(size: int):
    for candidate in ("/System/Library/Fonts/Menlo.ttc", "/System/Library/Fonts/SFNSMono.ttf"):
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


class _HighResolutionDraw:
    """Draw in logical pixels while retaining crisp, twice-size export text."""

    def __init__(self, image, scale=2):
        self.draw, self.scale = ImageDraw.Draw(image), scale

    def _points(self, points):
        if points and isinstance(points[0], (tuple, list)):
            return [tuple(v * self.scale for v in point) for point in points]
        return tuple(v * self.scale for v in points)

    def _font(self, font):
        return font.font_variant(size=font.size * self.scale) if hasattr(font, "font_variant") else font

    def textbbox(self, xy, text, font):
        return tuple(v / self.scale for v in self.draw.textbbox(self._points(xy), text, font=self._font(font)))

    def text(self, xy, text, font, fill):
        self.draw.text(self._points(xy), text, font=self._font(font), fill=fill)

    def line(self, xy, fill, width=1):
        self.draw.line(self._points(xy), fill=fill, width=width * self.scale)

    def rectangle(self, xy, fill, outline=None):
        self.draw.rectangle(self._points(xy), fill=fill, outline=outline, width=self.scale)

    def polygon(self, xy, fill, outline=None):
        self.draw.polygon(self._points(xy), fill=fill, outline=outline, width=self.scale)


def _junction_display_window(features: list[Element], length: int, gap_start: int, gap_length: int, core_bases: int = 12) -> tuple[int, int]:
    """Frame the actual flanking core annotations, not a fixed cut-site flank."""
    spans = []
    for feature in features:
        name = feature.get("name", "")
        if not (_FIXED_FRAGMENT_LABEL.fullmatch(name) or _VARIABLE_FRAGMENT_LABEL.fullmatch(name)):
            continue
        for segment in feature.findall("Segment"):
            a, b = map(int, segment.get("range").split("-"))
            spans.extend((a - 1 + shift, b + shift) for shift in (-length, 0, length))
    if not spans:
        return gap_start - core_bases, gap_start + gap_length + core_bases
    left = max((s for s in spans if s[1] <= gap_start), key=lambda s: s[1])
    right = min((s for s in spans if s[0] >= gap_start + gap_length), key=lambda s: s[0])
    return min(max(left[0], left[1] - core_bases), gap_start - 4), max(min(right[1], right[0] + core_bases), gap_start + gap_length + 4)


def _write_schematic_junction_image(path: Path, map_path: Path, gap_start: int, gap_length: int, clean: bool = False, joined: bool = False, fragment_ends: bool = False) -> None:
    """Reference-style sequence view, rendered from the exported map itself.

    This is a diagram of pre-assembly ends, not an order sequence or an
    assembled sequence. The separate *_assembled.png shows the latter.
    """
    clean = clean or joined
    sequence = str(SeqIO.read(map_path, "snapgene").seq)
    length = len(sequence)
    features = _snapgene_features_xml(map_path)
    first, last = _junction_display_window(features, length, gap_start, gap_length)
    windows = [(first, gap_start), (gap_start + gap_length, last)] if clean else [(first, last)]
    spacer_start, spacer_end = gap_start + (gap_length - 6) // 2, gap_start + (gap_length + 6) // 2
    if fragment_ends:
        windows = [(first, spacer_start), (spacer_end, last)]
    data, offset, strands = map_path.read_bytes(), 0, None
    while offset < len(data):
        kind, size = struct.unpack(">BI", data[offset:offset + 5])
        if kind == 0x14:
            strands = fromstring(data[offset + 5:offset + 5 + size])
        offset += 5 + size
    colors = {}
    for strand in ("TopStrand", "BottomStrand"):
        values = ["#a0a0a0"] * length
        for item in strands.findall(f"{strand}/ColorRange") if strands is not None else []:
            a, b = map(int, item.get("range").split(".."))
            values[a:b + 1] = [_display_color(item.get("colors", ""))] * (b - a + 1)
        colors[strand] = values
    def projected(position):
        return position - gap_length - 4 if joined and position >= gap_start + gap_length else position
    visible = []
    seen = set()
    for feature in features:
        segments = []
        aa_offset = 0
        for segment in feature.findall("Segment"):
            a, b = map(int, segment.get("range").split("-"))
            a -= 1
            for shift in (-length, 0, length):
                if any(a + shift < end and b + shift > start for start, end in windows):
                    segments.append((segment, a + shift, b + shift, aa_offset))
            aa_offset += b - a
        if segments:
            key = (feature.get("name"), tuple((projected(a), projected(b)) for _, a, b, _ in segments))
            if joined and key in seen:
                continue
            seen.add(key)
            visible.append((feature, segments))
    if fragment_ends:
        # Don't repeat a clipped annotation solely on a copied fusion when
        # its larger, intrinsic feature is already visible on the other end.
        visible = [(f, ss) for f, ss in visible if not (
            f.get("name") != "└┐" and not re.fullmatch("[ACGT]{4}", f.get("name", "")) and
            all((gap_start - 4 <= a < b <= gap_start) or
                (gap_start + gap_length <= a < b <= gap_start + gap_length + 4) for _, a, b, _ in ss) and
            any(other is not f and other.get("name") == f.get("name") and
                sum(d - c for _, c, d, _ in other_ss) > sum(b - a for _, a, b, _ in ss)
                for other, other_ss in visible))]
    if joined:
        # A feature carried on both sides of a fusion should not create two
        # annotation tracks when one clipped copy is contained in the other.
        filtered = []
        for feature, segments in visible:
            extents = [(projected(a), projected(b)) for _, a, b, _ in segments]
            if any(other.get("name") == feature.get("name") and
                   all(any(c <= a and b <= d for _, c0, d0, _ in other_segments
                           for c, d in [(projected(c0), projected(d0))]) for a, b in extents)
                   for other, other_segments in filtered):
                continue
            filtered = [(other, other_segments) for other, other_segments in filtered
                        if not (other.get("name") == feature.get("name") and
                                all(any(a <= projected(c) and projected(d) <= b for a, b in extents)
                                    for _, c, d, _ in other_segments))]
            filtered.append((feature, segments))
        visible = filtered
    coding_columns = {projected(p) for feature, segments in visible for segment, a, b, _ in segments
                      if _coding_segment(feature, segment) for p in range(a, b)}
    recognition_columns = set()
    if fragment_ends:
        for feature, segments in visible:
            if feature.get("name") != "└┐":
                continue
            for segment, a, b, _ in segments:
                label = segment.get("name", "")
                recognition = ENZYMES.get(label.strip("←→"), {}).get("reverse_site" if label.startswith("←") else "forward_site", "")
                offset = "".join(sequence[p % length] for p in range(a, b)).find(recognition) if recognition else -1
                if offset >= 0:
                    recognition_columns.update(range(a + offset, a + offset + len(recognition)))
    # Reserve annotation lanes using visible extents; labels remain separate
    # from fragment bars and from nearby protein translations.
    lanes, protein = [], []
    for feature, segments in visible:
        if feature.get("name") == "└┐" or _FIXED_FRAGMENT_LABEL.fullmatch(feature.get("name", "")) or _VARIABLE_FRAGMENT_LABEL.fullmatch(feature.get("name", "")) or re.fullmatch("[ACGT]{4}", feature.get("name", "")):
            continue
        a, b = min(a for _, a, _, _ in segments), max(b for _, _, b, _ in segments)
        a, b = max(first, a), min(last, b)
        a, b = projected(a), projected(b)
        lane = next((i for i, spans in enumerate(lanes) if all(b <= s or a >= e for s, e in spans)), len(lanes))
        if lane == len(lanes):
            lanes.append([])
        lanes[lane].append((a, b))
        protein.append((feature, segments, lane))
    cell, margin = 16, 42
    open_gap = 90 if fragment_ends else 240
    width = (last - first - gap_length) * cell + open_gap + margin * 2 if clean else (last - first) * cell + margin * 2
    if joined:
        width = (last - first - gap_length - 4) * cell + margin * 2
    if fragment_ends:
        width = (last - first - 6) * cell + open_gap + margin * 2
    protein_y = 184 if fragment_ends else 147
    height = protein_y + 6 + 72 * max(1, len(lanes))
    scale = 2
    image = Image.new("RGB", (width * scale, height * scale), "#ffffff" if clean else "#f2f2f2")
    draw = _HighResolutionDraw(image, scale)
    dna_font, small = _font(17), _font(12)
    def x(position):
        if joined:
            return margin + (projected(position) - first) * cell
        if fragment_ends and position >= spacer_end:
            return margin + (position - first - 6) * cell + open_gap
        if clean and position >= gap_start + gap_length:
            return margin + (gap_start - first) * cell + open_gap + (position - gap_start - gap_length) * cell
        return margin + (position - first) * cell
    def centered(label, a, b, y, font=small, fill="#111111"):
        span = draw.textbbox((0, 0), label, font=font)
        draw.text(((x(a) + x(b) - (span[2] - span[0])) / 2, y), label, font=font, fill=fill)
    for position in range(first, last):
        if fragment_ends and spacer_start <= position < spacer_end:
            continue
        if clean and gap_start <= position < gap_start + gap_length + (4 if joined else 0):
            continue
        base = sequence[position % length]
        complement = base.translate(str.maketrans("ACGT", "TGCA"))
        if (joined or fragment_ends) and projected(position) not in coding_columns and position not in recognition_columns:
            base, complement = base.lower(), complement.lower()
        for strand, text, y in (("TopStrand", base, 15), ("BottomStrand", complement, 52)):
            value = colors[strand][position % length]
            if joined and gap_start - 4 <= position < gap_start:
                value = colors["BottomStrand"][position % length]
            # Center the glyph in its base column so cut lines occupy the
            # spaces between bases, never the letters themselves.
            if joined or fragment_ends:
                centered(text, position, position + 1, y, font=dna_font, fill=value)
            else:
                draw.text((x(position), y), text, font=dna_font, fill=value)
        if not joined or position < gap_start - 4 or position >= gap_start:
            draw.line((x(position) + cell / 2, 39, x(position) + cell / 2, 47), fill="#666666", width=1)
    axis_windows = [(first, gap_start - 4), (gap_start + gap_length + 4, last)] if joined else windows
    for start, end in axis_windows:
        draw.line((x(start), 43, x(end), 43), fill="#72777d" if clean else "#666666", width=1 if clean else 2)
    if joined:
        draw.line([(x(gap_start - 4), 7), (x(gap_start - 4), 43), (x(gap_start), 43), (x(gap_start), 76)], fill="#444444", width=3)
    if fragment_ends:
        # Both fragments carry a 5′ four-base overhang. The gray top copy on
        # the left and gray bottom copy on the right are removed by digestion.
        for cut in (gap_start - 4, gap_start + gap_length):
            draw.line([(x(cut), 7), (x(cut), 43), (x(cut + 4), 43), (x(cut + 4), 76)], fill="#111111", width=3)
    if joined:
        # These bars describe physical additions, not an extension of the
        # intrinsic feature itself. Draw first so native core bars remain solid
        # wherever a copied fusion overlaps a neighbouring core annotation.
        for feature in features:
            color = _display_color(_qualifier(feature, "fragment_color"), "#777777")
            for a, b in json.loads(_qualifier(feature, "assigned_extensions", "[]")):
                for shift in (-length, 0, length):
                    for begin, finish in windows:
                        left, right = max(a + shift, begin), min(b + shift, finish)
                        if left < right:
                            _draw_hatched_region(draw, (x(left), 80, x(right), 99), color)
    for feature, segments in visible:
        name = feature.get("name", "")
        if name == "└┐":
            if clean:
                continue
            for segment, a, b, _ in segments:
                a, b = max(first, a), min(last, b)
                label = segment.get("name", "")
                if fragment_ends:
                    if label == "//":
                        continue
                    enzyme_name = label.strip("←→")
                    recognition = ENZYMES.get(enzyme_name, {}).get("reverse_site" if label.startswith("←") else "forward_site", "")
                    shown = "".join(sequence[p % length] for p in range(a, b))
                    offset = shown.find(recognition) if recognition else -1
                    if offset >= 0:
                        a, b = a + offset, a + offset + len(recognition)
                    xa, xb = x(a), x(b)
                    points = ([(xb, 87), (xa + 11, 87), (xa, 96), (xa + 11, 105), (xb, 105)] if label.startswith("←") else
                              [(xa, 87), (xb - 11, 87), (xb, 96), (xb - 11, 105), (xa, 105)])
                    draw.polygon(points, fill="#303030", outline="#111111")
                    centered(enzyme_name, a, b, 88, fill="#ffffff")
                    continue
                if label == "//":
                    centered("//", a, b, 78)
                else:
                    draw.line((x(a), 88, x(b), 88), fill="#111111", width=3)
                    # Leave a gray label backing so the line does not cross text.
                    text_width = draw.textbbox((0, 0), label, font=small)[2]
                    mid = (x(a) + x(b)) / 2
                    draw.rectangle((mid - text_width / 2 - 3, 79, mid + text_width / 2 + 3, 96), fill="#f2f2f2")
                    centered(label, a, b, 79)
            continue
        is_core = bool(_FIXED_FRAGMENT_LABEL.fullmatch(name) or _VARIABLE_FRAGMENT_LABEL.fullmatch(name))
        is_fusion = bool(re.fullmatch("[ACGT]{4}", name)) and sum(b - a for _, a, b, _ in segments) == 4
        if is_core or is_fusion:
            if clean and is_fusion:
                continue
            for segment, a, b, _ in segments:
                a, b = max(first, a), min(last, b)
                if a >= b:
                    continue
                if fragment_ends and is_fusion and a < gap_start:
                    continue  # Label the shared overhang once, on the right.
                draw.rectangle((x(a), 87 if fragment_ends else 80, x(b), 106 if fragment_ends else 99 if clean else 97), fill="#ffffff" if fragment_ends and is_fusion else segment.get("color", "#cccccc"), outline="#a8adb3" if clean else "#777777")
                centered(name, a, b, (88 if is_core else 110) if fragment_ends else 81 if is_core else 101)
    for feature, segments, lane in protein:
        y = protein_y + lane * 72
        codons = _visible_codons(feature, sequence, first, last)
        for segment, a, b, aa_offset in segments:
            left, right = max(first, a), min(last, b)
            fill = segment.get("color", "#cccccc")
            if fill == "noColor":
                fill = "#cccccc"
            reverse_arrow = feature.get("directionality") == "2" and a >= first
            arrow = feature.get("directionality") == "1" and b <= last
            tip = min(10, (x(right) - x(left)) / 2)
            points = ([(x(right), y), (x(left) + tip, y), (x(left), y + 9), (x(left) + tip, y + 18), (x(right), y + 18)] if reverse_arrow else
                      [(x(left), y), (x(right) - (tip if arrow else 0), y), (x(right), y + 9), (x(right) - (tip if arrow else 0), y + 18), (x(left), y + 18)])
            draw.polygon(points, fill=fill, outline="#777777")
            label = segment.get("name") or feature.get("name", "")
            label_width = draw.textbbox((0, 0), label, font=small)[2]
            inside = (clean or fragment_ends) and label_width + 16 < x(right) - x(left)
            # Dark feature colors need white type for legibility.
            label_color = "#ffffff" if inside and sum(int(fill[i:i + 2], 16) * weight for i, weight in ((1, .299), (3, .587), (5, .114))) < 150 else "#111111"
            centered(label, left, right, y + 1 if inside else y + 23, fill=label_color)
            if segment.get("name") and segment.get("name") != feature.get("name"):
                centered(feature.get("name", ""), left, right, y + 38)
            if _coding_segment(feature, segment):
                draw.line((x(left), y - 31, x(right), y - 31), fill="#888888")
                for codon in codons:
                    if left <= codon["start"] + 1 < right:
                        centered(codon["aa"], codon["start"], codon["end"], y - 23)
                        if codon["number"] == 1 or codon["number"] % 5 == 0:
                            centered(str(codon["number"]), codon["start"], codon["end"], y - 46)
    if not clean and not fragment_ends:
        draw.text((margin, height - 20), "Cloning schematic · gray adapters / N spacer removed during assembly · see *_assembled.png for final junction", font=_font(10), fill="#777777")
    image.save(path, "PNG")


def _write_junction_screencap(
    path: Path, junction_id: str, fusion: str, sequence: str, fusion_start: int,
    enzyme: str, left_name: str, right_name: str, left_color: str, right_color: str,
    features: list[Element] | None = None,
) -> None:
    """Render a sequence-detail PNG of the actual selected circular junction."""
    flank = 15
    context = "".join(sequence[(fusion_start + offset) % len(sequence)] for offset in range(-flank, 4 + flank))
    fusion_at = flank
    # Aligned lower strand is a complement (3′→5′), not a reverse complement.
    reverse = context.translate(str.maketrans("ACGT", "TGCA"))
    reverse_fusion_at = fusion_at
    first, last = fusion_start - flank, fusion_start + 4 + flank
    nearby = []
    for feature in features or []:
        name = feature.get("name", "")
        if _FIXED_FRAGMENT_LABEL.fullmatch(name) or _VARIABLE_FRAGMENT_LABEL.fullmatch(name) or name.startswith("junction_"):
            continue
        spans = []
        for segment in feature.findall("Segment"):
            a, b = map(int, segment.get("range").split("-"))
            for shift in (-len(sequence), 0, len(sequence)):
                left, right = max(first, a - 1 + shift), min(last, b + shift)
                if left < right:
                    spans.append((segment, left, right))
        if spans:
            nearby.append((feature, spans))
    font, small, title_font = _font(22), _font(16), _font(25)
    char_width, x0, width, height = 25, 52, max(1420, (len(context) + 4) * 25 + 100), 290
    height += 106 * len(nearby)
    image = Image.new("RGB", (width, height), "#f8faf8")
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle((18, 18, width - 18, height - 18), radius=14, fill="#ffffff", outline="#d9e1db", width=2)
    draw.text((x0, 38), f"{junction_id}  ·  selected fusion {fusion}", fill="#18382c", font=title_font)
    draw.text((x0, 73), f"{enzyme} scarless junction  ·  {left_name} → {right_name}", fill="#65766e", font=small)

    def draw_strand(text: str, y: int, start: int, reverse_row: bool = False) -> None:
        draw.text((x0 - 29, y + 3), "3′" if reverse_row else "5′", fill="#65766e", font=small)
        for index, base in enumerate(text):
            x = x0 + index * char_width
            if start <= index < start + 4:
                draw.rounded_rectangle((x - 2, y - 3, x + char_width - 3, y + 28), radius=4, fill="#fde6bc")
                color = "#a45e00"
            elif index < start:
                color = left_color
            else:
                color = right_color
            draw.text((x, y), base, fill=color, font=font)
        draw.line((x0, y + 35, x0 + len(text) * char_width - 5, y + 35), fill="#76847d", width=2)

    draw_strand(context, 116, fusion_at)
    draw_strand(reverse, 181, reverse_fusion_at, True)
    draw.text((x0 + fusion_at * char_width, 242), "selected 4 bp fusion", fill="#a45e00", font=small)
    def x(position):
        return x0 + (position - first) * char_width
    for index, (feature, spans) in enumerate(nearby):
        y = 282 + index * 106
        coding = any(_coding_segment(feature, s) for s, _, _ in spans)
        label = feature.get('name', '') + (' · amino acids' if coding else '')
        draw.text((x0, y), label, fill="#334b40", font=small)
        for segment, a, b in spans:
            color = segment.get("color", "#cccccc")
            if color == "noColor":
                color = "#eeeeee"
            tip = min(10, (x(b) - x(a)) / 2) if feature.get("directionality") in {"1", "2"} else 0
            reverse = feature.get("directionality") == "2"
            points = ([(x(b), y + 27), (x(a) + tip, y + 27), (x(a), y + 34), (x(a) + tip, y + 41), (x(b), y + 41)] if reverse else
                      [(x(a), y + 27), (x(b) - tip, y + 27), (x(b), y + 34), (x(b) - tip, y + 41), (x(a), y + 41)])
            draw.polygon(points, fill=color, outline="#829389")
        for codon in _visible_codons(feature, sequence, first, last):
            mid = (x(codon["start"]) + x(codon["end"])) / 2
            box = draw.textbbox((0, 0), codon["aa"], font=font)
            draw.text((mid - (box[2] - box[0]) / 2, y + 46), codon["aa"], fill="#243e32", font=font)
            if codon["number"] == 1 or codon["number"] % 5 == 0:
                number = str(codon["number"])
                box = draw.textbbox((0, 0), number, font=small)
                draw.text((mid - (box[2] - box[0]) / 2, y + 74), number, fill="#65766e", font=small)
    image.save(path, "PNG")


def _write_plasmid_junction_figure(path, rows, starts, retained, sequence, junction_rows, junction_locations, joined_views):
    """To-scale linearized plasmid with alternating sequence-detail callouts."""
    length = len(sequence)
    backbone = next((i for i, row in enumerate(rows) if row.get("is_backbone")), 0)
    rotation = (starts[backbone] + (len(retained[backbone]) - 4) // 2) % length
    items = sorted(zip(junction_rows, joined_views), key=lambda pair: (junction_locations[str(pair[0]["junction_id"])][0] - rotation) % length)
    images = [Image.open(filename).convert("RGB") for _, filename in items]
    panel_width = max(image.width for image in images)
    panel_height = max(image.height for image in images)
    slots = (len(items) + 1) // 2
    margin, gutter = 80, 90
    width = margin * 2 + slots * panel_width + (slots - 1) * gutter
    map_y = 120 + panel_height + 210
    bottom_y = map_y + 210
    height = bottom_y + panel_height + 140
    figure = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(figure)
    colors = _fragment_palette(rows)
    font, small = _font(27), _font(21)
    map_left, map_right = margin + 50, width - margin - 50
    def x(base):
        return map_left + base / length * (map_right - map_left)
    draw.text((margin, 24), f"Assembled plasmid · {length:,} bp", font=_font(35), fill="#18382c")
    draw.text((margin, 69), "Linearized within the backbone · junction details alternate above and below", font=small, fill="#647168")
    # Draw expansion wedges behind the map and the white detail panels.
    for rank, ((junction, filename), image) in enumerate(zip(items, images)):
        column = rank // 2
        panel_x = margin + column * (panel_width + gutter)
        top = rank % 2 == 0
        panel_y = 125 if top else bottom_y
        position, left_index, right_index = junction_locations[str(junction["junction_id"])]
        position = (position - rotation) % length
        a, b = x(max(0, position - 22)), x(min(length, position + 26))
        edge = panel_y + panel_height + 36 if top else panel_y - 16
        near = map_y - 27 if top else map_y + 27
        points = [(panel_x + 20, edge), (panel_x + panel_width - 20, edge), (b, near), (a, near)]
        draw.polygon(points, fill="#f4f6f8")
        draw.line((panel_x + 20, edge, a, near), fill="#cad2d8", width=2)
        draw.line((panel_x + panel_width - 20, edge, b, near), fill="#cad2d8", width=2)
        title = f"{rank + 1}  {junction['selected_fusion']} · {rows[left_index]['element_name']} → {rows[right_index]['element_name']}"
        draw.text((panel_x + 18, panel_y), title, font=font, fill="#253c30")
        figure.paste(image, (panel_x, panel_y + 37))
    for index, row in enumerate(rows):
        # Display each shared fusion with its downstream physical fragment;
        # the underlying assembled sequence and native annotations are intact.
        start = (starts[index] - 4 - rotation) % length
        end = start + len(retained[index]) - 4
        parts = [(start, min(length, end))] + ([(0, end - length)] if end > length else [])
        color = colors[index]
        for a, b in parts:
            draw.rectangle((x(a), map_y - 24, x(b), map_y + 24), fill=color)
        for a, b in _fragment_extensions(row, len(retained[index])):
            extension_start = (start + a) % length
            extension_end = extension_start + b - a
            spans = [(extension_start, min(length, extension_end))] + ([(0, extension_end - length)] if extension_end > length else [])
            for a, b in spans:
                _draw_hatched_region(draw, (x(a), map_y - 24, x(b), map_y + 24), "white", background=color, spacing=5)
        core_start, core_end = _fragment_core_bounds(row, len(retained[index]))
        label_start = (start + core_start) % length
        label_end = label_start + max(0, min(core_end, len(retained[index]) - 4) - core_start)
        label_parts = [(label_start, min(length, label_end))] + ([(0, label_end - length)] if label_end > length else [])
        for a, b in label_parts:
            label = row["element_name"]
            label_width = draw.textbbox((0, 0), label, font=font)[2]
            if x(b) - x(a) > label_width + 12:
                draw.text(((x(a) + x(b) - label_width) / 2, map_y - 17), label, font=font, fill="#18382c")
    for tick in range(6):
        base = round(length * tick / 5)
        draw.line((x(base), map_y + 28, x(base), map_y + 39), fill="#778079", width=2)
        label = f"{base:,}"
        span = draw.textbbox((0, 0), label, font=small)[2]
        draw.text((x(base) - span / 2, map_y + 43), label, font=small, fill="#647168")
    for rank, (junction, _) in enumerate(items):
        position = (junction_locations[str(junction["junction_id"])][0] + 2 - rotation) % length
        anchor = x(position)
        draw.line((anchor, map_y - 31, anchor, map_y + 31), fill="#283d32", width=3)
        badge_y = map_y - 66 if rank % 2 == 0 else map_y + 106
        draw.ellipse((anchor - 19, badge_y - 19, anchor + 19, badge_y + 19), fill="#253c30")
        label = str(rank + 1)
        span = draw.textbbox((0, 0), label, font=small)[2]
        draw.text((anchor - span / 2, badge_y - 14), label, font=small, fill="white")
    _draw_hatched_region(draw, (margin, height - 87, margin + 40, height - 62), _fragment_color(0))
    draw.text((margin + 53, height - 90), "Hatched: assigned DNA outside the annotated core. Solid: annotated core.", font=small, fill="#647168")
    draw.text((margin, height - 46), "Map is to scale; short extensions are clearer in the zoom-ins. DNA colors retain native core ownership.", font=small, fill="#647168")
    figure.save(path, "PNG")


def _assembled_circular_map(fragment_rows: list[dict[str, Any]]) -> tuple[str, list[int], list[str]]:
    """Return the scarless circular DNA, core starts, and retained fragment sequences."""
    retained = []
    for row in fragment_rows:
        retained.append(
            (row["left_fusion"] if row["left_fusion_owner"] == "extension" else "") + row["core_sequence"] +
            (row["right_fusion"] if row["right_fusion_owner"] == "extension" else "")
        )
    if not retained or not retained[0].startswith(fragment_rows[0]["left_fusion"]):
        raise ValueError("The first fragment must retain its left fusion for circular assembly")
    sequence, starts = retained[0][4:], [0]
    for index in range(1, len(retained)):
        if retained[index - 1][-4:] != retained[index][:4]:
            raise ValueError(f"Assembly order is incompatible at {fragment_rows[index - 1]['fragment_id']} → {fragment_rows[index]['fragment_id']}")
        starts.append(len(sequence))
        sequence += retained[index][4:]
    if retained[-1][-4:] != retained[0][:4]:
        raise ValueError("Final Golden Gate fusion does not close the circular assembly")
    return sequence, starts, retained


def _span_covers_boundary(start: int, end: int, boundary: int, length: int) -> bool:
    """Whether a non-wrapping feature span covers/touches a circular boundary."""
    return start <= boundary <= end or (end == length and boundary == 0)


def _is_puc_origin_feature(name: str) -> bool:
    return bool(_PUC_ORIGIN_LABEL.fullmatch(name.strip()))


def _spans_overlap(first: list[tuple[int, int]], second: list[tuple[int, int]]) -> bool:
    return any(start_a < end_b and start_b < end_a for start_a, end_a in first for start_b, end_b in second)


def _merge_adjacent_spans(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Normalize feature sublocations without changing circular feature order."""
    merged: list[tuple[int, int]] = []
    for start, end in spans:
        if merged and merged[-1][1] == start:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def plan_from_annotated_snapgene(path: str | Path, enzyme: str = "BsmBI") -> dict[str, Any]:
    """Build a Golden Gate plan from semantic SnapGene feature annotations.

    A feature whose entire label is ``[name]`` is a fixed fragment core, and
    one labelled ``{name}`` is a variable fragment core. The core overlapping
    a ``pUC ori`` (or ``ori pUC``) feature is the backbone. A feature whose
    entire label is ``(name)`` or exactly ``-`` is an allowed junction-placement region. Each
    region is paired to the core boundary it touches on the circular map.
    """
    path = Path(path).expanduser().resolve()
    if path.suffix.lower() != ".dna":
        raise ValueError("Annotated Golden Gate design requires a SnapGene .dna file")
    scaffold = Scaffold.from_path(path)
    components: list[dict[str, Any]] = []
    junction_regions: list[dict[str, Any]] = []
    origin_spans: list[tuple[int, int]] = []
    for name, spans in scaffold.features.items():
        if _is_puc_origin_feature(name):
            origin_spans.extend(spans)
    for name, spans in scaffold.features.items():
        label = name.strip()
        fixed = _FIXED_FRAGMENT_LABEL.fullmatch(label)
        variable = _VARIABLE_FRAGMENT_LABEL.fullmatch(label)
        junction = _JUNCTION_PLACEMENT_LABEL.fullmatch(label)
        if fixed or variable:
            if not spans:
                continue
            normalized_spans = _merge_adjacent_spans(spans)
            components.append({
                "feature": name,
                "name": (fixed or variable).group(1).strip(),
                "kind": "fixed" if fixed else "variable",
                "spans": normalized_spans,
                "start": normalized_spans[0][0],
                "end": normalized_spans[-1][1],
            })
        elif _is_puc_origin_feature(name):
            continue
        elif junction:
            for start, end in spans:
                junction_regions.append({"label": (junction.group(1) or "-").strip(), "start": start, "end": end})
    for index, component in enumerate(components):
        for other in components[index + 1:]:
            if _spans_overlap(component["spans"], other["spans"]):
                raise ValueError(
                    f"Fragment core annotations overlap: {component['feature']!r} and {other['feature']!r}. "
                    "Each DNA base may belong to only one bracketed fragment core."
                )
    if len(components) < 2:
        raise ValueError("Add at least two fragment cores labelled [fixed name] or {variable name}")
    backbones = [component for component in components if _spans_overlap(component["spans"], origin_spans)]
    if len(backbones) != 1:
        if not backbones:
            raise ValueError("Add a pUC ori feature overlapping exactly one fragment core to identify the backbone")
        names = ", ".join(component["feature"] for component in backbones)
        raise ValueError(f"The pUC ori feature overlaps multiple fragment cores ({names}); it must identify one backbone")
    backbone = backbones[0]
    for component in components:
        component["is_backbone"] = component is backbone
        component["source"] = "PCR" if component["is_backbone"] else "SYN"
    components.sort(key=lambda component: component["start"])
    if not junction_regions:
        raise ValueError("Add highlighted junction-placement features labelled (name) or -")
    assigned_junctions: list[dict[str, Any]] = []
    for index, component in enumerate(components):
        next_component = components[(index + 1) % len(components)]
        matches = [
            region for region in junction_regions
            if _span_covers_boundary(region["start"], region["end"], next_component["start"], len(scaffold.sequence))
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Between fragment cores {component['feature']!r} and {next_component['feature']!r}, "
                f"add exactly one junction-placement feature labelled (name) or - that touches their boundary (found {len(matches)})"
            )
        assigned_junctions.append(matches[0])
    if len({(region["start"], region["end"]) for region in assigned_junctions}) != len(junction_regions):
        raise ValueError("Every (name) or - junction-placement feature must touch one fragment-core boundary")
    plan = {
        "project": path.stem,
        "scaffold": str(path),
        "enzyme": enzyme,
        "coordinate_assembly": True,
        "junction_candidates": [
            {"id": f"junction_{index + 1}", "scaffold_range": [region["start"], region["end"]], "label": region["label"]}
            for index, region in enumerate(assigned_junctions)
        ],
        "fragments": [],
    }
    junction_ids = [item["id"] for item in plan["junction_candidates"]]
    for index, component in enumerate(components):
        plan["fragments"].append(
            {
                "id": f"fragment_{index + 1}",
                "name": component["name"],
                "annotation_kind": component["kind"],
                "is_backbone": component["is_backbone"],
                "source": component["source"],
                "core": {"scaffold_feature": component["feature"]},
                "left_junction": junction_ids[index - 1],
                "right_junction": junction_ids[index],
            }
        )
    return plan


def _coordinate_fragment_plan(plan: dict[str, Any], scaffold: Scaffold, chosen: dict[str, Candidate]) -> tuple[dict[str, Any], str | None]:
    """Cut native DNA at selected fusion coordinates; never add scar bases.

    Each retained fragment starts at its left fusion and ends after its right
    fusion. Neighboring fragments overlap by exactly four original bases.
    Variable entries replace their labelled core, retaining native flanks and
    validating any fusion that lies inside the supplied variant sequence.
    """
    if not plan.get("coordinate_assembly"):
        return plan, None
    resolved = copy.deepcopy(plan)
    length = len(scaffold.sequence)
    junction_specs = {str(item["id"]): item for item in plan["junction_candidates"]}
    cuts = {}
    for junction_id, candidate in chosen.items():
        start, _ = junction_specs[junction_id]["scaffold_range"]
        cuts[junction_id] = (int(start) + candidate.offset) % length
    if len(set(cuts.values())) != len(cuts):
        raise ValueError("Selected junctions must occupy distinct scaffold positions")
    expected_parts = []
    source_cores = []
    for fragment in plan["fragments"]:
        original, spans = scaffold.feature(fragment["core"]["scaffold_feature"])
        source_cores.append((original, spans[0][0], spans[-1][1] % length))
    for index, fragment in enumerate(resolved["fragments"]):
        original, core_start, core_end = source_cores[index]
        fixed_prefix, fixed_suffix = "", ""
        if fragment.get("annotation_kind") == "variable":
            fixed_positions = {i for i in range(len(original))
                               if any(spec.get("label") == "-" and int(spec["scaffold_range"][0]) <= (core_start + i) % length < int(spec["scaffold_range"][1])
                                      for spec in junction_specs.values())}
            editable = [i for i in range(len(original)) if i not in fixed_positions]
            if not editable or any(i in fixed_positions for i in range(editable[0], editable[-1] + 1)):
                raise ValueError(f"{fragment['name']}: '-' regions must leave one contiguous variable core")
            fixed_prefix, fixed_suffix = original[:editable[0]], original[editable[-1] + 1:]
        start = cuts[str(fragment["left_junction"])]
        right_start = cuts[str(fragment["right_junction"])]
        window_length = (right_start - start) % length + 4
        if window_length > length:
            raise ValueError(f"{fragment['name']}: overlapping or reversed junction windows")
        end = (right_start + 4) % length
        native_window = "".join(scaffold.sequence[(start + i) % length] for i in range(window_length))
        relative_start = (core_start - start) % length
        if relative_start >= window_length:
            relative_start -= length
        relative_end = relative_start + len(original)
        overlap_start, overlap_end = max(0, relative_start), min(window_length, relative_end)
        if overlap_end <= overlap_start:
            raise ValueError(f"{fragment['name']}: selected junction window does not cover its core")
        fragment["variant_window"] = {
            "fixed_core_prefix": fixed_prefix, "fixed_core_suffix": fixed_suffix,
            "prefix": native_window[:overlap_start], "suffix": native_window[overlap_end:],
            "trim_left": max(0, -relative_start), "trim_right": max(0, relative_end - window_length),
            "discarded_prefix": original[:max(0, -relative_start)],
            "discarded_suffix": original[len(original) - max(0, relative_end - window_length):],
        }
        fragment["core"] = {"scaffold_range": [start, end]}
        fragment["native_window"] = [start, end]
        fragment["native_core_feature"] = plan["fragments"][index]["core"]["scaffold_feature"]
        fragment["left_fusion_owner"] = fragment["right_fusion_owner"] = "core"
        # Independent target reconstruction uses full cores + intervening native
        # DNA, not the overlap-trimming assembly routine under test.
        variants = fragment.get("variants")
        expected_parts.append(fixed_prefix + _dna(str(variants[0]["sequence"])) + fixed_suffix if variants else original)
        next_start = source_cores[(index + 1) % len(source_cores)][1]
        expected_parts.append("" if core_end == next_start else _range(scaffold.sequence, core_end, next_start))
    return resolved, "".join(expected_parts)


def _run_plan_data(
    plan: dict[str, Any], scaffold_path: Path, output_root: str | Path,
    design_name: str, plan_reference: str,
) -> dict[str, Any]:
    """Run already-resolved plan data and write only below output_root."""
    scaffold_path = scaffold_path.expanduser().resolve()
    scaffold = Scaffold.from_path(scaffold_path)
    enzyme_input = plan.get("enzyme", "BsmBI")
    enzyme = {"name": enzyme_input, **ENZYMES[enzyme_input]} if isinstance(enzyme_input, str) and enzyme_input in ENZYMES else enzyme_input
    if not isinstance(enzyme, dict) or "forward_site" not in enzyme or "reverse_site" not in enzyme:
        raise ValueError("Use BsmBI/BsaI or provide an enzyme mapping with forward_site and reverse_site")
    forward_site, reverse_site = _dna(str(enzyme["forward_site"])), _dna(str(enzyme["reverse_site"]))
    defaults = plan.get("defaults", {})
    padding, spacer = _dna(str(defaults.get("padding", "GGCTAC"))), _dna(str(defaults.get("spacer", "A")))
    anneal_bp = int(defaults.get("pcr_anneal_bp", 22))
    chosen, candidate_rows = _select_junctions(plan, scaffold)
    plan, expected_target = _coordinate_fragment_plan(plan, scaffold, chosen)
    fragment_rows: list[dict[str, Any]] = []; primer_rows: list[dict[str, Any]] = []; fasta_records: list[str] = []; warnings: list[str] = []
    for fragment in _expanded_fragment_specs(plan):
        fragment_id, source = str(fragment["id"]), str(fragment.get("source", "SYN")).upper()
        if source not in {"SYN", "PCR"}:
            raise ValueError(f"{fragment_id}: source must be SYN or PCR")
        left, right = chosen[str(fragment["left_junction"])].fusion, chosen[str(fragment["right_junction"])].fusion
        core = _core(fragment["core"], scaffold, scaffold_path.parent)
        left_owner, right_owner = str(fragment.get("left_fusion_owner", "extension")).lower(), str(fragment.get("right_fusion_owner", "extension")).lower()
        if left_owner not in {"core", "extension"} or right_owner not in {"core", "extension"}:
            raise ValueError(f"{fragment_id}: fusion owner must be core or extension")
        if left_owner == "core" and not core.startswith(left):
            raise ValueError(f"{fragment_id}: core must start with its owned left fusion {left}")
        if right_owner == "core" and not core.endswith(right):
            raise ValueError(f"{fragment_id}: core must end with its owned right fusion {right}")
        if not fragment.get("allow_repeated_fusion", False):
            for side, owner, already_present in (("left", left_owner, core.startswith(left)), ("right", right_owner, core.endswith(right))):
                if owner == "extension" and already_present:
                    raise ValueError(f"{fragment_id}: {side} fusion already occurs at the core end; adding it again would duplicate 4 bp. Declare {side}_fusion_owner: core, or explicitly allow_repeated_fusion for an intentional repeat")
        five = padding + forward_site + spacer + (left if left_owner == "extension" else "")
        three = (right if right_owner == "extension" else "") + spacer + reverse_site + padding
        full = five + core + three
        internal = sorted({site for site in (forward_site, reverse_site) if site in core})
        if internal:
            warnings.append(f"{fragment_id}: internal {enzyme.get('name', 'Type IIS')} site(s): {', '.join(internal)}")
        row = {"fragment_id": fragment_id, "name": fragment.get("name", fragment_id), "source": source,
               "element_id": fragment["element_id"], "element_name": fragment["element_name"],
               "assembly_position": fragment["assembly_position"], "option_index": fragment["option_index"],
               "variant_name": fragment["variant_name"], "annotation_kind": fragment.get("annotation_kind", ""),
               "core_source": fragment["core"], "native_core_feature": fragment.get("native_core_feature"), "native_window": fragment.get("native_window"), "variant_window": fragment.get("variant_window"), "is_backbone": bool(fragment.get("is_backbone")),
               "left_junction": fragment["left_junction"], "left_fusion": left, "left_fusion_owner": left_owner,
               "right_junction": fragment["right_junction"], "right_fusion": right, "right_fusion_owner": right_owner,
               "core_bp": len(core), "five_prime_extension": five, "core_sequence": core, "three_prime_extension": three,
               "final_fragment_bp": len(full), "final_fragment_sequence": full, "internal_enzyme_sites": ";".join(internal)}
        fragment_rows.append(row)
        if source == "SYN":
            fasta_records.append(to_fasta(fragment_id, full))
        else:
            bases = int(fragment.get("pcr_anneal_bp", anneal_bp))
            if len(core) < bases:
                raise ValueError(f"{fragment_id}: core is shorter than its PCR annealing window")
            primer_rows.append({"fragment_id": fragment_id, "name": row["name"], "forward_primer_5_to_3": five + core[:bases],
                                "reverse_primer_5_to_3": reverse_complement(three) + reverse_complement(core[-bases:]),
                                "annealing_bp": bases, "expected_amplicon_bp": len(full)})
    # Validate reconstruction before creating ANY order or map outputs.
    representative_rows = [row for row in fragment_rows if row["option_index"] == 1]
    map_variants = [{key: row[key] for key in ("element_name", "variant_name", "fragment_id")} for row in representative_rows]
    assembled_sequence, component_starts, retained = _assembled_circular_map(representative_rows)
    if expected_target is not None and (len(assembled_sequence) != len(expected_target) or assembled_sequence not in expected_target + expected_target):
        raise ValueError("Assembly does not reproduce the intended core replacements and native junction DNA; no valid product can be exported")
    output_root = Path(output_root).expanduser().resolve()
    output_dir = output_root / f"{design_name}_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    junction_rows = [{"junction_id": spec["id"], "selected_fusion": chosen[str(spec["id"])].fusion,
                      "offset_in_allowed_region": chosen[str(spec["id"])].offset, "candidate_region_source": chosen[str(spec["id"])].source,
                      "score": f"{chosen[str(spec['id'])].score:.2f}"} for spec in plan.get("junction_candidates", [])]
    automatic_junctions = []
    for spec in plan.get("junction_candidates", []):
        if spec.get("label") != "-" or "scaffold_range" not in spec:
            continue
        start, end = spec["scaffold_range"]
        region = _range(scaffold.sequence, start, end)
        candidate = chosen[str(spec["id"])]
        left_fragment = next(row for row in representative_rows if row["right_junction"] == spec["id"])
        right_fragment = next(row for row in representative_rows if row["left_junction"] == spec["id"])
        automatic_junctions.append({"junction_id": spec["id"], "scaffold_region": region,
                                   "left_fragment": left_fragment["element_name"], "left_end": region[:candidate.offset + 4],
                                   "right_fragment": right_fragment["element_name"], "right_end": region[candidate.offset:],
                                   "shared_fusion": candidate.fusion})
    _tsv(output_dir / "automatic_junction_bases.tsv", automatic_junctions,
         ["junction_id", "scaffold_region", "left_fragment", "left_end", "right_fragment", "right_end", "shared_fusion"])
    _tsv(output_dir / "junctions.tsv", junction_rows); _tsv(output_dir / "candidate_options.tsv", candidate_rows); _tsv(output_dir / "fragments.tsv", fragment_rows)
    recipe_keys = ["assembly_position", "element_id", "element_name", "option_index", "variant_name", "fragment_id", "name", "source", "left_fusion", "right_fusion", "final_fragment_bp"]
    _tsv(output_dir / "assembly_recipe.tsv", [{key: row[key] for key in recipe_keys} for row in fragment_rows])
    order_keys = ["element_id", "element_name", "variant_name", "fragment_id", "name", "left_fusion", "right_fusion", "final_fragment_bp", "final_fragment_sequence"]
    _tsv(output_dir / "synthesis_order.tsv", [{key: row[key] for key in order_keys} for row in fragment_rows if row["source"] == "SYN"], order_keys)
    _tsv(output_dir / "pcr_primers.tsv", primer_rows, ["fragment_id", "name", "forward_primer_5_to_3", "reverse_primer_5_to_3", "annealing_bp", "expected_amplicon_bp"])
    (output_dir / "synthesis_order.fasta").write_text("".join(fasta_records), encoding="utf-8")
    schematic_path = output_dir / "assembly_schematic.dna"
    schematic_junctions = _write_assembly_schematic(schematic_path, representative_rows, fragment_rows, scaffold, enzyme, padding, spacer)
    map_features: list[dict[str, Any]] = []
    component_colors = _fragment_palette(representative_rows)
    for index, row in enumerate(representative_rows):
        color = component_colors[index]
        kind = row.get("annotation_kind")
        display_name = f"[{row['name']}]" if kind == "fixed" else f"{{{row['name']}}}" if kind == "variable" else row["name"]
        map_features.append({
            "name": display_name, "start": component_starts[index], "end": component_starts[index] + len(retained[index]) - 4,
            "color": color, "note": f"{row['source']} fragment; selected fusions {row['left_fusion']} → {row['right_fusion']}",
        })
    junction_locations: dict[str, tuple[int, int, int]] = {}
    for index, row in enumerate(representative_rows):
        position = len(assembled_sequence) - 4 if index == 0 else component_starts[index] - 4
        junction_locations[str(row["left_junction"])] = (position, (index - 1) % len(representative_rows), index)
        map_features.append({"name": f"{row['left_junction']} · {row['left_fusion']}", "start": position, "end": position + 4,
                             "color": _JUNCTION_COLOR, "note": "Selected Golden Gate fusion overhang"})
    map_path = output_dir / "assembled_product_map.dna"
    native_mapping = {}
    for index, row in enumerate(representative_rows):
        left_extension = 4 if row["left_fusion_owner"] == "extension" else 0
        native_mapping.update({p: target % len(assembled_sequence) for p, target in
                               _native_row_mapping(row, scaffold, component_starts[index] - 4 + left_extension).items()})
    map_features.extend(_project_annotations(_snapgene_features_xml(scaffold.path), scaffold.sequence, native_mapping, require_complete=True))
    _write_snapgene_map(map_path, assembled_sequence, map_features, f"{plan.get('project', design_name)} assembled Golden Gate product")
    # Read-back check catches malformed binary output before exposing it to the user.
    map_readback = SeqIO.read(map_path, "snapgene")
    if str(map_readback.seq) != assembled_sequence:
        raise ValueError("Generated SnapGene map did not round-trip correctly")
    screencap_dir = output_dir / "junction_screencaps"
    screencap_dir.mkdir(exist_ok=True)
    screencaps = []
    clean_views = []
    joined_views = []
    for row in junction_rows:
        junction_id = str(row["junction_id"])
        position, left_index, right_index = junction_locations[junction_id]
        filename = re.sub(r"[^A-Za-z0-9_.-]+", "_", junction_id).strip("._") or "junction"
        image_path = screencap_dir / f"{filename}.png"
        _write_junction_screencap(
            screencap_dir / f"{filename}_assembled.png", junction_id, str(row["selected_fusion"]), assembled_sequence, position,
            str(enzyme.get("name", "Type IIS")), str(representative_rows[left_index]["name"]), str(representative_rows[right_index]["name"]),
            component_colors[left_index], component_colors[right_index],
            features=_snapgene_features_xml(map_path),
        )
        _write_schematic_junction_image(image_path, schematic_path, schematic_junctions[junction_id], len(spacer) * 2 + len(padding) * 2 + len(enzyme["forward_site"]) + len(enzyme["reverse_site"]) + 6, fragment_ends=True)
        clean_path = screencap_dir / f"{filename}_clean.png"
        _write_schematic_junction_image(clean_path, schematic_path, schematic_junctions[junction_id], len(spacer) * 2 + len(padding) * 2 + len(enzyme["forward_site"]) + len(enzyme["reverse_site"]) + 6, clean=True)
        clean_views.append(str(clean_path))
        joined_path = screencap_dir / f"{filename}_joined.png"
        _write_schematic_junction_image(joined_path, schematic_path, schematic_junctions[junction_id], len(spacer) * 2 + len(padding) * 2 + len(enzyme["forward_site"]) + len(enzyme["reverse_site"]) + 6, joined=True)
        joined_views.append(str(joined_path))
        screencaps.append(str(image_path))
    # Overview begins with the backbone-to-insert junction, following the
    # circular assembly order. Individual images stay available for slides.
    backbone_index = next((i for i, row in enumerate(representative_rows) if row.get("is_backbone")), 0)
    ordered_views = sorted(zip(junction_rows, clean_views), key=lambda pair: (junction_locations[str(pair[0]["junction_id"])][1] - backbone_index) % len(representative_rows))
    panels = [Image.open(path).convert("RGB") for _, path in ordered_views]
    overview_path = output_dir / "junction_clean_overview.png"
    overview = Image.new("RGB", (max(panel.width for panel in panels), sum(panel.height for panel in panels) + 32 * (len(panels) - 1)), "#ffffff")
    panel_y = 0
    for panel in panels:
        overview.paste(panel, (0, panel_y))
        panel_y += panel.height + 32
    overview.save(overview_path, "PNG")
    joined_overview_path = output_dir / "junction_joined_overview.png"
    ordered_joined = sorted(zip(junction_rows, joined_views), key=lambda pair: (junction_locations[str(pair[0]["junction_id"])][1] - backbone_index) % len(representative_rows))
    joined_panels = [Image.open(path).convert("RGB") for _, path in ordered_joined]
    joined_overview = Image.new("RGB", (max(panel.width for panel in joined_panels), sum(panel.height for panel in joined_panels) + 32 * (len(joined_panels) - 1)), "#ffffff")
    panel_y = 0
    for panel in joined_panels:
        joined_overview.paste(panel, (0, panel_y))
        panel_y += panel.height + 32
    joined_overview.save(joined_overview_path, "PNG")
    plasmid_junction_figure = output_dir / "plasmid_with_junction_blowups.png"
    _write_plasmid_junction_figure(plasmid_junction_figure, representative_rows, component_starts, retained, assembled_sequence, junction_rows, junction_locations, joined_views)
    sequence_validation = "Exact circular match to uploaded plasmid with selected core replacements" if expected_target is not None else "Explicit-plan overlap assembly; not compared with the full scaffold"
    (output_dir / "design.json").write_text(json.dumps({"plan": plan_reference, "scaffold": str(scaffold.path), "enzyme": enzyme, "junctions": junction_rows, "fragments": fragment_rows, "assembledMap": str(map_path), "schematicMap": str(schematic_path), "sequenceValidation": sequence_validation, "mapVariants": map_variants, "junctionScreencaps": screencaps, "plasmidJunctionFigure": str(plasmid_junction_figure), "automaticJunctions": automatic_junctions, "junctionJoinedViews": joined_views, "junctionJoinedOverview": str(joined_overview_path), "junctionCleanViews": clean_views, "junctionCleanOverview": str(overview_path), "warnings": warnings}, indent=2) + "\n", encoding="utf-8")
    return {"ok": True, "outputDir": str(output_dir), "project": plan.get("project", design_name), "scaffoldBp": len(scaffold.sequence), "assembledBp": len(assembled_sequence), "enzyme": enzyme.get("name", "custom"), "junctions": junction_rows, "fragments": [{key: row[key] for key in ("fragment_id", "name", "element_name", "assembly_position", "variant_name", "source", "left_fusion", "right_fusion", "final_fragment_bp")} for row in fragment_rows], "assembledMap": str(map_path), "schematicMap": str(schematic_path), "sequenceValidation": sequence_validation, "mapVariants": map_variants, "junctionScreencaps": screencaps, "plasmidJunctionFigure": str(plasmid_junction_figure), "automaticJunctions": automatic_junctions, "junctionJoinedViews": joined_views, "junctionJoinedOverview": str(joined_overview_path), "junctionCleanViews": clean_views, "junctionCleanOverview": str(overview_path), "warnings": warnings}


def run_plan(plan_path: str | Path, output_root: str | Path) -> dict[str, Any]:
    """Run a YAML plan and return a compact UI-safe result."""
    plan_path = Path(plan_path).expanduser().resolve()
    plan = yaml.safe_load(plan_path.read_text(encoding="utf-8"))
    if not isinstance(plan, dict):
        raise ValueError("The plan's top-level YAML value must be a mapping")
    raw_scaffold = Path(str(plan["scaffold"])).expanduser()
    scaffold_path = raw_scaffold if raw_scaffold.is_absolute() else plan_path.parent / raw_scaffold
    return _run_plan_data(plan, scaffold_path, output_root, plan_path.stem, str(plan_path))


def run_annotated_snapgene_design(path: str | Path, output_root: str | Path, enzyme: str = "BsmBI", variable_texts: dict[str, str] | None = None) -> dict[str, Any]:
    """Run Golden Gate design directly from a conventionally annotated .dna map."""
    scaffold_path = Path(path).expanduser().resolve()
    plan = copy.deepcopy(plan_from_annotated_snapgene(scaffold_path, enzyme))
    variable_fragments = {item["id"]: item for item in plan["fragments"] if item.get("annotation_kind") == "variable"}
    for element_id, text in (variable_texts or {}).items():
        if element_id not in variable_fragments:
            raise ValueError(f"Unknown variable fragment: {element_id}")
        fragment = variable_fragments[element_id]
        variants = parse_variable_elements(text, fragment["name"])
        if variants:
            fragment["variants"] = variants
    return _run_plan_data(plan, scaffold_path, output_root, scaffold_path.stem, str(scaffold_path))
