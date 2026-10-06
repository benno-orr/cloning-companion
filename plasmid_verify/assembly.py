from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .fasta import clean_sequence, reverse_complement
from .models import AssemblyResult, BaseAnnotation


class AssemblyError(ValueError):
    pass


@dataclass(frozen=True)
class TypeIISEnzyme:
    name: str
    recognition: str
    top_offset: int
    bottom_offset: int

    @property
    def overhang_length(self) -> int:
        return self.bottom_offset - self.top_offset


TYPE_IIS_ENZYMES: Dict[str, TypeIISEnzyme] = {
    "BsaI": TypeIISEnzyme("BsaI", "GGTCTC", 1, 5),
    "BsmBI / Esp3I": TypeIISEnzyme("BsmBI / Esp3I", "CGTCTC", 1, 5),
    "BbsI": TypeIISEnzyme("BbsI", "GAAGAC", 2, 6),
    "SapI": TypeIISEnzyme("SapI", "GCTCTTC", 1, 4),
    "AarI": TypeIISEnzyme("AarI", "CACCTGC", 4, 8),
}


def _all_occurrences(sequence: str, motif: str) -> List[int]:
    hits: List[int] = []
    start = 0
    while True:
        hit = sequence.find(motif, start)
        if hit < 0:
            return hits
        hits.append(hit)
        start = hit + 1


def _circular_intervals(sequence: str, left: str, right: str) -> List[Tuple[int, int, bool]]:
    """Return left-start/right-end intervals, allowing a single origin crossing."""
    n = len(sequence)
    doubled = sequence + sequence
    candidates: List[Tuple[int, int, bool]] = []
    for left_start in _all_occurrences(sequence, left):
        search_from = left_start + len(left)
        right_start = doubled.find(right, search_from, left_start + n + len(right))
        if right_start >= 0:
            right_end = right_start + len(right)
            candidates.append((left_start, right_end % n, right_end > n))
    return candidates


def _choose_arm_interval(sequence: str, insert: str, arm_length: int) -> Tuple[int, int, bool, str]:
    if len(insert) < arm_length * 2:
        raise AssemblyError(f"Insert is shorter than two {arm_length}-nt homology arms")
    possibilities: List[Tuple[int, int, bool, str, int]] = []
    for orientation, candidate in (("forward", insert), ("reverse complement", reverse_complement(insert))):
        left, right = candidate[:arm_length], candidate[-arm_length:]
        for start, end, wraps in _circular_intervals(sequence, left, right):
            span = (end - start) % len(sequence)
            if span == 0:
                span = len(sequence)
            possibilities.append((start, end, wraps, orientation, span))
    if not possibilities:
        raise AssemblyError(
            f"Could not find both {arm_length}-nt terminal homology arms in the current plasmid"
        )
    possibilities.sort(key=lambda item: item[4])
    best = possibilities[0]
    if len(possibilities) > 1 and possibilities[1][4] == best[4] and possibilities[1][:2] != best[:2]:
        raise AssemblyError("Homology arms map ambiguously to more than one equally short interval")
    return best[0], best[1], best[2], best[3]


def _insert_annotations(
    name: str,
    length: int,
    junction_length: int,
    noncoding_5: int = 0,
    noncoding_3: int = 0,
) -> List[BaseAnnotation]:
    annotations: List[BaseAnnotation] = []
    feature_start = min(junction_length + max(0, noncoding_5), length)
    feature_end = max(feature_start, length - junction_length - max(0, noncoding_3))
    for index in range(length):
        if index < junction_length or index >= length - junction_length:
            region = "junction"
        elif index < feature_start or index >= feature_end:
            region = "non-coding"
        else:
            region = "coding"
        annotations.append(BaseAnnotation(region, name, feature_start, feature_end))
    return annotations


def _replace_circular(
    sequence: str,
    annotations: List[BaseAnnotation],
    start: int,
    end: int,
    wraps: bool,
    replacement: str,
    replacement_annotations: List[BaseAnnotation],
) -> Tuple[str, List[BaseAnnotation]]:
    if wraps:
        # Rotate at the left arm. This keeps the assembled locus contiguous and is
        # biologically equivalent for a circular plasmid.
        rotated_sequence = sequence[start:] + sequence[:start]
        rotated_annotations = annotations[start:] + annotations[:start]
        interval_length = (end - start) % len(sequence)
        return (
            replacement + rotated_sequence[interval_length:],
            replacement_annotations + rotated_annotations[interval_length:],
        )
    return (
        sequence[:start] + replacement + sequence[end:],
        annotations[:start] + replacement_annotations + annotations[end:],
    )


def _normalize_feature_coordinates(annotations: List[BaseAnnotation]) -> List[BaseAnnotation]:
    bounds: Dict[str, Tuple[int, int]] = {}
    for index, annotation in enumerate(annotations):
        if annotation.region != "coding":
            continue
        if annotation.source not in bounds:
            bounds[annotation.source] = (index, index + 1)
        else:
            bounds[annotation.source] = (bounds[annotation.source][0], index + 1)
    return [
        BaseAnnotation(
            annotation.region,
            annotation.source,
            bounds.get(annotation.source, (None, None))[0],
            bounds.get(annotation.source, (None, None))[1],
        )
        for annotation in annotations
    ]


def assemble_gibson(
    parent: str,
    inserts: Sequence[Tuple[str, str]],
    target_id: str,
    arm_length: int = 15,
    noncoding: Optional[Dict[str, Tuple[int, int]]] = None,
) -> AssemblyResult:
    if not inserts:
        raise AssemblyError("At least one insert is required")
    sequence = clean_sequence(parent)
    annotations = [BaseAnnotation("backbone", "parent") for _ in sequence]
    notes: List[str] = []
    for name, raw_insert in inserts:
        insert = clean_sequence(raw_insert)
        start, end, wraps, orientation = _choose_arm_interval(sequence, insert, arm_length)
        if orientation == "reverse complement":
            insert = reverse_complement(insert)
        noncoding_5, noncoding_3 = (noncoding or {}).get(name, (0, 0))
        replacement_annotations = _insert_annotations(
            name, len(insert), arm_length, noncoding_5, noncoding_3
        )
        sequence, annotations = _replace_circular(
            sequence, annotations, start, end, wraps, insert, replacement_annotations
        )
        notes.append(
            f"{name}: {orientation}; replaced the circular interval bracketed by "
            f"{arm_length}-nt terminal homology arms"
        )
    return AssemblyResult(
        target_id, sequence, _normalize_feature_coordinates(annotations), "Gibson", notes
    )


@dataclass(frozen=True)
class _GoldenGatePart:
    name: str
    clean_sequence: str
    left_overhang: str
    right_overhang: str
    core_start: int
    core_end: int
    orientation: str


def _extract_golden_gate_part(raw: str, name: str, enzyme: TypeIISEnzyme) -> _GoldenGatePart:
    sequence = clean_sequence(raw)
    motif = enzyme.recognition
    reverse_motif = reverse_complement(motif)
    candidates: List[_GoldenGatePart] = []
    for orientation, candidate in (("forward", sequence), ("reverse complement", reverse_complement(sequence))):
        for forward_start in _all_occurrences(candidate, motif):
            reverse_start = candidate.find(reverse_motif, forward_start + len(motif))
            if reverse_start < 0:
                continue
            clean_start = forward_start + len(motif) + enzyme.top_offset
            clean_end = reverse_start - enzyme.top_offset
            left_end = forward_start + len(motif) + enzyme.bottom_offset
            right_start = reverse_start - enzyme.bottom_offset
            if clean_start >= clean_end or left_end > clean_end or right_start < clean_start:
                continue
            left = candidate[clean_start:left_end]
            right = candidate[right_start:clean_end]
            clean = candidate[clean_start:clean_end]
            if len(left) == enzyme.overhang_length and len(right) == enzyme.overhang_length:
                candidates.append(
                    _GoldenGatePart(
                        name,
                        clean,
                        left,
                        right,
                        len(left),
                        len(clean) - len(right),
                        orientation,
                    )
                )
    if not candidates:
        raise AssemblyError(
            f"{name}: no inward-facing {enzyme.name} site pair was found. "
            "Include both Type IIS sites and their fusion overhangs in each insert FASTA."
        )
    candidates.sort(key=lambda part: len(part.clean_sequence), reverse=True)
    return candidates[0]


def _find_vector_cassette(parent: str, enzyme: TypeIISEnzyme) -> Tuple[str, int, int, str, str]:
    motif = enzyme.recognition
    reverse_motif = reverse_complement(motif)
    candidates: List[Tuple[int, str, int, int, str, str]] = []
    for forward_start in _all_occurrences(parent, motif):
        rotated = parent[forward_start:] + parent[:forward_start]
        reverse_start = rotated.find(reverse_motif, len(motif))
        if reverse_start < 0:
            continue
        clean_start = len(motif) + enzyme.top_offset
        clean_end = reverse_start - enzyme.top_offset
        left = rotated[clean_start : len(motif) + enzyme.bottom_offset]
        right = rotated[reverse_start - enzyme.bottom_offset : clean_end]
        if len(left) == enzyme.overhang_length and len(right) == enzyme.overhang_length:
            cassette_end = reverse_start + len(reverse_motif)
            candidates.append((cassette_end, rotated, 0, cassette_end, left, right))
    if not candidates:
        raise AssemblyError(
            f"Parent plasmid has no inward-facing {enzyme.name} site pair defining a destination cassette"
        )
    candidates.sort(key=lambda item: item[0])
    _, rotated, start, end, left, right = candidates[0]
    return rotated, start, end, left, right


def assemble_golden_gate(
    parent: str,
    inserts: Sequence[Tuple[str, str]],
    target_id: str,
    enzyme_name: str = "BsaI",
    noncoding: Optional[Dict[str, Tuple[int, int]]] = None,
) -> AssemblyResult:
    if not inserts:
        raise AssemblyError("At least one insert is required")
    if enzyme_name not in TYPE_IIS_ENZYMES:
        raise AssemblyError(f"Unsupported Type IIS enzyme: {enzyme_name}")
    enzyme = TYPE_IIS_ENZYMES[enzyme_name]
    parent = clean_sequence(parent)
    rotated_parent, cassette_start, cassette_end, vector_left, vector_right = _find_vector_cassette(
        parent, enzyme
    )
    parts = [_extract_golden_gate_part(sequence, name, enzyme) for name, sequence in inserts]
    if parts[0].left_overhang != vector_left:
        raise AssemblyError(
            f"First insert 5′ overhang {parts[0].left_overhang} does not match vector overhang {vector_left}"
        )
    for previous, current in zip(parts, parts[1:]):
        if previous.right_overhang != current.left_overhang:
            raise AssemblyError(
                f"Incompatible fusion overhangs: {previous.name} ends {previous.right_overhang}, "
                f"but {current.name} starts {current.left_overhang}"
            )
    if parts[-1].right_overhang != vector_right:
        raise AssemblyError(
            f"Last insert 3′ overhang {parts[-1].right_overhang} does not match vector overhang {vector_right}"
        )

    assembled = parts[0].clean_sequence
    first_noncoding_5, first_noncoding_3 = (noncoding or {}).get(parts[0].name, (0, 0))
    assembled_annotations = _insert_annotations(
        parts[0].name,
        len(parts[0].clean_sequence),
        enzyme.overhang_length,
        first_noncoding_5,
        first_noncoding_3,
    )
    for part in parts[1:]:
        # The fusion overhang is present on both design fragments but only once in
        # the ligated molecule.
        trim = enzyme.overhang_length
        noncoding_5, noncoding_3 = (noncoding or {}).get(part.name, (0, 0))
        part_annotations = _insert_annotations(
            part.name, len(part.clean_sequence), trim, noncoding_5, noncoding_3
        )
        assembled += part.clean_sequence[trim:]
        assembled_annotations += part_annotations[trim:]

    backbone = rotated_parent[cassette_end:]
    sequence = assembled + backbone
    annotations = assembled_annotations + [BaseAnnotation("backbone", "parent") for _ in backbone]
    notes = [
        f"{part.name}: {part.orientation}; {part.left_overhang} → {part.right_overhang}"
        for part in parts
    ]
    notes.append(
        f"Removed the {enzyme.name} destination cassette and both recognition sites; "
        f"assembled {len(parts)} insert(s) by compatible {enzyme.overhang_length}-nt fusion overhangs"
    )
    return AssemblyResult(
        target_id,
        sequence,
        _normalize_feature_coordinates(annotations),
        f"Golden Gate ({enzyme.name})",
        notes,
    )
