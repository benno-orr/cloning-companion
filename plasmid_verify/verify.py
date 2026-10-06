from __future__ import annotations

from collections import Counter
from typing import Dict, List, Sequence, Tuple

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from .fasta import clean_sequence, reverse_complement
from .models import AssemblyResult, BaseAnnotation, Mutation, VerificationResult


def _best_rotation(target: str, observed: str, k: int = 25) -> Tuple[str, str, int]:
    """Orient and rotate a circular consensus using exact anchors spread over target."""
    if not observed:
        return observed, "forward", 0
    anchor_length = min(k, max(8, len(target) // 10), len(target), len(observed))
    best = None
    for orientation, candidate in (("forward", observed), ("reverse complement", reverse_complement(observed))):
        doubled = candidate + candidate
        for target_pos in range(0, max(1, len(target) - anchor_length + 1), max(1, len(target) // 24)):
            anchor = target[target_pos : target_pos + anchor_length]
            observed_pos = doubled.find(anchor)
            if observed_pos < 0 or observed_pos >= len(candidate):
                continue
            rotation = (observed_pos - target_pos) % len(candidate)
            rotated = candidate[rotation:] + candidate[:rotation]
            # A sparse identity score disambiguates repeats without a full alignment.
            sample = min(len(target), len(rotated), 1000)
            score = sum(target[i] == rotated[i] for i in range(sample))
            item = (score, rotated, orientation, rotation)
            if best is None or item[0] > best[0]:
                best = item
    if best is not None:
        return best[1], best[2], best[3]

    # Fallback for very noisy sequences: find the best shared 8-mer.
    short = min(8, len(target), len(observed))
    for orientation, candidate in (("forward", observed), ("reverse complement", reverse_complement(observed))):
        for target_pos in range(max(1, len(target) - short + 1)):
            observed_pos = candidate.find(target[target_pos : target_pos + short])
            if observed_pos >= 0:
                rotation = (observed_pos - target_pos) % len(candidate)
                return candidate[rotation:] + candidate[:rotation], orientation, rotation
    return observed, "forward", 0


def _gapped_alignment(target: str, observed: str) -> Tuple[str, str]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -9
    aligner.extend_gap_score = -1
    alignment = aligner.align(target, observed)[0]
    coordinates = alignment.coordinates
    target_parts: List[str] = []
    observed_parts: List[str] = []
    for column in range(coordinates.shape[1] - 1):
        t0, t1 = int(coordinates[0, column]), int(coordinates[0, column + 1])
        o0, o1 = int(coordinates[1, column]), int(coordinates[1, column + 1])
        target_chunk = target[t0:t1]
        observed_chunk = observed[o0:o1]
        width = max(len(target_chunk), len(observed_chunk))
        target_parts.append(target_chunk.ljust(width, "-"))
        observed_parts.append(observed_chunk.ljust(width, "-"))
    return "".join(target_parts), "".join(observed_parts)


def _translation(codon: str) -> str:
    if len(codon) != 3 or "-" in codon or "N" in codon:
        return "?"
    return str(Seq(codon).translate())


def _coding_effect(
    event_type: str,
    target_start: int,
    target_bases: str,
    observed_bases: str,
    annotation: BaseAnnotation,
    target: str,
    target_to_observed: Dict[int, str],
) -> Tuple[str, str]:
    if event_type in {"insertion", "deletion"}:
        changed_length = len(observed_bases) if event_type == "insertion" else len(target_bases)
        if changed_length % 3:
            return "frameshift", f"frameshift at nt {target_start + 1}"
        return f"in-frame {event_type}", f"in-frame {event_type} ({changed_length} nt)"
    feature_start = annotation.feature_start or 0
    local_position = target_start - feature_start
    codon_start = target_start - (local_position % 3)
    expected_codon = target[codon_start : codon_start + 3]
    observed_codon = "".join(target_to_observed.get(i, "-") for i in range(codon_start, codon_start + 3))
    expected_aa = _translation(expected_codon)
    observed_aa = _translation(observed_codon)
    aa_position = (local_position // 3) + 1
    if expected_aa == observed_aa and expected_aa != "?":
        return "silent", f"{expected_aa}{aa_position}="
    if observed_aa == "*" and expected_aa != "*":
        return "nonsense", f"{expected_aa}{aa_position}*"
    if expected_aa == "*" and observed_aa != "*":
        return "stop loss", f"*{aa_position}{observed_aa}"
    return "missense", f"{expected_aa}{aa_position}{observed_aa}"


def _event_region(annotations: Sequence[BaseAnnotation], start: int, end: int) -> BaseAnnotation:
    if not annotations:
        return BaseAnnotation("unknown", "unknown")
    indices = range(max(0, start), min(len(annotations), max(start + 1, end)))
    choices = [annotations[index] for index in indices]
    if not choices:
        choices = [annotations[min(max(start - 1, 0), len(annotations) - 1)]]
    # Junction takes precedence for boundary-spanning events, then coding.
    for preferred in ("junction", "coding", "non-coding", "backbone"):
        for annotation in choices:
            if annotation.region == preferred:
                return annotation
    return choices[0]


def _mutations_from_alignment(
    target: str,
    aligned_target: str,
    aligned_observed: str,
    annotations: Sequence[BaseAnnotation],
) -> List[Mutation]:
    target_index = 0
    target_positions: List[int] = []
    target_to_observed: Dict[int, str] = {}
    for expected, observed in zip(aligned_target, aligned_observed):
        target_positions.append(target_index)
        if expected != "-":
            target_to_observed[target_index] = observed
            target_index += 1

    events: List[Mutation] = []
    index = 0
    while index < len(aligned_target):
        expected = aligned_target[index]
        observed = aligned_observed[index]
        if expected == observed:
            index += 1
            continue
        event_type = "substitution"
        if expected == "-":
            event_type = "insertion"
        elif observed == "-":
            event_type = "deletion"
        end = index + 1
        while end < len(aligned_target):
            next_expected, next_observed = aligned_target[end], aligned_observed[end]
            next_type = "substitution"
            if next_expected == "-":
                next_type = "insertion"
            elif next_observed == "-":
                next_type = "deletion"
            if next_expected == next_observed or next_type != event_type:
                break
            # Keep substitutions as individual SNVs so codon effects remain clear.
            if event_type == "substitution":
                break
            end += 1
        target_start = target_positions[index]
        target_bases = aligned_target[index:end].replace("-", "")
        observed_bases = aligned_observed[index:end].replace("-", "")
        target_end = target_start + len(target_bases)
        annotation = _event_region(annotations, target_start, target_end)
        category = annotation.region
        protein_change = ""
        if annotation.region == "coding":
            category, protein_change = _coding_effect(
                event_type,
                target_start,
                target_bases,
                observed_bases,
                annotation,
                target,
                target_to_observed,
            )
        context_start = max(0, target_start - 8)
        context_end = min(len(target), max(target_end, target_start + 1) + 8)
        events.append(
            Mutation(
                event_type,
                category,
                annotation.region,
                target_start + 1,
                target_bases or "–",
                observed_bases or "–",
                annotation.source,
                protein_change,
                target[context_start:context_end],
            )
        )
        index = end
    return events


def _verdict(mutations: Sequence[Mutation]) -> str:
    if not mutations:
        return "PASS"
    failing = {"frameshift", "nonsense", "stop loss", "missense", "in-frame insertion", "in-frame deletion", "junction"}
    if any(mutation.category in failing for mutation in mutations):
        return "FAIL"
    return "REVIEW"


def verify_consensus(assembly: AssemblyResult, observed: str) -> VerificationResult:
    observed = clean_sequence(observed)
    rotated, orientation, rotation = _best_rotation(assembly.sequence, observed)
    aligned_target, aligned_observed = _gapped_alignment(assembly.sequence, rotated)
    mutations = _mutations_from_alignment(
        assembly.sequence, aligned_target, aligned_observed, assembly.annotations
    )
    return VerificationResult(
        assembly.target_id,
        _verdict(mutations),
        len(assembly.sequence),
        len(observed),
        orientation,
        rotation,
        mutations,
        aligned_target,
        aligned_observed,
    )

