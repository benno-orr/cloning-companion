from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass(frozen=True)
class BaseAnnotation:
    region: str
    source: str
    feature_start: Optional[int] = None
    feature_end: Optional[int] = None


@dataclass
class AssemblyResult:
    target_id: str
    sequence: str
    annotations: List[BaseAnnotation]
    method: str
    notes: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.annotations):
            raise ValueError("Every target base must have exactly one annotation")


@dataclass
class Mutation:
    mutation_type: str
    category: str
    region: str
    target_position: int
    target_bases: str
    observed_bases: str
    source: str
    protein_change: str = ""
    context: str = ""

    def as_dict(self) -> Dict[str, object]:
        return {
            "position": self.target_position,
            "mutation": self.mutation_type,
            "category": self.category,
            "region": self.region,
            "source": self.source,
            "target": self.target_bases,
            "observed": self.observed_bases,
            "protein_change": self.protein_change,
            "context": self.context,
        }


@dataclass
class VerificationResult:
    target_id: str
    verdict: str
    expected_length: int
    observed_length: int
    orientation: str
    rotation: int
    mutations: List[Mutation]
    aligned_target: str
    aligned_observed: str

    @property
    def counts(self) -> Dict[str, int]:
        result: Dict[str, int] = {}
        for mutation in self.mutations:
            result[mutation.category] = result.get(mutation.category, 0) + 1
        return result

