from __future__ import annotations

import re
from io import BytesIO
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple, Union

from Bio import SeqIO


DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str


def clean_sequence(sequence: str) -> str:
    sequence = re.sub(r"\s+", "", sequence).upper().replace("U", "T")
    if not sequence:
        raise ValueError("Sequence is empty")
    if not DNA_RE.fullmatch(sequence):
        bad = sorted(set(sequence) - set("ACGTRYSWKMBDHVN"))
        raise ValueError(f"Unsupported sequence characters: {', '.join(bad)}")
    return sequence


def parse_fasta(data: Union[str, bytes], fallback_name: str = "sequence") -> List[FastaRecord]:
    if isinstance(data, bytes):
        text = data.decode("utf-8-sig")
    else:
        text = data
    records: List[FastaRecord] = []
    name = fallback_name
    parts: List[str] = []
    saw_header = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if parts:
                records.append(FastaRecord(name, clean_sequence("".join(parts))))
                parts = []
            name = line[1:].strip().split()[0] or fallback_name
            saw_header = True
        else:
            parts.append(line)
    if parts:
        records.append(FastaRecord(name, clean_sequence("".join(parts))))
    if not records:
        raise ValueError("No FASTA sequence found")
    if not saw_header and len(records) == 1:
        records[0] = FastaRecord(fallback_name, records[0].sequence)
    return records


def parse_sequence_file(data: bytes, filename: str) -> List[FastaRecord]:
    """Read FASTA or SnapGene .dna, retaining only the DNA sequence."""
    if Path(filename).suffix.lower() == ".dna":
        try:
            record = SeqIO.read(BytesIO(data), "snapgene")
        except (ValueError, UnicodeError, EOFError) as exc:
            raise ValueError(f"Invalid SnapGene .dna file: {exc}") from exc
        name = str(record.id or record.name or Path(filename).stem)
        if name == "<unknown id>":
            name = Path(filename).stem
        return [FastaRecord(name, clean_sequence(str(record.seq)))]
    return parse_fasta(data, fallback_name=Path(filename).stem)


def first_record(upload) -> FastaRecord:
    records = parse_sequence_file(upload.getvalue(), upload.name)
    return records[0]


def infer_target_id(filename: str) -> str:
    """Infer the target prefix from a convention like T01_consensus.fasta."""
    stem = Path(filename).stem
    stem = re.sub(
        r"(?i)(?:[_\-.](?:consensus|wps|plasmid|sequence|seq|insert|fragment|part))(?:[_\-.].*)?$",
        "",
        stem,
    )
    if "__" in stem:
        stem = stem.split("__", 1)[0]
    if "_" in stem:
        stem = stem.split("_", 1)[0]
    return stem.strip("_.- ") or Path(filename).stem


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return sequence.upper().translate(table)[::-1]


def to_fasta(name: str, sequence: str, width: int = 80) -> str:
    lines = [f">{name}"]
    lines.extend(sequence[i : i + width] for i in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"
