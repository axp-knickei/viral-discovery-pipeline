from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from loguru import logger


@dataclass
class QCResult:
    cleaned_reads: list[Path]
    report: Path


@dataclass
class ViralResult:
    contigs: Path
    report: Path


@dataclass
class AssemblyResult:
    contigs: Path
    report: Path


@dataclass
class HostRemovalResult:
    filtered_reads: list[Path]
    report: Path


@dataclass
class ClusteringResult:
    clusters: Path
    abundance_table: Path


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def validate_fastq(path: Path, min_read_length: int) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input FASTQ not found: {path}")
    if path.suffix not in {".fastq", ".fq", ".fastq.gz", ".fq.gz"}:
        raise ValueError(f"Unsupported FASTQ extension: {path}")
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for idx, line in enumerate(handle, start=1):
            if idx % 4 == 2:
                if len(line.strip()) < min_read_length:
                    raise ValueError(
                        f"Read length below minimum in {path} at record {idx // 4}"
                    )
            if idx > 40:
                break
    logger.debug("Validated FASTQ: {}", path)


def verify_outputs(outputs: Iterable[Path]) -> None:
    missing = [path for path in outputs if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing expected outputs: {missing}")
    logger.debug("Verified outputs: {}", outputs)
