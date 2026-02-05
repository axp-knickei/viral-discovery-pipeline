from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

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


def check_dependency(tool_name: str) -> Path:
    """Check if a tool is available in the system PATH."""
    path = shutil.which(tool_name)
    if not path:
        raise RuntimeError(f"Required tool '{tool_name}' not found in PATH.")
    return Path(path)


def group_paired_reads(files: List[Path]) -> Dict[str, Tuple[Path, Optional[Path]]]:
    """
    Groups reads into R1 and R2 pairs based on filename patterns.
    Returns a dictionary mapping sample ID to (R1, R2) tuple.
    R2 can be None for single-end reads.
    """
    groups: Dict[str, Dict[str, Path]] = {}

    for file_path in files:
        name = file_path.name
        # Heuristic matching for paired files
        if "_R1" in name:
            sample_id = name.split("_R1")[0]
            read_type = "R1"
        elif "_R2" in name:
            sample_id = name.split("_R2")[0]
            read_type = "R2"
        elif "_1." in name: # e.g. sample_1.fastq
            sample_id = name.split("_1.")[0]
            read_type = "R1"
        elif "_2." in name:
            sample_id = name.split("_2.")[0]
            read_type = "R2"
        else:
            # Assume single end or unidentified pattern
            sample_id = file_path.stem.split(".")[0] # clean name
            read_type = "R1" # Treat as R1 (single)

        if sample_id not in groups:
            groups[sample_id] = {}
        groups[sample_id][read_type] = file_path

    result = {}
    for sample_id, pairs in groups.items():
        r1 = pairs.get("R1")
        r2 = pairs.get("R2")
        if r1:
            result[sample_id] = (r1, r2)
        elif r2 and not r1:
            # Only R2 found? Treat as single end or error?
            # Assign R2 as primary if R1 missing is unlikely but for now let's skip or log
            logger.warning(f"Found R2 but no R1 for sample {sample_id}. Skipping.")

    return result
