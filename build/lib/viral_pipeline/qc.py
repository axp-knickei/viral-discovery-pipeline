from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

from loguru import logger

from viral_pipeline.config import AppConfig
from viral_pipeline.utils import QCResult, ensure_dir, validate_fastq


@dataclass
class QCInputs:
    input_files: List[Path]
    output_dir: Path


class QualityController:
    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def run_quality_control(self, input_files: List[Path], output_dir: Path) -> QCResult:
        ensure_dir(output_dir)
        for path in input_files:
            validate_fastq(path, self.config.quality.min_read_length)
        cleaned = [output_dir / f"{path.stem}.clean.fastq" for path in input_files]
        report = output_dir / "qc_report.json"
        logger.info("Simulating QC for {} inputs", len(input_files))
        for path in cleaned:
            path.write_text("", encoding="utf-8")
        report.write_text("{}", encoding="utf-8")
        return QCResult(cleaned_reads=cleaned, report=report)
