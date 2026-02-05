from __future__ import annotations

from pathlib import Path

from loguru import logger

from viral_pipeline.config import AppConfig
from viral_pipeline.utils import ViralResult, ensure_dir


class ViralDetector:
    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def identify_viral_sequences(self, contigs: Path, output_dir: Path) -> ViralResult:
        ensure_dir(output_dir)
        viral_contigs = output_dir / "viral_contigs.fasta"
        report = output_dir / "viral_detection_report.json"
        logger.info("Simulating viral detection using {}", self.config.database.genomad_db)
        viral_contigs.write_text(">viral_contig_1\nATGC\n", encoding="utf-8")
        report.write_text("{}", encoding="utf-8")
        return ViralResult(contigs=viral_contigs, report=report)
