from __future__ import annotations

from pathlib import Path
from typing import List

from loguru import logger

from viral_pipeline.config import AppConfig
from viral_pipeline.utils import AssemblyResult, ensure_dir


class Assembler:
    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def assemble(self, reads: List[Path], output_dir: Path) -> AssemblyResult:
        ensure_dir(output_dir)
        contigs = output_dir / "contigs.fasta"
        report = output_dir / "assembly_report.json"
        logger.info("Simulating assembly with {} reads", len(reads))
        contigs.write_text(">contig_1\nATGC\n", encoding="utf-8")
        report.write_text("{}", encoding="utf-8")
        return AssemblyResult(contigs=contigs, report=report)
