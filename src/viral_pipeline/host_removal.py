from __future__ import annotations

from pathlib import Path
from typing import List

from loguru import logger

from viral_pipeline.config import AppConfig
from viral_pipeline.utils import HostRemovalResult, ensure_dir


class HostRemoval:
    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def filter_host(self, input_files: List[Path], output_dir: Path) -> HostRemovalResult:
        ensure_dir(output_dir)
        filtered = [output_dir / f"{path.stem}.host_filtered.fastq" for path in input_files]
        report = output_dir / "host_removal_report.json"
        logger.info("Simulating host removal with database {}", self.config.database.host_genome)
        for path in filtered:
            path.write_text("", encoding="utf-8")
        report.write_text("{}", encoding="utf-8")
        return HostRemovalResult(filtered_reads=filtered, report=report)
