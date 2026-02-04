from __future__ import annotations

from pathlib import Path

from loguru import logger

from viral_pipeline.config import AppConfig
from viral_pipeline.utils import ClusteringResult, ensure_dir


class Clusterer:
    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def cluster(self, contigs: Path, output_dir: Path) -> ClusteringResult:
        ensure_dir(output_dir)
        clusters = output_dir / "clusters.tsv"
        abundance = output_dir / "abundance.tsv"
        logger.info(
            "Simulating clustering with ANI {}", self.config.clustering.ani_threshold
        )
        clusters.write_text("cluster\tcontig\ncluster1\tcontig_1\n", encoding="utf-8")
        abundance.write_text("cluster\tabundance\ncluster1\t1\n", encoding="utf-8")
        return ClusteringResult(clusters=clusters, abundance_table=abundance)
