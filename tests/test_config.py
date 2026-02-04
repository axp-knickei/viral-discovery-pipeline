from pathlib import Path

import pytest

from viral_pipeline.config import AppConfig, load_config


def test_load_config_default(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
        database:
          host_genome: "/data/db/hg38/hg38.fa.gz"
          dvf_models: "/data/db/DeepVirFinder/models"
          genomad_db: "/data/db/genomad_db"
          checkv_db: "/data/db/checkv-db-v1.4"
        tools:
          threads: 4
          memory_gb: 16
          min_contig_length: 5000
        quality:
          quality_threshold: 20
          min_read_length: 30
          adapter_removal: true
        clustering:
          ani_threshold: 0.95
          min_coverage: 0.85
        logging:
          level: "INFO"
          log_dir: "logs"
        """,
        encoding="utf-8",
    )
    config = load_config(config_path)
    assert isinstance(config, AppConfig)
    assert config.tools.threads == 4


def test_invalid_config(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
        database:
          host_genome: "/data/db/hg38/hg38.fa.gz"
          dvf_models: "/data/db/DeepVirFinder/models"
          genomad_db: "/data/db/genomad_db"
          checkv_db: "/data/db/checkv-db-v1.4"
        tools:
          threads: 0
          memory_gb: 16
          min_contig_length: 5000
        quality:
          quality_threshold: 20
          min_read_length: 30
          adapter_removal: true
        clustering:
          ani_threshold: 0.95
          min_coverage: 0.85
        """,
        encoding="utf-8",
    )
    with pytest.raises(Exception):
        load_config(config_path)
