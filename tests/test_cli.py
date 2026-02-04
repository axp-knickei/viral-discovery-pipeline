from pathlib import Path

from typer.testing import CliRunner

from viral_pipeline.cli import app


runner = CliRunner()


def test_cli_validate(tmp_path: Path) -> None:
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
        """,
        encoding="utf-8",
    )
    result = runner.invoke(app, ["validate", "--config-path", str(config_path)])
    assert result.exit_code == 0
    assert "Configuration is valid" in result.stdout
