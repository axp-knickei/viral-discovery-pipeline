from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from loguru import logger
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from viral_pipeline import load_config
from viral_pipeline.assembly import Assembler
from viral_pipeline.clustering import Clusterer
from viral_pipeline.config import ConfigPaths
from viral_pipeline.host_removal import HostRemoval
from viral_pipeline.qc import QualityController
from viral_pipeline.viral_detection import ViralDetector

app = typer.Typer(help="Modern viral discovery pipeline")
console = Console()


def _configure_logging(verbose: bool) -> None:
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(lambda msg: console.print(msg, end=""), level=level)


@app.command()
def run(
    input_dir: Path = typer.Option(..., help="Directory with FASTQ inputs"),
    output_dir: Path = typer.Option(Path("outputs"), help="Output directory"),
    config_path: Optional[Path] = typer.Option(None, help="Config YAML file"),
    env_config: Optional[Path] = typer.Option(None, help="Override config YAML"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show steps without running"),
    resume: bool = typer.Option(False, "--resume", help="Resume from existing outputs"),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose logging"),
) -> None:
    _configure_logging(verbose)
    config_paths = ConfigPaths.from_repo(Path.cwd())
    config = load_config(config_path or config_paths.default, env_config)

    if dry_run:
        console.print("[yellow]Dry run enabled: no outputs will be generated.[/yellow]")
        return

    input_files = sorted(input_dir.glob("*.fastq")) + sorted(input_dir.glob("*.fq"))
    if not input_files:
        raise typer.BadParameter("No FASTQ files found in input directory")

    qc = QualityController(config)
    host_removal = HostRemoval(config)
    assembler = Assembler(config)
    detector = ViralDetector(config)
    clusterer = Clusterer(config)

    with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as progress:
        task = progress.add_task("Running QC", start=True)
        qc_result = qc.run_quality_control(input_files, output_dir / "qc")
        progress.update(task, description="Removing host reads")
        host_result = host_removal.filter_host(qc_result.cleaned_reads, output_dir / "host")
        progress.update(task, description="Assembling contigs")
        assembly_result = assembler.assemble(host_result.filtered_reads, output_dir / "assembly")
        progress.update(task, description="Detecting viral sequences")
        viral_result = detector.identify_viral_sequences(
            assembly_result.contigs, output_dir / "viral"
        )
        progress.update(task, description="Clustering vOTUs")
        clusterer.cluster(viral_result.contigs, output_dir / "clustering")
        progress.update(task, description="Pipeline complete", completed=1)

    if resume:
        console.print("[green]Resume flag acknowledged (no-op for now).[/green]")


@app.command()
def init(
    output_path: Path = typer.Option(Path("config.yaml"), help="Where to write config"),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose logging"),
) -> None:
    _configure_logging(verbose)
    config_paths = ConfigPaths.from_repo(Path.cwd())
    template = config_paths.default.read_text(encoding="utf-8")
    output_path.write_text(template, encoding="utf-8")
    console.print(f"[green]Wrote configuration template to {output_path}[/green]")


@app.command()
def validate(
    config_path: Path = typer.Option(..., help="Path to config YAML"),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose logging"),
) -> None:
    _configure_logging(verbose)
    load_config(config_path)
    console.print(f"[green]Configuration is valid: {config_path}[/green]")


@app.command()
def config(
    config_path: Optional[Path] = typer.Option(None, help="Config YAML to load"),
    env_config: Optional[Path] = typer.Option(None, help="Optional override config"),
    validate: bool = typer.Option(False, "--validate", help="Validate configuration"),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose logging"),
) -> None:
    _configure_logging(verbose)
    config_paths = ConfigPaths.from_repo(Path.cwd())
    path = config_path or config_paths.default
    if validate:
        load_config(path, env_config)
        console.print(f"[green]Configuration is valid: {path}[/green]")
    else:
        config = load_config(path, env_config)
        console.print(config.model_dump_json(indent=2))


if __name__ == "__main__":
    app()
