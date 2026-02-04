from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

import yaml
from loguru import logger
from pydantic import BaseModel, Field, ValidationError


class DatabaseConfig(BaseModel):
    host_genome: Path
    dvf_models: Path
    genomad_db: Path
    checkv_db: Path


class ToolsConfig(BaseModel):
    threads: int = Field(ge=1)
    memory_gb: int = Field(ge=1)
    min_contig_length: int = Field(ge=1)


class QualityConfig(BaseModel):
    quality_threshold: int = Field(ge=0)
    min_read_length: int = Field(ge=0)
    adapter_removal: bool = True


class ClusteringConfig(BaseModel):
    ani_threshold: float = Field(ge=0.0, le=1.0)
    min_coverage: float = Field(ge=0.0, le=1.0)


class LoggingConfig(BaseModel):
    level: str = "INFO"
    log_dir: Path = Path("logs")


class AppConfig(BaseModel):
    database: DatabaseConfig
    tools: ToolsConfig
    quality: QualityConfig
    clustering: ClusteringConfig
    logging: LoggingConfig = LoggingConfig()


@dataclass
class ConfigPaths:
    root: Path
    default: Path
    development: Path
    production: Path

    @classmethod
    def from_repo(cls, repo_root: Path) -> "ConfigPaths":
        config_root = repo_root / "config"
        return cls(
            root=config_root,
            default=config_root / "default.yaml",
            development=config_root / "development.yaml",
            production=config_root / "production.yaml",
        )


def _read_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid config structure in {path}")
    return data


def _merge_dicts(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    result = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _merge_dicts(result[key], value)
        else:
            result[key] = value
    return result


def load_config(default_path: Path, env_path: Path | None = None) -> AppConfig:
    base = _read_yaml(default_path)
    if env_path:
        base = _merge_dicts(base, _read_yaml(env_path))
    logger.debug("Loaded configuration", config=base)
    return AppConfig.model_validate(base)


def validate_config(config_path: Path) -> AppConfig:
    try:
        return load_config(config_path)
    except ValidationError as exc:
        logger.error("Configuration validation failed: {}", exc)
        raise
