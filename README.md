# Viral Discovery Pipeline

A modern, configurable bioinformatics pipeline for viral discovery with typed configuration, a Python CLI,
container recipes, and workflow orchestration support.

## Quickstart

```bash
pip install -e .
viral-pipeline --help
viral-pipeline init --output-path config.yaml
viral-pipeline config --validate --config-path config.yaml
viral-pipeline run --input-dir data/ --output-dir outputs/
```

## Configuration

Configuration lives in `config/` and is validated with Pydantic. Use the CLI to validate or render the
merged configuration.

```bash
viral-pipeline config --validate --config-path config/default.yaml
```

## Development

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .[all]
pytest
```

## Documentation

- MkDocs user guides live in `docs/`
- Sphinx API documentation lives in `docs/api/`

## Containers

- `Dockerfile` for multi-stage builds
- `docker-compose.yml` for local dev
- `Singularity.def` for HPC environments
