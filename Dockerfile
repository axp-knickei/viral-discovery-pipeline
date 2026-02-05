# syntax=docker/dockerfile:1
FROM python:3.11-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

FROM base AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
 && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md /app/
COPY src/ /app/src/
RUN pip install --upgrade pip \
 && pip wheel --no-cache-dir --wheel-dir /wheels .

FROM base AS runtime

COPY --from=builder /wheels /wheels
COPY pyproject.toml README.md /app/
RUN pip install --no-cache-dir /wheels/*

COPY src/ /app/src/
RUN pip install -e /app

ENTRYPOINT ["viral-pipeline"]
CMD ["--help"]
