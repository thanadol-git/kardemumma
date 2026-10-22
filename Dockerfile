# syntax=docker/dockerfile:1.6
FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_ROOT_USER_ACTION=ignore \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/matplotlib

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        procps && \
    rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md LICENSE /opt/kardemumma/
COPY src /opt/kardemumma/src
RUN pip install --no-cache-dir /opt/kardemumma

WORKDIR /work

CMD ["kardemumma-preview", "--help"]
