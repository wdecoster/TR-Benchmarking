# ATaRVa from PyPI. Dependencies are resolved as they were on the release date (uv --exclude-newer),
# so the image reproduces the environment a user installing that release would have had.
FROM python:3.11-slim-bookworm
ARG ATARVA_VERSION
ARG EXCLUDE_NEWER
RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev procps \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir uv==0.9.5 \
    && uv pip install --system --no-cache --exclude-newer "${EXCLUDE_NEWER}" "ATaRVa==${ATARVA_VERSION}" \
    && uv pip freeze --system > /opt/atarva-${ATARVA_VERSION}.lock.txt \
    && atarva --version || atarva -v
