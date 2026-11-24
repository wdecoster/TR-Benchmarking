# ATaRVa from PyPI. Dependencies are resolved as they were on the release date (uv --exclude-newer),
# so the image reproduces the environment a user installing that release would have had.
FROM python:3.11-slim-bookworm
ARG ATARVA_VERSION
ARG EXCLUDE_NEWER
# pyabpoa (a dependency) has no wheel for this Python and is compiled here with -march=native, which on the CI runner
# (AVX-512) gives an image that dies with "Illegal instruction" on any CPU without AVX-512 (found on a laptop,
# 2026-11-24). AVX2=1 makes its setup.py build for AVX2 instead; results do not depend on the SIMD width.
ENV AVX2=1
RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev procps \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir uv==0.9.5 \
    && uv pip install --system --no-cache --exclude-newer "${EXCLUDE_NEWER}" "ATaRVa==${ATARVA_VERSION}" \
    && uv pip freeze --system > /opt/atarva-${ATARVA_VERSION}.lock.txt \
    && atarva --version || atarva -v
