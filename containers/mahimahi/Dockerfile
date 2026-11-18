# mahi-mahi (Oxford Nanopore) from its GitHub release wheel, checked against the release's SHA256SUMS.
# Dependencies are resolved as they were on the release date (uv --exclude-newer), as for ATaRVa.
# License: Oxford Nanopore Technologies PLC. Public License v1.0 (research use with nanopore data); the license
# text is included in the image at /opt/mahi-mahi/LICENSE.md.
FROM python:3.12-slim-bookworm
ARG MAHIMAHI_VERSION
ARG EXCLUDE_NEWER
RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev curl ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/mahi-mahi
RUN base=https://github.com/nanoporetech/mahi-mahi/releases/download/v${MAHIMAHI_VERSION} \
    && whl=mahi_mahi-${MAHIMAHI_VERSION}-py3-none-any.whl \
    && curl -fsSLO $base/$whl && curl -fsSLO $base/SHA256SUMS \
    && grep " $whl\$" SHA256SUMS | sha256sum -c - \
    && curl -fsSL -o LICENSE.md https://raw.githubusercontent.com/nanoporetech/mahi-mahi/v${MAHIMAHI_VERSION}/LICENSE.md \
    && pip install --no-cache-dir uv==0.9.5 \
    && uv pip install --system --no-cache --exclude-newer "${EXCLUDE_NEWER}" ./$whl \
    && uv pip freeze --system > mahi-mahi-${MAHIMAHI_VERSION}.lock.txt \
    && rm $whl \
    && mahi-mahi --version
WORKDIR /
