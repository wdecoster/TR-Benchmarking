# inquiSTR from the GitHub release binary (static musl build). Used for releases that have no bioconda image yet
# (0.28.2 was released on 2026-11-22). The download is checked against the SHA256 below, recorded when the release was
# first fetched. procps: Nextflow needs `ps` in the image.
FROM debian:bookworm-slim
ARG INQUISTR_VERSION
ARG INQUISTR_SHA256
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl procps \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL -o /usr/local/bin/inquiSTR "https://github.com/wdecoster/inquiSTR/releases/download/v${INQUISTR_VERSION}/inquiSTR-linux-musl" \
    && echo "${INQUISTR_SHA256}  /usr/local/bin/inquiSTR" | sha256sum -c - \
    && chmod +x /usr/local/bin/inquiSTR \
    && inquiSTR --version | grep -q "${INQUISTR_VERSION}" \
    && apt-get purge -y curl && apt-get autoremove -y
