# STRdust from the GitHub release binary (static musl build), for releases that have no bioconda image yet
# (1.4.2 was released on 2026-11-22). The download is checked against the SHA256 below, recorded when the release was
# first fetched. procps: Nextflow needs `ps` in the image; awk, sort and cut come with the base image.
FROM debian:bookworm-slim
ARG STRDUST_VERSION
ARG STRDUST_SHA256
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl procps \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL -o /usr/local/bin/STRdust "https://github.com/wdecoster/STRdust/releases/download/v${STRDUST_VERSION}/STRdust-linux-musl" \
    && echo "${STRDUST_SHA256}  /usr/local/bin/STRdust" | sha256sum -c - \
    && chmod +x /usr/local/bin/STRdust \
    && STRdust --version | grep -q "${STRDUST_VERSION}" \
    && apt-get purge -y curl && apt-get autoremove -y
