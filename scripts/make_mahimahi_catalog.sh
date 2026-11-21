#!/usr/bin/env bash
# mahi-mahi catalog from the medaka catalog (same loci and 0-based coordinates as the STRdust/medaka catalogs):
# chrom start end motif repeat_id variant_id, with repeat_id = variant_id = the locus ID (no composite TRs).
# Loci with a single-base motif are left out: mahi-mahi stops on them ("Single-base motifs are not supported", a
# limitation of its plotting library, checked before any genotyping), so a catalog with one cannot be run at all.
#   scripts/make_mahimahi_catalog.sh catalogs/benchmark-catalog-v2.medaka.bed > catalogs/benchmark-catalog-v2.mahimahi.bed
set -euo pipefail
awk -F'\t' -v OFS='\t' 'BEGIN { print "#chrom\tstart\tend\tmotif\trepeat_id\tvariant_id" } length($4) > 1 { print $1, $2, $3, $4, $5, $5 }' "$1"
