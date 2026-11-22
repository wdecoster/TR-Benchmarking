# Adaptations: what we changed to get each tool to run

Running a tool "as published" often needs small changes: a pinned dependency, a catalog converted to the tool's
format, a flag the tool requires. Each of those is a choice that can influence a result, and some of them are
choices the tool's developers might make differently. This page lists **every such change** in one place, with the
reason and the evidence, so that tool developers and readers can check them and say so.

**How to comment:** open an issue with the template "Comment on an adaptation" and give the ID. A change is removed
or replaced when its developers (or anyone, with evidence) show a better option; the benchmark is then rerun for
that tool and the change is recorded here. Maintainer's conflict of interest (STRdust, inquiSTR): see
[FORK.md](FORK.md); the same list applies to every tool.

Where a change affects results the effect is stated if known, otherwise "not measured".

## Input data (all tools)
| ID | What we did | Why | Evidence | Effect on results |
|---|---|---|---|---|
| A-input-1 | The input CRAMs have no `@SQ UR` tags. | Without it htslib cannot look for an obsolete path; the paper did the same for the same reason ("we had to re-header these samples"). | medaka [#574](https://github.com/nanoporetech/medaka/issues/574); the paper's methods | None expected. |
| A-input-2 | Tools that do not hand the FASTA to htslib (vamos, ATaRVa 0.7.1, mahi-mahi, medaka) read the CRAM via `REF_PATH`, an MD5-keyed cache built from the same FASTA. | Without it, vamos exited 0 with an empty VCF, ATaRVa failed, medaka skipped every locus. | Run log, 2026-11-17 and 2026-11-22 | None expected: same reference sequence. |
| A-input-3 | A task fails when htslib reports a decode error (`Failure to decode slice`, `Unable to fetch reference`). | Failing silently produced empty or partial output that looks like a result. | As above | Prevents invalid results. |

## Per tool
| ID | Tool | What we did | Why | Evidence | Effect on results |
|---|---|---|---|---|---|
| A-medaka-1 | medaka 2.1.1 (paper version) | Container with **pyabpoa 1.5.5**, and the version check in `medaka/tandem/tandem.py` changed from `1.5.1` to `1.5.5` (the only code change). Build: `containers/medaka`. | With pyabpoa 1.5.1 (the version 2.1.1 checks for) the run dies silently in the consensus stage on HG002 (two attempts here). Upstream attributes this to abpoa; 1.5.5 resolves it; medaka 2.2.1 ships with a newer one. | [#583](https://github.com/nanoporetech/medaka/issues/583) (reporters used the same patch) | Unknown. We could not confirm which pyabpoa the paper used. |
| A-medaka-2 | medaka 2.1.1 | Image based on Debian 12. | Newer glibc refuses the PyTorch build that 2.1.1 resolves to. | `containers/medaka` | None expected. |
| A-atarva-1 | ATaRVa 0.5.0 | Own image from PyPI 0.5.0 with dependencies as of the release date (`containers/atarva`). | No image of this version exists; the paper used a local conda environment. | `containers/atarva` | Not measured. |
| A-atarva-2 | ATaRVa ≥ 0.7 | `atarva genotype` instead of `atarva` as the command. | The command changed in 0.7. | tool help | None. |
| A-vamos-1 | vamos "3.0.5" | The upstream tags v3.0.4, v3.0.5 and v3.0.6 point to the same commit (d45c4d0), whose binary reports 3.0.4; "3.0.5" is that commit. | Tag situation upstream. | FORK.md | None: same code. |
| A-strdust-1 | STRdust 1.x | `--unphased` needs a strategy; `ward` is used. | The option requires one; `ward` is the documented general choice. | STRdust help | Not measured. |
| A-strdust-2 | STRdust (post-processing) | VCF sorted with `sort` without `-V`. | The 1.x biocontainer has busybox `sort` without `-V`. | Run log 2026-11-17 | None: same order for chr1–22, X, Y. |
| A-mahimahi-1 | mahi-mahi 1.0.0 | `region_merge_distance=25`. | The catalog has loci 26 bp apart; with the default (50) mahi-mahi stops instead of genotyping them separately. | FORK.md | Loci within 50 bp but ≥ 25 bp apart are genotyped separately rather than merged. |
| A-mahimahi-2 | mahi-mahi 1.0.0 | Per-locus plots off. | Plots are meant for a few dozen loci, not 43,009. | FORK.md | None for calls. |
| A-mahimahi-3 | mahi-mahi 1.0.0 | **1,000 loci with a 1 bp motif removed** from its catalog (42,009 loci); reported as "not supported", not as failed. | mahi-mahi refuses single-base motifs and one such locus stops the whole run. | Run log 2026-11-21 | Homopolymer loci are not scored for mahi-mahi. |
| A-mahimahi-4 | mahi-mahi 1.0.0 | Own catalog file converted from the medaka catalog (`scripts/make_mahimahi_catalog.sh`). | Different catalog format. | FORK.md | None expected. |
| A-mahimahi-5 | mahi-mahi 1.0.0 | In sensitive mode `min_region_coverage=1`, `min_read_support_per_allele=1`. | Closest equivalent of the other tools' "min reads 1". | FORK.md | mahi-mahi still keeps a minimum of 3 reads per cluster. |
| A-inquistr-1 | inquiSTR 0.28.1 | Same catalog file as STRdust; phased mode; `--haploid chrX,chrY` for XY; `--support 1` in sensitive mode; `--unphased` only without HP tags. | Mirrors the closest comparable tool. | FORK.md | Output is length-only; sequence-based metrics are not available. |
| A-inquistr-2 | inquiSTR 0.28.2 | Own image from the GitHub release binary (sha256 checked). | No bioconda image yet (released 2026-11-22). | `containers/inquistr` | None: the release binary. |
| A-inquistr-3 | inquiSTR 0.28.2 | Also run with `--unphased` on the haplotagged input (`inquistr_unphased`, reported next to the default). | Compares phased and unphased calling; the paper ran STRdust with `--unphased` on haplotagged data. | FORK.md | A second result set, not a replacement. |

## Evaluation and timing
| ID | What we did | Why |
|---|---|---|
| A-eval-1 | Timings from `timing.nf`: one chromosome (chr20), single-threaded, one run at a time in shuffled order, with an anchor tool to compare sessions. | A controlled measurement on a shared machine; genome-wide runs are single measurements under varying load. Fixed costs weigh more on a small subset. |
| A-eval-3 | Truth: the T2T HG002 v1.1 assembly (maternal and paternal haplotypes) aligned to GRCh38 with `minimap2 -ax asm5 --secondary=no`; allele lengths per catalog locus as in the paper's `assembly_concordance.py` (contigs with MAPQ >= 5 spanning the locus, all must agree; diploid loci need both haplotypes). | The paper only says "asm5"; the other options are not given. Details in `prepare/prepare_truth.nf`. |
| A-eval-4 | All metrics are computed on allele **lengths** for every tool. Where a tool gives no length field, it is derived: STRdust 0.16 `QUICKREF` records (GT 0\|0, no lengths) count as the reference length; vamos' `.fixed.vcf` has a REF one base shorter than the catalog interval, so the ALT sequence length is used; mahi-mahi writes one value for a homozygous call; a call with one allele missing counts as no call (as in the paper's script). | inquiSTR and straglr report lengths only, so one length-based path is needed for all. | 
| A-eval-2 | Results are reported by motif-length class (1–6 bp and > 6 bp). | The catalog is a stratified sample of isolated loci with about 1,000 loci per motif length; most loci have motifs > 6 bp, which not every tool is designed for. |

## Not changed
Straglr 1.5.5, LongTR 1.2, STRkit 0.24.2 and 0.25.0, ATaRVa 0.5.0/0.7.1 (apart from the above) run with the
commands of the original workflow. Straglr's memory use (peak 117 GB here, up to ~70 GB in the paper) is a
measurement, not an adaptation.
