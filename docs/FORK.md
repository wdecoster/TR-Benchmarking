# About this fork

This is a fork of [dashnowlab/TR-Benchmarking](https://github.com/dashnowlab/TR-Benchmarking), the workflow behind
Aliyev et al., *A comprehensive assessment of tandem repeat genotyping methods for Nanopore long-read genomes*,
Genome Biology 27:343 (2026). All credit for the benchmark design, catalog and original code goes to its authors.

Aims:
1. reproduce the published results at the published tool versions;
2. keep the benchmark current: re-run tools when they release a new major or minor version, add new tools,
   and publish the results openly.

Maintained by Wouter De Coster, a co-author of the paper and the developer of STRdust and inquiSTR, two of the
benchmarked tools. To keep that conflict of interest manageable, everything is public: the code, the exact
command for every tool, run logs and results. The same rules apply to every tool. **Tool developers are explicitly
invited to open issues or PRs** to improve how their tool is run (parameters, catalog conversion, output parsing).

Every change made to get a tool running is listed, with reasons, in **[ADAPTATIONS.md](ADAPTATIONS.md)**; developers
are invited to comment on them.

## What changed compared to upstream

| Upstream | This fork |
|---|---|
| Five workflow copies (`run_benchmarking_*.nf`) with paths on the authors' cluster | `main.nf`: one workflow. `--mode default` = commands of `run_benchmarking_hg.nf`; `--mode sensitive` = the "min reads" commands of `run_benchmarking_deveson.nf`. The originals are kept unchanged for reference. |
| Tools from local installs, conda envs and containers | Containers only. Version profiles `conf/versions_paper.config` and `conf/versions_current.config` pin the image of every tool. |
| Slurm (CU Boulder Alpine) only | `-profile local` (one machine, Apptainer); the original settings remain as `-profile alpine` |
| GIAB preparation not in the repo | `prepare/prepare_giab.nf` (stream from ONT open data, 30x, seed 42) |

Results go to `<outdir>/<mode>/<tool>/<version>/`, so versions never overwrite each other.

### Tool builds versus the paper
| Tool | Paper | Container used here | Same build as the paper? |
|---|---|---|---|
| ATaRVa | 0.5.0 (local conda env) | `containers/atarva` (PyPI 0.5.0, dependencies as of the release date) | same version |
| LongTR | 1.2 | the paper's Seqera container | yes |
| Medaka Tandem | 2.1.1 (conda env) | `containers/medaka` (bioconda 2.1.1 with pyabpoa 1.5.5, version check patched)† | same version, patched |
| Straglr | 1.5.5 (git checkout) | the 1.5.5 container declared in the original workflow | same version, presumably |
| STRdust | 0.16.0 (static binary) | biocontainers strdust 0.16.0 | same version |
| STRkit | 0.24.2 | the paper's ghcr container | yes |
| vamos | 3.0.5 and 2.1.7 (compiled) | biocontainers vamos 3.0.5 / 2.1.7 | same version* |

\* The upstream tags vamos-v3.0.4, v3.0.5 and v3.0.6 all point to the same commit (d45c4d0, 2026-02-09), whose
binary reports "Version: 3.0.4". "vamos 3.0.5" is therefore that commit.

† The bioconda recipes of medaka 2.1.1 do not pin pyabpoa, and `medaka tandem` 2.1.1 refuses to start with anything
but pyabpoa 1.5.1, so none of the published 2.1.1 images can run it. With pyabpoa 1.5.1 the run also fails on HG002
(medaka #583); the image used here has pyabpoa 1.5.5 and a patched version check, see
[ADAPTATIONS.md](ADAPTATIONS.md) (A-medaka-1).

Versions reported by every container: `tests/containers.nf` (all 21 images checked 2026-11-17).

### CRAM reference
The input CRAMs carry no `@SQ UR` tags (removed because of medaka issue #574). Tools that do not hand the FASTA to
htslib as the CRAM reference (vamos takes no FASTA; ATaRVa 0.7.1 does not use it for CRAM) find the reference through
`REF_PATH`, pointing to an MD5-keyed cache built from the same FASTA (process `ref_cache`). Without it, vamos exits
successfully with an empty VCF; such tasks now fail on htslib's decode errors.

### Command changes in current versions
- ATaRVa ≥ 0.7: `atarva genotype ...` (same options).
- STRdust 1.x: `--unphased` requires a strategy; `ward` is used (the documented general choice).

## Running

```bash
# 1. inputs: HG002 (+ parents for the Mendelian analysis), 30x, as in the paper
nextflow run prepare/prepare_giab.nf -profile local -c conf/local/<site>.config \
    --samples HG002,HG003,HG004 --basecall sup --outdir <results>/inputs
# 2. genotyping at the paper versions, single-threaded like the paper's timing runs
nextflow run main.nf -profile local,paper -c conf/local/<site>.config --list <results>/inputs/giab_sup.samples.list
# 3. current versions, only tools whose major/minor version changed
nextflow run main.nf -profile local,current -c conf/local/<site>.config --list ... --tools atarva,medaka,strdust,strkit,vamos
```
`conf/site.config.example` shows the site settings (output folder, reference, core budget). Launch from the
results folder. The trace in `pipeline_info/` records wall time, %CPU (CPU time = realtime × %cpu) and peak memory
per task; each task directory also records the machine load at its start and end (`.loadavg_start/_end`).

## Timing benchmark (`timing.nf`)
Runtimes from the genome-wide runs are single measurements, taken while other jobs share the machine. `timing.nf`
measures them more carefully on a small subset:
- input CRAM and catalog restricted to one chromosome (`--chrom chr20`: 992 of the 43,009 loci);
- every tool version of the chosen version sets (`--version_sets paper,current`; a version shared by both runs
  once) × `--reps` replicates, single-threaded, **one run at a time**, in shuffled order. Nextflow does not
  keep the input order when scheduling, so the actual order is recorded (start time per run), not fixed;
- the same commands as `main.nf`: both workflows take them from `lib/TrCommands.groovy`;
- output `timing/timing_chr20.tsv`: wall time of the tool command, load average at start and end, and a digest of
  the calls (to check that replicates are identical); CPU time and peak memory are in the Nextflow trace.

```bash
nextflow run timing.nf -profile local -c conf/local/<site>.config \
    --cram <results>/inputs/giab_sup/HG002.30x.haplotagged.cram --karyotype XY --reps 3
```
Incremental: like the genotyping, only new tool versions are timed (`--version_sets current --tools <changed>`);
earlier timings are reused. Every session writes its own table (`timing_chr20_<date>.tsv`) and always also times an
**anchor** (`--anchor longtr:paper`, a version that never changes): sessions are comparable only when the anchor's
time is stable, and its ratio between sessions shows how much faster or slower a session ran.

Caveats: fixed costs (start-up, model loading) weigh more on one chromosome than genome-wide; LongTR's
`--max-tr-len` comes from the catalog it is given, so on chr20 it is chr20's longest locus.

## Reproduction targets
- Pathogenic loci (paper Fig. 6): `scripts/pathogenic_results.{default,min_reads}.tsv` hold the per-sample,
  per-tool allele lengths and sequences, so this part can be reproduced call by call (data: SRA PRJNA786382).
- Timing (Fig. 1): relative between tools only, since the hardware differs.

## Re-run rule
A tool is re-run when its **major or minor** version changes; patch releases are not re-run, unless the patch changes
results and its developer says so (STRdust 1.4.2, ADAPTATIONS A-strdust-3).

## Tools added after the paper
### mahi-mahi 1.0.0 (Oxford Nanopore, version set `current`)
- Container: `containers/mahimahi`, built from the release wheel (checked against the release's SHA256SUMS) with
  dependencies as of the release date; CI runs ONT's install check (one HG002 locus, byte-identical VCF) before
  pushing. License: ONT Public License v1.0 (research use with nanopore data), included in the image.
- Catalog: `catalogs/benchmark-catalog-v2.mahimahi.bed` (`scripts/make_mahimahi_catalog.sh`, from the medaka
  catalog: same coordinates, repeat_id = variant_id = locus ID; without the 1 bp motifs, see below).
- Settings, all visible in `lib/TrCommands.groovy`:
  - `region_merge_distance=25`: the catalog has loci 26 bp apart; with the default (50) mahi-mahi stops with an
    error rather than genotype them separately;
  - per-locus plots off (meant for a few dozen loci, not 43,009);
  - sensitive mode: `min_region_coverage=1`, `min_read_support_per_allele=1` (closest equivalent of the other
    tools' "min reads 1"; mahi-mahi keeps a fixed minimum of 3 reads per cluster);
  - `--bypass-haplotag-check` only for input without haplotags (`--haplotagged false`, e.g. the pathogenic panel).
- Homopolymers: mahi-mahi refuses single-base motifs ("Single-base motifs are not supported", checked before any
  genotyping, one such locus stops the whole run). The 1,000 loci with a 1 bp motif are therefore left out of its
  catalog (42,009 loci) and count as "not supported" for mahi-mahi, not as failed calls or zeros.
- Scope: mahi-mahi is designed for STRs and states it is not suited to long-motif TRs; 36,940 of the 43,009
  catalog loci have a motif longer than 6 bp. It is run on the full catalog like every tool; results should be
  read per motif-length class (the catalog has ~1,000 loci for each motif length 1-6 bp).
- Published: VCF, log and resolved config; the per-read TSV and tagged BAM stay in the work directory.
- ONT developers are invited to check these choices.

### inquiSTR 0.28.2 (version set `current`; 0.28.1 was run first)
**Conflict of interest:** inquiSTR is developed by the maintainer of this fork. It is run under the same rules as
every other tool: its default settings, the same catalog loci, the same inputs and the same metrics. Its
command is chosen to mirror the closest comparable tool (STRdust, same developer, same catalog file), and
anything that departs from the defaults is listed here.
- Container: `containers/inquistr`, the release's static binary (SHA256 checked); the bioconda image does not exist for
  0.28.2 yet (0.28.1 ran with `quay.io/biocontainers/inquistr:0.28.1--hdcf344f_0`). 0.28.2 fixes the `NaN` genotypes that
  0.28.1 wrote as `1|2` (loci without a call).
- Two variants of the same binary: `inquistr` (phased mode, default) and `inquistr_unphased` (`--unphased` on the
  haplotagged input as well, as the paper ran STRdust).
- Catalog: `catalogs/benchmark-catalog-v2.strdust.bed`. It has the same BED layout as STRdust's: a 0-based start,
  the motif and the locus ID. inquiSTR keeps its default 10 bp padding around each locus.
- Settings, all visible in `lib/TrCommands.groovy`:
  - phased mode (HP tags), inquiSTR's default. With `--haplotagged false` (e.g. the pathogenic panel) it runs
    with `--unphased`, because inquiSTR stops on input without HP tags. The paper ran STRdust with `--unphased`
    even on haplotagged data; that choice is kept for STRdust to reproduce the paper;
  - `--haploid chrX,chrY` for XY samples, as for STRdust;
  - sensitive mode: `--support 1`, as for STRdust;
  - `--reference` for CRAM decoding (no ref cache needed); the decode check still applies.
- Output is length-only. Alleles are symbolic (`<STR>`); lengths are in `FORMAT/FRB` (full allele length in bp)
  and `FORMAT/RB` (relative to the reference). Metrics that need allele sequences (e.g. edit distance) do not
  apply to inquiSTR and are reported as not available, not as zero.

## Planned tools
None at the moment. Suggestions are welcome as issues.

## Known gaps in reproducing the paper
See the issues in this fork. In short: the paper does not state hac vs sup basecalls for GIAB or the exact
downsampling fractions; the HPRC sample metadata index used by `align_HPRC/` is no longer online; scripts that
aggregate per-sample outputs into the figure inputs are not in the repository.
