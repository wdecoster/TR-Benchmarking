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
| Medaka Tandem | 2.1.1 (conda env) | biocontainers medaka 2.1.1 (py310) | same version |
| Straglr | 1.5.5 (git checkout) | the 1.5.5 container declared in the original workflow | same version, presumably |
| STRdust | 0.16.0 (static binary) | biocontainers strdust 0.16.0 | same version |
| STRkit | 0.24.2 | the paper's ghcr container | yes |
| vamos | 3.0.5 and 2.1.7 (compiled) | biocontainers vamos 3.0.5 / 2.1.7 | same version* |

\* The upstream tags vamos-v3.0.4, v3.0.5 and v3.0.6 all point to the same commit (d45c4d0, 2026-02-09), whose
binary reports "Version: 3.0.4". "vamos 3.0.5" is therefore that commit.

Versions reported by every container: `tests/containers.nf` (all 21 images checked 2026-11-17).

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

## Reproduction targets
- Pathogenic loci (paper Fig. 6): `scripts/pathogenic_results.{default,min_reads}.tsv` hold the per-sample,
  per-tool allele lengths and sequences, so this part can be reproduced call by call (data: SRA PRJNA786382).
- Timing (Fig. 1): relative between tools only, since the hardware differs.

## Re-run rule
A tool is re-run when its **major or minor** version changes; patch releases are not re-run.

## Planned tools
- [mahi-mahi](https://github.com/nanoporetech/mahi-mahi) (ONT, v1.0.0)
- [inquiSTR](https://github.com/wdecoster/inquiSTR)

## Known gaps in reproducing the paper
See the issues in this fork. In short: the paper does not state hac vs sup basecalls for GIAB or the exact
downsampling fractions; the HPRC sample metadata index used by `align_HPRC/` is no longer online; scripts that
aggregate per-sample outputs into the figure inputs are not in the repository.
