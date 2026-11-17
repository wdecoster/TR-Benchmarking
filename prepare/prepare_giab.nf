#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Prepare the GIAB inputs as described in the paper (Methods, "Sequencing data preparation"):
 * EPI2ME GIAB 2025.01 aligned CRAMs, both flow cells of a sample pooled and downsampled to ~30x with
 * `samtools view -s 42.<fraction>`.
 *
 * Differences from the paper's procedure, none of which changes which reads are kept:
 *  - each flow cell is subsampled while streaming it from S3, then the two are merged. samtools -s
 *    decides per read from a hash of the read name and the seed, so subsample-then-merge keeps exactly
 *    the reads that merge-then-subsample keeps;
 *  - the fraction is computed from the EPI2ME mosdepth summaries: target / (cov1 + cov2). The paper
 *    does not give its fractions, so this may differ in the last digits;
 *  - @SQ UR tags (paths on the EPI2ME system) are removed from the header, as the paper did for medaka
 *    (medaka issue #574);
 *  - optional --slice_pad N keeps only reads overlapping the catalog loci +/- N bp (lighter; must be
 *    shown to give the same calls before use). NOTE 2026-11-17: indexed region reads from the S3 CRAMs fail
 *    ("Container header CRC32 failure" for every region; the .crai seems not to match), so slicing needs a
 *    local copy first.
 *
 * Open question: hac or sup basecalls? The paper does not say; --basecall chooses (default sup).
 *
 * Usage:
 *   nextflow run prepare/prepare_giab.nf -profile local -c conf/local/<site>.config \
 *       --samples HG002,HG003,HG004 --basecall sup --outdir <results>/inputs
 */

params.samples    = 'HG002'
params.table      = "${projectDir}/giab_2025.01.tsv"
params.basecall   = 'sup'
params.target_cov = 30
params.seed       = 42
params.decimals   = 3
params.slice_pad  = null
params.catalog    = "${projectDir}/../catalogs/benchmark-catalog-v2.strkit.bed"
params.ref        = null
params.outdir     = null
params.threads    = 4
params.samtools_container = 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

if (!params.ref || !params.outdir) error("--ref and --outdir are required")
if (!(params.basecall in ['hac', 'sup'])) error("--basecall must be hac or sup")

def baseUrl(sample) { "https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/wf-human-variation/${params.basecall}/${sample}" }

workflow {
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)
    wanted = params.samples.tokenize(',')*.trim()

    rows = Channel.fromPath(params.table, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .filter { it.sample in wanted }
        .map { row ->
            // coverage per flow cell from the EPI2ME mosdepth summaries (small text files, read here because the
            // samtools container has no curl); fraction = target / sum, as the seed.fraction argument of samtools -s
            def fcs = row.flowcells.tokenize(',')
            def covs = fcs.collect { fc ->
                def line = new URL("${baseUrl(row.sample)}/${fc}/output/SAMPLE.mosdepth.summary.txt").text.readLines().find { it.startsWith('total\t') }
                if (!line) error("no total coverage for ${row.sample} ${fc}")
                line.split('\t')[3] as BigDecimal
            }
            def frac = params.target_cov / covs.sum()
            if (frac >= 1) error("${row.sample}: coverage ${covs.sum()} below target ${params.target_cov}")
            def s = "${params.seed}" + String.format("%.${params.decimals}f", frac).substring(1)
            tuple(row.sample, row.karyotype, fcs, covs, s)
        }

    regions = params.slice_pad != null ? slice_regions(file(params.catalog, checkIfExists: true)) : file('NO_REGIONS')
    prepared = giab_30x(rows, ref, fai, regions).cram

    prepared
        .map { sample, karyotype, cram, crai -> "${params.outdir}/giab_${params.basecall}/${cram.name}\t${karyotype}" }
        .collectFile(name: "giab_${params.basecall}.samples.list", storeDir: params.outdir, newLine: true, sort: true)
}

process slice_regions {
    container params.samtools_container
    input:
    path catalog
    output:
    path 'slice_regions.bed'
    script:
    """
    awk -v OFS='\\t' -v p=${params.slice_pad} '{s=\$2-p; if (s<0) s=0; print \$1, s, \$3+p}' ${catalog} \\
        | sort -k1,1 -k2,2n \\
        | awk -v OFS='\\t' '\$1!=c || \$2>e {if (c!="") print c,s,e; c=\$1; s=\$2; e=\$3; next} {if (\$3>e) e=\$3} END {print c,s,e}' \\
        > slice_regions.bed
    """
}

process giab_30x {
    tag "${sample}"
    container params.samtools_container
    cpus params.threads
    memory 8.GB
    publishDir "${params.outdir}/giab_${params.basecall}", mode: 'copy'

    input:
    tuple val(sample), val(karyotype), val(flowcells), val(covs), val(subsample)
    path ref
    path fai
    path regions

    output:
    tuple val(sample), val(karyotype), path("${sample}.30x.haplotagged.cram"), path("${sample}.30x.haplotagged.cram.crai"), emit: cram
    path "${sample}.coverage.tsv"

    script:
    def base = baseUrl(sample)
    def region_opt = regions.name == 'NO_REGIONS' ? '' : "-M -L ${regions}"
    def cov_lines = [flowcells, covs].transpose().collect { fc, c -> "${sample}\\t${fc}\\t${c}" }.join('\\n')
    """
    set -euo pipefail
    printf '${cov_lines}\\n${sample}\\ttotal\\t${covs.sum()}\\tsamtools -s ${subsample}\\n' > ${sample}.coverage.tsv

    # the flow cells are streamed in parallel: a single HTTPS stream from S3 reached only ~6 MB/s (2026-11-17)
    pids=()
    for fc in ${flowcells.join(' ')}; do
        samtools view -@ ${task.cpus} -T ${ref} -s ${subsample} -C -o \$fc.cram ${region_opt} ${base}/\$fc/output/SAMPLE.haplotagged.cram &
        pids+=(\$!)
    done
    for p in "\${pids[@]}"; do wait \$p; done
    samtools merge -@ ${task.cpus} --reference ${ref} -O CRAM -o merged.cram ${flowcells.collect { it + '.cram' }.join(' ')}
    samtools view -H merged.cram | sed -E '/^@SQ/s/\\tUR:[^\\t]*//' > header.sam
    samtools reheader -i header.sam merged.cram
    mv merged.cram ${sample}.30x.haplotagged.cram
    samtools index ${sample}.30x.haplotagged.cram
    rm -f ${flowcells.collect { it + '.cram' }.join(' ')}
    """
}
