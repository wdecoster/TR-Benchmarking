#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Pathogenic-expansion panel (paper Fig. 6): adaptive-sampling ONT data of Stevanovski et al. 2022
 * (Sci. Adv., SRA PRJNA786382), downloaded from ENA and aligned to GRCh38.
 *
 * prepare/pathogenic_runs.tsv lists the 27 samples in scripts/pathogenic_results.*.tsv with their ENA runs,
 * FASTQ URLs and MD5s (generated 2026-11-17 from the ENA portal API). MB382 is the run titled
 * 'NOTCH2NLC_STR_Readfish' (inferred: no run carries the ID, and MB382 is the NOTCH2NLC 15/127 GGC sample).
 *
 * The paper only says "aligned to GRCh38 using minimap2": preset, version and options are not given. Used here:
 * minimap2 2.30 -ax map-ont (the default preset family of the repo's HPRC alignment), all reads (no length
 * filter), one read group per sample, sorted CRAM; @SQ UR tags removed (medaka issue #574).
 *
 *   nextflow run prepare/prepare_pathogenic.nf -profile local -c conf/local/<site>.config --outdir <results>/inputs
 */

params.table   = "${projectDir}/pathogenic_runs.tsv"
params.samples = null          // comma-separated subset; default all
params.ref     = null
params.outdir  = null
params.align_threads = 16
params.minimap2_container = 'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce'   // minimap2 2.30, samtools 1.23.1 (nf-core module image)

if (!params.ref || !params.outdir) error("--ref and --outdir are required")

workflow {
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)
    wanted = params.samples ? params.samples.tokenize(',')*.trim() : null

    rows = Channel.fromPath(params.table, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .filter { wanted == null || it.sample in wanted }

    fastqs = rows.flatMap { r ->
        def urls = r.fastq_urls.tokenize(';')
        def md5s = r.fastq_md5.tokenize(';')
        if (urls.size() != md5s.size()) error("${r.sample}: ${urls.size()} URLs but ${md5s.size()} MD5s")
        [urls, md5s].transpose().collect { u, m -> tuple(r.sample, r.karyotype, u, m) }   // downloaded by process 'download' (Nextflow's own URL staging does not resume truncated transfers)
    }

    // groupKey tells Nextflow how many FASTQs a sample has, so a sample is aligned as soon as its own files are in
    // (a plain groupTuple would wait for every download of every sample)
    def n_fastq = [:]
    file(params.table).splitCsv(header: true, sep: '\t').each { n_fastq[it.sample] = it.fastq_urls.tokenize(';').size() }
    downloaded = download(fastqs)
        .map { sample, karyotype, fq -> tuple(groupKey(sample, n_fastq[sample]), karyotype, fq) }
        .groupTuple(by: [0, 1])

    aligned = align(downloaded, ref, fai)

    aligned
        .map { sample, karyotype, cram, crai -> "${params.outdir}/pathogenic/${cram.name}\t${karyotype}" }
        .collectFile(name: 'pathogenic.samples.list', storeDir: params.outdir, newLine: true, sort: true)
}

process download {
    // curl on the host (no container): resumes a dropped connection (-C -) and retries; the MD5 must match
    tag "${sample} ${url.tokenize('/')[-1]}"
    cpus 1
    maxForks 4        // polite to ENA

    input:
    tuple val(sample), val(karyotype), val(url), val(md5)

    output:
    tuple val(sample), val(karyotype), path("*.fastq.gz")

    script:
    def fq = url.tokenize('/')[-1]
    """
    ok=0
    for i in \$(seq 1 30); do
        if curl -sS -L -C - --connect-timeout 30 --retry 3 --retry-delay 10 -o ${fq} ${url}; then ok=1; break; fi
        echo "download attempt \$i failed, resuming" >&2; sleep 30
    done
    [ \$ok = 1 ] || { echo "download of ${url} failed after 30 attempts" >&2; exit 1; }
    echo "${md5}  ${fq}" | md5sum -c -
    """
}

process align {
    tag "${sample}"
    container params.minimap2_container
    cpus params.align_threads
    memory 48.GB
    maxForks 2        // two alignments at once (2 x 16 threads); the site config must allow it (executor.cpus)
    publishDir "${params.outdir}/pathogenic", mode: 'copy'

    input:
    tuple val(sample), val(karyotype), path(fastqs)
    path ref
    path fai

    output:
    tuple val(sample), val(karyotype), path("${sample}.cram"), path("${sample}.cram.crai")

    script:
    """
    set -euo pipefail
    cat ${fastqs.join(' ')} \\
        | minimap2 -t ${task.cpus} -ax map-ont -R '@RG\\tID:${sample}\\tSM:${sample}' ${ref} - \\
        | samtools sort -@ 4 -m 2G --reference ${ref} -O CRAM -o ${sample}.cram -
    samtools view -H ${sample}.cram | sed -E '/^@SQ/s/\\tUR:[^\\t]*//' > header.sam
    samtools reheader -i header.sam ${sample}.cram
    samtools index ${sample}.cram
    """
}
