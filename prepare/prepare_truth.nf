#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Truth for the accuracy metrics: the haplotypes of the T2T HG002 assembly (v1.1, Q100) aligned to GRCh38, and the
 * allele length of each haplotype at every catalog locus.
 *
 *   maternal: hg002v1.1.mat.fasta.gz (chr1-22, chrX)     paternal: hg002v1.1.pat.fasta.gz (chr1-22, chrY)
 *   both from https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/ (public, ~1.8 GB in total)
 *
 * Alignment: minimap2 -ax asm5 --secondary=no. The paper only says assemblies were aligned with `-x asm5`; the
 * other options are not given (see docs/ADAPTATIONS.md, A-truth-1). Truth lengths: evaluate/truth_from_assembly.py,
 * the method of the paper's assembly_concordance.py.
 *
 *   nextflow run prepare/prepare_truth.nf -profile local -c conf/local/<site>.config --outdir <results>/inputs
 * Output: <outdir>/truth/<asm>.{mat,pat}.bam(.bai), their sha256 (of the downloaded FASTA), and HG002.truth.tsv.
 */

params.ref      = null
params.outdir   = null
params.asm_url  = 'https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies'
params.asm      = 'hg002v1.1'
params.catalog  = "${projectDir}/../catalogs/benchmark-catalog-v2.strdust.bed"
params.karyotype = 'XY'
params.align_threads = 16
params.minimap2_container = 'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce'   // minimap2 2.30, samtools 1.23.1
params.pysam_container = 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'

if (!params.ref || !params.outdir) error("--ref and --outdir are required")

workflow {
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)
    asms = Channel.of('mat', 'pat').map { hap -> tuple(hap, file("${params.asm_url}/${params.asm}.${hap}.fasta.gz")) }   // Nextflow downloads the URL

    verified = verify(asms)
    aligned = align(verified.fa, ref, fai)
    bams = aligned.map { hap, bam, bai -> tuple(hap, bam, bai) }.toSortedList { a, b -> a[0] <=> b[0] }
    // sorted by haplotype: mat first, pat second
    truth(bams.map { l -> l[0][1] }, bams.map { l -> l[0][2] }, bams.map { l -> l[1][1] }, bams.map { l -> l[1][2] },
          file(params.catalog, checkIfExists: true), file("${projectDir}/../evaluate"), file("${projectDir}/../scripts"))
}

process verify {
    tag "${hap}"
    container params.minimap2_container
    cpus 1
    memory 1.GB
    publishDir "${params.outdir}/truth", pattern: '*.sha256', mode: 'copy'

    input:
    tuple val(hap), path(fasta)

    output:
    tuple val(hap), path(fasta), emit: fa
    path "${params.asm}.${hap}.fasta.gz.sha256"

    script:
    """
    gzip -t ${fasta}
    sha256sum ${fasta} > ${params.asm}.${hap}.fasta.gz.sha256
    """

    stub:
    """
    echo stub > ${params.asm}.${hap}.fasta.gz.sha256
    """
}

process align {
    tag "${hap}"
    container params.minimap2_container
    cpus params.align_threads
    memory 32.GB
    publishDir "${params.outdir}/truth", mode: 'copy'

    input:
    tuple val(hap), path(fasta)
    path ref
    path fai

    output:
    tuple val(hap), path("${params.asm}.${hap}.bam"), path("${params.asm}.${hap}.bam.bai")

    script:
    """
    minimap2 -t ${task.cpus} -ax asm5 --secondary=no ${ref} ${fasta} \\
        | samtools sort -@ 4 -m 2G -O BAM -o ${params.asm}.${hap}.bam -
    samtools index ${params.asm}.${hap}.bam
    """

    stub:
    """
    touch ${params.asm}.${hap}.bam ${params.asm}.${hap}.bam.bai
    """
}

process truth {
    container params.pysam_container
    cpus 1
    memory 8.GB
    publishDir "${params.outdir}/truth", mode: 'copy'

    input:
    path mat
    path mat_bai
    path pat
    path pat_bai
    path catalog
    path evaluate
    path scripts

    output:
    path 'HG002.truth.tsv'

    script:
    """
    python3 evaluate/truth_from_assembly.py --mat ${mat} --pat ${pat} --catalog ${catalog} --karyotype ${params.karyotype} --out HG002.truth.tsv
    """

    stub:
    """
    echo stub > HG002.truth.tsv
    """
}
