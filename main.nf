#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Portable version of the benchmark workflow (fork of dashnowlab/TR-Benchmarking).
 *
 * Tool commands are copied from the original workflows:
 *   mode 'default'   = run_benchmarking_hg.nf          (GIAB runs; paper Fig. 1, 3, 4, 5)
 *   mode 'sensitive' = run_benchmarking_deveson.nf     ("min reads" runs; paper Fig. 6)
 * What changed: no hardcoded cluster paths (catalogs from the repo, output in --outdir),
 * tool versions and containers from a version profile (conf/versions_*.config), and
 * per-version CLI differences in that profile instead of in the code.
 * The original *.nf files are kept unchanged for provenance.
 *
 * Usage:
 *   nextflow run main.nf -profile local,paper -c <site.config> \
 *       --list samples.list --ref GRCh38_no_alt.fa --outdir results [--mode sensitive] [--tools strdust,vamos]
 *   samples.list: <cram or bam><TAB><karyotype XX|XY>, sample name = file name without .cram/.bam
 */

params.list       = null
params.ref        = null
params.outdir     = null
params.mode       = 'default'
params.tools      = 'atarva,longtr,medaka,straglr,strdust,strkit,vamos'
params.threads    = 1          // the paper timed every tool single-threaded
params.catalog_dir = "${projectDir}/catalogs"
params.catalog    = 'benchmark-catalog-v2'
params.mosdepth   = false
params.haplotagged = true     // input reads carry HP/PS tags (GIAB EPI2ME CRAMs); false for the pathogenic panel (mahi-mahi, inquiSTR)

def die(msg) { error(msg) }

if (!params.list || !params.ref || !params.outdir) die("--list, --ref and --outdir are required")
if (!(params.mode in ['default', 'sensitive'])) die("--mode must be 'default' or 'sensitive'")
if (!params.version_set) die("choose a version profile: -profile ...,paper or ...,current")

def selected = params.tools.tokenize(',')*.trim()
def known = ['atarva', 'inquistr', 'inquistr_unphased', 'longtr', 'mahimahi', 'medaka', 'straglr', 'strdust', 'strkit', 'vamos', 'vamos_old']
selected.each { if (!(it in known)) die("unknown tool '${it}' (known: ${known.join(',')})") }
selected.each { if (!params.tool_versions.containsKey(it)) die("tool '${it}' has no entry in version set '${params.version_set}'") }

// one output folder per tool, version and mode, so versions never overwrite each other
def outFor = { String tool -> "${params.outdir}/${params.mode}/${tool}/${params.tool_versions[tool].version}" }
def cat = { String suffix -> "${params.catalog_dir}/${params.catalog}.${suffix}" }
def sensitive = params.mode == 'sensitive'

include { ref_cache } from './modules/ref_cache.nf'

def sampleName = { aln ->
    aln.name.replaceFirst(/\.(cram|bam)$/, '').replaceFirst(/\.merged$/, '')
}

workflow {
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)

    samples = Channel.fromPath(params.list, checkIfExists: true)
        .splitText()
        .map { it.replaceAll(/\r?\n$/, '') }
        .filter { it.trim() && !it.trim().startsWith('#') }
        .map { line ->
            def parts = line.split('\t', -1)
            if (parts.size() < 2) error("Expected <bam|cram><TAB><karyotype>, got: '${line}'")
            def aln = file(parts[0].trim(), checkIfExists: true)
            def karyo = parts[1].trim()
            if (!(karyo in ['XX', 'XY'])) error("karyotype must be XX or XY, got '${karyo}' for ${aln}")
            def idx
            if (aln.name.endsWith('.cram')) idx = file("${aln}.crai", checkIfExists: true)
            else if (aln.name.endsWith('.bam')) {
                idx = file("${aln}.bai").exists() ? file("${aln}.bai") : file(aln.toString().replaceFirst(/\.bam$/, '.bai'), checkIfExists: true)
            } else error("Unsupported file type: ${aln}")
            tuple(sampleName(aln), aln, idx, karyo)
        }

    // Tools that do not pass the FASTA to htslib as the CRAM reference (vamos takes no FASTA; ATaRVa 0.7) find the
    // reference through REF_PATH: an MD5-keyed cache of the same FASTA (the input CRAMs carry no UR tags).
    refcache = (selected.any { it in TrCommands.NEEDS_REF_CACHE }) ? ref_cache(ref, fai) : null

    if (params.mosdepth) mosdepth(samples, ref, fai)
    if ('atarva' in selected)    atarva(samples, ref, fai, file(cat('atarva.bed.gz')), file(cat('atarva.bed.gz.tbi')), refcache)
    def inqTools = ['inquistr', 'inquistr_unphased'].findAll { it in selected }
    if (inqTools) inquistr(samples.combine(Channel.fromList(inqTools)), ref, fai, file(cat('strdust.bed')))   // one task per sample and variant
    if ('longtr' in selected)    longtr(samples, ref, fai, file(cat('longtr.bed')))
    if ('mahimahi' in selected)  mahimahi(samples, ref, fai, file(cat('mahimahi.bed')), refcache)
    if ('medaka' in selected)    medaka(samples, ref, fai, file(cat('medaka.bed')), refcache)
    if ('straglr' in selected)   straglr(samples, ref, fai, file(cat('strglr.bed')))
    if ('strdust' in selected)   strdust(samples, ref, fai, file(cat('strdust.bed')))
    if ('strkit' in selected)    strkit(samples, ref, fai, file(cat('strkit.bed')))
    if ('vamos' in selected) {
        vamos(samples, ref, fai, file(cat('vamos.bed')), 'vamos', refcache)
        vamos_fix(vamos.out.raw, ref, fai, file("${projectDir}/scripts/${params.tool_versions.vamos.fix_script}"), 'vamos')
    }
    if ('vamos_old' in selected) {   // vamos 2.1.7, used next to 3.0.5 for the pathogenic loci
        vamos_old(samples, ref, fai, file(cat('vamos.bed')), 'vamos_old', refcache)
        vamos_old_fix(vamos_old.out.raw, ref, fai, file("${projectDir}/scripts/${params.tool_versions.vamos_old.fix_script}"), 'vamos_old')
    }
}

/* -------------------------------------------------------------------------- */

process mosdepth {
    container params.tool_versions.mosdepth.container
    cpus 4
    memory 8.GB
    publishDir "${params.outdir}/mosdepth", mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.mosdepth.summary.txt"

    script:
    """
    mosdepth -n -x -t ${task.cpus} -f ${ref} ${sample} ${aln}
    """
}

process atarva {
    tag "${sample}"
    container params.tool_versions.atarva.container
    cpus params.threads
    memory { 32.GB * task.attempt }
    publishDir outFor('atarva'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog
    path catalog_tbi
    path refcache

    output:
    path "${sample}.atarva.vcf"

    script:
    TrCommands.script('atarva', params.tool_versions.atarva, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

process longtr {
    tag "${sample}"
    container params.tool_versions.longtr.container
    cpus 1          // LongTR is single-threaded
    memory { 32.GB * task.attempt }
    publishDir outFor('longtr'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog

    output:
    path "${sample}.longTR.vcf.gz"

    script:
    TrCommands.script('longtr', params.tool_versions.longtr, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

process inquistr {
    tag "${sample}"
    container { params.tool_versions[tool].container }
    cpus params.threads
    memory { 16.GB * task.attempt }
    publishDir { outFor(tool) }, mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype), val(tool)
    path ref
    path fai
    path catalog

    output:
    path "${sample}.inquistr.vcf"

    script:
    TrCommands.script(tool, params.tool_versions[tool], [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive, haplotagged: params.haplotagged as boolean])
}

process mahimahi {
    tag "${sample}"
    container params.tool_versions.mahimahi.container
    cpus 1          // mahi-mahi is single-threaded
    memory { 32.GB * task.attempt }
    publishDir outFor('mahimahi'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog
    path refcache

    output:
    // the VCF, run log and resolved config; the per-read TSV and tagged BAM stay in the work directory
    path "${sample}.mahimahi/${sample}.{vcf,log,config.yaml}"

    script:
    TrCommands.script('mahimahi', params.tool_versions.mahimahi, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive, haplotagged: params.haplotagged as boolean])
}

process medaka {
    tag "${sample}"
    container params.tool_versions.medaka.container
    cpus params.threads
    memory { 32.GB * task.attempt }
    publishDir outFor('medaka'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog
    path refcache

    output:
    path "${sample}.medaka.vcf"     // a directory: medaka tandem writes several files, incl. medaka_to_ref.TR.vcf

    script:
    TrCommands.script('medaka', params.tool_versions.medaka, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

process straglr {
    tag "${sample}"
    container params.tool_versions.straglr.container
    cpus params.threads
    memory { 80.GB * task.attempt }   // paper: mean 37 GB, up to 70 GB
    publishDir outFor('straglr'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog

    output:
    path "${sample}.vcf"

    script:
    TrCommands.script('straglr', params.tool_versions.straglr, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

process strdust {
    tag "${sample}"
    container params.tool_versions.strdust.container
    cpus params.threads
    memory { 32.GB * task.attempt }
    publishDir outFor('strdust'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog

    output:
    path "${sample}.strdust.vcf"
    path "${sample}.strdust.sorted.vcf"

    script:
    TrCommands.script('strdust', params.tool_versions.strdust, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

process strkit {
    tag "${sample}"
    container params.tool_versions.strkit.container
    cpus params.threads
    memory { 32.GB * task.attempt }
    publishDir outFor('strkit'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog

    output:
    path "${sample}.strkit*.vcf"

    script:
    TrCommands.script('strkit', params.tool_versions.strkit, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

// vamos 3.x (process 'vamos') and 2.1.7 (process 'vamos_old'): same command apart from -Z
process vamos {
    tag "${sample}"
    container params.tool_versions.vamos.container
    cpus params.threads
    memory { 32.GB * task.attempt }
    publishDir outFor('vamos'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog
    val tool
    path refcache

    output:
    tuple val(sample), path("${sample}.vamos.vcf"), emit: raw

    script:
    TrCommands.script('vamos', params.tool_versions[tool], [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

process vamos_old {
    tag "${sample}"
    container params.tool_versions.vamos_old.container
    cpus params.threads
    memory { 32.GB * task.attempt }
    publishDir outFor('vamos_old'), mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai
    path catalog
    val tool
    path refcache

    output:
    tuple val(sample), path("${sample}.vamos.vcf"), emit: raw

    script:
    TrCommands.script('vamos_old', params.tool_versions[tool], [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: karyotype, cpus: task.cpus, sensitive: sensitive])
}

// fix-vcf*.py needs pysam, which the vamos containers lack; split out so vamos timing stays clean
process vamos_fix {
    tag "${sample}"
    container params.tool_versions.pysam.container
    cpus 1
    memory 4.GB
    publishDir outFor('vamos'), mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path ref
    path fai
    path fix_script
    val tool

    output:
    path "${sample}.vamos.fixed.vcf"

    script:
    """
    python3 ${fix_script} ${vcf} ${sample}.vamos.fixed.vcf --ref ${ref}
    """
}

process vamos_old_fix {
    tag "${sample}"
    container params.tool_versions.pysam.container
    cpus 1
    memory 4.GB
    publishDir outFor('vamos_old'), mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path ref
    path fai
    path fix_script
    val tool

    output:
    path "${sample}.vamos.fixed.vcf"

    script:
    """
    python3 ${fix_script} ${vcf} ${sample}.vamos.fixed.vcf --ref ${ref}
    """
}
