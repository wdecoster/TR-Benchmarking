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
params.samtools_container = 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

def die(msg) { error(msg) }

if (!params.list || !params.ref || !params.outdir) die("--list, --ref and --outdir are required")
if (!(params.mode in ['default', 'sensitive'])) die("--mode must be 'default' or 'sensitive'")
if (!params.version_set) die("choose a version profile: -profile ...,paper or ...,current")

def selected = params.tools.tokenize(',')*.trim()
def known = ['atarva', 'longtr', 'medaka', 'straglr', 'strdust', 'strkit', 'vamos', 'vamos_old']
selected.each { if (!(it in known)) die("unknown tool '${it}' (known: ${known.join(',')})") }
selected.each { if (!params.tool_versions.containsKey(it)) die("tool '${it}' has no entry in version set '${params.version_set}'") }

// one output folder per tool, version and mode, so versions never overwrite each other
def outFor = { String tool -> "${params.outdir}/${params.mode}/${tool}/${params.tool_versions[tool].version}" }
def cat = { String suffix -> "${params.catalog_dir}/${params.catalog}.${suffix}" }
def sensitive = params.mode == 'sensitive'

// exported in tool scripts that need the cache; fails the task if htslib still could not decode CRAM slices,
// because some tools (vamos) exit 0 with an empty result in that case
def refEnv = 'export REF_PATH="$PWD/ref_cache/%2s/%2s/%s"; export REF_CACHE="$REF_PATH"'
// (the tool's stderr goes to tool.err first, then to the task log, so the check does not race the log writer)
def toolErr = '2> tool.err || { cat tool.err >&2; exit 1; }'
def decodeCheck = 'cat tool.err >&2; if grep -q "Failure to decode slice\\|Unable to fetch reference" tool.err; then echo "CRAM decoding failed (reference not found)" >&2; exit 1; fi'

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
    refcache = (selected.any { it in ['atarva', 'vamos', 'vamos_old'] }) ? ref_cache(ref, fai) : null

    if (params.mosdepth) mosdepth(samples, ref, fai)
    if ('atarva' in selected)    atarva(samples, ref, fai, file(cat('atarva.bed.gz')), file(cat('atarva.bed.gz.tbi')), refcache)
    if ('longtr' in selected)    longtr(samples, ref, fai, file(cat('longtr.bed')))
    if ('medaka' in selected)    medaka(samples, ref, fai, file(cat('medaka.bed')))
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

// htslib reference cache (REF_PATH layout %2s/%2s/%s): MD5 of each upper-cased sequence, as in the @SQ M5 tags
process ref_cache {
    container params.samtools_container
    cpus 1
    memory 4.GB
    storeDir "${params.outdir}/ref_cache"

    input:
    path ref
    path fai

    output:
    path 'ref_cache'

    script:
    """
    set -euo pipefail
    mkdir -p tmp_cache
    cut -f1 ${fai} | while read -r ctg; do
        samtools faidx ${ref} "\$ctg" | tail -n +2 | tr -d '\\n' | tr a-z A-Z > seq
        m=\$(md5sum seq | cut -c1-32)
        mkdir -p tmp_cache/\${m:0:2}/\${m:2:2}
        mv seq tmp_cache/\${m:0:2}/\${m:2:2}/\${m:4}
    done
    mv tmp_cache ref_cache
    """
}

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
    def extra = sensitive ? '--min-reads 1' : ''
    """
    ${refEnv}
    ${params.tool_versions.atarva.cmd} -t ${task.cpus} -f ${ref} -b ${aln} -r ${catalog} --format cram ${extra} --karyotype ${karyotype} -o ${sample}.atarva.vcf ${toolErr}
    ${decodeCheck}
    """
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
    def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
    def haploid_args = (karyotype == 'XY') ? '--haploid-chrs chrX,chrY' : ''
    def extra = sensitive ? "--min-reads 1 --min-mean-qual 1 --bam-samps ${sample} --bam-libs ${sample}" : ''
    """
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${catalog} | sort -n | tail -n 1)"
    LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --max-tr-len \$MAX_TR_LEN \\
        --regions ${catalog} \\
        ${haploid_args} \\
        ${extra} \\
        --bams ${aln} \\
        --tr-vcf ${sample}.longTR.vcf.gz
    """
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

    output:
    path "${sample}.medaka.vcf"     // a directory: medaka tandem writes several files, incl. medaka_to_ref.TR.vcf

    script:
    def sex = (karyotype == 'XX') ? 'female' : 'male'
    def extra = sensitive ? '--min_depth 1' : ''
    """
    medaka tandem ${aln} ${ref} ${catalog} ${sex} ${sample}.medaka.vcf --workers ${task.cpus} --sample_name ${sample} ${extra} --debug
    """
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
    def sex = (karyotype == 'XX') ? 'f' : 'm'
    def extra = sensitive ? '--min_support 1 --min_cluster_size 1 --max_num_clusters 2' : ''
    """
    ${params.tool_versions.straglr.cmd} ${aln} ${ref} ${sample} --loci ${catalog} --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y --sex ${sex} ${extra} --genotype_in_size --nprocs ${task.cpus}
    """
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
    def haploid = (karyotype == 'XY') ? '--haploid chrX,chrY' : ''
    def extra = sensitive ? '--support 1' : ''
    """
    ${params.tool_versions.strdust.cmd} \\
      -R ${catalog} ${params.tool_versions.strdust.unphased} ${extra} -t ${task.cpus} --sample ${sample} \\
      ${haploid} ${ref} ${aln} > ${sample}.strdust.vcf

    # same order as the original "sort -k1,1V -k2,2n" for chr1-22/X/Y/M (all the catalog has), written without
    # -V because the STRdust 1.x container only has busybox sort
    { grep '^#' ${sample}.strdust.vcf; \\
      grep -v '^#' ${sample}.strdust.vcf \\
        | awk -F '\\t' -v OFS='\\t' '{ c = \$1; sub(/^chr/, "", c); k = (c ~ /^[0-9]+\$/) ? sprintf("0%09d", c) : "1" c; print k, \$0 }' \\
        | LC_ALL=C sort -t "\$(printf '\\t')" -k1,1 -k3,3n | cut -f2-; \\
    } > ${sample}.strdust.sorted.vcf
    """
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
    if (sensitive)
        """
        strkit call ${aln} --min-reads 1 --ploidy ${karyotype} --sex-chr ${karyotype} --realign --ref ${ref} --loc ${catalog} --vcf ${sample}.strkit-min_read_1.vcf --processes ${task.cpus} --min-read-align-score 0.1
        """
    else
        """
        strkit call ${aln} --ploidy ${karyotype} --realign --sex-chr ${karyotype} --ref ${ref} --loc ${catalog} --vcf ${sample}.strkit.vcf --processes ${task.cpus}
        """
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
    """
    ${refEnv}
    vamos --read -b ${aln} -r ${catalog} -s ${sample} -o ${sample}.vamos.vcf -S ${params.tool_versions[tool].extra} -t ${task.cpus} ${toolErr}
    ${decodeCheck}
    """
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
    """
    ${refEnv}
    vamos --read -b ${aln} -r ${catalog} -s ${sample} -o ${sample}.vamos.vcf -S ${params.tool_versions[tool].extra} -t ${task.cpus} ${toolErr}
    ${decodeCheck}
    """
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
