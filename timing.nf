#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Runtime benchmark on one chromosome: every tool version of the chosen version sets, N replicates each, run one
 * at a time in a random order. One task at a time keeps the tools from competing with each other for CPU, memory
 * bandwidth and I/O; the random order spreads slow drifts in machine load (other jobs, thermal state, page cache)
 * over all tools instead of favouring whichever runs first. The runs are shuffled with a fixed seed, but Nextflow
 * does not keep that order when it schedules them (it picks among the waiting tasks), so the actual order is not
 * reproducible: it is recorded instead (start time per run in the output table).
 *
 * The commands are those of main.nf (lib/TrCommands.groovy), single-threaded like the paper's timing runs, on the
 * input CRAM and catalog restricted to --chrom. Each task first reads its inputs once (page cache) so the timing is
 * not dominated by cold reads. Recorded per run: wall time of the tool command itself (date, inside the
 * container), and in the Nextflow trace realtime, %cpu and peak RSS; the load average at start and end; a digest of
 * the calls (VCF records without the header) to check that replicates give identical output.
 *
 * Caveat: LongTR's --max-tr-len is taken from the catalog it is given, as in main.nf, so on one chromosome it is
 * that chromosome's longest locus.
 *
 *   nextflow run timing.nf -profile local -c conf/local/<site>.config \
 *       --cram <results>/inputs/giab_sup/HG002.30x.haplotagged.cram --karyotype XY --reps 3
 */

params.cram         = null
params.karyotype    = null
params.ref          = null
params.outdir       = null
params.chrom        = 'chr20'
params.reps         = 3
params.seed         = 42
params.version_sets = 'paper,current'
params.tools        = 'atarva,longtr,medaka,straglr,strdust,strkit,vamos'
params.threads      = 1
params.catalog_dir  = "${projectDir}/catalogs"
params.catalog      = 'benchmark-catalog-v2'

include { ref_cache } from './modules/ref_cache.nf'

if (!params.cram || !params.karyotype || !params.ref || !params.outdir) error("--cram, --karyotype, --ref and --outdir are required")
if (!(params.karyotype in ['XX', 'XY'])) error("--karyotype must be XX or XY")
if ((params.reps as int) < 1) error("--reps must be >= 1")

// tool versions of each set, read from the same files as the -profile paper/current of main.nf
def loadSet(String name) {
    def f = file("${projectDir}/conf/versions_${name}.config", checkIfExists: true)
    def cfg = new ConfigSlurper().parse(f.text)
    if (!cfg.params.tool_versions) error("no tool_versions in ${f}")
    cfg.params.tool_versions
}

// (names of script-level variables must differ from process inputs: a script-level 'sets' would shadow val(sets))
def toolNames = params.tools.tokenize(',')*.trim()
def setNames = params.version_sets.tokenize(',')*.trim()

// one configuration per distinct (tool, version, container); a version shared by several sets runs once
def configs = [:]
setNames.each { s ->
    def tv = loadSet(s)
    toolNames.each { t ->
        def v = tv[t]
        if (!v) error("tool '${t}' has no entry in version set '${s}'")
        def key = "${t}\t${v.version}\t${v.container}".toString()   // String, not GString: GString keys do not match in a Map
        if (configs.containsKey(key)) configs[key].sets << s
        else configs[key] = [tool: t, v: v as Map, sets: [s]]
    }
}

def runs = []
configs.values().each { c -> (1..(params.reps as int)).each { r -> runs << c + [rep: r] } }
Collections.shuffle(runs, new Random(params.seed as long))
runs = runs.withIndex().collect { c, i -> tuple(i + 1, c.tool, c.v.version, c.sets.join('+'), c.rep, c.v.container, c.v) }

workflow {
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)
    cram = file(params.cram, checkIfExists: true)
    crai = file("${params.cram}.crai", checkIfExists: true)

    log.info "timing: ${configs.size()} tool versions x ${params.reps} replicates = ${runs.size()} runs on ${params.chrom}, order (seed ${params.seed}):"
    runs.each { log.info "  ${it[0]}\t${it[1]} ${it[2]} (${it[3]}) rep ${it[4]}" }

    sub = subset_cram(cram, crai, ref, fai)
    catalogs = subset_catalogs(toolNames.collect { file("${params.catalog_dir}/${params.catalog}.${TrCommands.catalogSuffix(it)}", checkIfExists: true) })
    refcache = ref_cache(ref, fai)

    timed(Channel.fromList(runs), sub, catalogs, ref, fai, refcache)
        .collectFile(name: "timing_${params.chrom}.tsv", storeDir: "${params.outdir}/timing", sort: { it.text.tokenize('\t')[0] as int },
                     seed: "order\ttool\tversion\tsets\trep\tstart_epoch\twall_s\tload_start\tload_end\tcalls_md5\n")
}

process subset_cram {
    container params.samtools_container
    cpus 4
    memory 4.GB
    storeDir "${params.outdir}/timing/inputs"

    input:
    path cram
    path crai
    path ref
    path fai

    output:
    tuple path("${cram.baseName}.${params.chrom}.cram"), path("${cram.baseName}.${params.chrom}.cram.crai")

    script:
    def out = "${cram.baseName}.${params.chrom}.cram"
    """
    set -euo pipefail
    samtools view -@ ${task.cpus} -T ${ref} -C -o tmp.cram ${cram} ${params.chrom}
    samtools view -H tmp.cram | sed -E '/^@SQ/s/\\tUR:[^\\t]*//' > header.sam
    samtools reheader -i header.sam tmp.cram
    mv tmp.cram ${out}
    samtools index ${out}
    """
}

process subset_catalogs {
    container params.samtools_container
    cpus 1
    memory 2.GB
    storeDir "${params.outdir}/timing/inputs"

    input:
    path cats, stageAs: 'full/*'

    output:
    path "catalogs_${params.chrom}"

    script:
    """
    set -euo pipefail
    mkdir -p tmp
    for f in full/*; do
        n=\$(basename \$f)
        case \$n in
            *.bed.gz) zcat \$f | awk -v c=${params.chrom} '\$1 == c' | bgzip > tmp/\$n; tabix -p bed tmp/\$n ;;
            *)        awk -v c=${params.chrom} '\$1 == c' \$f > tmp/\$n ;;
        esac
        [ -s tmp/\$n ] || { echo "no ${params.chrom} loci in \$n" >&2; exit 1; }
    done
    mv tmp catalogs_${params.chrom}
    """
}

process timed {
    tag "${order} ${tool} ${version} rep ${rep}"
    container "${image}"
    cpus params.threads
    memory { tool == 'straglr' ? 80.GB : 32.GB }
    maxForks 1                  // one run at a time
    errorStrategy 'terminate'   // no retries: a retry would break the order and double-count load

    input:
    tuple val(order), val(tool), val(version), val(in_sets), val(rep), val(image), val(v)
    tuple path(aln), path(idx)
    path catalogs
    path ref
    path fai
    path refcache

    output:
    path 'timing.tsv'

    script:
    def sample = aln.baseName
    def catalog = "${catalogs}/${params.catalog}.${TrCommands.catalogSuffix(tool)}"
    def out = TrCommands.output(tool, sample, false)
    def cmd = TrCommands.script(tool, v, [sample: sample, aln: aln, ref: ref, catalog: catalog, karyotype: params.karyotype,
                                          cpus: task.cpus, sensitive: false])
    """
    set -euo pipefail
    cat ${aln} ${idx} ${catalog} > /dev/null    # page cache: equal footing for every run
    load_start=\$(cut -d' ' -f1-3 /proc/loadavg)
    t0=\$(date +%s.%N)
    ${cmd}
    t1=\$(date +%s.%N)
    load_end=\$(cut -d' ' -f1-3 /proc/loadavg)
    calls() { if [ -d "\$1" ]; then cat "\$1"/*.vcf; elif [ "\${1%.gz}" != "\$1" ]; then zcat "\$1"; else cat "\$1"; fi; }
    md5=\$(calls ${out} | grep -v '^#' | md5sum | cut -c1-32)
    wall=\$(awk -v a=\$t0 -v b=\$t1 'BEGIN { printf "%.1f", b - a }')
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' ${order} ${tool} ${version} ${in_sets} ${rep} \$t0 \$wall "\$load_start" "\$load_end" \$md5 > timing.tsv
    """
}
