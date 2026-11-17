#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Smoke test of the containers in a version profile: pulls every image (into the regular cache, so later runs
 * reuse it) and records what each tool reports as its version.
 *   nextflow run tests/containers.nf -profile local,paper -c conf/local/<site>.config
 */

params.samtools_container = 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

def checks = [
    atarva  : { v -> "${v.cmd.tokenize(' ')[0]} --version" },
    longtr  : { v -> 'LongTR --help 2>&1 | head -5' },
    medaka  : { v -> 'medaka --version' },
    straglr : { v -> "${v.cmd} --help 2>&1 | head -3; python3 -c 'import importlib.metadata as m; print(\"straglr\", m.version(\"straglr\"))'" },
    strdust : { v -> "${v.cmd} --version" },
    strkit  : { v -> 'strkit --version' },
    vamos   : { v -> 'vamos --version 2>&1 || vamos -h 2>&1 | head -5' },
    vamos_old: { v -> 'vamos --version 2>&1 || vamos -h 2>&1 | head -5' },
    mosdepth: { v -> 'mosdepth --version' },
    pysam   : { v -> 'python3 -c "import pysam; print(\'pysam\', pysam.__version__)"' },
]

workflow {
    items = params.tool_versions.collect { tool, v -> tuple(tool, v.version, v.container, checks[tool](v)) }
    items << tuple('samtools', '1.22.1', params.samtools_container, 'samtools --version | head -1')
    check(Channel.fromList(items))
        .collectFile(name: "containers_${params.version_set}.txt", storeDir: "${params.outdir}/pipeline_info", sort: true)
}

process check {
    tag "${tool} ${version}"
    container "${image}"
    cpus 1
    memory 2.GB

    input:
    tuple val(tool), val(version), val(image), val(cmd)

    output:
    path 'report.txt'

    script:
    """
    { echo "== ${tool} (expected ${version}) ${image}"; ( ${cmd} ) 2>&1 | head -8; echo; } > report.txt
    """
}
