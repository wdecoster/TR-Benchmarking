#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Accuracy metrics of the tool outputs of one sample against the assembly truth (length level, and sequence level when
 * the truth has sequences): call rate and concordance per motif-length class, all-loci and no-homopolymer analyses.
 * One run evaluates every tool and version found under <results>/<mode>/ ; the JSON is what the dashboard export reads.
 *
 *   nextflow run evaluate/evaluate.nf -profile local -c conf/local/<site>.config \
 *       --results <outdir> --mode default --truth <outdir>/inputs/truth/HG002.truth.tsv \
 *       --truth_seq <outdir>/inputs/truth/HG002.truth.seq.tsv --evaldir <outdir>/evaluation
 * Output: <evaldir>/evaluation.<mode>.<label>.{tsv,json}   (label = --label, default the sample)
 * The step is a pure function of the tool outputs and the truth: rerun it for every new tool version (-resume is safe here
 * because the inputs are the files themselves, but never across a change of --mode of a tool run).
 */

params.results   = null
params.mode      = 'default'
params.sample    = 'HG002.30x.haplotagged'
params.karyotype = 'XY'
params.catalog   = "${projectDir}/../catalogs/benchmark-catalog-v2.strdust.bed"
params.truth     = null
params.truth_seq = null
params.evaldir   = null
params.label     = null
params.eval_container = 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'    // only python3 is needed

def label = params.label ?: params.sample

if (!params.results || !params.truth || !params.evaldir) error("--results, --truth and --evaldir are required")

workflow {
    evaluate(file("${params.results}/${params.mode}", checkIfExists: true),
             file(params.catalog, checkIfExists: true),
             file(params.truth, checkIfExists: true),
             params.truth_seq ? file(params.truth_seq, checkIfExists: true) : file('NO_FILE'),
             file(projectDir))
}

process evaluate {
    container params.eval_container
    cpus 1
    memory 8.GB
    publishDir params.evaldir, mode: 'copy'

    input:
    path results, stageAs: 'results'
    path catalog
    path truth
    path truth_seq
    path evalcode, stageAs: 'evalcode'

    output:
    path "evaluation.${params.mode}.${label}.tsv"
    path "evaluation.${params.mode}.${label}.json"

    script:
    def seq = params.truth_seq ? "--truth-seq ${truth_seq}" : ''
    """
    python3 evalcode/evaluate.py --catalog ${catalog} --results results --sample ${params.sample} --karyotype ${params.karyotype} \\
        --truth ${truth} ${seq} --json evaluation.${params.mode}.${label}.json > evaluation.${params.mode}.${label}.tsv
    """
}
