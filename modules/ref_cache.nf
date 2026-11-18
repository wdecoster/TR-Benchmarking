// Shared by main.nf and timing.nf
params.samtools_container = 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

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
