/*
 * The command of every tool, in one place. Used by main.nf (one process per tool) and timing.nf (all tools in a
 * single process), so the benchmark and the timing runs execute exactly the same commands.
 *
 *   script(tool, v, a)
 *     tool: atarva, inquistr, inquistr_unphased, longtr, mahimahi, medaka, straglr, strdust, strkit, vamos, vamos_old
 *     v:    the tool's entry in params.tool_versions (version, container, cmd, unphased, extra, ...)
 *     a:    sample, aln, ref, catalog, karyotype (XX/XY), cpus, sensitive (boolean),
 *           haplotagged (boolean, default true; mahi-mahi and inquiSTR refuse reads without HP tags unless told otherwise)
 *
 * Commands of ATaRVa, mahi-mahi and vamos read the CRAM through an htslib reference cache in ./ref_cache (REF_PATH) and fail
 * when htslib still could not decode CRAM slices: vamos exits 0 with an empty VCF in that case.
 */
class TrCommands {

    static final String REF_ENV = 'export REF_PATH="$PWD/ref_cache/%2s/%2s/%s"; export REF_CACHE="$REF_PATH"'
    // the tool's stderr goes to tool.err first, then to the task log, so the check does not race the log writer
    static final String TOOL_ERR = '2> tool.err || { cat tool.err >&2; exit 1; }'
    static final String DECODE_CHECK = 'cat tool.err >&2; if grep -q "Failure to decode slice\\|Unable to fetch reference" tool.err; then echo "CRAM decoding failed (reference not found)" >&2; exit 1; fi'

    // which tools need ./ref_cache staged in the task directory
    static final List<String> NEEDS_REF_CACHE = ['atarva', 'mahimahi', 'medaka', 'vamos', 'vamos_old']

    // main output file (or directory) per tool, relative to the task directory
    static String output(String tool, String sample, boolean sensitive) {
        switch (tool) {
            case 'atarva':  return "${sample}.atarva.vcf"
            case 'inquistr':
            case 'inquistr_unphased': return "${sample}.inquistr.vcf"
            case 'longtr':  return "${sample}.longTR.vcf.gz"
            case 'mahimahi': return "${sample}.mahimahi/${sample}.vcf"
            case 'medaka':  return "${sample}.medaka.vcf"
            case 'straglr': return "${sample}.vcf"
            case 'strdust': return "${sample}.strdust.sorted.vcf"
            case 'strkit':  return sensitive ? "${sample}.strkit-min_read_1.vcf" : "${sample}.strkit.vcf"
            case 'vamos':
            case 'vamos_old': return "${sample}.vamos.vcf"
        }
        throw new IllegalArgumentException("unknown tool ${tool}")
    }

    // catalog file suffix per tool (catalogs/<catalog>.<suffix>)
    static String catalogSuffix(String tool) {
        switch (tool) {
            case 'atarva':  return 'atarva.bed.gz'
            case 'inquistr':
            case 'inquistr_unphased': return 'strdust.bed'    // same BED layout (0-based start, motif, ID) as STRdust
            case 'longtr':  return 'longtr.bed'
            case 'mahimahi': return 'mahimahi.bed'
            case 'medaka':  return 'medaka.bed'
            case 'straglr': return 'strglr.bed'
            case 'strdust': return 'strdust.bed'
            case 'strkit':  return 'strkit.bed'
            case 'vamos':
            case 'vamos_old': return 'vamos.bed'
        }
        throw new IllegalArgumentException("unknown tool ${tool}")
    }

    static String script(String tool, Map v, Map a) {
        def sample = a.sample, aln = a.aln, ref = a.ref, catalog = a.catalog, karyotype = a.karyotype
        def cpus = a.cpus
        boolean sensitive = a.sensitive as boolean
        switch (tool) {

        case 'atarva':
            def extra = sensitive ? '--min-reads 1' : ''
            return """
            ${REF_ENV}
            ${v.cmd} -t ${cpus} -f ${ref} -b ${aln} -r ${catalog} --format cram ${extra} --karyotype ${karyotype} -o ${sample}.atarva.vcf ${TOOL_ERR}
            ${DECODE_CHECK}
            """.stripIndent()

        case 'inquistr':
        case 'inquistr_unphased':
            // - phased mode (HP tags) is inquiSTR's default; it stops on input without HP tags, hence -u then
            // - output is length-only: symbolic <STR> alleles, allele lengths in FORMAT/FRB (full length, bp) and
            //   FORMAT/RB (relative to the reference); no allele sequences
            // - reads the CRAM with --reference (no ./ref_cache needed); the decode check still applies
            def haploid = (karyotype == 'XY') ? '--haploid chrX,chrY' : ''
            def extra = sensitive ? '--support 1' : ''
            // v.minlen: explicit --minlen (0.30.0 made the bound inclusive, so the default 1 now counts 1 bp indels; 2 = old behaviour)
            if (v.minlen) extra += " --minlen ${v.minlen}"
            // inquistr_unphased: the same version with --unphased on haplotagged input too (as the paper ran STRdust)
            def unphased = (tool == 'inquistr_unphased' || (a.haplotagged != null && !a.haplotagged)) ? '--unphased' : ''
            return """
            inquiSTR call ${aln} -R ${catalog} --reference ${ref} -t ${cpus} --sample-name ${sample} ${haploid} ${extra} ${unphased} --vcf \\
                > ${sample}.inquistr.vcf ${TOOL_ERR}
            ${DECODE_CHECK}
            """.stripIndent()

        case 'longtr':
            def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
            def haploid_args = (karyotype == 'XY') ? '--haploid-chrs chrX,chrY' : ''
            def extra = sensitive ? "--min-reads 1 --min-mean-qual 1 --bam-samps ${sample} --bam-libs ${sample}" : ''
            return """
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
            """.stripIndent()

        case 'mahimahi':
            // - region_merge_distance 25: the catalog has loci 26 bp apart; with the default (50) mahi-mahi stops
            //   ("loci within region_merge_distance but different repeat IDs") instead of genotyping them separately
            // - no plots: one PDF per locus is meant for a few dozen loci, not tens of thousands
            // - sensitive mode: the closest equivalent of the other tools' "min reads 1" (mahi-mahi keeps a fixed
            //   minimum of 3 reads per cluster)
            def extra = sensitive ? '--set min_region_coverage=1 --set min_read_support_per_allele=1' : ''
            def haplotag = (a.haplotagged == null || a.haplotagged) ? '' : '--bypass-haplotag-check'
            return """
            ${REF_ENV}
            mahi-mahi --cram ${aln} --ref ${ref} --loci ${catalog} --sample ${sample} --sex ${karyotype} -o ${sample}.mahimahi \\
                --set region_merge_distance=25 --set generate_clustering_plots=false --set generate_motif_plots=false \\
                --set generate_methylation_plots=false ${extra} ${haplotag} ${TOOL_ERR}
            ${DECODE_CHECK}
            """.stripIndent()

        case 'medaka':
            def sex = (karyotype == 'XX') ? 'female' : 'male'
            def extra = sensitive ? '--min_depth 1' : ''
            // medaka 2.1.1 (own image) cannot find the reference of a UR-less CRAM by itself: it then skips every locus
            // ("Retrieved too few reads (0 < 3)") and still exits 0. Hence the reference cache and the decode check.
            return """
            ${REF_ENV}
            medaka tandem ${aln} ${ref} ${catalog} ${sex} ${sample}.medaka.vcf --workers ${cpus} --sample_name ${sample} ${extra} --debug ${TOOL_ERR}
            ${DECODE_CHECK}
            """.stripIndent()

        case 'straglr':
            def sex = (karyotype == 'XX') ? 'f' : 'm'
            def extra = sensitive ? '--min_support 1 --min_cluster_size 1 --max_num_clusters 2' : ''
            return """
            ${v.cmd} ${aln} ${ref} ${sample} --loci ${catalog} --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y --sex ${sex} ${extra} --genotype_in_size --nprocs ${cpus}
            """.stripIndent()

        case 'strdust':
            def haploid = (karyotype == 'XY') ? '--haploid chrX,chrY' : ''
            def extra = sensitive ? '--support 1' : ''
            return """
            ${v.cmd} \\
              -R ${catalog} ${v.unphased} ${extra} -t ${cpus} --sample ${sample} \\
              ${haploid} ${ref} ${aln} > ${sample}.strdust.vcf

            # same order as the original "sort -k1,1V -k2,2n" for chr1-22/X/Y/M (all the catalog has), written without
            # -V because the STRdust containers only have busybox sort
            { grep '^#' ${sample}.strdust.vcf; \\
              grep -v '^#' ${sample}.strdust.vcf \\
                | awk -F '\\t' -v OFS='\\t' '{ c = \$1; sub(/^chr/, "", c); k = (c ~ /^[0-9]+\$/) ? sprintf("0%09d", c) : "1" c; print k, \$0 }' \\
                | LC_ALL=C sort -t "\$(printf '\\t')" -k1,1 -k3,3n | cut -f2-; \\
            } > ${sample}.strdust.sorted.vcf
            """.stripIndent()

        case 'strkit':
            if (sensitive)
                return """
                strkit call ${aln} --min-reads 1 --ploidy ${karyotype} --sex-chr ${karyotype} --realign --ref ${ref} --loc ${catalog} --vcf ${sample}.strkit-min_read_1.vcf --processes ${cpus} --min-read-align-score 0.1
                """.stripIndent()
            return """
            strkit call ${aln} --ploidy ${karyotype} --realign --sex-chr ${karyotype} --ref ${ref} --loc ${catalog} --vcf ${sample}.strkit.vcf --processes ${cpus}
            """.stripIndent()

        case 'vamos':
        case 'vamos_old':
            return """
            ${REF_ENV}
            vamos --read -b ${aln} -r ${catalog} -s ${sample} -o ${sample}.vamos.vcf -S ${v.extra} -t ${cpus} ${TOOL_ERR}
            ${DECODE_CHECK}
            """.stripIndent()
        }
        throw new IllegalArgumentException("unknown tool ${tool}")
    }
}
