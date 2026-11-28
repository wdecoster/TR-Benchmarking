#!/usr/bin/env python3
"""Metrics for the tool outputs of one sample: call rate per motif class and, with a truth file, concordance with it.

usage: evaluate.py --catalog benchmark-catalog-v2.strdust.bed --results <outdir>/default [--truth truth.tsv] [--truth-seq truth.seq.tsv] [--json out.json]

Output files are found as written by main.nf: <results>/<tool>/<version>/<sample>.<suffix>. The truth file has one
line per catalog locus with allele lengths of the assembly haplotypes: locus ID, then length_h1[, length_h2]
(tab separated, haploid loci have one length). It is made by `truth_from_assembly.py` from the aligned assemblies, which also
writes the allele sequences (`--seq-out`, same layout) for the sequence-level metrics (`--truth-seq`; tools without sequences,
inquiSTR, mahi-mahi and straglr, have none).
"""
import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calls as C
import metrics as M

# <results>/<tool>/<version>/<this> as written by main.nf (see TrCommands.output)
OUTPUTS = {
    'atarva': '{sample}.atarva.vcf', 'inquistr': '{sample}.inquistr.vcf', 'inquistr_unphased': '{sample}.inquistr.vcf',
    'longtr': '{sample}.longTR.vcf.gz', 'mahimahi': '{sample}.mahimahi/{sample}.vcf',
    'medaka': '{sample}.medaka.vcf/medaka_to_ref.TR.vcf', 'straglr': '{sample}.vcf', 'strdust': '{sample}.strdust.sorted.vcf',
    'strkit': '{sample}.strkit.vcf', 'vamos': '{sample}.vamos.fixed.vcf',
}


# when the output of a process is a directory (medaka, mahi-mahi), the VCF inside it
INNER = {'medaka': 'medaka_to_ref.TR.vcf', 'mahimahi': '{sample}.vcf'}


def resolve(tool, path, sample):
    """The VCF of one tool output given as the workflow's output path (a file, or a directory for medaka and mahi-mahi)."""
    if os.path.isdir(path):
        inner = INNER.get(tool)
        if not inner:
            sys.exit(f'{path} is a directory but {tool} writes a single file')
        path = os.path.join(path, inner.format(sample=sample))
    return path


def find_outputs(results, sample):
    for tool, pattern in OUTPUTS.items():
        for vdir in sorted(glob.glob(os.path.join(results, tool, '*'))):
            path = os.path.join(vdir, pattern.format(sample=sample))
            if os.path.exists(path):
                yield tool, os.path.basename(vdir), path


def read_truth(path):
    truth = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[0].startswith('#') or len(f) < 2:
                continue
            truth[f[0]] = [int(x) for x in f[1:] if x not in ('', 'NA')]
    return truth


def read_truth_seq(path):
    truth = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[0].startswith('#') or len(f) < 2:
                continue
            truth[f[0]] = [x for x in f[1:]]
    return truth


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--results', help='folder <results>/<tool>/<version>/<sample>.* as written by main.nf')
    ap.add_argument('--tool', help='with --vcf: tool of this output (see calls.TOOLS)')
    ap.add_argument('--version', default='', help='with --vcf: tool version, for the labels')
    ap.add_argument('--vcf', help='one tool output (file, or directory for medaka and mahi-mahi) instead of --results')
    ap.add_argument('--sample', default='HG002.30x.haplotagged')
    ap.add_argument('--karyotype', default='XY', choices=['XX', 'XY'], help='haploid chrX/chrY for XY')
    ap.add_argument('--truth')
    ap.add_argument('--truth-seq', help='assembly allele sequences (truth_from_assembly.py --seq-out): adds the sequence-level metrics')
    ap.add_argument('--json')
    a = ap.parse_args()
    catalog = C.read_catalog(a.catalog)
    truth = read_truth(a.truth) if a.truth else {}
    seq_truth = read_truth_seq(a.truth_seq) if a.truth_seq else None
    out = {}
    if a.vcf:
        if not a.tool:
            sys.exit('--vcf needs --tool')
        found = [(a.tool, a.version, resolve(a.tool, a.vcf, a.sample))]
    elif a.results:
        found = list(find_outputs(a.results, a.sample))
    else:
        sys.exit('give --results or --vcf')
    if not found:
        sys.exit(f'no tool outputs for sample {a.sample} under {a.results}')
    for tool, version, path in found:
        calls = C.read_calls(tool, path, catalog, a.karyotype)
        out[f'{tool} {version}'] = {an: M.concordance(calls, truth, catalog, C.not_supported(tool, catalog), an) for an in M.ANALYSES}
        if seq_truth is not None:
            for an in M.ANALYSES:
                sq = M.sequence_concordance(calls, seq_truth, catalog, an)
                if sq:
                    for cls, r in sq.items():
                        out[f'{tool} {version}'][an][cls].update(r)
    cols = ['n_loci', 'n_not_supported', 'n_called', 'call_rate', 'call_rate_unflagged', 'n_eval', 'MATCH', 'ONE-OFF', 'MOTIF-OFF', 'MISMATCH',
            'exact_genotype', 'exact_of_truth_loci', 'mean_abs_diff']
    if seq_truth is not None:
        cols += ['n_truth_seq', 'n_seq_eval', 'seq_exact_alleles', 'seq_exact_genotype', 'seq_exact_of_truth_loci',
                 'mean_edit_distance', 'mean_norm_edit_distance']
    print('tool version\tanalysis\tclass\t' + '\t'.join(cols))
    f = lambda x: '' if x is None else (f'{x:.4f}' if isinstance(x, float) else str(x))
    for k, by_analysis in out.items():
        for an, by_class in by_analysis.items():
            for cls, r in by_class.items():
                print(f'{k}\t{an}\t{cls}\t' + '\t'.join(f(r.get(c)) for c in cols))
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
