#!/usr/bin/env python3
"""Allele lengths of the assembly haplotypes at every catalog locus: the truth for evaluate.py.

Same method as the paper's scripts/assembly_concordance.py (it calls `cooper` of scripts/aam_mod.py): every
assembly contig that aligns with MAPQ >= 5 and spans the whole locus gives an allele length from its CIGAR; a
haplotype has a truth length only when all its spanning contigs agree. A diploid locus needs both haplotypes; on chrX
and chrY of an XY sample the single haplotype counts (chrX from the maternal, chrY from the paternal assembly).

usage: truth_from_assembly.py --mat mat.bam --pat pat.bam --catalog benchmark-catalog-v2.strdust.bed --karyotype XY --out truth.tsv
Output: locus ID <TAB> length (<TAB> length), shorter allele first; loci without truth are absent.
"""
import argparse
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from aam_mod import cooper      # noqa: E402


def truth(mat_bam, pat_bam, catalog_bed, karyotype):
    with tempfile.TemporaryDirectory() as d:
        sorted_bed = os.path.join(d, 'catalog.bed')
        with open(catalog_bed) as fh, open(sorted_bed, 'w') as out:
            rows = [l.split('\t')[:3] for l in fh if l.strip() and not l.startswith('#')]
            for r in sorted(rows, key=lambda r: (r[0], int(r[1]), int(r[2]))):
                out.write('\t'.join(r) + '\n')
        pysam.tabix_compress(sorted_bed, sorted_bed + '.gz', force=True)
        pysam.tabix_index(sorted_bed + '.gz', preset='bed', force=True)
        tbx = pysam.TabixFile(sorted_bed + '.gz')
        res = {}
        for name, path in (('mat', mat_bam), ('pat', pat_bam)):
            with pysam.AlignmentFile(path, 'rb') as bam:
                res[name] = {k: v[1] for k, v in cooper(bam, tbx, 5, 0).items()}
        tbx.close()
    out = {}
    for k in set(res['mat']) | set(res['pat']):
        chrom = k.split('-')[0]
        m, p = res['mat'].get(k), res['pat'].get(k)
        if karyotype == 'XY' and chrom == 'chrX':
            lens = [m] if m is not None else []
        elif karyotype == 'XY' and chrom == 'chrY':
            lens = [p] if p is not None else []
        else:
            lens = [m, p] if m is not None and p is not None else []
        if lens:
            out[k] = sorted(lens)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mat', required=True)
    ap.add_argument('--pat', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--karyotype', required=True, choices=['XX', 'XY'])
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    t = truth(a.mat, a.pat, a.catalog, a.karyotype)
    with open(a.out, 'w') as fh:
        fh.write('#locus\tlength_h1\tlength_h2\n')
        for k in sorted(t):
            fh.write(k + '\t' + '\t'.join(map(str, t[k])) + '\n')
    n = sum(1 for l in open(a.catalog) if l.strip() and not l.startswith('#'))
    print(f'{len(t)} of {n} catalog loci have a truth length', file=sys.stderr)


if __name__ == '__main__':
    main()
