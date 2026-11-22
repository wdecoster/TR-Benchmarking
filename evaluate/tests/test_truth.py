"""Needs pysam (skipped otherwise): python3 evaluate/tests/test_truth.py"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
try:
    import pysam
except ImportError:
    pysam = None


def make_bam(path, reads):
    hdr = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': c, 'LN': 1000} for c in ('chr1', 'chrX', 'chrY')]}
    with pysam.AlignmentFile(path, 'wb', header=hdr) as out:
        for chrom, start, cigar, n in reads:
            a = pysam.AlignedSegment()
            a.query_name = f'c{chrom}{start}'
            a.reference_id = ['chr1', 'chrX', 'chrY'].index(chrom)
            a.reference_start, a.mapping_quality, a.cigarstring = start, 60, cigar
            a.query_sequence = 'A' * n
            out.write(a)
    pysam.sort('-o', path, path)
    pysam.index(path)


def test_truth_lengths_and_sex_chromosomes():
    from truth_from_assembly import truth
    with tempfile.TemporaryDirectory() as d:
        make_bam(f'{d}/mat.bam', [('chr1', 50, '60M3I90M', 153), ('chrX', 150, '150M', 150)])
        make_bam(f'{d}/pat.bam', [('chr1', 50, '150M', 150), ('chrY', 250, '150M', 150)])
        with open(f'{d}/cat.bed', 'w') as fh:
            fh.write('chr1\t100\t130\tCAG\nchrX\t200\t220\tAC\nchrY\t300\t320\tAC\n')
        t = truth(f'{d}/mat.bam', f'{d}/pat.bam', f'{d}/cat.bed', 'XY')
    assert t == {'chr1-100-130': [30, 33], 'chrX-200-220': [20], 'chrY-300-320': [20]}


if __name__ == '__main__':
    if pysam is None:
        print('pysam not installed: skipped')
    else:
        test_truth_lengths_and_sex_chromosomes()
        print('1 test passed')
