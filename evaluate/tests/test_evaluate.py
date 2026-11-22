import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import calls as C
import metrics as M

CAT = {
    'chr1-100-130': C.Locus('chr1', 100, 130, 'CAG'),         # 30 bp, STR
    'chr1-200-260': C.Locus('chr1', 200, 260, 'CAGCAGCAGCAGCAG'),   # VNTR
    'chr2-50-60': C.Locus('chr2', 50, 60, 'A'),               # homopolymer, mahi-mahi cannot take it
}


def vcf(tmp_path, lines, name='x.vcf'):
    p = tmp_path / name
    p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n' + '\n'.join(lines) + '\n')
    return str(p)


def test_status():
    assert [M.status(d, 3) for d in (0, 1, -1, 2, 3, -3, 4)] == ['MATCH', 'ONE-OFF', 'ONE-OFF', 'MOTIF-OFF', 'MOTIF-OFF', 'MOTIF-OFF', 'MISMATCH']


def test_pairing_ignores_haplotype_order():
    assert M.pair_alleles([30, 45], [45, 30]) == [(30, 30), (45, 45)] or sorted(M.pair_alleles([30, 45], [45, 30])) == [(30, 30), (45, 45)]
    assert M.pair_alleles([40], [30, 50]) == [(40, 30), (40, 50)]       # one called allele vs two truth alleles
    assert M.pair_alleles([30, 60], [31]) == [(30, 31)]                   # haploid truth: closest called allele


def test_strdust_uses_frb(tmp_path):
    p = vcf(tmp_path, ['chr1\t101\tchr1-100-130\tA\tT\t.\t.\tEND=130\tGT:RB:FRB\t0/1:0,3:30,33'])
    assert C.read_calls('strdust', p, CAT) == {'chr1-100-130': C.Call([30, 33], '')}


def test_strdust_quickref_record_without_lengths_is_reference(tmp_path):
    p = vcf(tmp_path, ['chr1\t101\t.\t' + 'A' * 30 + '\t.\t.\t.\tQUICKREF;END=130\tGT:SUP\t0|0:.,.'])
    assert C.read_calls('strdust', p, CAT)['chr1-100-130'].lengths == [30, 30]


def test_vcf_without_genotype_is_a_no_call(tmp_path):
    p = vcf(tmp_path, ['chr1\t101\tchr1-100-130\tA\t.\t.\t.\tEND=130\tGT:RB:FRB\t./.:.:.'])
    assert C.read_calls('strdust', p, CAT)['chr1-100-130'].lengths is None


def test_ref_alt_length_is_independent_of_padding(tmp_path):
    # strkit anchors one base before the locus: REF 31 bp, ALT 3 bp longer than REF -> 33 bp allele
    p = vcf(tmp_path, ['chr1\t100\tlocus1\t' + 'A' * 31 + '\t' + 'A' * 34 + '\t.\t.\tBED_START=100;BED_END=130\tGT\t1/0'])
    assert C.read_calls('strkit', p, CAT)['chr1-100-130'].lengths == [34 - 31 + 30, 30]


def test_mahimahi_homozygous_has_one_value_and_homopolymers_are_not_supported(tmp_path):
    p = vcf(tmp_path, ['chr1\t100\t.\tA\t<STR3>\t.\tMotifBias\tEND=130;REPID=chr1-100-130\tGT:BPCN\t1/1:30'])
    c = C.read_calls('mahimahi', p, CAT)['chr1-100-130']
    assert c.lengths == [30, 30] and c.flags == 'MotifBias'
    assert C.not_supported('mahimahi', CAT) == ['chr2-50-60']
    assert C.not_supported('strdust', CAT) == []


def test_concordance_by_class_and_not_supported():
    calls = {'chr1-100-130': C.Call([30, 33]), 'chr1-200-260': C.Call(None)}
    truth = {'chr1-100-130': [33, 30], 'chr1-200-260': [60, 60], 'chr2-50-60': [10, 10]}
    r = M.concordance(calls, truth, CAT, not_supported=['chr2-50-60'])
    assert r['STR']['n_loci'] == 2 and r['STR']['n_not_supported'] == 1       # the homopolymer
    assert r['STR']['call_rate'] == 1.0 and r['STR']['MATCH'] == 1.0 and r['STR']['exact_genotype'] == 1.0
    assert r['VNTR']['n_record'] == 1 and r['VNTR']['n_called'] == 0 and r['VNTR']['call_rate'] == 0.0
    assert r['VNTR']['MATCH'] is None                                         # nothing to compare
    assert r['all']['n_loci'] == 3 and r['all']['n_called'] == 1


if __name__ == '__main__':      # runs without pytest: python3 tests/test_evaluate.py
    import inspect
    import pathlib
    import tempfile
    n = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            with tempfile.TemporaryDirectory() as d:
                fn(pathlib.Path(d)) if inspect.signature(fn).parameters else fn()
            n += 1
    print(f'{n} tests passed')
