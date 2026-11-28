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


CATX = dict(CAT, **{'chrX-500-530': C.Locus('chrX', 500, 530, 'CAG')})


def vcf(tmp_path, lines, name='x.vcf'):
    p = tmp_path / name
    p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n' + '\n'.join(lines) + '\n')
    return str(p)


def test_inquistr_haploid_nan_second_value_is_a_valid_call_on_chrx_only(tmp_path):
    p = vcf(tmp_path, ['chrX\t501\t.\tN\t<STR>\t.\t.\tEND=530\tGT:RB:FRB\t1|2:3,NaN:33,NaN',
                       'chr1\t101\t.\tN\t<STR>\t.\t.\tEND=130\tGT:RB:FRB\t1|2:3,NaN:33,NaN'])
    c = C.read_calls('inquistr', p, CATX)
    assert c['chrX-500-530'].lengths == [33]                      # haploid: only the first value counts
    assert c['chr1-100-130'].lengths is None                      # diploid with a missing allele: no call
    assert C.read_calls('inquistr', p, CATX, karyotype='XX')['chrX-500-530'].lengths is None


def test_two_allele_tool_is_truncated_to_one_allele_on_chrx(tmp_path):
    p = vcf(tmp_path, ['chrX\t501\t.\t' + 'A' * 30 + '\t' + 'A' * 33 + '\t.\t.\tEND=530\tGT\t0/1'])
    assert C.read_calls('longtr', p, CATX)['chrX-500-530'].lengths == [30]


def test_deleted_repeat_array_is_an_empty_allele(tmp_path):
    p = vcf(tmp_path, ['chr1\t101\tchr1-100-130\t' + 'A' * 30 + '\t<DEL>\t.\t.\tEND=130\tGT\t0/1'])
    assert C.read_calls('longtr', p, CAT)['chr1-100-130'].lengths == [30, 0]


def test_straglr_gt_forms(tmp_path):
    p = vcf(tmp_path, ['chr1\t100\t.\tA\t<CNV:TR>\t.\tPASS\tEND=130;RB=45\tGT\t1',        # homozygous, one index
                       'chr1\t200\t.\tA\t\t.\tPASS\tEND=260\tGT\t0'])                       # reference, no RB
    c = C.read_calls('straglr', p, CAT)
    assert c['chr1-100-130'].lengths == [45, 45] and c['chr1-200-260'].lengths == [60, 60]


def test_wrong_key_scheme_is_an_error_not_zero_calls(tmp_path):
    p = vcf(tmp_path, ['chr1\t999\tx\tA\tT\t.\t.\tEND=1000\tGT:FRB\t0/1:1,2'] * 3)
    try:
        C.read_calls('inquistr', p, CAT)
    except ValueError:
        return
    raise AssertionError('expected a ValueError')


def test_replicate_file_names_with_underscored_tools():
    import replicate_agreement as R
    m = R.NAME.search('calls_2026-11-22_13-28_3_inquistr_unphased_0.28.2_rep2.tsv.gz')
    assert (m['tool'], m['version'], m['rep']) == ('inquistr_unphased', '0.28.2', '2')
    assert R.NAME.search('calls_s_1_vamos_old_2.1.7_rep1.tsv.gz')['tool'] == 'vamos_old'


def test_status():
    assert [M.status(d, 3) for d in (0, 1, -1, 2, 3, -3, 4)] == ['MATCH', 'ONE-OFF', 'ONE-OFF', 'MOTIF-OFF', 'MOTIF-OFF', 'MOTIF-OFF', 'MISMATCH']


def test_pairing_ignores_haplotype_order():
    assert M.pair_alleles([30, 45], [45, 30]) == [(30, 30), (45, 45)] or sorted(M.pair_alleles([30, 45], [45, 30])) == [(30, 30), (45, 45)]
    assert M.pair_alleles([40], [30, 50]) == [(40, 30), (40, 50)]       # one called allele vs two truth alleles
    assert M.pair_alleles([30, 60], [31]) == [(30, 31)]                   # haploid truth: closest called allele


def test_strdust_uses_frb(tmp_path):
    p = vcf(tmp_path, ['chr1\t101\tchr1-100-130\tA\tT\t.\t.\tEND=130\tGT:RB:FRB\t0/1:0,3:30,33'])
    c = C.read_calls('strdust', p, CAT)['chr1-100-130']
    assert c.lengths == [30, 33] and c.flags == ''      # the sequences (here the one-base REF/ALT) are not used for the length


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


def test_concordance_all_loci_counts_unsupported_as_not_called_and_homopolymer_free_analysis():
    calls = {'chr1-100-130': C.Call([30, 33]), 'chr1-200-260': C.Call(None)}
    truth = {'chr1-100-130': [33, 30], 'chr1-200-260': [60, 60], 'chr2-50-60': [10, 10]}
    r = M.concordance(calls, truth, CAT, not_supported=['chr2-50-60'])                 # main analysis: every locus
    assert r['STR']['n_loci'] == 2 and r['STR']['n_not_supported'] == 1                  # the homopolymer stays in
    assert r['STR']['call_rate'] == 0.5 and r['STR']['MATCH'] == 1.0 and r['STR']['exact_genotype'] == 1.0
    assert r['STR']['n_truth'] == 2 and r['STR']['exact_of_truth_loci'] == 0.5           # the unsupported locus counts as wrong
    assert r['VNTR']['n_record'] == 1 and r['VNTR']['n_called'] == 0 and r['VNTR']['call_rate'] == 0.0
    assert r['VNTR']['MATCH'] is None                                                    # nothing to compare
    assert r['all']['n_loci'] == 3 and r['all']['n_called'] == 1
    h = M.concordance(calls, truth, CAT, not_supported=['chr2-50-60'], analysis='no_homopolymers')
    assert h['STR']['n_loci'] == 1 and h['STR']['n_not_supported'] == 0 and h['STR']['call_rate'] == 1.0
    assert h['all']['n_loci'] == 2 and h['all']['exact_of_truth_loci'] == 0.5

def test_edit_distance_matches_a_reference_implementation():
    import random

    def dp(a, b):
        prev = list(range(len(b) + 1))
        for i, x in enumerate(a, 1):
            cur = [i]
            for j, y in enumerate(b, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
            prev = cur
        return prev[-1]
    rnd = random.Random(1)
    for _ in range(500):
        al = rnd.choice(['AC', 'ACGT', 'A'])
        a = ''.join(rnd.choice(al) for _ in range(rnd.randint(0, 30)))
        b = ''.join(rnd.choice(al) for _ in range(rnd.randint(0, 30)))
        assert M.edit_distance(a, b) == dp(a, b)


def test_allele_sequence_is_cut_to_the_catalog_interval(tmp_path):
    # catalog chr1:100-130 (30 bp). The record has 5 flank bases on both sides in REF and ALT (LongTR style)
    flank_l, flank_r = 'TTTTT', 'GGGGG'
    ref = 'CAG' * 10
    alt = 'CAG' * 12
    p = vcf(tmp_path, [f'chr1\t96\tchr1-100-130\t{flank_l}{ref}{flank_r}\t{flank_l}{alt}{flank_r}\t.\t.\tEND=130\tGT\t1|0'])
    c = C.read_calls('longtr', p, CAT)['chr1-100-130']
    assert c.lengths == [36, 30] and c.seqs == [alt, ref]


def test_length_only_tool_has_no_sequences_and_vamos_alt_is_taken_whole(tmp_path):
    p = vcf(tmp_path, ['chr1\t101\t.\tN\t<STR>\t.\t.\tEND=130\tGT:RB:FRB\t1|1:3,3:33,33'])
    assert C.read_calls('inquistr', p, CAT)['chr1-100-130'].seqs is None
    alt = 'CAG' * 11
    p = vcf(tmp_path, [f'chr1\t101\t.\tA\t{alt}\t.\t.\tEND=130\tGT\t1/1'], 'v.vcf')
    c = C.read_calls('vamos', p, CAT)['chr1-100-130']
    assert c.seqs == [alt, alt]


def test_sequence_concordance_counts_a_missing_call_as_wrong_and_pairs_by_edit_distance():
    ref = 'CAG' * 10
    other = 'CAG' * 9 + 'CAA'
    cat = {'chr1-100-130': C.Locus('chr1', 100, 130, 'CAG'), 'chr1-200-260': C.Locus('chr1', 200, 260, 'CAGCAGCAGCAGCAG')}
    calls = {'chr1-100-130': C.Call([30, 30], '', [other, ref])}               # haplotype order swapped
    truth = {'chr1-100-130': [ref, other], 'chr1-200-260': ['A' * 60, 'A' * 60]}
    r = M.sequence_concordance(calls, truth, cat)
    assert r['all']['n_truth_seq'] == 2 and r['all']['n_seq_eval'] == 1
    assert r['all']['seq_exact_genotype'] == 1.0 and r['all']['seq_exact_of_truth_loci'] == 0.5
    assert M.sequence_concordance({'x': C.Call([1], '', None)}, truth, cat) is None


def test_sequence_edit_distance_is_reported():
    cat = {'chr1-100-130': C.Locus('chr1', 100, 130, 'CAG')}
    ref = 'CAG' * 10
    calls = {'chr1-100-130': C.Call([33, 30], '', [ref + 'CAG', ref])}
    r = M.sequence_concordance(calls, {'chr1-100-130': [ref, ref]}, cat)
    assert r['STR']['mean_edit_distance'] == 1.5 and r['STR']['seq_exact_alleles'] == 0.5


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
