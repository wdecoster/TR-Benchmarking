"""Metrics on length-only calls (see calls.py). All functions are pure; the scripts do the file handling."""
from collections import Counter
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from calls import Call, Locus

# Outcome of one allele against the assembly, as in the paper's assembly_concordance.py (|difference| in bp):
# 0 MATCH, 1 ONE-OFF, <= the motif length MOTIF-OFF (a whole repeat unit), larger MISMATCH.
STATUSES = ('MATCH', 'ONE-OFF', 'MOTIF-OFF', 'MISMATCH')


def status(diff: int, motif_len: int) -> str:
    d = abs(diff)
    if d == 0:
        return 'MATCH'
    if d == 1:
        return 'ONE-OFF'
    if d <= motif_len:
        return 'MOTIF-OFF'
    return 'MISMATCH'


def pair_alleles(tool: Sequence[int], truth: Sequence[int]) -> List[Tuple[int, int]]:
    """Pair the called alleles with the assembly alleles so that the total length difference is smallest.

    Diploid vs diploid: the better of the two pairings (haplotype order is not trusted). Anything else (haploid
    chrX/chrY of XY samples, or a tool that gave one value where the truth has two): each truth allele is paired
    with the closest called allele.
    """
    t, a = list(tool), list(truth)
    if len(t) == 2 and len(a) == 2:
        p1 = [(t[0], a[0]), (t[1], a[1])]
        p2 = [(t[0], a[1]), (t[1], a[0])]
        return p1 if sum(abs(x - y) for x, y in p1) <= sum(abs(x - y) for x, y in p2) else p2
    return [(min(t, key=lambda x: abs(x - y)), y) for y in a]


ANALYSES = {
    'all_loci': lambda l: True,                          # the main analysis: every catalog locus, nothing left out
    'no_homopolymers': lambda l: len(l.motif) > 1,       # the same without 1 bp motifs, for every tool alike
}


def concordance(calls: Dict[str, Call], truth: Dict[str, List[int]], catalog: Dict[str, Locus],
                not_supported: Iterable[str] = (), analysis: str = 'all_loci') -> Dict[str, dict]:
    """Per motif class ('STR', 'VNTR', 'all'): how many catalog loci the tool covers and how its alleles compare.

    analysis          'all_loci' (default, the main analysis) or 'no_homopolymers' (loci with a 1 bp motif left out
                      for every tool); loci a tool cannot take count as not called in 'all_loci'
    n_loci            loci in the analysis
    n_not_supported   of those, loci the tool cannot take at all (mahi-mahi and 1 bp motifs); informational, they are
                      not called and stay in the denominators
    n_record, n_called  loci with a VCF record / with a genotype
    call_rate         n_called / n_loci
    n_flagged         called loci whose record carries a flag (STRdust CLUSTERFAILURE, a mahi-mahi filter):
                      the genotype is still used; call_rate_unflagged leaves such loci out of the numerator
    n_eval            called loci that also have an assembly genotype; their alleles are compared
    <status>          fraction of compared alleles per status; also 'exact_genotype' (every allele MATCH, of the
                      evaluated loci) and 'mean_abs_diff' (bp)
    n_truth, exact_of_truth_loci
                      loci with an assembly genotype, and the fraction of them with an exact genotype: loci without
                      a call count as wrong, so a tool cannot gain by skipping hard loci
    """
    skip = set(not_supported)
    keep = ANALYSES[analysis]
    out = {}
    for cls in ('STR', 'VNTR', 'all'):
        loci = [l for l in catalog.values() if keep(l) and (cls == 'all' or l.motif_class == cls)]
        n_rec = n_called = n_flagged = n_eval = n_alleles = n_exact = n_truth = 0
        counts: Counter = Counter()
        abs_sum = 0
        for l in loci:
            t = truth.get(l.id)
            n_truth += bool(t)
            c = calls.get(l.id)
            if c is None:
                continue
            n_rec += 1
            if c.lengths is None:
                continue
            n_called += 1
            n_flagged += bool(c.flags)
            if not t:
                continue
            n_eval += 1
            sts = []
            for x, y in pair_alleles(c.lengths, t):
                s = status(x - y, len(l.motif))
                counts[s] += 1
                abs_sum += abs(x - y)
                n_alleles += 1
                sts.append(s)
            n_exact += all(s == 'MATCH' for s in sts)
        n = len(loci)
        r = dict(n_loci=n, n_not_supported=sum(l.id in skip for l in loci), n_record=n_rec, n_called=n_called,
                 call_rate=n_called / n if n else None, n_flagged=n_flagged,
                 call_rate_unflagged=(n_called - n_flagged) / n if n else None, n_eval=n_eval, n_truth=n_truth)
        for s in STATUSES:
            r[s] = counts[s] / n_alleles if n_alleles else None
        r['exact_genotype'] = n_exact / n_eval if n_eval else None
        r['exact_of_truth_loci'] = n_exact / n_truth if n_truth else None
        r['mean_abs_diff'] = abs_sum / n_alleles if n_alleles else None
        out[cls] = r
    return out


# ---------------------------------------------------------------------------------------------------------------
# Sequence level: edit distance between the called and the assembly allele sequence over the catalog interval.

def edit_distance(a: str, b: str) -> int:
    """Levenshtein distance (bit-parallel algorithm of Myers 1999 in Hyyro's form, on Python integers: no dependency)."""
    if a == b:
        return 0
    if not a or not b:
        return len(a) + len(b)
    if len(a) < len(b):
        a, b = b, a                       # the bit vectors run over the shorter string
    m = len(b)
    mask = (1 << m) - 1
    high = 1 << (m - 1)
    peq: Dict[str, int] = {}
    for i, c in enumerate(b):
        peq[c] = peq.get(c, 0) | (1 << i)
    pv, mv, score = mask, 0, m
    for c in a:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) & mask) ^ pv) | eq
        ph = (mv | ~(xh | pv)) & mask
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph = ((ph << 1) | 1) & mask
        mh = (mh << 1) & mask
        pv = (mh | ~(xv | ph)) & mask
        mv = ph & xv
    return score


def pair_sequences(tool: Sequence[str], truth: Sequence[str]) -> List[Tuple[int, int, int]]:
    """Pair called and assembly sequences like pair_alleles, but on edit distance. Returns (edit distance, length of the
    called allele, length of the assembly allele) per assembly allele."""
    t, a = list(tool), list(truth)

    def d(x, y):
        return edit_distance(x, y)
    if len(t) == 2 and len(a) == 2:
        p1 = [(t[0], a[0]), (t[1], a[1])]
        p2 = [(t[0], a[1]), (t[1], a[0])]
        d1 = [d(x, y) for x, y in p1]
        d2 = [d(x, y) for x, y in p2]
        pairs, ds = (p1, d1) if sum(d1) <= sum(d2) else (p2, d2)
    else:
        pairs, ds = [], []
        for y in a:
            best = min(t, key=lambda x: d(x, y))
            pairs.append((best, y))
            ds.append(d(best, y))
    return [(e, len(x), len(y)) for e, (x, y) in zip(ds, pairs)]


def sequence_concordance(calls: Dict[str, Call], seq_truth: Dict[str, List[str]], catalog: Dict[str, Locus],
                         analysis: str = 'all_loci') -> Optional[Dict[str, dict]]:
    """Per motif class: how well the called allele sequences equal the assembly sequences. None for a tool without sequences.

    n_truth_seq            loci with assembly sequences in the analysis
    n_seq_eval             of those, loci called with sequences (the alleles that are compared)
    seq_exact_alleles      fraction of compared alleles with edit distance 0
    seq_exact_genotype     fraction of the compared loci with every allele exact
    seq_exact_of_truth_loci  exact loci / n_truth_seq: a locus without a call counts as wrong
    mean_edit_distance     mean edit distance per compared allele (bp)
    mean_norm_edit_distance  mean of edit distance / max(called length, assembly length) per compared allele
    """
    if not any(c.seqs for c in calls.values()):
        return None
    keep = ANALYSES[analysis]
    out = {}
    for cls in ('STR', 'VNTR', 'all'):
        n_truth = n_eval = n_alleles = n_exact_alleles = n_exact_loci = 0
        ed_sum = 0
        norm_sum = 0.0
        for l in catalog.values():
            if not keep(l) or (cls != 'all' and l.motif_class != cls):
                continue
            t = seq_truth.get(l.id)
            if not t:
                continue
            n_truth += 1
            c = calls.get(l.id)
            if c is None or not c.seqs:
                continue
            n_eval += 1
            res = pair_sequences(c.seqs, t)
            for e, lx, ly in res:
                n_alleles += 1
                n_exact_alleles += e == 0
                ed_sum += e
                norm_sum += e / max(lx, ly, 1)
            n_exact_loci += all(e == 0 for e, _, _ in res)
        out[cls] = dict(n_truth_seq=n_truth, n_seq_eval=n_eval,
                        seq_exact_alleles=n_exact_alleles / n_alleles if n_alleles else None,
                        seq_exact_genotype=n_exact_loci / n_eval if n_eval else None,
                        seq_exact_of_truth_loci=n_exact_loci / n_truth if n_truth else None,
                        mean_edit_distance=ed_sum / n_alleles if n_alleles else None,
                        mean_norm_edit_distance=norm_sum / n_alleles if n_alleles else None)
    return out
