#!/usr/bin/env python3
"""Replicate agreement from the calls that timing.nf keeps (timing/calls/calls_<session>_<order>_<tool>_<version>_rep<N>.tsv.gz).

For every tool version with at least two replicates of the same input:
  n_loci            loci with a record in at least one replicate
  record_identical  loci whose VCF record is byte-identical in all replicates (and present in all)
  length_identical  loci called in every replicate with identical allele lengths (haplotype order ignored)
  callset_stable    loci called in all replicates or in none; the rest are called in some replicates only
  max_len_diff, median_len_diff   over loci called in every replicate whose lengths differ: the largest
                    allele-length difference between any two replicates, in bp (blank if no locus differs)
With three replicates this is descriptive; no confidence intervals are given.
"""
import argparse
import glob
import gzip
import json
import os
import re
import statistics
import sys
from collections import defaultdict
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calls as C

NAME = re.compile(r'calls_(?P<session>.+?)_(?P<order>\d+)_(?P<tool>[a-z]+)_(?P<version>.+)_rep(?P<rep>\d+)\.tsv\.gz$')


def records(path):
    """The VCF record lines of a calls file."""
    with gzip.open(path, 'rt') as fh:
        return [l.rstrip('\n') for l in fh if not l.startswith('#')]


def agreement(tool, files, catalog):
    reps = []
    for f in files:
        # the calls file is a header-less VCF, which read_calls reads as it is
        reps.append((records(f), C.read_calls(tool, f, catalog)))
    keys = sorted(set().union(*[set(c) for _, c in reps]))
    # raw record per locus (key -> line) for the byte-identity check
    raw = []
    for lines, c in reps:
        d = {}
        for line in lines:
            rec = line.split('\t')
            k = C._key(tool, rec, C._info(rec[7]))
            if k in catalog:
                d[k] = line
        raw.append(d)
    n = len(reps)
    rec_ident = len_ident = stable = 0
    diffs = []
    for k in keys:
        called = [c[k].lengths if k in c and c[k].lengths else None for _, c in reps]
        if all(k in r for r in raw) and len({r[k] for r in raw}) == 1:
            rec_ident += 1
        n_called = sum(x is not None for x in called)
        if n_called in (0, n):
            stable += 1
        if n_called == n:
            sets = [tuple(sorted(x)) for x in called]
            if len(set(sets)) == 1:
                len_ident += 1
            else:
                # largest allele-length difference between two replicates (a different ploidy counts as 0 here)
                diffs.append(max((max(abs(a - b) for a, b in zip(s1, s2))
                                  for s1, s2 in combinations(sets, 2) if len(s1) == len(s2)), default=0))
    return dict(n_reps=n, n_loci=len(keys), record_identical=rec_ident / len(keys) if keys else None,
                length_identical=len_ident / len(keys) if keys else None,
                callset_stable=stable / len(keys) if keys else None,
                max_len_diff=max(diffs) if diffs else None, median_len_diff=statistics.median(diffs) if diffs else None,
                n_loci_length_differs=len(diffs))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('calls_dir', help='timing/calls')
    ap.add_argument('--catalog', required=True, help='full catalog (benchmark-catalog-v2.strdust.bed) or its chromosome subset')
    ap.add_argument('--json', help='write the result as JSON')
    ap.add_argument('--min-reps', type=int, default=2)
    a = ap.parse_args()
    catalog = C.read_catalog(a.catalog)
    groups = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(a.calls_dir, 'calls_*.tsv.gz'))):
        m = NAME.search(os.path.basename(f))
        if m:
            groups[(m['tool'], m['version'])].append(f)
    result = {}
    for (tool, version), files in sorted(groups.items()):
        if len(files) >= a.min_reps:
            result[f'{tool} {version}'] = agreement(tool, files, catalog)
    print('tool version\treps\tloci\trecord_identical\tlength_identical\tcallset_stable\tloci_len_differs\tmedian_diff\tmax_diff')
    for k, r in result.items():
        f = lambda x: '' if x is None else (f'{x:.4f}' if isinstance(x, float) else str(x))
        print('\t'.join([k, str(r['n_reps']), str(r['n_loci']), f(r['record_identical']), f(r['length_identical']),
                         f(r['callset_stable']), str(r['n_loci_length_differs']), f(r['median_len_diff']), f(r['max_len_diff'])]))
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(result, fh, indent=1)


if __name__ == '__main__':
    main()
