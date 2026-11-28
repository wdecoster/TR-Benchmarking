#!/usr/bin/env python3
"""Summary report for the public dashboard (docs/dashboard.md of the planning repo; schema_version 1).

Reads what the pipeline computed and writes ONE json with aggregates only: no sample IDs besides the public reference sample,
no per-locus rows, no calls, no server paths.

  scripts/export_report.py --report-id 2026-11-28 --commit <sha> \
      --eval default=eval_truth.json default=eval_medaka.json sensitive=eval_inq_sens.json \
      --timing timing_chr20_<session>.tsv --determinism replicate_agreement.json --out report_2026-11-28.json

--eval MODE=FILE      evaluate.py --json output (tool version -> analysis -> class -> metrics); repeat for several files
--timing FILE         timing table of timing.nf for one session (rows: tool, version, rep, wall_s); the first --timing is the reference session
--determinism FILE    evaluate/replicate_agreement.py --json output
"""
import argparse
import json
import statistics
import sys
from collections import defaultdict

# tool -> (label, output kind, developed by the maintainer of this benchmark, legacy)
TOOLS = {
    'atarva': ('ATaRVa', 'sequence', False, False), 'inquistr': ('inquiSTR (phased)', 'length-only', True, False),
    'inquistr_unphased': ('inquiSTR (unphased)', 'length-only', True, False), 'longtr': ('LongTR', 'sequence', False, False),
    'mahimahi': ('mahi-mahi', 'length-only', False, False), 'medaka': ('medaka tandem', 'sequence', False, False),
    'straglr': ('Straglr', 'length-only', False, True), 'strdust': ('STRdust', 'sequence', True, False),
    'strkit': ('STRkit', 'sequence', False, False), 'vamos': ('vamos', 'sequence', False, False),
}
# version sets of conf/versions_paper.config: the versions of the benchmark paper
PAPER = {'atarva': '0.5.0', 'longtr': '1.2', 'medaka': '2.1.1', 'straglr': '1.5.5', 'strdust': '0.16.0', 'strkit': '0.24.2', 'vamos': '3.0.5'}

CLASSES = {'STR': 'str', 'VNTR': 'long', 'all': 'all'}
# metric id -> (label, description, higher_is_better, how to compute from the evaluate.py numbers)
METRICS = {
    'call_rate': ('Call rate', 'Share of catalog loci with a genotype. Loci a tool cannot take count as not called.', True,
                  lambda r: r['call_rate']),
    'exact_of_truth_loci': ('Exact genotype, all truth loci', 'Share of the loci with an assembly genotype where every called allele has exactly the assembly length. A locus without a call counts as wrong, so skipping hard loci cannot help.', True,
                            lambda r: r['exact_of_truth_loci']),
    'allele_exact': ('Alleles exact', 'Of the alleles of called loci, the share with exactly the assembly length.', True, lambda r: r['MATCH']),
    'allele_within_1bp': ('Alleles within 1 bp', 'Share of called alleles whose length is within 1 bp of the assembly allele.', True,
                          lambda r: r['MATCH'] + r['ONE-OFF']),
    'allele_gross_error': ('Alleles off by more than a repeat unit', 'Share of called alleles more than one motif length away from the assembly allele (MISMATCH of the paper).', False,
                           lambda r: r['MISMATCH']),
}


def tool_key(tool, version):
    return f'{tool}@{version}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--report-id', required=True)
    ap.add_argument('--commit', required=True)
    ap.add_argument('--eval', nargs='+', required=True, metavar='MODE=FILE')
    ap.add_argument('--timing', nargs='*', default=[])
    ap.add_argument('--determinism')
    ap.add_argument('--notes', nargs='*', default=[])
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    tools, metrics = {}, []
    for spec in a.eval:
        mode, _, path = spec.partition('=')
        if mode not in ('default', 'sensitive'):
            sys.exit(f'--eval needs default=FILE or sensitive=FILE, got {spec}')
        for label, by_analysis in json.load(open(path)).items():
            tool, version = label.split(' ', 1)
            key = tool_key(tool, version)
            name, out, coi, legacy = TOOLS[tool]
            t = tools.setdefault(key, dict(key=key, tool=tool, label=name, version=version, output=out, coi=coi, legacy=legacy,
                                           paper_version=PAPER.get(tool) == version, modes=[]))
            if mode not in t['modes']:
                t['modes'].append(mode)
            for analysis, by_class in by_analysis.items():
                for cls, r in by_class.items():
                    for mid, (_, _, _, f) in METRICS.items():
                        v = f(r)
                        if v is None:
                            continue                        # not applicable / not computable: a missing row, never a 0
                        n = r['n_truth'] if mid == 'exact_of_truth_loci' else (r['n_loci'] if mid == 'call_rate' else r['n_eval'])
                        metrics.append(dict(tool=key, mode=mode, analysis=analysis, **{'class': CLASSES[cls]}, metric=mid,
                                            value=round(v, 5), n=n))

    timing = None
    if a.timing:
        rows = defaultdict(list)
        sessions = []
        for path in a.timing:
            with open(path) as fh:
                head = fh.readline().rstrip('\n').split('\t')
                for line in fh:
                    r = dict(zip(head, line.rstrip('\n').split('\t')))
                    rows[(r['tool'], r['version'])].append(float(r['wall_s']))
            sessions.append(path)
        runs = [dict(tool=tool_key(t, v), reps=len(w), wall_s=w, median_s=statistics.median(w), spread=round(max(w) / min(w), 3))
                for (t, v), w in sorted(rows.items())]
        timing = dict(scope='chr20 of HG002, one thread, three repetitions, one laptop session', runs=runs)

    determinism = None
    if a.determinism:
        determinism = []
        for label, r in json.load(open(a.determinism)).items():
            tool, version = label.split(' ', 1)
            determinism.append(dict(tool=tool_key(tool, version), reps=r['n_reps'], loci=r['n_loci'],
                                    record_identical=r['record_identical'], loci_length_differs=r['n_loci_length_differs'],
                                    max_len_diff=r['max_len_diff'], lengths_assessed=not tool.startswith('vamos')))

    report = dict(
        schema_version=1, report_id=a.report_id, generated=a.report_id,
        benchmark=dict(repo='https://github.com/wdecoster/TR-Benchmarking', commit=a.commit,
                       catalog=dict(id='benchmark-catalog-v2', loci=43009,
                                    description='the benchmark paper\'s catalog: isolated loci sampled from TRExplorer (up to 1,000 per motif size) plus STRchive disease loci; not genome-wide'),
                       reference='GRCh38 no-alt',
                       input=dict(sample='HG002', data='GIAB EPI2ME 2025.01, sup basecalls, about 30x', truth='T2T HG002 v1.1 assembly'),
                       classes=dict(str='motif 1-6 bp', long='motif > 6 bp')),
        tools=sorted(tools.values(), key=lambda t: (t['tool'], t['version'])),
        metric_defs=[dict(id=k, label=v[0], description=v[1], higher_is_better=v[2]) for k, v in METRICS.items()],
        metrics=metrics, timing=timing, determinism=determinism, notes=a.notes)
    with open(a.out, 'w') as fh:
        json.dump(report, fh, indent=1)
    print(f'{len(tools)} tool versions, {len(metrics)} metric rows -> {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
