# Evaluation

Length-only metrics on the tool outputs, one code path for every tool (including inquiSTR and straglr, which report
no sequences). Run the tests with `python3 evaluate/tests/test_evaluate.py`.

| File | What |
|---|---|
| `calls.py` | reads any tool's VCF into allele lengths per catalog locus; where each length comes from is documented per tool at the top of the file |
| `metrics.py` | call rate and concordance with a truth, per motif-length class (STR 1-6 bp, VNTR > 6 bp); outcome per allele as in the paper (MATCH, ONE-OFF, MOTIF-OFF, MISMATCH) |
| `evaluate.py` | all tool outputs of a results folder → table / JSON; with `--truth` also concordance |
| `replicate_agreement.py` | agreement of the replicate runs of `timing.nf` (`timing/calls/`) |

| `truth_from_assembly.py` | allele lengths of the assembly haplotypes at every catalog locus (the truth for the concordance); made on the server by `prepare/prepare_truth.nf` |

Not yet here: the export of the dashboard JSON. Sequence-based metrics (edit distance) remain in
`scripts/assembly_concordance.py`; they do not apply to length-only tools.
Judgement calls in reading the outputs (for example vamos, QUICKREF records, partial calls) are listed in
`docs/ADAPTATIONS.md` as they are decided.
