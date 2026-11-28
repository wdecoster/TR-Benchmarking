# Evaluation

Length-only metrics on the tool outputs, one code path for every tool (including inquiSTR and straglr, which report
no sequences). Run the tests with `python3 evaluate/tests/test_evaluate.py`.

| File | What |
|---|---|
| `calls.py` | reads any tool's VCF into allele lengths per catalog locus; where each length comes from is documented per tool at the top of the file |
| `metrics.py` | call rate and concordance with a truth, per motif-length class (STR 1-6 bp, VNTR > 6 bp), for two analyses: `all_loci` (the main one, no locus left out; loci a tool cannot take count as not called) and `no_homopolymers` (1 bp motifs left out for every tool); outcome per allele as in the paper (MATCH, ONE-OFF, MOTIF-OFF, MISMATCH) |
| `evaluate.py` | all tool outputs of a results folder → table / JSON; with `--truth` also concordance |
| `replicate_agreement.py` | agreement of the replicate runs of `timing.nf` (`timing/calls/`) |

| `truth_from_assembly.py` | allele lengths of the assembly haplotypes at every catalog locus (the truth for the concordance); made on the server by `prepare/prepare_truth.nf` |

Not yet here: the export of the dashboard JSON. The paper's own sequence concordance is still in
`scripts/assembly_concordance.py` (see the sequence-level section below for ours).
Judgement calls in reading the outputs (for example vamos, QUICKREF records, partial calls) are listed in
`docs/ADAPTATIONS.md` as they are decided.

## Sequence level and the workflow step
- `truth_from_assembly.py --seq-out` also writes the assembly allele sequences over the catalog interval (`HG002.truth.seq.tsv`, made by
  `prepare/prepare_truth.nf`). `evaluate.py --truth-seq` adds, for the tools whose VCF carries sequences (atarva, longtr, medaka, strdust,
  strkit, vamos; `calls.SEQ_TOOLS`), the edit distance of each called allele to its assembly allele: `seq_exact_alleles`,
  `seq_exact_genotype`, `seq_exact_of_truth_loci` (a locus without a call counts as wrong, as for the length metric), `mean_edit_distance`
  and `mean_norm_edit_distance` (edit distance / longer allele). inquiSTR, mahi-mahi and straglr write lengths only and get no sequence
  columns (shown as "not applicable", never as zero).
- How a called allele is cut to the catalog interval, and what that assumes, is documented in `calls._allele_seq`; edit distances use
  `metrics.edit_distance` (Myers bit-parallel, no dependency, tested against a reference implementation).
- **Part of the workflow:** with `--truth` (and `--truth_seq`) `main.nf` evaluates every tool output of the run in a containerised step
  (`evaluate_calls`, one task per tool and sample) and publishes `<outdir>/evaluation/<mode>/<sample>.<tool>.<version>.eval.{json,tsv}`,
  so a new tool version is evaluated by the same run that produces its calls. `evaluate/evaluate.nf` re-evaluates existing results folders
  with the same `evaluate.py` (for runs made before the step existed, or after a change of the metrics).
