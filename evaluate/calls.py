"""Read the VCF of any benchmarked tool into one common, length-only form.

A *call* is the pair (or single value, on haploid chrX/chrY of XY samples) of allele lengths of one catalog locus:
the length in bp of the repeat allele over the catalog interval, which is the number that every tool can give and
that is comparable between tools. Sequences are optional (`Call.seqs`, used by the sequence-level metrics in metrics.py): the tools that
write REF/ALT sequences (SEQ_TOOLS) give the sequence of each called allele over the catalog interval; inquiSTR,
mahi-mahi and straglr report lengths only and have no sequences.

Per tool, the length comes from the most direct field of its VCF, not from REF/ALT, wherever the tool has one:
  strdust, inquistr   FORMAT/FRB  (full allele length)
  atarva              FORMAT/AL   (allele length)
  straglr             INFO/RB     (allele sizes; GT indexes [reference, RB...])
  mahimahi            FORMAT/BPCN (allele length in bp per haplotype); filtered records (FILTER not PASS) are
                      kept and marked in `flags`
For the others the length is the catalog length plus len(ALT) - len(REF) of the called allele: this does not depend
on how much flanking sequence a tool puts in REF (strkit and vamos shift the anchor by one base):
  longtr, medaka, strkit
  vamos (the `.fixed.vcf` with sequence alleles) uses len(ALT), see read_calls

Records without a usable genotype are returned as `None` (a record exists, the tool made no call). Loci without a
record are simply absent from the result.
"""
import gzip
from typing import Dict, List, NamedTuple, Optional

SEQ_TOOLS = ('atarva', 'longtr', 'medaka', 'strdust', 'strkit', 'vamos')     # tools whose VCF carries the allele sequences
TOOLS = ('atarva', 'inquistr', 'inquistr_unphased', 'longtr', 'mahimahi', 'medaka', 'straglr', 'strdust', 'strkit', 'vamos')


class Locus(NamedTuple):
    chrom: str
    start: int      # 0-based, as in the catalog
    end: int
    motif: str

    @property
    def id(self) -> str:
        return f'{self.chrom}-{self.start}-{self.end}'

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def motif_class(self) -> str:
        """'STR' for motifs of 1-6 bp, 'VNTR' for longer ones (the split the paper's catalog design asks for)."""
        return 'STR' if 1 <= len(self.motif) <= 6 else 'VNTR'


# Loci a tool cannot take at all. They are reported as "not supported", not as failed calls or zeros.
def not_supported(tool: str, catalog: Dict[str, Locus]) -> List[str]:
    if tool == 'mahimahi':                       # refuses single-base motifs (see docs/ADAPTATIONS.md, A-mahimahi-3)
        return [k for k, l in catalog.items() if len(l.motif) == 1]
    return []


class Call(NamedTuple):
    lengths: Optional[List[int]]    # one value per haplotype; None = the tool made no call at this locus
    flags: str = ''                 # tool-specific remarks (for example STRdust's CLUSTERFAILURE), not interpreted
    seqs: Optional[List[str]] = None    # allele sequences over the catalog interval, one per length; None for length-only tools


def read_catalog(path: str) -> Dict[str, Locus]:
    """chrom, start (0-based), end, motif, ... as in benchmark-catalog-v2.strdust.bed."""
    loci = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            locus = Locus(f[0], int(f[1]), int(f[2]), f[3])
            loci[locus.id] = locus
    return loci


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _info(field: str) -> Dict[str, str]:
    out = {}
    for item in field.split(';'):
        k, _, v = item.partition('=')
        out[k] = v
    return out


def _gt(gt: str) -> Optional[List[Optional[int]]]:
    """Allele indices of a GT string; None for a record without genotype at all."""
    if gt in ('.', './.', '.|.', ''):
        return None
    return [None if a == '.' else int(a) for a in gt.replace('|', '/').split('/')]


def _ints(text: str) -> List[Optional[int]]:
    """Integers of a comma-separated field; anything else ('NaN', '.') is None."""
    out = []
    for x in text.split(','):
        try:
            out.append(int(x))
        except ValueError:
            out.append(None)
    return out


def _key(tool: str, rec, info: Dict[str, str]) -> str:
    """Catalog locus ID (chrom-start0-end) of a record."""
    chrom, pos, id_ = rec[0], int(rec[1]), rec[2]
    if tool in ('strdust', 'longtr') and id_.count('-') == 2:
        return '-'.join(id_.split('-')[:3])
    if tool == 'mahimahi':
        return info['REPID']                        # the catalog's own locus ID
    if tool == 'medaka':
        return '_'.join(id_.split('_')[:3]).replace('_', '-')
    if tool == 'atarva':
        return f'{chrom}-{info["START"]}-{info["END"]}'
    if tool == 'strkit':
        return f'{chrom}-{info["BED_START"]}-{info["BED_END"]}'
    if tool == 'straglr':
        return f'{chrom}-{pos}-{info["END"]}'       # straglr writes the 0-based start as POS
    return f'{chrom}-{pos - 1}-{info["END"]}'       # inquistr, vamos: 1-based POS, END as in the catalog


def _from_alleles(rec, gt, cat_len, alt_is_full=False) -> List[Optional[int]]:
    """Catalog length + (len(ALT) - len(REF)) of each called allele; other symbolic or missing alleles give None; <DEL> and * give 0.
    With alt_is_full the ALT sequence length is the allele length (vamos, see read_calls)."""
    ref, alts = rec[3], rec[4].split(',') if rec[4] not in ('.', '') else []
    lengths = []
    for a in gt:
        if a is None:
            lengths.append(None)
        elif a == 0:
            lengths.append(cat_len)
        elif a <= len(alts) and alts[a - 1] in ('<DEL>', '*'):
            lengths.append(0)                       # deletion of the whole repeat array: an empty allele, as in the paper
        elif a <= len(alts) and not alts[a - 1].startswith('<') and alts[a - 1] != '.':
            lengths.append(len(alts[a - 1]) if alt_is_full else cat_len + len(alts[a - 1]) - len(ref))
        else:
            lengths.append(None)
    return lengths


def _allele_seq(rec, a, locus, full=False) -> Optional[str]:
    """Sequence of allele `a` (0 = REF) over the catalog interval, or None if the record has no sequence for it.

    The record's REF may reach beyond the catalog interval (flanking bases: LongTR, STRkit, ATaRVa and medaka write
    different amounts). The same number of bases is cut from both ends of the ALT sequence, so the allele covers
    the catalog interval exactly; this assumes the tool did not change its flanks. A REF that does not span the interval
    (vamos, which is one base shorter) cannot be trimmed: with `full` the ALT sequence is taken as the allele."""
    ref, alts = rec[3], rec[4].split(',') if rec[4] not in ('.', '') else []
    if a == 0:
        s = ref
    elif a <= len(alts) and alts[a - 1] in ('<DEL>', '*'):
        return ''
    elif a <= len(alts) and not alts[a - 1].startswith('<') and alts[a - 1] != '.':
        s = alts[a - 1]
    else:
        return None
    if full:
        return s.upper()
    pos = int(rec[1])
    left, right = locus.start + 1 - pos, pos + len(ref) - 1 - locus.end
    if left < 0 or right < 0 or len(s) < left + right:
        return s.upper()
    return s[left:len(s) - right].upper()


HAPLOID_XY = ('chrX', 'chrY')


def read_calls(tool: str, path: str, catalog: Dict[str, Locus], karyotype: str = 'XY') -> Dict[str, Call]:
    """Calls of one VCF, keyed by catalog locus ID. Records of loci that are not in the catalog are ignored, but
    at least half of the records must be in it (a different key scheme would otherwise read as "no calls").

    On chrX and chrY of an XY sample only the first allele counts, for every tool, as in the paper's script
    (`v.genotypes[0][:ploidy]`): tools that write two alleles there get no second chance. A call is a no-call when
    a value that counts is missing (so inquiSTR's `200,NaN` is a valid haploid call, `NaN,NaN` and `200,NaN` on an
    autosome are not).
    """
    if tool not in TOOLS:
        raise ValueError(f'unknown tool {tool}')
    tool = 'inquistr' if tool == 'inquistr_unphased' else tool      # same output format
    calls: Dict[str, Call] = {}
    n_rec = 0
    with _open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            n_rec += 1
            rec = line.rstrip('\n').split('\t')
            info = _info(rec[7])
            key = _key(tool, rec, info)
            locus = catalog.get(key)
            if locus is None:
                continue
            fmt = dict(zip(rec[8].split(':'), rec[9].split(':')))
            gt = _gt(fmt.get('GT', '.'))
            flags = []
            if 'CLUSTERFAILURE' in rec[7].split(';'):       # QUICKREF is not a flag: it is a shortcut to a 0/0 call
                flags.append('CLUSTERFAILURE')
            if tool == 'mahimahi' and rec[6] not in ('PASS', '.'):
                flags.append(rec[6])
            flags = ','.join(flags)
            if gt is None:
                calls[key] = Call(None, flags)
                continue
            lengths: Optional[List[Optional[int]]]
            field = {'strdust': 'FRB', 'inquistr': 'FRB', 'atarva': 'AL', 'mahimahi': 'BPCN'}.get(tool)
            if field and field in fmt:
                lengths = _ints(fmt[field])
            elif field:
                # no length field: STRdust 0.16 writes QUICKREF records (INFO flag) as `GT:SUP` with GT 0|0, i.e. the
                # reference allele, for loci it did not genotype in detail
                lengths = _from_alleles(rec, gt, locus.length)
            elif tool == 'straglr':
                rb = [x for x in _ints(info['RB']) if x is not None] if 'RB' in info else []
                table = [locus.length] + rb
                idx = gt * 2 if len(gt) == 1 else gt        # straglr writes a single index for a homozygous call
                lengths = [table[i] if i is not None and i < len(table) else None for i in idx]
            else:
                # vamos' REF is always one base shorter than the catalog interval, and its ALT sequences are in the
                # frame of the catalog interval: 96% of its lengths equal inquiSTR's with len(ALT), and none with
                # the REF-relative length. To be confirmed against the assembly truth (a constant +-1 would show).
                lengths = _from_alleles(rec, gt, locus.length, alt_is_full=(tool == 'vamos'))
            if tool == 'mahimahi' and len(lengths) == 1 and len(gt) == 2:
                lengths = lengths * 2                   # mahi-mahi writes one value for a homozygous call
            seqs = None
            if tool in SEQ_TOOLS:
                seqs = [_allele_seq(rec, a, locus, full=(tool == 'vamos' and a != 0)) for a in gt]
                if any(x is None for x in seqs):
                    seqs = None
            if karyotype == 'XY' and locus.chrom in HAPLOID_XY:
                lengths = lengths[:1]
                seqs = seqs[:1] if seqs else seqs
            if not lengths or any(x is None for x in lengths):
                calls[key] = Call(None, flags)
            else:
                calls[key] = Call(list(lengths), flags, seqs)
    if n_rec and len(calls) < n_rec / 2:
        raise ValueError(f'{path}: only {len(calls)} of {n_rec} records are catalog loci for {tool}: wrong catalog or key scheme?')
    return calls
