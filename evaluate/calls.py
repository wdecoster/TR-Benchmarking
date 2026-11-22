"""Read the VCF of any benchmarked tool into one common, length-only form.

A *call* is the pair (or single value, on haploid chrX/chrY of XY samples) of allele lengths of one catalog locus:
the length in bp of the repeat allele over the catalog interval, which is the number that every tool can give and
that is comparable between tools. Sequences are not used here (inquiSTR and straglr report lengths only); the
sequence-based metrics of the paper stay with `scripts/assembly_concordance.py`.

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

TOOLS = ('atarva', 'inquistr', 'longtr', 'mahimahi', 'medaka', 'straglr', 'strdust', 'strkit', 'vamos')


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


def _ints(text: str) -> Optional[List[int]]:
    try:
        return [int(x) for x in text.split(',')]
    except ValueError:
        return None


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
    """Catalog length + (len(ALT) - len(REF)) of each called allele; symbolic or missing alleles give None.
    With alt_is_full the ALT sequence length is the allele length (vamos, see read_calls)."""
    ref, alts = rec[3], rec[4].split(',') if rec[4] not in ('.', '') else []
    lengths = []
    for a in gt:
        if a is None:
            lengths.append(None)
        elif a == 0:
            lengths.append(cat_len)
        elif a <= len(alts) and not alts[a - 1].startswith('<') and alts[a - 1] not in ('*', '.'):
            lengths.append(len(alts[a - 1]) if alt_is_full else cat_len + len(alts[a - 1]) - len(ref))
        else:
            lengths.append(None)
    return lengths


def read_calls(tool: str, path: str, catalog: Dict[str, Locus]) -> Dict[str, Call]:
    """Calls of one VCF, keyed by catalog locus ID. Records of loci that are not in the catalog are ignored."""
    if tool not in TOOLS:
        raise ValueError(f'unknown tool {tool}')
    calls: Dict[str, Call] = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            rec = line.rstrip('\n').split('\t')
            info = _info(rec[7])
            key = _key(tool, rec, info)
            locus = catalog.get(key)
            if locus is None:
                continue
            fmt = dict(zip(rec[8].split(':'), rec[9].split(':')))
            gt = _gt(fmt.get('GT', '.'))
            flags = ''
            if 'CLUSTERFAILURE' in rec[7]:
                flags = 'CLUSTERFAILURE'
            if gt is None:
                calls[key] = Call(None, flags)
                continue
            lengths: Optional[List[Optional[int]]]
            field = {'strdust': 'FRB', 'inquistr': 'FRB', 'atarva': 'AL', 'mahimahi': 'BPCN'}.get(tool)
            if tool == 'mahimahi' and rec[6] not in ('PASS', '.'):
                flags = rec[6]
            if field and field in fmt:
                lengths = _ints(fmt[field])
            elif field:
                # no length field: STRdust 0.16 writes QUICKREF records (INFO flag) as `GT:SUP` with GT 0|0, i.e. the
                # reference allele, for loci it did not genotype in detail
                lengths = _from_alleles(rec, gt, locus.length)
            elif tool == 'straglr':
                rb = _ints(info['RB']) if 'RB' in info else []
                table = [locus.length] + (rb or [])
                idx = gt * 2 if len(gt) == 1 else gt        # straglr writes a single index for a homozygous call
                lengths = [table[i] if i is not None and i < len(table) else None for i in idx]
            else:
                # vamos' REF is always one base shorter than the catalog interval, and its ALT sequences are in the
                # frame of the catalog interval: 96% of its lengths equal inquiSTR's with len(ALT), and none with
                # the REF-relative length. To be confirmed against the assembly truth (a constant +-1 would show).
                lengths = _from_alleles(rec, gt, locus.length, alt_is_full=(tool == 'vamos'))
            if lengths is not None and len(lengths) == 1 and len(gt) == 2:
                lengths = lengths * 2                   # mahi-mahi writes one value for a homozygous call
            if lengths is None or any(x is None for x in lengths):
                calls[key] = Call(None, flags)
            else:
                calls[key] = Call(list(lengths), flags)
    return calls
