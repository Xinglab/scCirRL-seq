"""
Step 2 of the split-parallel-merge barcode calling workflow.

Collect candidate reference cell barcodes from one or more reads TSVs produced
by scCirRL_extract_reads.  Mirrors collect_candidate_bc.sp_collect_cand_ref_bc
but operates on TSVs instead of a BAM file, so all chunks can be processed in a
single fast pass without re-opening the BAM.

Outputs:
  ref_bc.tsv           – one barcode per line (the reference barcode list)
  high_qual_bc.rank.tsv – BC rank table (BC / Rank / UMICount)
"""
import argparse
import gzip
import sys
import os
from collections import defaultdict as dd

from .collect_candidate_bc import (
    get_cumulative_cnt,
    get_putative_knee_point,
    write_perfect_bu_cnt,
    _max_cand_bc_for_knee_calling,
)
from .utils import err_log_format_time
from .parameter import _bc_len, _umi_len


def _open(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def collect_ref_bc_from_tsvs(tsv_fns, bc_len, umi_len,
                              expect_cell_count, high_qual_bu_fn, log_fn=''):
    """Scan *tsv_fns* for perfect-match reads; detect knee point; return ref BCs.

    Parameters
    ----------
    tsv_fns : list[str]
        Paths to reads TSVs produced by scCirRL_extract_reads.
    bc_len, umi_len : int
    expect_cell_count : int
        If > 0 the top-N barcodes are returned directly, bypassing knee calling.
    high_qual_bu_fn : str
        Path for the BC-rank output table.
    log_fn : str
        Path for log messages (empty = stderr only).

    Returns
    -------
    list[str]  – reference barcodes
    """
    err_log_format_time(log_fn, str='Collecting candidate cell barcodes from {} TSV file(s) ...'.format(len(tsv_fns)))

    perfect_bc_umi_to_read = dd(lambda: dd(lambda: []))  # {bc: {umi: [qnames]}}

    for tsv_fn in tsv_fns:
        with _open(tsv_fn, 'r') as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                ele = line.rstrip('\n').split('\t')
                if len(ele) < 9:
                    continue
                qname, chrom, start, end, bc, umi, bu, is_perfect, ts_tag = ele[:9]
                if is_perfect != '1':
                    continue
                if len(bc) != bc_len or len(umi) != umi_len:
                    continue
                perfect_bc_umi_to_read[bc][umi].append(qname)

    n_perfect_reads = sum(
        len(reads)
        for umi_reads in perfect_bc_umi_to_read.values()
        for reads in umi_reads.values())

    perfect_bc_to_umi_cnt = {bc: len(umis) for bc, umis in perfect_bc_umi_to_read.items()}
    perfect_bc_to_umi_cnt = dict(sorted(perfect_bc_to_umi_cnt.items(), key=lambda d: -d[1]))

    write_perfect_bu_cnt(perfect_bc_to_umi_cnt, high_qual_bu_fn)

    if expect_cell_count > 0:
        ref_bcs = list(perfect_bc_to_umi_cnt.keys())[:expect_cell_count]
        err_log_format_time(log_fn, str='Top {} (-c {}) barcodes picked.'.format(len(ref_bcs), expect_cell_count))
    else:
        cumulative = get_cumulative_cnt(
            perfect_bc_to_umi_cnt.values(),
            max_n=_max_cand_bc_for_knee_calling,
            min_cnt=2)
        knee_rank = get_putative_knee_point(log_fn, cumulative)
        err_log_format_time(log_fn,
            str='Knee point barcode rank: {} (from {} barcodes, {} high-quality reads)'.format(
                knee_rank, len(cumulative), n_perfect_reads))
        if knee_rank is None:
            return []
        ref_bcs = list(perfect_bc_to_umi_cnt.keys())[:knee_rank]

    err_log_format_time(log_fn,
        str='Collecting candidate cell barcodes done! ({})'.format(len(ref_bcs)))
    return ref_bcs


def parser_argv():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='scCirRL collect-ref-bc: reads TSVs -> reference barcode list '
                    '(step 2 of the split-parallel-merge workflow)')
    parser.add_argument('tsv_files', metavar='reads.tsv', nargs='+',
                        help='One or more reads TSVs from scCirRL_extract_reads')
    parser.add_argument('out_ref_bc', metavar='ref_bc.tsv',
                        help='Output reference barcode list (one BC per line)')
    parser.add_argument('--high-qual-bc-rank', metavar='high_qual_bc.rank.tsv',
                        default='high_qual_bc.rank.tsv',
                        help='Output BC rank table')
    parser.add_argument('-c', '--cell-count', type=int, default=-1,
                        help='Force top-N cell count; -1 = auto knee-point detection')
    parser.add_argument('-b', '--bc-len',  type=int, default=_bc_len)
    parser.add_argument('-u', '--umi-len', type=int, default=_umi_len)
    parser.add_argument('-l', '--bc-list', type=str, default='',
                        help='Use a pre-existing reference barcode list instead of '
                             'calling barcodes from the TSVs (skips knee point step)')
    return parser.parse_args()


def _load_bc_list(bc_list_fn, bc_len):
    import re
    ref_bcs = []
    open_fn = gzip.open if bc_list_fn.endswith('.gz') else open
    mode = 'rt' if bc_list_fn.endswith('.gz') else 'r'
    with open_fn(bc_list_fn, mode) as fp:
        for line in fp:
            ele = re.split(r'-| |\t', line.rstrip())
            if len(ele[0]) == bc_len and ele[0] not in ref_bcs:
                ref_bcs.append(ele[0])
    return ref_bcs


def main():
    args = parser_argv()

    if args.bc_list:
        ref_bcs = _load_bc_list(args.bc_list, args.bc_len)
        sys.stderr.write('Loaded {} reference barcodes from {}\n'.format(
            len(ref_bcs), args.bc_list))
    else:
        ref_bcs = collect_ref_bc_from_tsvs(
            tsv_fns=args.tsv_files,
            bc_len=args.bc_len,
            umi_len=args.umi_len,
            expect_cell_count=args.cell_count,
            high_qual_bu_fn=args.high_qual_bc_rank,
            log_fn='')

    if not ref_bcs:
        sys.stderr.write('ERROR: no reference barcodes detected.\n')
        sys.exit(1)

    with open(args.out_ref_bc, 'w') as fp:
        for bc in ref_bcs:
            fp.write('{}\n'.format(bc))
    sys.stderr.write('Wrote {} reference barcodes to {}\n'.format(
        len(ref_bcs), args.out_ref_bc))


if __name__ == '__main__':
    main()
