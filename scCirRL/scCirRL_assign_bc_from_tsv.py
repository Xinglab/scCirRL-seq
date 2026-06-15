"""
Step 3 of the split-parallel-merge barcode calling workflow.

Assign reference barcodes and cluster UMIs for one reads TSV chunk produced by
scCirRL_extract_reads, using the reference barcode list from scCirRL_collect_ref_bc.
Run one instance per chunk in parallel; outputs are merged by scCirRL_merge_bc_umi.

Outputs
-------
out_bc_umi.tsv   – same column layout as the main bc_umi.tsv
out_chunk.bam    – (optional) tagged BAM chunk with CB/UB/BE/TC/TI/GI/GN tags

The BAM pass re-reads the original BAM for the same region used during extraction
so that only assigned reads (those present in umi_clu_res_dict) are emitted and
every written record carries the full set of scCirRL tags.
"""
import argparse
import gzip
import sys
import os
from collections import defaultdict as dd

import pysam as ps

from .assign_bc import get_uniq_match_ref_bc
from .cluster_umi import umi_clustering
from .parse_quant_file import collect_trans_to_gene, collect_read_to_trans, collect_trans_to_chrom_span
from .seq_utils import check_for_skip
from .parameter import scrl_para_class, _bc_len, _umi_len, _five_ada, _bc_max_ed, _umi_max_ed
from .utils import err_log_format_time
from .__init__ import __version__, __program__, __cmd__


BC_UMI_HEADER = '#CellBarcode\tUMI\tReadCount\tTranscriptID\tGeneID\tGeneName\tReadNames\tTransCate\tSpliceTag\n'


def _open(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def _load_ref_bcs(ref_bc_fn):
    ref_bcs = []
    with open(ref_bc_fn) as fp:
        for line in fp:
            bc = line.strip()
            if bc:
                ref_bcs.append(bc)
    return ref_bcs


def _parse_region(region):
    """Return (chrom,) or (chrom, start, end) suitable for pysam fetch()."""
    if not region:
        return ()
    if ':' in region:
        chrom, span = region.split(':', 1)
        start, end = span.split('-')
        return (chrom, int(start), int(end))
    return (region,)


def _region_filter(region):
    """Parse region string into (chrom, start, end) or (chrom, None, None).

    Returns (None, None, None) when no region is specified.
    """
    if not region:
        return None, None, None
    if ':' in region:
        chrom, span = region.split(':', 1)
        start, end = span.split('-')
        return chrom, int(start), int(end)
    return region, None, None


def _read_in_region(read_chrom, read_start, read_end, reg_chrom, reg_start, reg_end):
    """True if the read overlaps [reg_start, reg_end) on reg_chrom."""
    if read_chrom != reg_chrom:
        return False
    if reg_start is None:  # whole-chromosome region
        return True
    return int(read_start) < reg_end and int(read_end) > reg_start


def _search_ref_bc_from_row(qname, bc, umi, bu, is_perfect,
                             ref_bcs, cand_ref_bc_seq,
                             bc_len, bc_max_ed, umi_len):
    """BC assignment from pre-extracted (bc, umi, bu) strings."""
    if not bc or not umi:
        return '', '', -1, bool(is_perfect), False

    is_exact = False
    if bc in ref_bcs:
        is_exact = True
        return bc, umi, 0, bool(is_perfect), is_exact

    ref_bc, ref_umi, bc_ed = get_uniq_match_ref_bc(
        qname, bc, bu, ref_bcs, cand_ref_bc_seq, bc_len, bc_max_ed, umi_len)
    if ref_bc and ref_umi:
        return ref_bc, ref_umi, bc_ed, bool(is_perfect), is_exact

    return '', '', -1, bool(is_perfect), is_exact


def _write_tagged_bam(in_bam_fn, out_bam_fn, region,
                      umi_clu_res_dict, bc_eds, read_to_cate,
                      only_primary, skip_chimeric, log_fn):
    """Second pass: re-read the original BAM and emit assigned reads with tags.

    Only reads present in *umi_clu_res_dict* are written.  The same
    check_for_skip filter used during extraction is applied so that
    unmapped / secondary / supplementary reads are excluded.
    """
    para = scrl_para_class()
    para.only_primary  = only_primary
    para.skip_chimeric = skip_chimeric
    para.log_fn        = log_fn

    fetch_args = _parse_region(region)

    err_log_format_time(log_fn, str='Writing tagged BAM chunk: {}'.format(out_bam_fn))
    with ps.AlignmentFile(in_bam_fn) as in_bam:
        hdr = in_bam.header.to_dict()
        pg_entry = {'ID': __program__, 'PN': __program__,
                    'VN': __version__, 'CL': __cmd__}
        if 'PG' in hdr:
            existing_pns = [pg.get('PN') for pg in hdr['PG']]
            if existing_pns.count(__program__) > 0:
                pg_entry['ID'] = '{}.{}'.format(
                    __program__, existing_pns.count(__program__) + 1)
            hdr['PG'].append(pg_entry)
        else:
            hdr['PG'] = [pg_entry]

        with ps.AlignmentFile(out_bam_fn, 'wb', header=hdr) as out_bam:
            fetch_iter = in_bam.fetch(*fetch_args) if fetch_args else in_bam.fetch()
            for r in fetch_iter:
                skip, _ = check_for_skip(r, para)
                if skip:
                    continue
                qname = r.query_name
                if qname not in umi_clu_res_dict:
                    continue
                bc, umi, cmpt_trans, cmpt_gene_id, cmpt_gene_names = umi_clu_res_dict[qname]
                r.set_tag('CB', bc,                      'Z')
                r.set_tag('UB', umi,                     'Z')
                r.set_tag('BE', bc_eds[qname],           'i')
                r.set_tag('TC', read_to_cate[qname],     'Z')
                r.set_tag('TI', ','.join(cmpt_trans),    'Z')
                r.set_tag('GI', ','.join(cmpt_gene_id),  'Z')
                r.set_tag('GN', ','.join(cmpt_gene_names), 'Z')
                out_bam.write(r)

    err_log_format_time(log_fn, str='Indexing {}'.format(out_bam_fn))
    ps.index(out_bam_fn)


def assign_bc_from_tsv(tsv_fn, ref_bcs, cand_ref_bc_seq,
                        trans_to_gene_id_name, read_to_trans, read_to_cate,
                        out_bc_umi_fn,
                        bc_len, bc_max_ed, umi_len, umi_max_ed,
                        in_bam_fn=None, out_bam_fn=None, region=None,
                        only_primary=True, skip_chimeric=False,
                        log_fn=''):
    """Assign barcodes and cluster UMIs for reads in *tsv_fn*.

    Pass 1 (TSV): assign BCs, collect per-read tag data, cluster UMIs,
    write bc_umi TSV.

    Pass 2 (BAM, optional): re-read *in_bam_fn* for *region* and emit
    assigned reads with CB/UB/BE/TC/TI/GI/GN tags into *out_bam_fn*.

    UMI clustering is done per barcode across the whole chunk.  Reads from
    different genomic loci sharing the same BC+UMI map to different
    transcripts, so the UMI network graph separates them correctly without
    needing per-region sub-clustering (see cluster_umi.get_umi_trans).
    """
    bu_res         = []   # [qname, ref_bc, umi, is_perfect]
    bc_eds         = {}   # {qname: edit_distance}
    read_to_splice = {}   # {qname: 'Y'|'N'}

    n_perfect_in_ref, n_perfect_uniq_to_ref = 0, 0
    n_imperfect_in_ref, n_imperfect_uniq_to_ref = 0, 0
    n_total = 0

    reg_chrom, reg_start, reg_end = _region_filter(region)

    # ── Pass 1: TSV ──────────────────────────────────────────────────────────
    with _open(tsv_fn, 'r') as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rstrip('\n').split('\t')
            if len(ele) < 9:
                continue
            qname, chrom, start, end, bc, umi, bu, is_perfect_str, ts_tag = ele[:9]
            if reg_chrom and not _read_in_region(chrom, start, end, reg_chrom, reg_start, reg_end):
                continue
            is_perfect = is_perfect_str == '1'
            n_total += 1

            ref_bc, ref_umi, bc_ed, _is_perfect, is_exact = _search_ref_bc_from_row(
                qname, bc, umi, bu, is_perfect,
                ref_bcs, cand_ref_bc_seq, bc_len, bc_max_ed, umi_len)

            if ref_bc and ref_umi:
                if _is_perfect:
                    if is_exact: n_perfect_in_ref     += 1
                    else:        n_perfect_uniq_to_ref += 1
                else:
                    if is_exact: n_imperfect_in_ref     += 1
                    else:        n_imperfect_uniq_to_ref += 1
                bu_res.append([qname, ref_bc, ref_umi, _is_perfect])
                bc_eds[qname]         = bc_ed
                read_to_splice[qname] = 'Y' if ts_tag != '.' else 'N'

    umi_clu_res_dict, umi_clu_res_list = umi_clustering(
        log_fn, read_to_trans, trans_to_gene_id_name, bu_res, umi_max_ed)

    n_assigned = (n_perfect_in_ref + n_perfect_uniq_to_ref
                  + n_imperfect_in_ref + n_imperfect_uniq_to_ref)
    err_log_format_time(log_fn,
        str='[{}] total={}, assigned={} ({:.1f}%)'.format(
            os.path.basename(tsv_fn), n_total, n_assigned,
            100.0 * n_assigned / n_total if n_total else 0.0))

    with open(out_bc_umi_fn, 'w') as out_fp:
        out_fp.write(BC_UMI_HEADER)
        for bc, umi, reads, cmpt_trans, cmpt_gene_id, cmpt_gene_names in umi_clu_res_list:
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                bc, umi, len(reads),
                ','.join(cmpt_trans),
                ','.join(cmpt_gene_id),
                ','.join(cmpt_gene_names),
                ','.join(reads),
                ','.join([read_to_cate[r] for r in reads]),
                ','.join([read_to_splice.get(r, 'N') for r in reads])))

    # ── Pass 2: BAM (optional) ───────────────────────────────────────────────
    if in_bam_fn and out_bam_fn:
        _write_tagged_bam(
            in_bam_fn, out_bam_fn, region,
            umi_clu_res_dict, bc_eds, read_to_cate,
            only_primary, skip_chimeric, log_fn)

    return n_assigned, n_total


def _make_scrl_para(args):
    para = scrl_para_class()
    para.bc_len      = args.bc_len
    para.umi_len     = args.umi_len
    para.bc_max_ed   = args.bc_ed
    para.umi_max_ed  = args.umi_ed
    para.five_ada    = args.five_ada
    para.five_max_ed = int(len(args.five_ada) * 0.3)
    para.anno_gtf    = args.anno_gtf
    para.updated_gtf = args.updated_gtf
    para.cmpt_tsv    = args.cmpt_tsv
    para.log_fn      = ''
    return para


def parser_argv():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='scCirRL assign-bc-from-tsv: reads TSV + ref barcodes -> '
                    'bc_umi TSV (and optionally tagged BAM) chunk '
                    '(step 3 of the split-parallel-merge workflow)')

    parser.add_argument('reads_tsv',  metavar='reads.tsv[.gz]',
                        help='Reads TSV from scCirRL_extract_reads')
    parser.add_argument('ref_bc_tsv', metavar='ref_bc.tsv',
                        help='Reference barcode list from scCirRL_collect_ref_bc')
    parser.add_argument('out_bc_umi', metavar='out_bc_umi.tsv',
                        help='Output bc_umi TSV chunk')

    parser.add_argument('--region', type=str, default='',
                        help='Genomic region to process, e.g. chr1 or chr1:10000-50000. '
                             'Filters reads in the TSV by coordinate and restricts the '
                             'BAM pass to the same region. Default: process all reads.')

    bam_grp = parser.add_argument_group('BAM output (optional)')
    bam_grp.add_argument('--in-bam',  type=str, default='',
                         help='Original sorted BAM (required for BAM output)')
    bam_grp.add_argument('--out-bam', type=str, default='',
                         help='Output tagged BAM chunk path (e.g. chunk_chr1.bc_umi.bam)')

    input_grp = parser.add_argument_group('input files (optional)')
    input_grp.add_argument('-g', '--anno-gtf',    type=str, default='')
    input_grp.add_argument('-t', '--updated-gtf', type=str, default='')
    input_grp.add_argument('-m', '--cmpt-tsv',    type=str, default='')

    bc_grp = parser.add_argument_group('barcode and UMI options')
    bc_grp.add_argument('-b', '--bc-len',       type=int,  default=_bc_len)
    bc_grp.add_argument('-e', '--bc-ed',        type=int,  default=_bc_max_ed)
    bc_grp.add_argument('-u', '--umi-len',      type=int,  default=_umi_len)
    bc_grp.add_argument('-d', '--umi-ed',       type=int,  default=_umi_max_ed)
    bc_grp.add_argument('-5', '--five-ada',     type=str,  default=_five_ada)
    bc_grp.add_argument('-s', '--skip-chimeric', action='store_true', default=False)

    return parser.parse_args()


def main():
    args = parser_argv()

    if args.out_bam and not args.in_bam:
        sys.stderr.write('ERROR: --out-bam requires --in-bam\n')
        sys.exit(1)

    para = _make_scrl_para(args)

    ref_bcs = _load_ref_bcs(args.ref_bc_tsv)
    if not ref_bcs:
        sys.stderr.write('ERROR: no reference barcodes in {}\n'.format(args.ref_bc_tsv))
        sys.exit(1)
    cand_ref_bc_seq = ('N' * (para.bc_max_ed + 1)).join(ref_bcs)

    trans_to_gene_id_name = collect_trans_to_gene(para)
    region = args.region or None
    trans_to_span = collect_trans_to_chrom_span(para.updated_gtf) if region else None
    read_to_trans, read_to_cate = collect_read_to_trans(
        trans_to_gene_id_name, para, region=region, trans_to_span=trans_to_span)

    assign_bc_from_tsv(
        tsv_fn          = args.reads_tsv,
        ref_bcs         = ref_bcs,
        cand_ref_bc_seq = cand_ref_bc_seq,
        trans_to_gene_id_name = trans_to_gene_id_name,
        read_to_trans   = read_to_trans,
        read_to_cate    = read_to_cate,
        out_bc_umi_fn   = args.out_bc_umi,
        bc_len          = para.bc_len,
        bc_max_ed       = para.bc_max_ed,
        umi_len         = para.umi_len,
        umi_max_ed      = para.umi_max_ed,
        in_bam_fn       = args.in_bam   or None,
        out_bam_fn      = args.out_bam  or None,
        region          = region,
        only_primary    = True,
        skip_chimeric   = args.skip_chimeric,
        log_fn          = '')


if __name__ == '__main__':
    main()
