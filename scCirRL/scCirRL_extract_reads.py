"""
Step 1 of the split-parallel-merge barcode calling workflow.

Extract per-read BC/UMI candidate info from a BAM file (or a specific genomic
region) into a lightweight TSV.  Run one instance per BAM chunk in parallel;
the resulting TSVs feed into scCirRL_collect_ref_bc and scCirRL_assign_bc_from_tsv.

Output TSV columns (tab-separated, optionally gzipped):
  qname  chrom  start  end  bc  umi  bu  is_perfect  ts_tag
"""
import argparse
import gzip
import sys
import os
import pysam as ps

from .seq_utils import check_for_skip, get_bu_from_seq
from .parameter import scrl_para_class, _bc_len, _umi_len, _five_ada, _bc_max_ed

READS_TSV_HEADER = '#qname\tchrom\tstart\tend\tbc\tumi\tbu\tis_perfect\tts_tag\n'


def _open(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def _make_para(bc_len, umi_len, five_ada, only_primary, skip_chimeric):
    para = scrl_para_class()
    para.bc_len = bc_len
    para.umi_len = umi_len
    para.five_ada = five_ada
    para.five_max_ed = int(len(five_ada) * 0.3)
    para.only_primary = only_primary
    para.skip_chimeric = skip_chimeric
    para.log_fn = ''
    return para


def bam_to_reads_tsv(in_bam_fn, out_tsv_fn, bc_len, umi_len, five_ada,
                     only_primary=True, skip_chimeric=False, region=None):
    """Extract BC/UMI candidate info from *in_bam_fn* into *out_tsv_fn*.

    *region* is an optional string in samtools format, e.g. ``"chr1"`` or
    ``"chr1:10000-20000"``.  When omitted all mapped reads are processed.

    Returns the number of reads written.
    """
    para = _make_para(bc_len, umi_len, five_ada, only_primary, skip_chimeric)
    five_max_ed = para.five_max_ed
    n_written = 0

    out_dir = os.path.dirname(os.path.abspath(out_tsv_fn))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with ps.AlignmentFile(in_bam_fn) as in_bam, _open(out_tsv_fn, 'w') as out_fp:
        out_fp.write(READS_TSV_HEADER)

        if region:
            if ':' in region:
                chrom, span = region.split(':', 1)
                start, end = span.split('-')
                fetch_iter = in_bam.fetch(chrom, int(start), int(end))
            else:
                fetch_iter = in_bam.fetch(region)
        else:
            fetch_iter = in_bam.fetch()

        n_total = 0
        for r in fetch_iter:
            skip, _ = check_for_skip(r, para)
            if skip:
                continue
            bc, umi, bu, is_perfect = get_bu_from_seq(
                r, five_ada, five_max_ed, bc_len, umi_len)
            ts_tag = r.get_tag('ts') if r.has_tag('ts') else '.'
            out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                r.query_name,
                r.reference_name,
                r.reference_start,
                r.reference_end,
                bc, umi, bu,
                1 if is_perfect else 0,
                ts_tag))
            n_written += 1
            n_total += 1
            if n_total % 1_000_000 == 0:
                sys.stderr.write('  {} M reads processed\r'.format(n_total // 1_000_000))
                sys.stderr.flush()

        sys.stderr.write('  {} reads processed\n'.format(n_total))

    return n_written


def parser_argv():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='scCirRL extract-reads: BAM -> per-read BC/UMI candidate TSV '
                    '(step 1 of the split-parallel-merge workflow)')
    parser.add_argument('in_bam', metavar='in.sorted.bam',
                        help='Sorted, indexed BAM file (whole file or a region)')
    parser.add_argument('out_tsv', metavar='out.reads.tsv[.gz]',
                        help='Output TSV (append .gz for gzip compression)')

    parser.add_argument('-r', '--region', type=str, default='',
                        help='Genomic region to extract, e.g. chr1 or chr1:1000-5000. '
                             'Default: all reads')
    parser.add_argument('-b', '--bc-len',    type=int, default=_bc_len)
    parser.add_argument('-u', '--umi-len',   type=int, default=_umi_len)
    parser.add_argument('-5', '--five-ada',  type=str, default=_five_ada)
    parser.add_argument('-s', '--skip-chimeric', action='store_true', default=False,
                        help='Skip chimeric reads (reads with SA tag)')
    return parser.parse_args()


def main():
    args = parser_argv()
    region = args.region if args.region else None
    n = bam_to_reads_tsv(
        in_bam_fn=args.in_bam,
        out_tsv_fn=args.out_tsv,
        bc_len=args.bc_len,
        umi_len=args.umi_len,
        five_ada=args.five_ada,
        only_primary=True,
        skip_chimeric=args.skip_chimeric,
        region=region)
    sys.stderr.write('Wrote {} reads to {}\n'.format(n, args.out_tsv))


if __name__ == '__main__':
    main()
