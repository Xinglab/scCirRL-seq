"""
Step 4 of the split-parallel-merge barcode calling workflow.

Merge per-chunk bc_umi TSV files (and optionally tagged BAM chunks) produced by
scCirRL_assign_bc_from_tsv into the final bc_umi.tsv and bc_umi.sorted.bam.

TSV merge
---------
Simple header-deduplicating concatenation.  Each input TSV covers a distinct
genomic region, so there is no cross-chunk UMI duplication to resolve.

BAM merge
---------
pysam.merge is used to merge position-sorted BAM chunks into one sorted BAM,
which is then indexed.  Input chunks must each be sorted by position (which is
guaranteed when each chunk was produced from a position-sorted input BAM).
"""
import argparse
import gzip
import sys
import os

import pysam as ps

BC_UMI_HEADER = '#CellBarcode\tUMI\tReadCount\tTranscriptID\tGeneID\tGeneName\tReadNames\tTransCate\tSpliceTag\n'


def _open(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def merge_bc_umi_tsvs(in_fns, out_fn):
    """Concatenate *in_fns* into *out_fn*, emitting the header exactly once.

    Returns (n_chunks_merged, n_lines_written).
    """
    n_lines = 0
    with _open(out_fn, 'w') as out_fp:
        out_fp.write(BC_UMI_HEADER)
        for fn in in_fns:
            with _open(fn, 'r') as in_fp:
                for line in in_fp:
                    if line.startswith('#'):
                        continue
                    out_fp.write(line)
                    n_lines += 1
    return len(in_fns), n_lines


def merge_bam_chunks(in_bam_fns, out_bam_fn):
    """Merge position-sorted BAM chunks into one sorted, indexed BAM.

    Uses pysam.merge (samtools merge) which preserves sort order when all
    inputs are coordinate-sorted.  The output is indexed with pysam.index.
    """
    # pysam.merge("-f", out, in1, in2, ...) overwrites out if it exists
    ps.merge('-f', out_bam_fn, *in_bam_fns)
    ps.index(out_bam_fn)


def parser_argv():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='scCirRL merge-bc-umi: merge per-chunk bc_umi TSVs (and BAMs) '
                    'into final outputs (step 4 of the split-parallel-merge workflow)')

    parser.add_argument('out_bc_umi', metavar='bc_umi.tsv',
                        help='Merged output bc_umi TSV')

    tsv_grp = parser.add_argument_group('TSV inputs')
    tsv_grp.add_argument('--tsv-chunks', metavar='chunk.bc_umi.tsv',
                         nargs='+', required=True,
                         help='Per-chunk bc_umi TSV files from scCirRL_assign_bc_from_tsv')

    bam_grp = parser.add_argument_group('BAM inputs/output (optional)')
    bam_grp.add_argument('--bam-chunks', metavar='chunk.bc_umi.bam',
                         nargs='+', default=[],
                         help='Per-chunk tagged BAM files from scCirRL_assign_bc_from_tsv')
    bam_grp.add_argument('--out-bam', metavar='bc_umi.sorted.bam', type=str, default='',
                         help='Merged, sorted, indexed output BAM '
                              '(required when --bam-chunks is provided)')

    return parser.parse_args()


def main():
    args = parser_argv()

    if args.bam_chunks and not args.out_bam:
        sys.stderr.write('ERROR: --bam-chunks requires --out-bam\n')
        sys.exit(1)

    # ── merge TSVs ────────────────────────────────────────────────────────────
    n_chunks, n_lines = merge_bc_umi_tsvs(args.tsv_chunks, args.out_bc_umi)
    sys.stderr.write('Merged {} TSV chunks ({} UMI records) -> {}\n'.format(
        n_chunks, n_lines, args.out_bc_umi))

    # ── merge BAMs ────────────────────────────────────────────────────────────
    if args.bam_chunks:
        sys.stderr.write('Merging {} BAM chunks -> {} ...\n'.format(
            len(args.bam_chunks), args.out_bam))
        merge_bam_chunks(args.bam_chunks, args.out_bam)
        sys.stderr.write('Merged BAM written and indexed: {}\n'.format(args.out_bam))


if __name__ == '__main__':
    main()
