"""
Split a read-isoform compatible TSV (from ESPRESSO or Bambu) into per-chromosome
files, using transcript-to-chromosome information from an updated GTF.

A read is assigned to a chromosome if at least one of its compatible transcripts
is annotated on that chromosome.  Reads whose transcripts are all absent from the
GTF are written to a catch-all ``unassigned`` file.

Output files are written to *out_dir* with names:
  <prefix>.<chrom>.tsv   – one file per chromosome
  <prefix>.unassigned.tsv – reads with no GTF match
"""
import argparse
import os
import sys

from .parse_quant_file import collect_trans_to_chrom_span


def _get_trans_ids(ele):
    """Return (qname, list_of_trans_ids, raw_format) from a parsed TSV row.

    Supports Bambu (3 columns) and ESPRESSO (4 columns) formats.
    Returns (None, [], fmt) when the row is a header or unrecognised.
    """
    if len(ele) == 3:  # Bambu: readId, equalMatches, compatibleMatches
        qname, equal_trans, cmpt_trans = ele
        if qname == 'readId' or equal_trans == 'equalMatches' or cmpt_trans == 'compatibleMatches':
            return None, [], 'bambu'
        trans_ids = []
        if equal_trans != 'NA':
            trans_ids.extend(equal_trans.split(','))
        if cmpt_trans != 'NA':
            trans_ids.extend(cmpt_trans.split(','))
        return qname, [t for t in trans_ids if t], 'bambu'

    if len(ele) == 4:  # ESPRESSO: readID, sample/bam, category, transcriptID(s)
        qname, _, _, cmpt_trans = ele
        if cmpt_trans == 'NA':
            return qname, [], 'espresso'
        trans_ids = [t for t in cmpt_trans.split(',') if t]
        return qname, trans_ids, 'espresso'

    return None, [], 'unknown'


def split_cmpt_tsv(cmpt_tsv_fn, updated_gtf, out_dir, prefix):
    """Split *cmpt_tsv_fn* into per-chromosome files in *out_dir*.

    Parameters
    ----------
    cmpt_tsv_fn : str
        Read-isoform compatible TSV (Bambu or ESPRESSO format).
    updated_gtf : str
        Updated GTF file; transcript lines provide chromosome assignments.
    out_dir : str
        Output directory (created if absent).
    prefix : str
        Filename prefix for output files.

    Returns
    -------
    dict[str, int]  – {chrom: read_count} including 'unassigned'
    """
    sys.stderr.write('Loading transcript-to-chromosome map from {} ...\n'.format(updated_gtf))
    trans_to_span = collect_trans_to_chrom_span(updated_gtf)
    if not trans_to_span:
        sys.stderr.write('ERROR: no transcript entries found in {}\n'.format(updated_gtf))
        sys.exit(1)
    sys.stderr.write('  {} transcripts loaded\n'.format(len(trans_to_span)))

    os.makedirs(out_dir, exist_ok=True)

    out_fps = {}     # {chrom: file_object}
    counts  = {}     # {chrom: int}
    header_line = None

    def _get_fp(chrom):
        if chrom not in out_fps:
            path = os.path.join(out_dir, '{}.{}.tsv'.format(prefix, chrom))
            out_fps[chrom] = open(path, 'w')
            counts[chrom] = 0
            if header_line:
                out_fps[chrom].write(header_line)
        return out_fps[chrom]

    sys.stderr.write('Splitting {} ...\n'.format(cmpt_tsv_fn))
    n_total = 0
    with open(cmpt_tsv_fn) as fp:
        for line in fp:
            if line.startswith('#'):
                header_line = line
                for fp2 in out_fps.values():
                    fp2.write(line)
                continue

            ele = line.rstrip('\n').split('\t')
            qname, trans_ids, fmt = _get_trans_ids(ele)

            if fmt == 'unknown':
                sys.stderr.write('WARNING: unrecognised format ({}  columns), skipping line\n'.format(len(ele)))
                continue
            if qname is None:
                # header row embedded in file body (Bambu)
                header_line = line
                for fp2 in out_fps.values():
                    fp2.write(line)
                continue

            n_total += 1

            # find the chromosome(s) for this read's transcripts
            chroms = set()
            for tid in trans_ids:
                span = trans_to_span.get(tid)
                if span:
                    chroms.add(span[0])

            if not chroms:
                _get_fp('unassigned').write(line)
                counts['unassigned'] += 1
            else:
                for chrom in chroms:
                    _get_fp(chrom).write(line)
                    counts[chrom] += 1

            if n_total % 1_000_000 == 0:
                sys.stderr.write('  {} M reads processed\r'.format(n_total // 1_000_000))
                sys.stderr.flush()

    sys.stderr.write('  {} reads processed\n'.format(n_total))

    for fp2 in out_fps.values():
        fp2.close()

    return counts


def parser_argv():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Split a read-isoform compatible TSV (ESPRESSO/Bambu) into '
                    'per-chromosome files using transcript coordinates from an updated GTF.')
    parser.add_argument('cmpt_tsv',    metavar='read_isoform_compatible.tsv',
                        help='Read-isoform compatible TSV (ESPRESSO or Bambu format)')
    parser.add_argument('updated_gtf', metavar='updated.gtf',
                        help='Updated GTF file used to map transcripts to chromosomes')
    parser.add_argument('out_dir',     metavar='out_dir',
                        help='Output directory for per-chromosome TSV files')
    parser.add_argument('-p', '--prefix', type=str, default='cmpt',
                        help='Filename prefix for output files '
                             '(output: <prefix>.<chrom>.tsv)')
    return parser.parse_args()


def main():
    args = parser_argv()
    counts = split_cmpt_tsv(
        cmpt_tsv_fn = args.cmpt_tsv,
        updated_gtf = args.updated_gtf,
        out_dir     = args.out_dir,
        prefix      = args.prefix)

    sys.stderr.write('\nOutput summary:\n')
    for chrom, n in sorted(counts.items(), key=lambda x: (x[0] == 'unassigned', x[0])):
        sys.stderr.write('  {:20s}  {:>10,d} reads\n'.format(chrom, n))
    sys.stderr.write('Done.\n')


if __name__ == '__main__':
    main()
