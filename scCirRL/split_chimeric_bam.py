import sys
import os
import pysam as ps
from collections import defaultdict as dd
import argparse
from .__init__ import __version__, __program__

def split_chimeric_bam(in_bam, out_bam, cmd):
    with ps.AlignmentFile(in_bam) as in_fp:
        bam_header_dict = in_fp.header.to_dict()
        append_pg_dict = {'ID': __program__, 'PN': __program__, 'VN': __version__, 'CL': cmd}
        if 'PG' in bam_header_dict:
            bam_header_dict['PG'].append(append_pg_dict)
        else:
            bam_header_dict['PG'] = [append_pg_dict]
        with ps.AlignmentFile(out_bam, mode='wb', header=bam_header_dict) as out_fp:
            r = 0
            for read in in_fp:
                if read.is_unmapped:
                    continue
                if read.is_supplementary:
                    read.query_name += '_r{}'.format(hex(r))
                    # read.query_name += '_r0x{}'.format(f'{r:x}')
                    if read.flag == 2048:
                        read.flag = 0
                    elif read.flag == 2064:
                        read.flag = 16
                out_fp.write(read)
                r += 1
    
def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: split chimeric reads, i.e., mapped to multiple genomic regions, in BAM file into multiple subreads".format(os.path.basename(__file__)))
    parser.add_argument('in_bam', metavar='input.sam/bam', type=str, help='Input alignment file')
    parser.add_argument('out_bam', metavar='output.bam', type=str, help='Output alignment file')

    return parser.parse_args()

def main():
    args = parser_argv()
    cmd = ' '.join(sys.argv)
    in_bam, out_bam = args.in_bam, args.out_bam
    split_chimeric_bam(in_bam, out_bam, cmd)

if __name__ == '__main__':
    main()
