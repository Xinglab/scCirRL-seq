import sys
import os
import pysam as ps
import argparse
from .__init__ import __version__, __program__


def split_bam(in_bam_fn, bc_to_cell, out_dir):
    all_cells = set(bc_to_cell.values())
    out_bam_dict = {}

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    with ps.AlignmentFile(in_bam_fn) as in_bam:
        for cell in all_cells:
            out_bam = ps.AlignmentFile(out_dir + '/{}.bam'.format(cell), 'wb', template=in_bam)
            out_bam_dict[cell] = out_bam
        for r in in_bam:
            if r.has_tag('CB'):
                bc = r.get_tag('CB')
                if bc in bc_to_cell:
                    cell = bc_to_cell[bc]
                    out_bam = out_bam_dict[cell]
                    out_bam.write(r)

        for cell in all_cells:
            out_bam_dict[cell].close()
    for cell in all_cells:
        ps.index(out_dir + '/{}.bam'.format(cell))


def parse_bc_cell_list(in_fn):
    bc_to_cell = dict()
    with open(in_fn) as fp:
        for line in fp:
            ele = line.rstrip().rsplit('\t')
            if len(ele) == 2:
                bc = ele[0].rsplit('-')[0]
                clu = ele[1].replace('/', ',')
                clu = clu.replace(' ', '_')
                bc_to_cell[bc] = clu
    return bc_to_cell


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: split bam file by cell clusters".format(os.path.basename(__file__)))
    parser.add_argument('in_bam', metavar='input.sam/bam', type=str, help='Input alignment file')
    parser.add_argument('bc_to_cell_type', metavar='in.bc_to_cell_type.list', type=str, help='Input barcode to cell type list')
    parser.add_argument('out_dir', metavar='output_dir', type=str, help='Output directory')

    return parser.parse_args()

def main():
    args = parser_argv()
    in_bam_fn, bc_cell_list, out_dir = args.in_bam, args.bc_to_cell_type, args.out_dir
    bc_to_cell = parse_bc_cell_list(bc_cell_list)
    split_bam(in_bam_fn, bc_to_cell, out_dir)

if __name__ == '__main__':
    main()