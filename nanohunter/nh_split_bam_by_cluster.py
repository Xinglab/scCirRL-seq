import sys
import os
import pysam as ps


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


def main():
    if len(sys.argv) != 4:
        print('{} in.sorted.bam(with CB tag) in.bc_cell.list out_dir'.format(sys.argv[0]))
        sys.exit(1)

    in_bam_fn, bc_cell_list, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    bc_to_cell = parse_bc_cell_list(bc_cell_list)
    split_bam(in_bam_fn, bc_to_cell, out_dir)

if __name__ == '__main__':
    main()