import sys
import os
import pysam as ps
import edlib as ed
from collections import defaultdict as dd
import argparse

_5_ada = 'CTACACGACGCTCTTCCGATCT' # length: 22
_3_ada = 'CTCTGCGTTGATACCACTGCTT' # length: 22
comp_seq_dict = dd(lambda: 'N')
comp_seq_dict.update({'A':'T', 'a':'T', 'C':'G', 'c':'G', 'G':'C', 'g':'C', 'T':'A', 't':'A', 'U':'A', 'u':'A'})
_ed = 'editDistance'
_loc = 'locations'
task = "locations"
globa = "NW"
infix = "HW"
prefix = "SHW"

def revcomp(seq):
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([comp_seq_dict[r] for r in rev_seq])
    return rev_comp_seq

def get_five_cand_seq(r=ps.AlignedSegment(), cand_len=200):
    extra_len = 50
    cand_seq = ''
    clip_len = r.cigar[0][1] if r.cigar[0][0] == ps.CSOFT_CLIP else 0
    start = max(0, clip_len-cand_len)
    end = min(clip_len+extra_len, r.query_length)
    cand_seq = r.query_sequence[start:end]
    return cand_seq

def get_three_cand_seq(r=ps.AlignedSegment(), cand_len=200):
    extra_len = 50
    cand_seq = ''
    clip_len = r.cigar[-1][1] if r.cigar[-1][0] == ps.CSOFT_CLIP else 0
    start = max(-clip_len-extra_len, -r.query_length)
    end = -clip_len+cand_len if -clip_len+cand_len < 0 else r.query_length
    cand_seq = r.query_sequence[start:end]
    return cand_seq

def classify_read(read, five_ada=_5_ada, three_ada=_3_ada, five_max_ed=int(len(_5_ada)*0.3), three_max_ed=int(len(_3_ada)*0.3)):
    five_cand_seq = get_five_cand_seq(read)
    three_cand_seq = get_three_cand_seq(read)
    _5_res = ed.align(five_ada, five_cand_seq, task=task, mode=infix, k=five_max_ed)
    _3_res = ed.align(three_ada, three_cand_seq, task=task, mode=infix, k=five_max_ed)
    _rc_5_res = ed.align(five_ada, revcomp(three_cand_seq), task=task, mode=infix, k=five_max_ed)
    _rc_3_res = ed.align(three_ada, revcomp(five_cand_seq), task=task, mode=infix, k=three_max_ed)
    if (_5_res[_ed] != -1 and _3_res[_ed] != -1) or (_rc_5_res[_ed] != -1 and _rc_3_res[_ed] != -1):
        return 1, 1 # five, three
    elif _5_res[_ed] != -1 and _rc_5_res[_ed] != -1:
        return 2, 0
    elif _3_res[_ed] != -1 and _rc_3_res[_ed] != -1:
        return 0, 2
    elif _5_res[_ed] != -1 or _rc_5_res[_ed] != -1:
        return 1, 0
    elif _3_res[_ed] != -1 or _rc_3_res[_ed] != -1:
        return 0, 1
    else:
        return 0, 0

def classify_type(five, three):
    if five == 1:
        if three == 0:
            return 'RT'
        else: # three == 1:
            return 'Proper'
    elif five == 2:
        return 'RT-RT'
    else: # five == 0
        if three == 1:
            return 'TSO'
        elif three == 2:
            return 'TSO-TSO'
        else:
            return 'No-RT-TSO'

def classify_artifact(in_bam, out_fn, five_ada=_5_ada, three_ada=_3_ada):
    # return 0,0
    five_max_ed=int(len(_5_ada)*0.3)
    three_max_ed=int(len(_3_ada)*0.3)
    with ps.AlignmentFile(in_bam) as in_fp, open(out_fn, 'w') as out_fp:
        for read in in_fp:
            if read.is_unmapped:
                continue
            five, three = classify_read(read, five_ada, three_ada, five_max_ed, three_max_ed)
            read_type = classify_type(five, three)
            out_fp.write('{}\t{}\n'.format(read.query_name, read_type))

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: split chimeric reads, i.e., mapped to multiple genomic regions, in BAM file into multiple subreads".format(os.path.basename(__file__)))
    parser.add_argument('in_bam', metavar='input.sam/bam', type=str, help='Input alignment file')
    parser.add_argument('out_tsv', metavar='output.tsv', type=str, help='Output file, type for each read in the alignment')

    # optional arguments: 5'/3' adapter sequences
    parser.add_argument('-5', '--five-ada', type=str, default=_5_ada, help='5\' adapter sequence')
    parser.add_argument('-3', '--three-ada', type=str, default=_3_ada, help='3\' adapter sequence')
    
    return parser.parse_args()

def main():
    args = parser_argv()
    # cmd = ' '.join(sys.argv)
    in_bam, out_fn = args.in_bam, args.out_tsv
    five_ada, three_ada = args.five_ada, args.three_ada
    classify_artifact(in_bam, out_fn, five_ada, three_ada)

if __name__ == '__main__':
    main()
