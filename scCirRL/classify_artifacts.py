import sys
import os
import pysam as ps
import edlib as ed
from collections import defaultdict as dd
import argparse

_5_ada = 'CTACACGACGCTCTTCCGATCT' # length: 22
_3_ada = 'CTCTGCGTTGATACCACTGCTT' # length: 22
five_ada = 'CTACACGACGCTCTTCCGATCT' # length: 22
three_ada = 'CTCTGCGTTGATACCACTGCTT' # length: 22
polyA = 'A'*10
polyT = 'T'*10
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
    _poly = ed.align(polyT, five_cand_seq, task=task, mode=infix, k=3)
    _3_res = ed.align(three_ada, three_cand_seq, task=task, mode=infix, k=five_max_ed)
    _rc_5_res = ed.align(five_ada, revcomp(three_cand_seq), task=task, mode=infix, k=five_max_ed)
    _rc_3_res = ed.align(three_ada, revcomp(five_cand_seq), task=task, mode=infix, k=three_max_ed)
    _rc_poly = ed.align(polyT, revcomp(three_cand_seq), task=task, mode=infix, k=3)
    _5_cnt = 0
    _3_cnt = 0
    if (_5_res[_ed] != -1 and _3_res[_ed] != -1) or \
       (_rc_5_res[_ed] != -1 and _rc_5_res[_ed] != -1 and _rc_3_res[_ed] != -1):
        _5_cnt, _3_cnt = 1, 1
    elif _5_res[_ed] != -1 and _rc_5_res[_ed] != -1:
        if _5_res[_ed] < _rc_5_res[_ed]:
            if _rc_poly[_ed] != -1:
                _5_cnt, _3_cnt = 2, 0
            else:
                _5_cnt, _3_cnt = 1, 0
        elif _5_res[_ed] > _rc_5_res[_ed]:
            if _poly[_ed] != -1:
                _5_cnt, _3_cnt = 2, 0
            else:
                _5_cnt, _3_cnt = 1, 0
        else:
            _5_cnt, _3_cnt = 2, 0
    elif (_5_res[_ed] != -1) or (_rc_5_res[_ed] != -1):
        _5_cnt, _3_cnt = 1, 0
    elif _3_res[_ed] != -1 and _rc_3_res[_ed] != -1:
        _5_cnt, _3_cnt = 0, 2
    elif _3_res[_ed] != -1 or _rc_3_res[_ed] != -1:
         _5_cnt, _3_cnt =  0, 1
    else:
        _5_cnt, _3_cnt = 0, 0
    return _5_cnt, _3_cnt, _5_res[_ed], _3_res[_ed], _poly[_ed], _rc_5_res[_ed], _rc_3_res[_ed], _rc_poly[_ed]

def get_read_to_bc(bc_umi_tsv):
    sys.stderr.write('Reading barcode/UMI calling result file: {}\n'.format(bc_umi_tsv))
    read_to_bc = dd(lambda: False)
    with open(bc_umi_tsv) as in_fp:
        header = True
        idx = dict()
        for line in in_fp:
            ele = line.strip().split('\t')
            if header:
                header = False
                idx = {ele[i]:i for i in range(len(ele))}
                continue
            read_names = ele[idx['ReadNames']].split(',')
            for read_name in read_names:
                read_to_bc[read_name] = True
    sys.stderr.write('Done\n')
    return read_to_bc

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

def classify_artifact(in_bam, sample, read_to_bc, out_prefix, five_ada=_5_ada, three_ada=_3_ada):
    # return 0,0
    sys.stderr.write('Classifying artifacts for sample: {}\n'.format(sample))
    five_max_ed=int(len(_5_ada)*0.3)
    three_max_ed=int(len(_3_ada)*0.3)
    read_type_to_cnt = dd(lambda: 0)
    with ps.AlignmentFile(in_bam) as in_fp, open(out_prefix+'_detail.tsv', 'w') as detail_fp, \
         open(out_prefix+'_summary.tsv', 'w') as summary_fp:
        for read in in_fp:
            if read.is_unmapped:
                continue
            five, three, _5_ed, _3_ed, _poly_ed, _5_rc_ed, _3_rc_ed, _rc_poly_ed = classify_read(read, five_ada, three_ada, five_max_ed, three_max_ed)
            read_type = classify_type(five, three)
            if read.query_name in read_to_bc:
                read_type += '-with BC'
            else:
                read_type += '-no BC'
            read_type_to_cnt[read_type] += 1
            detail_fp.write('{}\t{}\t{}\t{}\n'.format(sample, read.query_name, read_type, ','.join([str(e) for e in [_5_ed, _poly_ed, _3_ed, _5_rc_ed, _rc_poly_ed, _3_rc_ed]])))
        for read_type, cnt in read_type_to_cnt.items():
            summary_fp.write('{}\t{}\t{}\n'.format(sample, read_type, cnt))
    sys.stderr.write('Done\n')

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: classify mapped reads into different categories, i.e., full-length, TSO/RT artifacts".format(os.path.basename(__file__)))
    parser.add_argument('in_bam', metavar='input.sam/bam', type=str, help='Input alignment file')
    parser.add_argument('bc_umi_tsv', metavar='bc_umi.tsv', type=str, help='Barcode/UMI calling result file')
    parser.add_argument('out_prefix', metavar='output_prefix', type=str, help='Output file, type for each read in the alignment')

    # optional arguments: 5'/3' adapter sequences
    parser.add_argument('-s', '--sample', type=str, default='', help='Sample name')
    parser.add_argument('-5', '--five-ada', type=str, default=_5_ada, help='5\' adapter sequence')
    parser.add_argument('-3', '--three-ada', type=str, default=_3_ada, help='3\' adapter sequence')
    
    return parser.parse_args()

def main():
    args = parser_argv()
    # cmd = ' '.join(sys.argv)
    in_bam, bc_umi_tsv, out_prefix = args.in_bam, args.bc_umi_tsv, args.out_prefix
    five_ada, three_ada = args.five_ada, args.three_ada
    sample = in_bam
    if '/' in sample:
        sample = sample.split('/')[-1]
    if args.sample != '':
        sample = args.sample
    read_to_bc = get_read_to_bc(bc_umi_tsv)
    classify_artifact(in_bam, sample, read_to_bc, out_prefix, five_ada, three_ada)

if __name__ == '__main__':
    main()
