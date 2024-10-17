import pysam as ps
import edlib as ed
import sys
import re
import gzip
from collections import defaultdict as dd

from . import utils as ut
from .parameter import scrl_para_class

# lib structure:
# 5' ada | BC16 | UMI12 | polyT-cDNA | 3' ada
_5_ada = 'CTACACGACGCTCTTCCGATCT' # length: 22
_3_ada = 'CTCTGCGTTGATACCACTGCTT' # length: 22
# _5_ada = 'CTACACGACGCTCTTCCGATCT'[-12:] # length: 22
# _3_ada = 'CTCTGCGTTGATACCACTGCTT'[:12] # length: 22


_ed = 'editDistance'
_loc = 'locations'
task = "locations"
globa = "NW"
infix = "HW"
prefix = "SHW"

comp_seq_dict = dd(lambda: 'N')
comp_seq_dict.update({'A':'T', 'a':'T', 'C':'G', 'c':'G', 'G':'C', 'g':'C', 'T':'A', 't':'A', 'U':'A', 'u':'A'})

def revcomp(seq):
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([comp_seq_dict[r] for r in rev_seq])
    return rev_comp_seq

# def scrl_init():
#     global n_perfect_reads
#     global n_perfect_in_ref_reads
#     global n_perfect_uniq_to_ref_reads
#     global n_imperfect_in_ref_reads
#     global n_imperfect_uniq_to_ref_reads
#     n_perfect_reads = 0
#     n_perfect_in_ref_reads = 0
#     n_perfect_uniq_to_ref_reads = 0
#     n_imperfect_in_ref_reads = 0
#     n_imperfect_uniq_to_ref_reads = 0

#     global scrl_ref_bcs
#     global scrl_cand_ref_bc_seq
#     # global scrl_read_to_ref_bc_umi
#     # global scrl_perfect_bc_to_umi
#     # global scrl_assigned_reads
#     scrl_ref_bcs = []
#     scrl_cand_ref_bc_seq = ''
#     # scrl_read_to_ref_bc_umi = dict()
#     # scrl_perfect_bc_to_umi = dict()
#     # scrl_assigned_reads = dict()

# skip if:
# . is_unmapped
# . is_supplementary or is_secondary, if only_primary
#   has supplementary (SA), if not skip_chimeric
def check_for_skip(r=ps.AlignedSegment(), scrl_para=scrl_para_class()):
    skip = False
    tag = 'proper'
    if r.is_unmapped or r.is_secondary:
        skip = True
        tag = 'unmapped'
    if r.is_supplementary:
        tag = 'non_primary'
        if scrl_para.only_primary:
            skip = True
    if r.has_tag('SA'):
        tag = 'chimeric'
        if scrl_para.skip_chimeric:
            skip = True
    return skip, tag

def collect_mp_fetch_set(scrl_para=scrl_para_class()):
    in_sam_fn = scrl_para.long_bam
    ut.err_log_format_time(scrl_para.log_fn, str="Pre-processing/scanning BAM ...")
    fetch_region_set = []
    n_total_reads = 0
    with ps.AlignmentFile(in_sam_fn) as in_sam:
        last_chrom, last_start, last_max_end = '', -1, -1
        first = True
        for r in in_sam:
            skip, tag = check_for_skip(r, scrl_para)
            if skip:
                continue
            qname, chrom, start, end = r.query_name, r.reference_name, r.reference_start, r.reference_end
            n_total_reads += 1
            # print('{}\t{}\t{}\t{}'.format(qname, chrom, start, end))
            if first:
                first = False
                last_chrom, last_start, last_max_end = chrom, start, end
            if chrom != last_chrom or start > last_max_end:
                # print('{}:{}-{}'.format(last_chrom, last_start, last_max_end))
                fetch_region_set.append((last_chrom, last_start, last_max_end))
                last_chrom, last_start, last_max_end = chrom, start, end
            if end > last_max_end:
                last_max_end = end
        fetch_region_set.append((last_chrom, last_start, last_max_end))
        # print('{}:{}-{}'.format(last_chrom, last_start, last_max_end))
    ut.err_log_format_time(scrl_para.log_fn, str="Pre-processing/scanning BAM done! ({} regions, {} total mapped reads to process)".format(len(fetch_region_set), n_total_reads))
    return fetch_region_set, n_total_reads

def is_rev_strand(r=ps.AlignedSegment()):
    if r.has_tag('ts'):
        ts = 1 if r.get_tag('ts') == '+' else -1
        if (r.is_forward and ts == 1) or (r.is_reverse and ts == -1):
            return False
        else:
            return True
    # no ts tag
    seq = r.query_sequence
    polyT, polyA = 'T'*15, 'A'*15

    # left
    cigar_op, cigar_len = r.cigar[0][0], r.cigar[0][1]
    if cigar_op == ps.CSOFT_CLIP:
        ed_res = ed.align(polyT, seq[:cigar_len+50], task="locations", mode="HW")
        if ed_res[_ed] < 2 and ed_res[_ed] != -1:
            return True
    # right
    cigar_op, cigar_len = r.cigar[-1][0], r.cigar[-1][1]
    if cigar_op == ps.CSOFT_CLIP:
        ed_res = ed.align(polyA, seq[-cigar_len-50:], task="locations", mode="HW")
        if ed_res[_ed] < 2 and ed_res[_ed] != -1:
            return False
    return None

# |-cand_seq-polyT-cDNA-|
# |-cDNA-polyA-cand_seq-|
# cand_seq: up to 200 bp, contains barcode/UMI/adapter
def get_cand_seq(r=ps.AlignedSegment(), cand_len=200):
    extra_len = 50
    cand_seq = ''
    is_rev = is_rev_strand(r)
    if is_rev is None:
        return ''
    elif is_rev:  # polyT
        clip_len = r.cigar[0][1] if r.cigar[0][0] == ps.CSOFT_CLIP else 0
        start = max(0, clip_len-cand_len)
        end = min(clip_len+extra_len, r.query_length)
        cand_seq = r.query_sequence[start:end]
    else:  # polyA
        clip_len = r.cigar[-1][1] if r.cigar[-1][0] == ps.CSOFT_CLIP else 0
        start = max(-clip_len-extra_len, -r.query_length)
        end = -clip_len+cand_len if -clip_len+cand_len < 0 else r.query_length
        cand_seq = revcomp(r.query_sequence[start:end])
    return cand_seq

# 1. 5' ada exists and no mis/indel
# 2. polyA/polyT & ts & strand match
# 3. BC+UMI = 28 bp (16+12)
def get_perfect_bu_from_seq(cand_seq, five_ada, bc_len, umi_len):
    bc, umi = '', ''
    if five_ada in cand_seq:
        start = cand_seq.index(five_ada) + len(five_ada)
        end = start + bc_len + umi_len
        if cand_seq[end:end+10] == 'T'*10:
            bc = cand_seq[start:start+bc_len]
            umi = cand_seq[end-umi_len:end]
    return bc, umi


# if prime match, break,
# else, try secondary
def get_perfect_bu(r, five_ada, bc_len, umi_len):
    cand_seq = get_cand_seq(r)
    bc, umi = get_perfect_bu_from_seq(cand_seq, five_ada, bc_len, umi_len)
    return bc, umi

def get_bu_from_seq(r, five_ada, five_max_ed, bc_len, umi_len):
    bc, umi, bu = '', '', ''
    is_perfect = False

    cand_seq = get_cand_seq(r)
    # check if perfect
    if five_ada in cand_seq:
        start = cand_seq.index(five_ada) + len(five_ada)
        end = start + bc_len + umi_len
        if cand_seq[end:end+10] == 'T'*10:
            is_perfect = True
    else:
        _5_res = ed.align(five_ada, cand_seq, task=task, mode=infix, k=five_max_ed)
        # _poly = ed.align('T'*10, cand_seq, task=task, mode=infix, k=3)
        if _5_res[_ed] == -1: # or _poly[_ed] == -1:
            return '', '', '', False
        start = _5_res[_loc][0][1] + 1

    bc = cand_seq[start:start+bc_len]
    umi = cand_seq[start+bc_len:start+bc_len+umi_len]
    bu = cand_seq[start:start+bc_len+umi_len+5]
    if len(bc) != bc_len or len(umi) != umi_len:
        return '', '', '', False
    else:
        return bc, umi, bu, is_perfect

def get_matched_bu(cand_bcs, bu, locs, bc_len, bc_max_ed, umi_len):
    bcs = dict()
    bus = []
    for loc in locs:
        # sys.stderr.write('{}\t{}\t{}\n'.format(loc, bc_len, max_ed))
        bc = cand_bcs[int(loc[1] / (bc_len + bc_max_ed + 1))]
        if bc not in bcs:
            bcs[bc] = 1
            res = ed.align(bc, bu, task="locations", mode="HW", k=bc_max_ed)
            if res['editDistance'] != -1:
                umi = bu[res['locations'][0][1]+1:res['locations'][0][1]+1+umi_len]
                if len(umi) == umi_len:
                    bus.append((bc, umi))
    return bus

def collect_ref_barcodes(scrl_para=scrl_para_class()):
    list_fn, bc_len = scrl_para.bc_list, scrl_para.bc_len
    ref_bcs = []
    if list_fn.endswith('.gz'):
        fp = gzip.open(list_fn, 'rt')
    else:
        fp = open(list_fn)
    for line in fp:
        ele = re.split('-| |\t', line.rstrip())
        if len(ele[0]) != bc_len:
            continue
        if ele[0] not in ref_bcs:
            ref_bcs.append(ele[0])

    ut.err_log_format_time(scrl_para.log_fn, str="Collected {} reference barcodes from {}".format(len(ref_bcs), list_fn))
    return ref_bcs
