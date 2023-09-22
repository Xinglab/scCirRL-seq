import pysam as ps
import edlib as ed
import sys
import re
import gzip
from collections import defaultdict as dd

from . import utils as ut

# stru:
# 5' ada | BC16 | UMI12 | polyT-cDNA | 3' ada
_5_ada = 'CTACACGACGCTCTTCCGATCT'
_3_ada = 'CCCATGTACTCTGCGTTGATACCACTGCTT'
_bc_len = 16
_umi_len = 12
_bu_len = 28
_ed_ratio = 0.3

_5_k = int(len(_5_ada) * _ed_ratio)
_ed = 'editDistance'
_loc = 'locations'
task = "locations"
globa = "NW"
infix = "HW"
prefix = "SHW"

_bc_max_ed = 2
_umi_max_ed = 1
_max_clip_len = 200

# import pyabpoa as pa
# poa_aln = pa.msa_aligner()


# mp_chunk_region_size = 10000


comp_seq_dict = dd(lambda: 'N')
comp_seq_dict.update({'A':'T', 'a':'T', 'C':'G', 'c':'G', 'G':'C', 'g':'C', 'T':'A', 't':'A', 'U':'A', 'u':'A'})


def revcomp(seq):
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([comp_seq_dict[r] for r in rev_seq])
    return rev_comp_seq


def nh_init():
    global n_perfect_reads
    global n_perfect_in_ref_reads
    global n_perfect_uniq_to_ref_reads
    global n_imperfect_in_ref_reads
    global n_imperfect_uniq_to_ref_reads
    n_perfect_reads = 0
    n_perfect_in_ref_reads = 0
    n_perfect_uniq_to_ref_reads = 0
    n_imperfect_in_ref_reads = 0
    n_imperfect_uniq_to_ref_reads = 0

    global nh_ref_bcs
    global nh_cand_ref_bc_seq
    # global nh_read_to_ref_bc_umi
    # global nh_perfect_bc_to_umi
    # global nh_assigned_reads
    nh_ref_bcs = []
    nh_cand_ref_bc_seq = ''
    # nh_read_to_ref_bc_umi = dict()
    # nh_perfect_bc_to_umi = dict()
    # nh_assigned_reads = dict()


# def get_mp_fetch_set(in_sam_fn, chunk_region_size):
#     fetch_set = []
#     with ps.AlignmentFile(in_sam_fn) as in_sam:
#         header = in_sam.header
#         for ref in header.references:
#             tot_len = header.get_reference_length(ref)
#             for i in range(int((tot_len-1)/chunk_region_size)+1):
#                 start = i * chunk_region_size
#                 end = (i+1) * chunk_region_size
#                 if end > tot_len:
#                     end = tot_len
#                 fetch_set.append((ref, start, end))
#     return fetch_set


def get_mp_fetch_set(in_sam_fn):
    ut.err_format_time(__name__, "Collecting BAM fetch regions ...")
    fetch_region_set = []
    with ps.AlignmentFile(in_sam_fn) as in_sam:
        last_chrom, last_start, last_max_end = '', -1, -1
        first = True
        for r in in_sam:
            if r.is_unmapped or r.is_supplementary or r.is_secondary:
                continue
            # filter out chimeric reads
            if (r.cigar[0][0] == ps.CSOFT_CLIP and r.cigar[0][1]) > _max_clip_len or (r.cigar[-1][0] == ps.CSOFT_CLIP and r.cigar[-1][1] > _max_clip_len):
                continue
            qname, chrom, start, end = r.query_name, r.reference_name, r.reference_start, r.reference_end
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
    ut.err_format_time(__name__, "Collecting BAM fetch regions done! ({} regions)".format(len(fetch_region_set)))
    return fetch_region_set


def is_rev_strand(r=ps.AlignedSegment()):
    if r.has_tag('ts'):
        ts = 1 if r.get_tag('ts') == '+' else -1
        if (r.flag == 0 and ts == 1) or (r.flag == 16 and ts == -1):
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
        if ed_res[_ed] < 2:
            return True
    # right
    cigar_op, cigar_len = r.cigar[-1][0], r.cigar[-1][1]
    if cigar_op == ps.CSOFT_CLIP:
        ed_res = ed.align(polyA, seq[-cigar_len-50:], task="locations", mode="HW")
        if ed_res[_ed] < 2:
            return False
    return None


def get_cDNAseq(r, bc, umi):
    extra_len = 50
    cand_seq = ''
    is_rev = is_rev_strand(r)
    if is_rev is None:
        return ''
    if is_rev:  # polyT
        whole_seq = r.query_sequence
        cand_seq = r.query_sequence[:r.cigar[0][1]+extra_len]
    else:  # polyA
        whole_seq = revcomp(r.query_sequence)
        cand_seq = revcomp(r.query_sequence[-r.cigar[-1][1]-extra_len:])

    res = ed.align(_5_ada+bc+umi, cand_seq, task=task, mode=infix, k=int(len(_5_ada+bc+umi)*_ed_ratio))
    if res[_ed] == -1:
        return ''
    pT_start = res[_loc][0][1] + 1  # polyT start
    cDNA_seq = whole_seq[pT_start+10:]
    return cDNA_seq


# def get_cons_seq(seqs):
#     if len(seqs) < 1:
#         return ''
#     elif len(seqs) <= 2:
#         return seqs[0]
#     else:
#         res = poa_aln.msa(seqs, out_cons=True, out_msa=False)
#         return res.cons_seq[0]


def get_cand_seq(r=ps.AlignedSegment()):
    extra_len = 50
    cand_seq = ''
    is_rev = is_rev_strand(r)
    if is_rev is None:
        return ''
    elif is_rev:  # polyT
        cand_seq = r.query_sequence[:r.cigar[0][1]+extra_len]
    else:  # polyA
        cand_seq = revcomp(r.query_sequence[-r.cigar[-1][1]-extra_len:])
    return cand_seq


# 1. 5' ada exists and no mis/indel
# 2. polyA/polyT & ts & strand match
# 3. BC+UMI = 28 bp (16+12)
def get_perfect_bu_from_seq(cand_seq, bc_len, umi_len):
    bc, umi = '', ''
    if _5_ada in cand_seq:
        start = cand_seq.index(_5_ada) + len(_5_ada)
        end = start + _bu_len
        if cand_seq[end:end+10] == 'T'*10:
            bc = cand_seq[start:start+bc_len]
            umi = cand_seq[end-umi_len:end]
    return bc, umi


# if prime match, break,
# else, try secondary
def get_perfect_bu(r, bc_len, umi_len):
    cand_seq = get_cand_seq(r)
    bc, umi = get_perfect_bu_from_seq(cand_seq, bc_len, umi_len)
    return bc, umi


def get_confident_cDNA_seq(r):
    extra_len = 50
    cand_seq = ''
    is_rev = is_rev_strand(r)
    if is_rev is None:
        return ''
    if is_rev:
        whole_seq = r.query_sequence
        cand_seq = r.query_sequence[:r.cigar[0][1]+extra_len]
    else:
        whole_seq = revcomp(r.query_sequence)
        cand_seq = revcomp(r.query_sequence[-r.cigar[-1][1]-extra_len:])
    cDNA_seq = ''
    if _5_ada in cand_seq and 'T'*10 in cand_seq:
        start = cand_seq.index(_5_ada) + len(_5_ada)
        end = cand_seq.index('T'*10)
        if end - start == _bu_len:
            cDNA_seq = whole_seq[end+10:]
    return cDNA_seq


def get_bu_from_seq(r, bc_len, umi_len):
    bc, umi, bu = '', '', ''
    is_perfect = False

    cand_seq = get_cand_seq(r)
    # check if perfect
    if _5_ada in cand_seq:
        start = cand_seq.index(_5_ada) + len(_5_ada)
        end = start + _bu_len
        if cand_seq[end:end+10] == 'T'*10:
            is_perfect = True
    else:
        _5_res = ed.align(_5_ada, cand_seq, task=task, mode=infix, k=_5_k)
        if _5_res[_ed] == -1:
            return '', '', '', False
        start = _5_res[_loc][0][1] + 1

    bc = cand_seq[start:start+bc_len]
    umi = cand_seq[start+bc_len:start+bc_len+umi_len]
    bu = cand_seq[start:start+bc_len+umi_len+5]
    if len(bc) != bc_len or len(umi) != umi_len:
        return '', '', '', False
    else:
        return bc, umi, bu, is_perfect


def get_min_ed(seq, seq_set):
    min_ed = len(seq)
    for s in seq_set:
        res = ed.align(seq, s, task="locations")
        if res['editDistance'] < min_ed:
            min_ed = res['editDistance']
    return min_ed


def get_matched_bu(cand_bcs, bu, locs, bc_len, bc_max_ed, umi_len):
    bcs = dict()
    bus = []
    for loc in locs:
        # sys.stderr.write('{}\t{}\t{}\n'.format(loc, bc_len, max_ed))
        bc = cand_bcs[int(loc[1] / (bc_len + bc_max_ed + 1))]
        if bc not in bcs:
            bcs[bc] = 1
            res = ed.align(bc, bu, task="locations", mode="HW")
            umi = bu[res['locations'][0][1]+1:res['locations'][0][1]+1+umi_len]
            if len(umi) == umi_len:
                bus.append((bc, umi))
    return bus


def get_ref_barcodes(list_fn, bc_len):
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

    ut.err_format_time(__name__, "Collected {} reference barcodes from {}.".format(len(ref_bcs), list_fn))
    return ref_bcs
