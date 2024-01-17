import sys, os
from collections import defaultdict as dd
import pysam as ps
from kneed import KneeLocator
import edlib as ed

from . import seq_utils as su
from . import utils as ut
from .parameter import nanohunter_para


def get_knee_point(umi_cnts):
    if len(umi_cnts) < 10:
        ut.err_fatal_format_time(__name__, 'Not able to detect knee point barcode: too few ({}) candidate barcodes.'.format(len(umi_cnts)))
    bc_ranks = [i for i in range(1, len(umi_cnts)+1)]
    bc_kneedle = KneeLocator(bc_ranks, umi_cnts, S=1.0, curve="concave", direction="increasing")
    return bc_kneedle.knee

def get_cumulative_cnt(cnt, max_n=10000, min_cnt=2):
    cumu_cnt = []
    cumu = 0
    n = 0
    for c in cnt:
        if c < min_cnt:
            break
        if n >= max_n:
            break
        cumu_cnt.append(c+cumu)
        n += 1
        cumu += c
    return cumu_cnt


def get_merged_umi_cnt(bc_umi_to_read, ref_bc, bc):
    umi1 = set(bc_umi_to_read[ref_bc].keys())
    umi2 = set(bc_umi_to_read[bc].keys())
    return len(umi1.union(umi2))


def has_ovlp(pos1, pos2):
    ref1, start1, end1 = pos1
    ref2, start2, end2 = pos2
    if ref1 != ref2 or end1 < start2 or end2 < start1:
        return False
    return True


# NO indel allowed here
def cand_ref_bc_merge1(bc, bc_len, bc_max_ed, cand_bcs, perfect_bc_umi_read, read_to_map_pos):
    cand_bc_seq = ('N' * (bc_max_ed+1)).join(cand_bcs)  # use joined bc or do multiple times of ed.align?
    ed_res = ed.align(bc, cand_bc_seq, mode="HW", task="locations", k=bc_max_ed)
    if ed_res['editDistance'] == -1:
        return ''
    elif ed_res['editDistance'] == 0:
        sys.stderr.write('Duplciated barcodes.\n')
        sys.exit(1)

    cand_matched_ref_bcs = []
    for loc in ed_res['locations']:
        bc1 = cand_bcs[int(loc[1] / (bc_len + bc_max_ed + 1))]
        if bc1 not in cand_matched_ref_bcs:
            cand_matched_ref_bcs.append(bc1)

    umis = perfect_bc_umi_read[bc].keys()
    matched_ref_bc = dict()
    for ref_bc1 in cand_matched_ref_bcs:
        shared_umis = list(set(umis).intersection(set(perfect_bc_umi_read[ref_bc1].keys())))
        for shared_umi in shared_umis:
            for ref_read in perfect_bc_umi_read[ref_bc1][shared_umi]:
                for read in perfect_bc_umi_read[bc][shared_umi]:
                    if has_ovlp(read_to_map_pos[ref_read], read_to_map_pos[read]):
                        matched_ref_bc[ref_bc1] = 1
                        break
    if len(matched_ref_bc) == 1:
        return list(matched_ref_bc.keys())[0]
    else:
        return ''


#@func: cand_ref_bc_merge
#       merge barcodes if they share same UMI
#       try to find matched bc for 2nd half of bc in the 1st half
#@return:
#       perfect_bc_to_umi_cnt:
def cand_ref_bc_merge(bc_len, bc_max_ed, perfect_bc_to_umi_cnt, perfect_knee_bc_rank, perfect_bc_umi_to_read, read_to_map_pos):
    ref_bcs = []

    cand_ref_bc_rank = int(perfect_knee_bc_rank * 1.5)
    all_perfect_bcs = list(perfect_bc_to_umi_cnt.keys())
    half_knee = int(perfect_knee_bc_rank/2)

    # for first half-knee cand bc, no merge check
    for bc in all_perfect_bcs[:half_knee]:
        ref_bcs.append(bc)

    # only first cand_ref_bc_rank cand bc can be added to ref_bc
    for ii, bc in enumerate(all_perfect_bcs[half_knee:]):
        i = ii + half_knee
        ref_bc = cand_ref_bc_merge1(bc, bc_len, bc_max_ed, ref_bcs, perfect_bc_umi_to_read, read_to_map_pos)
        if ref_bc:
            perfect_bc_to_umi_cnt[ref_bc] = get_merged_umi_cnt(perfect_bc_umi_to_read, ref_bc, bc)
            # perfect_bc_to_umi_cnt[ref_bc] += perfect_bc_to_umi_cnt[bc]
            del perfect_bc_to_umi_cnt[bc]
        else:
            if i < cand_ref_bc_rank:
                ref_bcs.append(bc)
    return perfect_bc_to_umi_cnt


def write_perfect_bu_cnt(bc_to_umi_cnt, bu_fn):
    with open(bu_fn, 'w') as fp:
        fp.write('BC\tRank\tUMICount\n')
        rank = 0
        for bc in bc_to_umi_cnt:
            rank += 1
            fp.write('{}\t{}\t{}\n'.format(bc, rank, bc_to_umi_cnt[bc]))

#       read_to_map_pos: {read: (ref, start, end)}
def collect_perfect_bc(in_sam_fn, five_ada, bc_len, umi_len, max_clip_len, use_all_alignments):
    perfect_bc_umi_to_read = dd(lambda: dd(lambda: []))  # {BC: {UMI: [reads]}}
    n_perfect_reads = 0

    with ps.AlignmentFile(in_sam_fn) as in_sam:  # one read may be processed multiple times in mp
        for r in in_sam:
            if r.is_unmapped:
                continue
            if not use_all_alignments:
                if r.is_supplementary or r.is_secondary:
                    continue
                if (r.cigar[0][0] == ps.CSOFT_CLIP and r.cigar[0][1]) > max_clip_len or (r.cigar[-1][0] == ps.CSOFT_CLIP and r.cigar[-1][1] > max_clip_len):
                    continue
            qname = r.query_name
            # sys.stderr.write('qname\t{}\n'.format(qname))
            # if qname == 'bcb94ec6-8ff4-4f55-a173-ff4c63a51d63_rep0_2.6':
                # print('ok')
            bc, umi = su.get_perfect_bu(r, five_ada, bc_len, umi_len)
            if bc and umi:
                perfect_bc_umi_to_read[bc][umi].append(qname)
                n_perfect_reads += 1
    return perfect_bc_umi_to_read, n_perfect_reads


def sp_collect_cand_ref_bc(nh_para=nanohunter_para()):
    in_sam_fn, five_ada, bc_len, umi_len, = nh_para.long_bam, nh_para.five_ada, nh_para.bc_len, nh_para.umi_len
    max_clip_len, use_all_alignments, high_qual_bu_fn = nh_para.max_clip_len, nh_para.all_alignments, nh_para.high_qual_bu_fn

    ut.err_log_format_time(nh_para.log_fn, __name__, "Collecting candidate cell barcodes ...")
    
    perfect_bc_umi_to_read, n_perfect_reads = collect_perfect_bc(in_sam_fn, five_ada, bc_len, umi_len, max_clip_len, use_all_alignments)
    perfect_bc_to_umi_cnt = dd(lambda: 0)
    for bc in perfect_bc_umi_to_read:
        umis = list(perfect_bc_umi_to_read[bc].keys())
        perfect_bc_to_umi_cnt[bc] = len(umis)  # UMI count

    # 1st knee point
    perfect_bc_to_umi_cnt = dict(sorted(perfect_bc_to_umi_cnt.items(), key=lambda d: -d[1]))
    # perfect_bc_to_umi_cnt = dict(sorted(perfect_bc_to_umi_cnt.items(), key=lambda d: (-d[1], d[0])))
    perfect_bc_to_cumulative_umi_cnt = get_cumulative_cnt(perfect_bc_to_umi_cnt.values(), max_n=10000, min_cnt=2)  # TODO max_n
    # first knee rank
    perfect_knee_bc_rank = get_knee_point(perfect_bc_to_cumulative_umi_cnt)  # multi by 1.1?
    ut.err_log_format_time(nh_para.log_fn, __name__, 'Knee point barcode rank: {} (detected from {} barcodes, {} high-quality reads)'.format(perfect_knee_bc_rank, len(perfect_bc_to_cumulative_umi_cnt), n_perfect_reads))
    
    if perfect_knee_bc_rank is None:
        return []
    write_perfect_bu_cnt(perfect_bc_to_umi_cnt, high_qual_bu_fn)
    ref_bcs = list(perfect_bc_to_umi_cnt.keys())[:perfect_knee_bc_rank]
    # merged_perfect_bc_to_umi_count = cand_ref_bc_merge(perfect_bc_to_umi_cnt, perfect_knee_bc_rank, perfect_bc_umi_to_read, read_to_map_pos)

    # 2nd knee point
    # merged_perfect_bc_to_umi_count = dict(sorted(merged_perfect_bc_to_umi_count.items(), key=lambda d: -d[1]))
    # merged_perfect_bc_to_cumulative_umi_cnt = get_cumulative_cnt(merged_perfect_bc_to_umi_count.values(), max_n=10000, min_cnt=2)
    # ref_knee_bc_rank = get_knee_point(merged_perfect_bc_to_cumulative_umi_cnt)  # multi by 1.1?
    # extended_ref_bc_rank = int(ref_knee_bc_rank * 1.0)
    # ut.err_format_time(__name__, 'Second knee rank: {} ({})'.format(extended_ref_bc_rank, ref_knee_bc_rank))
    # write_perfect_bu_cnt(merged_perfect_bc_to_umi_count, conf_bu_fn)
    # ref_bcs = list(merged_perfect_bc_to_umi_count.keys())[:extended_ref_bc_rank]
    ut.err_log_format_time(nh_para.log_fn, __name__, "Collecting candidate cell barcodes done! ({})".format(len(ref_bcs)))
    return ref_bcs


def collect_cand_ref_bc(nh_para=nanohunter_para()):
    # bc_len, umi_len, max_clip_len, use_all_alignments, in_sam_fn, high_qual_bu_fn):
    ref_bcs = sp_collect_cand_ref_bc(nh_para) #in_sam_fn, bc_len, umi_len, max_clip_len, use_all_alignments, high_qual_bu_fn)
    return ref_bcs
