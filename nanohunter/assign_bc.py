import edlib as ed
import sys
import pysam as ps
from collections import defaultdict as dd

from . import utils as ut
from . import seq_utils as su
from . import cluster_umi as cu
from .parameter import nanohunter_para
from .__init__ import __version__
from .__init__ import __program__
from .__init__ import __cmd__


debug_qname = 'f8c30655-ec79-4d65-81f3-f51908c74a88_18_799'

def get_uniq_match_ref_bc(qname, bc, bu, cand_ref_bcs, cand_ref_bc_seq, bc_len, bc_max_ed, umi_len):
    ed_res = ed.align(bc, cand_ref_bc_seq, mode="HW", task="locations", k=bc_max_ed)
    if ed_res['editDistance'] >= 0:
        ref_bus = su.get_matched_bu(cand_ref_bcs, bu, ed_res['locations'], bc_len, bc_max_ed, umi_len)
        if len(ref_bus) == 1:  # one best matched ref bc
            return ref_bus[0][0], ref_bus[0][1], ed_res['editDistance']
    return '', '', ed_res['editDistance']


# Four categories of barcode sequence:
# 1: Perfect In Reference
# 2: Perfect Unique to Reference
# 3: Imperfect In Reference
# 4: Imperfect Unique to Reference
#    during search, perfect/imperfect information in previous step (collect_candidate.mp_perfect_bc) are not used,
#    so all the reads will be processed again.
# allow indels in putative bc/umi
def search_for_matched_ref_bc(r, cand_ref_bcs, cand_ref_bc_seq, five_ada, five_max_ed, bc_len, bc_max_ed, umi_len):
    # check if there is an exact match of any ref_bc, no ed.align needed
    is_perfect, is_exact = False, False
    bc, umi, bu, is_perfect = su.get_bu_from_seq(r, five_ada, five_max_ed, bc_len, umi_len)
    if bc and umi:
        if bc in cand_ref_bcs:  # 1/3: PerfectInRef/ImperfectInRef
            is_exact = True
            return bc, umi, 0, is_perfect, is_exact
        else:  # search with ed.align
            is_exact = False
            ref_bc, umi, bc_ed = get_uniq_match_ref_bc(r.query_name, bc, bu, cand_ref_bcs, cand_ref_bc_seq, bc_len, bc_max_ed, umi_len)
            if ref_bc and umi:  # 2/4: PerfectUniqToRef/ImperfectUniqToRef
                return ref_bc, umi, bc_ed, is_perfect, is_exact
    return '', '', -1, is_perfect, is_exact


def get_cmpt_genes(cmpt_trans, trans_to_gene_id_name):
    if not cmpt_trans or not trans_to_gene_id_name: # no compatible transcript
        return {'NA'}, {'NA'}
    gene_ids, gene_names = set(), set()
    for trans in cmpt_trans:
        if trans not in trans_to_gene_id_name:
            continue
        gene_ids.add(trans_to_gene_id_name[trans]['id'])
        gene_names.add(trans_to_gene_id_name[trans]['name'])
    if gene_ids == set():
        gene_ids = {'NA'}
    if gene_names == set():
        gene_names = {'NA'}
    return gene_ids, gene_names

def assign_ref_bc1(mp_fetch_set1, in_bam, nh_ref_bcs, nh_cand_ref_bc_seq, read_to_trans, nh_para=nanohunter_para()): #, bc_len, bc_max_ed, umi_len, umi_max_ed):
    bu_res = []
    umi_cluster_res = []
    contig, start, end = mp_fetch_set1
    n_perfect_in_ref_reads, n_perfect_uniq_to_ref_reads, n_imperfect_in_ref_reads, n_imperfect_uniq_to_ref_reads = 0, 0, 0, 0
    # with ps.AlignmentFile(in_sam_fn) as in_sam,  open('{}.txt'.format(mp.current_process().name), 'a') as err_fp:
    bu_reads = set()
    bc_eds = dict()
    out_rs = {}
    n_processed_reads = 0
    for r in in_bam.fetch(contig, start, end):
        if r.is_unmapped:
            continue
        if not nh_para.all_alignments:
            if r.is_supplementary or r.is_secondary:
                continue
            if (r.cigar[0][0] == ps.CSOFT_CLIP and r.cigar[0][1]) > nh_para.max_clip_len or (r.cigar[-1][0] == ps.CSOFT_CLIP and r.cigar[-1][1] > nh_para.max_clip_len):
                continue
        qname = r.query_name
        n_processed_reads += 1
        # if qname == debug_qname:
            # print('OK')
        # assign bc/umi for each read
        ref_bc, umi, bc_ed, is_perfect, is_exact = search_for_matched_ref_bc(r, nh_ref_bcs, nh_cand_ref_bc_seq, 
                                                                             nh_para.five_ada, nh_para.five_max_ed,
                                                                             nh_para.bc_len, nh_para.bc_max_ed, nh_para.umi_len)
        if ref_bc and umi:
            if is_perfect:
                if is_exact: n_perfect_in_ref_reads += 1
                else: n_perfect_uniq_to_ref_reads += 1
            else:
                if is_exact: n_imperfect_in_ref_reads += 1
                else: n_imperfect_uniq_to_ref_reads += 1
            bu_res.append([qname, ref_bc, umi, is_perfect])
            bu_reads.add(qname)
            bc_eds[qname] = bc_ed
            out_rs[qname] = r
    # UMI clustering within each bc
    # return:
    # (bc, umi, reads, cmpt_trans)
    umi_cluster_res = cu.umi_clustering(read_to_trans, bu_res, nh_para.umi_max_ed)
    # mp_q_res.put([bu_res, sub_n_perfect_in_ref_reads, sub_n_perfect_uniq_to_ref_reads, sub_n_imperfect_in_ref_reads, sub_n_imperfect_uniq_to_ref_reads])
    res = [bu_reads, bc_eds, out_rs, umi_cluster_res, n_perfect_in_ref_reads, n_perfect_uniq_to_ref_reads, n_imperfect_in_ref_reads, n_imperfect_uniq_to_ref_reads]
    return res, n_processed_reads


# def sp_assign_ref_bc(mp_fetch_set, nh_ref_bcs, nh_cand_ref_bc_seq, bc_len, bc_max_ed, umi_len, umi_max_ed, in_sam_fn, read_to_trans, trans_to_gene_id_name, out_bu_fn, out_bu_bam):
def sp_assign_ref_bc(mp_fetch_set, n_total_reads, nh_ref_bcs, nh_cand_ref_bc_seq, trans_to_gene_id_name, read_to_trans, nh_para=nanohunter_para()):
    ut.err_log_format_time(nh_para.log_fn, __name__, "Assigning barcode & UMI ...")
    n_perfect_in_ref_reads, n_perfect_uniq_to_ref_reads, n_imperfect_in_ref_reads, n_imperfect_uniq_to_ref_reads = 0, 0, 0, 0
    bc_ed_count_dict = dd(lambda: 0) # {bc_ed: count}
    n_processed_reads = 0
    percentage_processed_reads = 10 # start from 5%
    with open(nh_para.out_bu_fn, 'w') as out_fp, ps.AlignmentFile(nh_para.long_bam) as in_bam:
        out_fp.write('#BC\tUMI\tReadCount\tTranscript\tGeneID\tGeneName\tReadNames\n')
        bam_header_dict = in_bam.header.to_dict()
        nhs_append_pg_dict = {'ID': __program__, 'PN': __program__, 'VN': __version__, 'CL': __cmd__}
        if 'PG' in bam_header_dict:
            bam_header_dict['PG'].append(nhs_append_pg_dict)
        else:
            bam_header_dict['PG'] = [nhs_append_pg_dict]
        with ps.AlignmentFile(nh_para.out_bu_bam, 'wb', header=bam_header_dict) as out_bam:
            for mp_fetch_set1 in mp_fetch_set:
                out_res, _n_processed_reads = assign_ref_bc1(mp_fetch_set1, in_bam, nh_ref_bcs, nh_cand_ref_bc_seq, read_to_trans, nh_para)
                n_processed_reads += _n_processed_reads
                if n_processed_reads / n_total_reads >= percentage_processed_reads / 100:
                    ut.err_log_format_time(nh_para.log_fn, __name__, '%d reads processed (%.1f%%)' % (n_processed_reads, n_processed_reads / n_total_reads * 100))
                    percentage_processed_reads = int(n_processed_reads / n_total_reads * 100) + 10
                bu_reads, bc_eds, out_rs, umi_cluster_res, sub_n_perfect_in_ref_reads, sub_n_perfect_uniq_to_ref_reads, sub_n_imperfect_in_ref_reads, sub_n_imperfect_uniq_to_ref_reads = out_res
                out_reads = []
                for bc, umi, reads, cmpt_trans in umi_cluster_res:
                    # write BC/UMI bu
                    cmpt_gene_id, cmpt_gene_names = get_cmpt_genes(cmpt_trans, trans_to_gene_id_name)
                    out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(bc, umi, len(reads), ','.join(cmpt_trans), ','.join(cmpt_gene_id), ','.join(cmpt_gene_names), ','.join(reads)))
                    out_reads.extend(reads)

                    # write BAM with BC/UMI tag
                    for read in reads:
                        r = out_rs[read]
                        bc_ed = bc_eds[read]
                        bc_ed_count_dict[bc_ed] += 1
                        r.set_tag('CB', bc, 'Z')
                        r.set_tag('UB', umi, 'Z')
                        r.set_tag('BE', bc_ed, 'i')
                        r.set_tag('TI', ','.join(cmpt_trans), 'Z')
                        r.set_tag('GI', ','.join(cmpt_gene_id), 'Z')
                        r.set_tag('GN', ','.join(cmpt_gene_names), 'Z')
                        out_bam.write(r)

                c_not_o = bu_reads - set(out_reads)
                o_not_c = set(out_reads) - bu_reads
                if c_not_o:
                    sys.stderr.write('Count but Not output:\t{}\n'.format(','.join(c_not_o)))
                if o_not_c:
                    sys.stderr.write('Ooutput but Not count:\t{}\n'.format(','.join(o_not_c)))
                # write BAM records with BC/UMI tags
                n_perfect_in_ref_reads += sub_n_perfect_in_ref_reads
                n_perfect_uniq_to_ref_reads += sub_n_perfect_uniq_to_ref_reads
                n_imperfect_in_ref_reads += sub_n_imperfect_in_ref_reads
                n_imperfect_uniq_to_ref_reads += sub_n_imperfect_uniq_to_ref_reads
    n_total_assigned_reads = n_perfect_in_ref_reads + n_perfect_uniq_to_ref_reads + n_imperfect_in_ref_reads + n_imperfect_uniq_to_ref_reads
    assign_ratio = '(%.1f%%)' % (n_total_assigned_reads / n_total_reads * 100)
    bc_ed_count_str = '\n'.join(['%21d : %d (%.1f%%)' % (k, v, v/n_total_reads*100) for k, v in sorted(bc_ed_count_dict.items())])
    
    ut.err_log_format_time(nh_para.log_fn, __name__, "Assigning barcode & UMI done!")
    ut.err_log_format_time(nh_para.log_fn, __name__, "Barcode calling summary:\n" +
f'''Total mapped reads    : {n_total_reads}
Total cell barcodes   : {len(nh_ref_bcs)}
Barcode-called reads  : {n_total_assigned_reads} {assign_ratio}
Barcode edit distance : read count\n{bc_ed_count_str}
''')
# Perfect In Reference:\t\t{n_perfect_in_ref_reads}
# Perfect Unique to Reference:\t\t{n_perfect_uniq_to_ref_reads}
# Imperfect In Reference:\t\t{n_imperfect_in_ref_reads}
# Imperfect Unique to Reference:\t\t{n_imperfect_uniq_to_ref_reads}\n
# ''')
    return


# 2nd round: assign bc and UMI clustering within each bc
# def assign_ref_bc(mp_fetch_set, nh_ref_bcs, nh_cand_ref_bc_seq, bc_len, bc_max_ed, umi_len, umi_max_ed, in_sam_fn, read_to_trans, trans_to_gene_id_name, out_bu_fn, out_bu_bam):
def assign_ref_bc(mp_fetch_set, n_total_reads, nh_ref_bcs, nh_cand_ref_bc_seq, trans_to_gene_id_name, read_to_trans, nh_para):
    sp_assign_ref_bc(mp_fetch_set, n_total_reads, nh_ref_bcs, nh_cand_ref_bc_seq, trans_to_gene_id_name, read_to_trans, nh_para)
