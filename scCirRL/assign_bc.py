import edlib as ed
import sys
import pysam as ps
from collections import defaultdict as dd

from . import utils as ut
from . import seq_utils as su
from . import cluster_umi as cu
from .parameter import scrl_para_class
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

def assign_ref_bc1(mp_fetch_set1, in_bam, scrl_ref_bcs, scrl_cand_ref_bc_seq, read_to_trans, trans_to_gene_id_name, scrl_para=scrl_para_class()): #, bc_len, bc_max_ed, umi_len, umi_max_ed):
    bu_res = []
    contig, start, end = mp_fetch_set1
    n_perfect_in_ref_reads, n_perfect_uniq_to_ref_reads, n_imperfect_in_ref_reads, n_imperfect_uniq_to_ref_reads = 0, 0, 0, 0
    # with ps.AlignmentFile(in_sam_fn) as in_sam,  open('{}.txt'.format(mp.current_process().name), 'a') as err_fp:
    bu_reads = []
    bc_eds = dict()
    out_rs = {}
    n_processed_reads = 0
    for r in in_bam.fetch(contig, start, end):
        qname = r.query_name
        skip, tag = su.check_for_skip(r, scrl_para)
        if skip:
            continue
        n_processed_reads += 1
        # if qname == debug_qname:
            # print('OK')
        # assign bc/umi for each read
        ref_bc, umi, bc_ed, is_perfect, is_exact = search_for_matched_ref_bc(r, scrl_ref_bcs, scrl_cand_ref_bc_seq, 
                                                                             scrl_para.five_ada, scrl_para.five_max_ed,
                                                                             scrl_para.bc_len, scrl_para.bc_max_ed, scrl_para.umi_len)
        if ref_bc and umi:
            if is_perfect:
                if is_exact: n_perfect_in_ref_reads += 1
                else: n_perfect_uniq_to_ref_reads += 1
            else:
                if is_exact: n_imperfect_in_ref_reads += 1
                else: n_imperfect_uniq_to_ref_reads += 1
            bu_res.append([qname, ref_bc, umi, is_perfect])
            bu_reads.append(qname)
            bc_eds[qname] = bc_ed
            out_rs[qname] = r
    # UMI clustering within each bc
    # return:
    # (bc, umi, reads, cmpt_trans)
    umi_clu_res_dict, umi_clu_res_list = cu.umi_clustering(scrl_para.log_fn, read_to_trans, trans_to_gene_id_name, bu_res, scrl_para.umi_max_ed)
    res = [bu_reads, bc_eds, out_rs, umi_clu_res_dict, umi_clu_res_list, n_perfect_in_ref_reads, n_perfect_uniq_to_ref_reads, n_imperfect_in_ref_reads, n_imperfect_uniq_to_ref_reads]
    return res, n_processed_reads


# def sp_assign_ref_bc(mp_fetch_set, scrl_ref_bcs, scrl_cand_ref_bc_seq, bc_len, bc_max_ed, umi_len, umi_max_ed, in_sam_fn, read_to_trans, trans_to_gene_id_name, out_bu_fn, out_bu_bam):
def assign_ref_bc(mp_fetch_set, n_total_reads, scrl_ref_bcs, scrl_cand_ref_bc_seq, 
                  trans_to_gene_id_name, read_to_trans, read_to_cate, scrl_para=scrl_para_class()):
    ut.err_log_format_time(scrl_para.log_fn, str="Assigning barcode & UMI ...")
    ut.err_log_progress_bar(scrl_para.log_fn)
    n_perfect_in_ref_reads, n_perfect_uniq_to_ref_reads, n_imperfect_in_ref_reads, n_imperfect_uniq_to_ref_reads = 0, 0, 0, 0
    bc_ed_count_dict = dd(lambda: 0) # {bc_ed: count}
    n_processed_reads = 0
    n_existing_stars = 0
    # percentage_processed_reads = 10 # start from 5%
    with open(scrl_para.out_bu_fn, 'w') as out_fp, ps.AlignmentFile(scrl_para.long_bam) as in_bam:
        out_fp.write('#CellBarcode\tUMI\tReadCount\tTranscriptID\tGeneID\tGeneName\tReadNames\tTransCate\tSpliceTag\n')
        bam_header_dict = in_bam.header.to_dict()
        scrl_append_pg_dict = {'ID': __program__, 'PN': __program__, 'VN': __version__, 'CL': __cmd__}
        if 'PG' in bam_header_dict:
            pg_PNs = [pg['PN'] for pg in bam_header_dict['PG']]
            if pg_PNs.count(__program__) > 0:
                scrl_append_pg_dict['ID'] = '{}.{}'.format(__program__, pg_PNs.count(__program__) + 1)
            bam_header_dict['PG'].append(scrl_append_pg_dict)
        else:
            bam_header_dict['PG'] = [scrl_append_pg_dict]
        with ps.AlignmentFile(scrl_para.out_bu_bam, 'wb', header=bam_header_dict) as out_bam:
            for mp_fetch_set1 in mp_fetch_set:
                out_res, _n_processed_reads = assign_ref_bc1(mp_fetch_set1, in_bam, 
                                                             scrl_ref_bcs, scrl_cand_ref_bc_seq, 
                                                             read_to_trans, trans_to_gene_id_name,
                                                             scrl_para)
                n_processed_reads += _n_processed_reads
                n_existing_stars = ut.err_log_progress_star(scrl_para.log_fn, n_total_reads, n_processed_reads, n_existing_stars)
                bu_reads, bc_eds, out_rs, umi_clu_res_dict, umi_clu_res_list, sub_n_perfect_in_ref_reads, sub_n_perfect_uniq_to_ref_reads, sub_n_imperfect_in_ref_reads, sub_n_imperfect_uniq_to_ref_reads = out_res
                
                splice_tags = dd(lambda: 'N') # N or Y
                for read in bu_reads:
                    if read in umi_clu_res_dict:
                        bc, umi, cmpt_trans, cmpt_gene_id, cmpt_gene_names = umi_clu_res_dict[read]
                        cate = read_to_cate[read]
                        r = out_rs[read]
                        bc_ed = bc_eds[read]
                        bc_ed_count_dict[bc_ed] += 1
                        r.set_tag('CB', bc, 'Z')
                        r.set_tag('UB', umi, 'Z')
                        r.set_tag('BE', bc_ed, 'i')
                        r.set_tag('TC', cate, 'Z')
                        r.set_tag('TI', ','.join(cmpt_trans), 'Z')
                        r.set_tag('GI', ','.join(cmpt_gene_id), 'Z')
                        r.set_tag('GN', ','.join(cmpt_gene_names), 'Z')
                        if r.has_tag('ts'):
                            splice_tags[read] = 'Y'
                        out_bam.write(r)

                for bc, umi, reads, cmpt_trans, cmpt_gene_id, cmpt_gene_names in umi_clu_res_list:
                    out_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(bc, umi, len(reads), ','.join(cmpt_trans), ','.join(cmpt_gene_id), ','.join(cmpt_gene_names), 
                                                                           ','.join(reads), ','.join([read_to_cate[read] for read in reads]),
                                                                           ','.join(splice_tags[read] for read in reads)))

                n_perfect_in_ref_reads += sub_n_perfect_in_ref_reads
                n_perfect_uniq_to_ref_reads += sub_n_perfect_uniq_to_ref_reads
                n_imperfect_in_ref_reads += sub_n_imperfect_in_ref_reads
                n_imperfect_uniq_to_ref_reads += sub_n_imperfect_uniq_to_ref_reads
    # build index for output bam
    ut.err_log_format_time(scrl_para.log_fn, str="Building index for {} ...".format(scrl_para.out_bu_bam))
    ps.index(scrl_para.out_bu_bam)
    ut.err_log_format_time(scrl_para.log_fn, str="Building index done!")
    n_total_assigned_reads = n_perfect_in_ref_reads + n_perfect_uniq_to_ref_reads + n_imperfect_in_ref_reads + n_imperfect_uniq_to_ref_reads
    assign_ratio = '(%.1f%%)' % (n_total_assigned_reads / n_total_reads * 100)
    bc_ed_count_str = '\n'.join(['%21d : %d (%.1f%%)' % (k, v, v/n_total_reads*100) for k, v in sorted(bc_ed_count_dict.items())])
    
    ut.err_log_format_time(scrl_para.log_fn, str="Assigning barcode & UMI done!")
    ut.err_log_format_time(scrl_para.log_fn, str="Barcode calling summary:\n" +
f'''Total mapped reads    : {n_total_reads}
Total cell barcodes   : {len(scrl_ref_bcs)}
Barcode-called reads  : {n_total_assigned_reads} {assign_ratio}
Barcode edit distance : read count\n{bc_ed_count_str}
''')
# Perfect In Reference:\t\t{n_perfect_in_ref_reads}
# Perfect Unique to Reference:\t\t{n_perfect_uniq_to_ref_reads}
# Imperfect In Reference:\t\t{n_imperfect_in_ref_reads}
# Imperfect Unique to Reference:\t\t{n_imperfect_uniq_to_ref_reads}\n
# ''')
    return