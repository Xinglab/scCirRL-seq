import sys
import os
import argparse
from collections import defaultdict as dd
import math

# R stats
from rpy2 import robjects as rb
from rpy2.robjects import r
from rpy2.robjects.packages import importr
stats = importr('stats')

fdr_thres = 0.05   # <= 0.05
p_thres = 0.05     # <= 0.05
delta_ratio = 0.05 # > 0.05

_trans_phased_ratio = 0.7
_trans_phased_cnt_per_allele = 5
_gene_total_phased_cnt = 20
_gene_phased_ratio = 0.7
_phased_trans_cnt_per_gene = 2
# filtered_out_gene_prefix = 'HLA-'


ALL_HAPS = ['H1', 'H2']  # 'none'

# gene wise allele specific info: 
#   as_gene -> cell types
#   as_transcript -> cell types
class gene_wise_allele_specific_info:
    def __init__(self, gene, gene_name):
        # 5 categories
        # LowExp: low expression, skipped in allele-specific splicing analysis
        # Unphased: overall phased ratio < 0.8, skipped in allele-specific splicing analysis
        # NoSplice: only one haplotype-resolved transcript, skipped in allele-specific splicing analysis
        # True: allele-specific splicing
        # False: non-allele-specific splicing
        self.as_gene_cate = dd(lambda: "Unknown") # {cell_type: cate}
        self.as_trans_cate = dd(lambda: dd(lambda: "None")) # {transcript: {cell_type: hap}}
        self.gene = gene
        self.gene_name = gene_name
        self.chrom = ''

    def add_as_cell_type(self, cell_type):
        self.as_gene_cate[cell_type] = 'True'
    def add_non_as_cell_type(self, cell_type):
        self.as_gene_cate[cell_type] = 'False'
    # def add_unphased_cell_type(self, cell_type): # NA: not enough data to determine if allele specific, none: no data
        # self.as_gene_cate[cell_type] = 'Unphased'
    # def add_nosplice_cell_type(self, cell_type):
        # self.as_gene_cate[cell_type] = 'NoSplice'
        
    def add_as_transcript_cell_type(self, cell_type, trans, hap):
        self.as_trans_cate[trans][cell_type] = hap
    def add_non_as_transcript_cell_type(self, cell_type, trans):
        self.as_trans_cate[trans][cell_type] = 'False'
    def set_chrom(self, chrom):
        self.chrom = chrom

    def get_as_gene_to_cell_type(self):
        return self.as_gene_cate
    def get_as_trans_to_cell_type(self):
        return self.as_trans_cate
    def get_all_as_trans_hap(self):
        trans_hap = set()
        for trans in self.as_trans:
            for cell_type, hap in self.as_trans[trans].items():
                trans_hap.add((trans, hap))
        return trans_hap

def get_chisq_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    cnt1 = [i for i in cnt1]
    cnt2 = [i for i in cnt2]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.chisq_test(df)
    return list(out[2])[0]

def get_fisher_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    cnt1 = [i for i in cnt1]
    cnt2 = [i for i in cnt2]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.fisher_test(df)
    return list(out[0])[0]

def get_fdr(ps=[]):
    rp = rb.FloatVector(ps)
    return list(stats.p_adjust(rp, method='fdr'))

def parse_bc_cell_list(in_fn):
    bc_to_cell = dict()
    with open(in_fn) as fp:
        for line in fp:
            ele = line.rstrip().rsplit('\t')
            if len(ele) == 2:
                bc_to_cell[ele[0].rsplit('-')[0]] = ele[1]
    return bc_to_cell

def get_reads_to_hap(hap_list):
    read_to_hap = dd(lambda: '')
    with open(hap_list) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            qname, hap, phase_set, chrome = ele
            if hap != 'none':
                read_to_hap[qname] = hap
        return read_to_hap

def get_hap(qnames, read_to_hap, trans_id, gene_name):
    hap_to_read_cnt = dd(lambda: 0)
    for qname in qnames:
        if qname in read_to_hap:
            hap_to_read_cnt[read_to_hap[qname]] += 1

    hap_to_read_cnt = dict(sorted(hap_to_read_cnt.items(), key=lambda d: -d[1]))
    if len(hap_to_read_cnt) == 0:  # no hap
        return 'none'
    elif len(hap_to_read_cnt) > 1:  # 2 haps
        if list(hap_to_read_cnt.values())[0] > list(hap_to_read_cnt.values())[1]:
            return list(hap_to_read_cnt.keys())[0]
        else:
            return 'none'
    else:  # 1 hap
        return list(hap_to_read_cnt.keys())[0]

def filter_by_read_type(_reads, _read_cates, read_type):
    if read_type == []:
        return _reads, _read_cates
    reads, read_cates = [], []
    for qname, cate in zip(_reads, _read_cates):
        if cate in read_type:
            reads.append(qname)
            read_cates.append(cate)
    return reads, read_cates

def get_cluster_wise_gene_to_trans_reads(scrl_bu_fn, use_read_cnt, read_type, bc_to_cluster, read_to_hap):
    cluster_wise_gene_to_trans_reads = dd(lambda: dd(lambda: dd(lambda: dd(lambda: 0))))  # {bc: {gene: {trans: {hap: cnt}}}}
    gene_id_to_name = dict()
    with open(scrl_bu_fn) as fp:
        header_idx = dict()
        for line in fp:
            ele = line.rsplit()
            if line.startswith('#'):
                ele[0] = ele[0][1:]
                header_idx = {ele[i]: i for i in range(len(ele))}
                continue
            # bc, umi, read_cnt, trans_id, gene_id, gene_name, read_names, read_cates, split_cates = ele
            bc = ele[header_idx['CellBarcode']]
            trans_id, gene_id, gene_name = ele[header_idx['TranscriptID']], ele[header_idx['GeneID']], ele[header_idx['GeneName']]
            read_names, read_cates = ele[header_idx['ReadNames']].rsplit(','), ele[header_idx['TransCate']].rsplit(',')
            read_names, read_cates = filter_by_read_type(read_names, read_cates, read_type)
            if trans_id == 'NA' or ',' in trans_id or gene_id == 'NA' or ',' in gene_id:  # only keep single-match reads
                continue
            # if gene_name.startswith(filtered_out_gene_prefix):
                # continue
            if bc not in bc_to_cluster:
                continue
            cluster = bc_to_cluster[bc]
            if use_read_cnt:
                for read in read_names:
                    hap = read_to_hap[read]
                    cluster_wise_gene_to_trans_reads[cluster][gene_id][trans_id][hap] += 1
            else:
                hap = get_hap(read_names, read_to_hap, trans_id, gene_name)  # H1/H2/none
                cluster_wise_gene_to_trans_reads[cluster][gene_id][trans_id][hap] += 1
            gene_id_to_name[gene_id] = gene_name
    return cluster_wise_gene_to_trans_reads, gene_id_to_name

# True/False: gene_total_count >= 10, n_trans >= 2, phased_count_per_trans >= 10, phased_ratio_per_trans >= 0.8
# Unable to phase (Unphased):  gene_phased_ratio < 0.8
# No alternative splicing (NoSplice): n_phased_transcript < 0.8
# Low expression (None): gene_total_count < 20
def set_as_gene(gene_to_as_info, gene, gene_name, cluster, as_type):
    if gene not in gene_to_as_info:
        gene_to_as_info[gene] = gene_wise_allele_specific_info(gene, gene_name)
    if as_type == 'True':
        gene_to_as_info[gene].add_as_cell_type(cluster)
    elif as_type == 'False':
        gene_to_as_info[gene].add_non_as_cell_type(cluster)
    # elif as_type == 'Unphased':
        # gene_to_as_info[gene].add_unphased_cell_type(cluster)
    # elif as_type == 'NoSplice':
        # gene_to_as_info[gene].add_nosplice_cell_type(cluster)

def as_gene_fdr_corr_and_asts(gene_list, p_list, 
                              gene_to_max_ratio, min_trans_cnt, 
                              cluster_wise_gene_to_trans_reads, gene_id_to_name,
                              asg_detailed_fp, asts_detailed_fp, 
                              inc_non_sign_gene, inc_non_sign_trans, gene_to_as_info):
    gene_fdrs = get_fdr(p_list)
    for gene_, fdr, p in zip(gene_list, gene_fdrs, p_list):
        cluster, gene_id = gene_[0], gene_[1]
        gene_name = gene_id_to_name[gene_id]
        gene_to_trans_reads = cluster_wise_gene_to_trans_reads[cluster]
        transs = gene_to_trans_reads[gene_id].keys()
        hap1_cnts = [gene_to_trans_reads[gene_id][trans][ALL_HAPS[0]] for trans in transs]
        hap2_cnts = [gene_to_trans_reads[gene_id][trans][ALL_HAPS[1]] for trans in transs]
        hap1_ratios = [cnt / sum(hap1_cnts) if sum(hap1_cnts) > 0 else 0 for cnt in hap1_cnts]
        hap2_ratios = [cnt / sum(hap2_cnts) if sum(hap2_cnts) > 0 else 0 for cnt in hap2_cnts]
        out_str = [cluster, gene_id, gene_name, fdr, 
                   ','.join(transs), 
                   ','.join(list(map(str, hap1_cnts))), 
                   ','.join(list(map(str, hap2_cnts))), 
                   ','.join(list(map(str, hap1_ratios))), 
                   ','.join(list(map(str, hap2_ratios)))]

        is_as_gene = True
        if math.isnan(fdr):
            # XXX write nan fdr, where gene is allele-specific-expressed
            is_as_gene = False
            if inc_non_sign_gene:
                asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            elif not inc_non_sign_trans:
                continue
        if len(gene_to_trans_reads[gene_id]) < min_trans_cnt:   # genes with < 2 haplotype-resolved transcripts
            print(f'{gene_id}:\tNoSplice')
            is_as_gene = False
            if inc_non_sign_gene:
                asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            # set_as_gene(gene_to_as_info, gene_id, gene_name, cluster, 'NoSplice')
            elif not inc_non_sign_trans:
                continue
        elif fdr > fdr_thres or gene_to_max_ratio[gene_] <= delta_ratio:
            is_as_gene = False
            set_as_gene(gene_to_as_info, gene_id, gene_name, cluster, 'False')
            if inc_non_sign_gene:
                asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            elif not inc_non_sign_trans:
                continue

        # write significant allele-specific spliced genes to output file
        if is_as_gene:
            set_as_gene(gene_to_as_info, gene_id, gene_name, cluster, 'True')
            asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))

        if is_as_gene or inc_non_sign_trans:
            # perform fisher exact test to identify allele-specific transcripts
            tot_cnt = dd(lambda: 0.0)  # {hap: tot_cnt of all trans}
            for hap in ALL_HAPS:
                for trans in gene_to_trans_reads[gene_id]:
                    tot_cnt[hap] += gene_to_trans_reads[gene_id][trans][hap]
            for trans in gene_to_trans_reads[gene_id]:
                cnts, ratios = [], []
                for hap in ALL_HAPS:
                    cnt = [gene_to_trans_reads[gene_id][trans][hap], tot_cnt[hap] - gene_to_trans_reads[gene_id][trans][hap]]
                    cnts.append(cnt)
                    if tot_cnt[hap] == 0:
                        ratios.append(0.0)
                    else:
                        ratios.append(gene_to_trans_reads[gene_id][trans][hap] / tot_cnt[hap])
                p = get_fisher_test_p_value(cnts)
                out_str = [cluster, trans, p, gene_id, gene_name, fdr, 
                        gene_to_trans_reads[gene_id][trans][ALL_HAPS[0]], 
                        gene_to_trans_reads[gene_id][trans][ALL_HAPS[1]], ratios[0], ratios[1]]
                as_hap = ALL_HAPS[0] if ratios[0] > ratios[1] else ALL_HAPS[1]

                if abs(ratios[0] - ratios[1]) <= delta_ratio or math.isnan(p) or p > p_thres:
                    if gene_id in gene_to_as_info:
                        gene_to_as_info[gene_id].add_non_as_transcript_cell_type(cluster, trans)
                    if inc_non_sign_trans:
                        asts_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
                    continue
                if gene_id in gene_to_as_info:
                    gene_to_as_info[gene_id].add_as_transcript_cell_type(cluster, trans, as_hap)
                # write significant allele-specific spliced transcripts to output file
                asts_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))

# phased count filtering criteria for allele-specific splicing:
# 1. For gene: at least 2 haplotype-resolved transcripts for each gene
# 2. For transcript: at least 10 phased reads for at least one haplotype for each transcripts
# 3. For transcript: at least 80% of reads are phased for each transcripts
# 4. For gene: overall used-to-total ratio of each gene is at least 80%, 
#    here `used` means the reads used for haplotype-resolved transcripts,
# .  excluding the transcripts not having >= 10 phased reads in at least one haplotype
def cluster_wise_allele_specific_splicing(cluster_wise_gene_to_trans_reads, gene_id_to_name,
                                          gene_detailed_out, trans_detailed_out,
                                          inc_non_sign_gene, inc_non_sign_trans,
                                          # ASS
                                          min_trans_phased_ratio=0.8, min_read_cnt=5,
                                          min_gene_phased_ratio=0.8, min_gene_total_phase_cnt=20, min_trans_cnt=2):
    gene_to_as_info = dd(lambda: None)
    with open(gene_detailed_out, 'w') as gene_detailed_fp, open(trans_detailed_out, 'w') as trans_detailed_fp:
        gene_detailed_fp.write('\t'.join(['cell_type', 'gene_id', 'gene_name', 'gene_fdr', 'transcripts', 
                                          'hap1_counts', 'hap2_counts', 'hap1_ratios', 'hap2_ratios']) + '\n')
        trans_detailed_fp.write('\t'.join(['cell_type', 'transcript_id', 'transcript_pval', 'gene_id', 'gene_name', 'gene_fdr', 
                                           'hap1_count', 'hap2_count', 'hap1_ratio', 'hap2_ratio']) + '\n')
        # chi-square test on genes for allele-specific splicing
        as_gene_list, as_p_list, as_gene_to_max_ratio = [], [], dict()
        # collect valid trans/gene for each cluster
        for cluster, gene_to_trans_reads in cluster_wise_gene_to_trans_reads.items():
            # only use isos with read count >= min_read_cnt in at least one allele
            # only use genes with trans >= min_trans_cnt
            genes = list(gene_to_trans_reads.keys())
            for gene in genes:
                trans_to_hap_cnt = gene_to_trans_reads[gene]
                transs = list(trans_to_hap_cnt.keys())
                if len(transs) < min_trans_cnt: # gene with <2 transcripts
                    del gene_to_trans_reads[gene]
                    continue
                gene_total_cnt = 0 # gene with gene_total_cnt < 20: None (low expression)
                gene_total_phased_cnt = 0 # total phased cnt of all trans that have >= 10 phased reads in at least one haplotype
                used_gene_phased_cnt = 0
                for trans in transs:
                    hap_to_cnt = trans_to_hap_cnt[trans]
                    filtered_trans = True
                    phased_cnt = 0
                    for hap, hap_cnt in hap_to_cnt.items():
                        gene_total_cnt += hap_cnt
                        if hap not in ALL_HAPS:
                            continue
                        phased_cnt += hap_cnt
                        if hap_cnt >= min_read_cnt:
                            filtered_trans = False
                    gene_total_phased_cnt += phased_cnt
                    # min_phased_trans_per_hap < 10 or transcript-wise phased_ratio < 0.8
                    if filtered_trans or phased_cnt / sum(hap_to_cnt.values()) < min_trans_phased_ratio:
                        del trans_to_hap_cnt[trans]
                    else:
                        used_gene_phased_cnt += phased_cnt
                if len(trans_to_hap_cnt) < min_trans_cnt or \
                   gene_total_phased_cnt / gene_total_cnt < min_gene_phased_ratio or \
                   used_gene_phased_cnt < min_gene_total_phase_cnt:
                    del gene_to_trans_reads[gene]
            for gene, trans_to_hap_cnt in gene_to_trans_reads.items():
                gene_ = (cluster, gene)
                cnts = []
                for hap in ALL_HAPS:
                    cnt = []
                    for trans, hap_to_cnt in trans_to_hap_cnt.items():
                        cnt.append(hap_to_cnt[hap])
                    cnts.append(cnt)
                # ratios of each iso
                rat1 = [x / sum(cnts[0]) if sum(cnts[0]) > 0 else 0 for x in cnts[0]]
                rat2 = [x / sum(cnts[1]) if sum(cnts[1]) > 0 else 0 for x in cnts[1]]
                as_gene_to_max_ratio[gene_] = max([abs(x-y) for x, y in zip(rat1, rat2)])
                p = get_chisq_test_p_value(cnts)
                as_gene_list.append(gene_)
                as_p_list.append(p)
            # within cluster FDR correction
            # comment following lines if across cluster FDR is applied
            # as_gene_fdr_corr_and_asts(as_gene_list, as_p_list, as_gene_to_max_ratio, min_trans_cnt, 
            #                           cluster_wise_gene_to_trans_reads, gene_id_to_name, 
            #                           gene_detailed_fp, trans_detailed_fp,
            #                           inc_non_sign_gene, inc_non_sign_trans, gene_to_as_info)
            # as_gene_list, as_p_list, as_gene_to_max_ratio = [], [], dict()
        # across cluster FDR correction
        # comment following lines if within cluster FDR is applied
        as_gene_fdr_corr_and_asts(as_gene_list, as_p_list, as_gene_to_max_ratio, min_trans_cnt, 
                                  cluster_wise_gene_to_trans_reads, gene_id_to_name,
                                  gene_detailed_fp, trans_detailed_fp, 
                                  inc_non_sign_gene, inc_non_sign_trans, gene_to_as_info)
    return gene_to_as_info

def output_allele_specific_splicing(all_cell_types, gene_to_allele_specific, assg_stat_out, asst_stat_out):
    gene_header = ['gene_name', 'gene_id']
    gene_header.extend(all_cell_types)
    trans_header = ['transcript_id', 'gene_name', 'gene_id']
    trans_header.extend(all_cell_types)
    with open(assg_stat_out, 'w') as gene_fp, open(asst_stat_out, 'w') as trans_fp:
        gene_fp.write('\t'.join(gene_header) + '\n')
        trans_fp.write('\t'.join(trans_header) + '\n')
        for gene in gene_to_allele_specific:
            # gene
            gene_as_info = gene_to_allele_specific[gene]
            gene_name = gene_as_info.gene_name
            gene_fp.write('\t'.join([gene_name, gene, '\t'.join([gene_as_info.as_gene_cate[cell_type] for cell_type in all_cell_types])]) + '\n')
            # transcript
            for trans, trans_as_info in gene_as_info.as_trans_cate.items():
                if ALL_HAPS[0] in trans_as_info.values() or ALL_HAPS[1] in trans_as_info.values():
                    trans_fp.write('\t'.join([trans, gene_name, gene, '\t'.join([trans_as_info[cell_type] for cell_type in all_cell_types])]) + '\n')

def get_gene_list(gene_list_file):
    gene_list = open(gene_list_file).read().splitlines()
    return gene_list

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: allele-specific splicing analysis".format(os.path.basename(__file__)))
    parser.add_argument('scrl_bu', metavar='bc_umi.tsv', type=str, help='Barcoded/UMI file, generated by NanoHunter')
    parser.add_argument('bc_cluster', metavar='bc_to_cluster.tsv', type=str, help='Barcode cell type/cluster file')
    parser.add_argument('hap_list', metavar='hap_list.tsv', type=str, help='Haplotype assignment of each read, generated by WhatsHap')
    parser.add_argument('out_pre', metavar='out_pre', type=str, help='Output prefix. Only significant allele-specific splicing genes/transcrpts will be output by default')

    # add optional arguments
    # group: phasing parameters
    phase_par = parser.add_argument_group('phasing options')
    phase_par.add_argument('-r', '--use-read-cnt', action='store_true', default=False, help='Use read count instead of UMI count')
    phase_par.add_argument('--read-type', type=str, default='', help='Read type to consider for phasing, i.e., FSM. Use comma to separate mutliple types. Default: all read types')
    phase_par.add_argument('--trans-phased-ratio', type=float, default=_trans_phased_ratio, 
                           help='Min. phase ratio of transcript')
    phase_par.add_argument('--trans-phased-cnt-per-allele', type=int, default=_trans_phased_cnt_per_allele, 
                           help='Min. phased UMI count of transcript in at least one allele')
    phase_par.add_argument('--gene-total-phased-cnt', type=int, default=_gene_total_phased_cnt, 
                           help='Min. total phased UMI count of gene')
    phase_par.add_argument('--gene-phased-ratio', type=float, default=_gene_phased_ratio, 
                           help='Min. total phase ratio of gene')
    phase_par.add_argument('--phased-trans-cnt-per-gene', type=int, default=_phased_trans_cnt_per_gene, 
                           help='Min. number of phased transcripts of gene')
    # output 
    out_par = parser.add_argument_group('output options')
    out_par.add_argument('-s', '--stat-out', action='store_true', default=False, help='Output additional allele-specific splicing result for all genes/transcripts')
    out_par.add_argument('-N', '--inc-non-sign-gene', action='store_true', default=False, help='Output both significant and non-significant allele-specific spliced genes')
    out_par.add_argument('-n', '--inc-non-sign-trans', action='store_true', default=False, help='Output both significant and non-significant allele-specific spliced transcripts for each significant gene')
    
    return parser.parse_args()

def main():
    args = parser_argv()

    scrl_bu_fn, bc_to_clu_tsv, hap_list, out_pre = args.scrl_bu, args.bc_cluster, args.hap_list, args.out_pre
    use_read_cnt = args.use_read_cnt
    read_types = args.read_type.rsplit(',')
    if '' in read_types:
        read_types.remove('')
    trans_phased_ratio, trans_phased_cnt_per_allele = args.trans_phased_ratio, args.trans_phased_cnt_per_allele
    gene_total_phased_cnt, gene_phased_ratio, phased_trans_cnt_per_gene = args.gene_total_phased_cnt, args.gene_phased_ratio, args.phased_trans_cnt_per_gene
    stat_out, inc_non_sign_gene, inc_non_sign_trans = args.stat_out, args.inc_non_sign_gene, args.inc_non_sign_trans

    assg_detailed_out = out_pre + '_ASS_genes.tsv'
    asst_detailed_out = out_pre + '_ASS_transcripts.tsv'

    bc_to_cluster = parse_bc_cell_list(bc_to_clu_tsv)
    all_cell_types = set(bc_to_cluster.values())
    reads_to_hap = get_reads_to_hap(hap_list)
    cluster_wise_gene_to_trans_reads, gene_id_to_name = get_cluster_wise_gene_to_trans_reads(scrl_bu_fn, use_read_cnt, read_types, 
                                                                                             bc_to_cluster, reads_to_hap)
    # allele-specific spliced gene/transcript
    gene_to_allele_specific = cluster_wise_allele_specific_splicing(cluster_wise_gene_to_trans_reads, gene_id_to_name, # input
                                                                    assg_detailed_out, asst_detailed_out, # output
                                                                    inc_non_sign_gene, inc_non_sign_trans, # output
                                                                    trans_phased_ratio, trans_phased_cnt_per_allele, # parameters
                                                                    gene_phased_ratio, gene_total_phased_cnt, phased_trans_cnt_per_gene)

    # write stat of all genes/transcripts to file
    if stat_out:
        assg_stat_out = out_pre + '_all_genes_stat.tsv'
        asst_stat_out = out_pre + '_all_transcripts_stat.tsv'
        output_allele_specific_splicing(all_cell_types, gene_to_allele_specific, assg_stat_out, asst_stat_out)

if __name__ == '__main__':
    main()