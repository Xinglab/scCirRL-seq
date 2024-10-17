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

# filtered_out_gene_prefix = 'HLA-'
_ase_gene_total_cnt = 20
_ase_gene_phased_ratio = 0.7
# _ase_gene_allele_ratio_thres = 0.7

ALL_HAPS = ['H1', 'H2']  # 'none'

def get_binomial_test_p_value(x, n, p=0.5, alternative='greater'):
    """
    Perform a binomial test using R's binom.test function via rpy2.

    Parameters:
        x: int
            The number of successes.
        n: int
            The number of trials.
        p: float, optional (default=0.5)
            The hypothesized probability of success.
        alternative: {'two.sided', 'less', 'greater'}, optional (default='two.sided')
            The alternative hypothesis. 'two.sided' (default) tests if the probability of success is different from p.
            'less' tests if the probability of success is less than p.
            'greater' tests if the probability of success is greater than p.

    Returns:
        p_value: float
            The p-value of the binomial test.
    """
    # Convert Python variables to R objects
    x_r = rb.IntVector([x])
    n_r = rb.IntVector([n])
    p_r = rb.FloatVector([p])

    # Import R's binom.test function
    r_binom_test = r['binom.test']
    # Perform the binomial test
    result = r_binom_test(x=x_r, n=n_r, p=p_r, alternative=alternative)
    # Extract the p-value from the result
    p_value = result.rx2('p.value')[0]
    return p_value

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

def get_gene_to_ae_p(trans_to_hap_cnt, ase_gene_total_cnt, ase_gene_phased_ratio):
    hap1_cnt, hap2_cnt, total_cnt = 0, 0, 0
    for trans, hap_to_cnt in trans_to_hap_cnt.items():
        cnt1, cnt2 = hap_to_cnt[ALL_HAPS[0]], hap_to_cnt[ALL_HAPS[1]]
        hap1_cnt += cnt1
        hap2_cnt += cnt2
        total_cnt += sum(hap_to_cnt.values())

    # only consider genes with phased ratio ≥ 0.8
    if total_cnt >= ase_gene_total_cnt and (hap1_cnt+hap2_cnt)/total_cnt >= ase_gene_phased_ratio:
        return get_binomial_test_p_value(max(hap1_cnt, hap2_cnt), hap1_cnt+hap2_cnt)
    else:
        return None

def ae_gene_fdr_corr(ae_gene_list, ae_p_list,
                     gene_to_trans_reads, gene_id_to_name, aseg_basic_fp):
    gene_fdrs = get_fdr(ae_p_list)
    for cluster_gene, fdr in zip(ae_gene_list, gene_fdrs):
        cluster, gene_id = cluster_gene[0], cluster_gene[1]
        if fdr <= fdr_thres: # and max_ratio >= ae_gene_ratio_thres:
            gene_name = gene_id_to_name[gene_id]
            trans_to_hap_cnt = gene_to_trans_reads[gene_id]
            hap1_cnts = [trans_to_hap_cnt[trans][ALL_HAPS[0]] for trans in trans_to_hap_cnt]
            hap2_cnts = [trans_to_hap_cnt[trans][ALL_HAPS[1]] for trans in trans_to_hap_cnt]
            non_hap_cnts = [trans_to_hap_cnt[trans]['none'] for trans in trans_to_hap_cnt] 
            hap1_ratios = [cnt / sum(hap1_cnts) if sum(hap1_cnts) > 0 else 0 for cnt in hap1_cnts]
            hap2_ratios = [cnt / sum(hap2_cnts) if sum(hap2_cnts) > 0 else 0 for cnt in hap2_cnts]
            non_hap_ratios = [cnt / sum(non_hap_cnts) if sum(non_hap_cnts) > 0 else 0 for cnt in non_hap_cnts]
            total_cnts = [sum(cnts) for cnts in [hap1_cnts, hap2_cnts, non_hap_cnts]]
            total_ratios = [cnt / sum(total_cnts) if sum(total_cnts) > 0 else 0 for cnt in total_cnts]
            aseg_basic_fp.write('\t'.join([cluster, gene_id, gene_name, str(fdr), ','.join(trans_to_hap_cnt.keys()),
                                ','.join(list(map(str, hap1_cnts))),
                                ','.join(list(map(str, hap2_cnts))),
                                ','.join(list(map(str, non_hap_cnts))),
                                ','.join(list(map(str, hap1_ratios))),
                                ','.join(list(map(str, hap2_ratios))),
                                ','.join(list(map(str, non_hap_ratios))),
                                ','.join(list(map(str, total_cnts))),
                                ','.join(list(map(str, total_ratios)))]) + '\n')

# phased count filtering criteria for allele-specific splicing:
# 1. For gene: at least 2 haplotype-resolved transcripts for each gene
# 2. For transcript: at least 10 phased reads for at least one haplotype for each transcripts
# 3. For transcript: at least 80% of reads are phased for each transcripts
# 4. For gene: overall used-to-total ratio of each gene is at least 80%, 
#    here `used` means the reads used for haplotype-resolved transcripts,
# .  excluding the transcripts not having >= 10 phased reads in at least one haplotype
def cluster_wise_allele_specific_expression(cluster_wise_gene_to_trans_reads, gene_id_to_name,
                                            aseg_detailed_out,
                                            ase_gene_total_cnt=_ase_gene_total_cnt, ase_gene_phased_ratio=_ase_gene_phased_ratio):
    with open(aseg_detailed_out, 'w') as aseg_detailed_fp:
        aseg_detailed_fp.write('\t'.join(['cell_type', 'gene_id', 'gene_name', 'gene_fdr', 'transcripts',
                                       'hap1_counts', 'hap2_counts', 'non_hap_counts',
                                       'hap1_ratios', 'hap2_ratios', 'non_hap_ratios',
                                       'total_cnts', 'total_ratios']) + '\n')
        # binomial test on genes for allele-specific expression
        ae_gene_list, ae_p_list = [], []
        # collect valid trans/gene for each cluster
        for cluster, gene_to_trans_reads in cluster_wise_gene_to_trans_reads.items():
            # only use isos with read count >= min_read_cnt in at least one allele
            # only use genes with trans >= min_trans_cnt
            genes = list(gene_to_trans_reads.keys())
            for gene in genes:
                trans_to_hap_cnt = gene_to_trans_reads[gene]
                ae_p = get_gene_to_ae_p(trans_to_hap_cnt, ase_gene_total_cnt, ase_gene_phased_ratio)
                if ae_p != None:
                    ae_gene_list.append((cluster, gene))
                    ae_p_list.append(ae_p)
            # within cluster FDR correction
            # comment following lines if across cluster FDR is applied
            # ae_gene_fdr_corr(ae_gene_list, ae_p_list, 
                            #  gene_to_trans_reads, gene_id_to_name, aseg_detailed_fp)
            # ae_gene_list, ae_p_list = [], []
        # across cluster FDR correction
        # comment following lines if within cluster FDR is applied
        ae_gene_fdr_corr(ae_gene_list, ae_p_list, 
                         gene_to_trans_reads, gene_id_to_name, aseg_detailed_fp)
    return

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
    phase_par.add_argument('--ase-gene-total-cnt', type=int, default=_ase_gene_total_cnt,
                              help='Min. total UMI count of gene to consider for allele-specific expression')
    phase_par.add_argument('--ase-gene-phased-ratio', type=float, default=_ase_gene_phased_ratio,
                              help='Min. phased ratio of gene to consider for allele-specific expression')
    return parser.parse_args()

def main():
    args = parser_argv()

    scrl_bu_fn, bc_to_clu_tsv, hap_list, out_pre = args.scrl_bu, args.bc_cluster, args.hap_list, args.out_pre
    use_read_cnt = args.use_read_cnt
    read_types = args.read_type.rsplit(',')
    if '' in read_types:
        read_types.remove('')
    ase_gene_total_cnt, ase_gene_phased_ratio = args.ase_gene_total_cnt, args.ase_gene_phased_ratio

    aseg_detailed_out = out_pre + '_ASE_genes.tsv'

    bc_to_cluster = parse_bc_cell_list(bc_to_clu_tsv)
    reads_to_hap = get_reads_to_hap(hap_list)
    cluster_wise_gene_to_trans_reads, gene_id_to_name = get_cluster_wise_gene_to_trans_reads(scrl_bu_fn, use_read_cnt, read_types, 
                                                                                             bc_to_cluster, reads_to_hap)
    # allele-specific expressed gene
    cluster_wise_allele_specific_expression(cluster_wise_gene_to_trans_reads, gene_id_to_name,  # input
                                            aseg_detailed_out, # output
                                            ase_gene_total_cnt, ase_gene_phased_ratio) # parameters
if __name__ == '__main__':
    main()