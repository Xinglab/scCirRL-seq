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

_trans_phased_ratio = 0.8
_trans_phased_cnt_per_allele = 10
_gene_phased_ratio = 0.8
_phased_trans_cnt_per_gene = 2
# filtered_out_gene_prefix = 'HLA-'
_ase_gene_total_cnt = 10
_ase_gene_phased_ratio = 0.8
# _ase_gene_allele_ratio_thres = 0.7

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

folds=1
def get_chisq_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    cnt1 = [i*folds for i in cnt1]
    cnt2 = [i*folds for i in cnt2]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.chisq_test(df)
    return list(out[2])[0]

def get_fisher_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    cnt1 = [i*folds for i in cnt1]
    cnt2 = [i*folds for i in cnt2]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.fisher_test(df)
    return list(out[0])[0]

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
    for qname in qnames.rsplit(','):
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

def get_cluster_wise_gene_to_trans_reads(scrl_bu_fn, bc_to_cluster, read_to_hap):
    cluster_wise_gene_to_trans_reads = dd(lambda: dd(lambda: dd(lambda: dd(lambda: 0))))  # {bc: {gene: {trans: {hap: cnt}}}}
    gene_id_to_name = dict()
    cluster_wise_stats = dd(lambda: dd(lambda: dd(lambda: 0)))
    with open(scrl_bu_fn) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            bc, umi, read_cnt, trans_id, gene_id, gene_name, read_names = ele
            if trans_id == 'NA' or ',' in trans_id or gene_id == 'NA' or ',' in gene_id:  # only keep single-match reads
                continue
            # if gene_name.startswith(filtered_out_gene_prefix):
                # continue
            if bc not in bc_to_cluster:
                continue
            cluster = bc_to_cluster[bc]
            hap = get_hap(read_names, read_to_hap, trans_id, gene_name)  # H1/H2/none
            cluster_wise_gene_to_trans_reads[cluster][gene_id][trans_id][hap] += 1
            gene_id_to_name[gene_id] = gene_name

            cluster_wise_stats[cluster]['cells'][bc] = 1
            cluster_wise_stats[cluster]['UMI']['total'] += 1
            if hap != 'none':
                cluster_wise_stats[cluster]['UMI']['phased'] += 1
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

        # allele-specific-expressed genes
        if math.isnan(fdr):
            # XXX write nan fdr, where gene is allele-specific-expressed
            if inc_non_sign_gene:
                asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            continue
        if len(gene_to_trans_reads[gene_id]) < min_trans_cnt:   # genes with < 2 haplotype-resolved transcripts
            print(f'{gene_id}:\tNoSplice')
            if inc_non_sign_gene:
                asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            # set_as_gene(gene_to_as_info, gene_id, gene_name, cluster, 'NoSplice')
            continue
        elif fdr > fdr_thres or gene_to_max_ratio[gene_] <= delta_ratio:
            set_as_gene(gene_to_as_info, gene_id, gene_name, cluster, 'False')
            if inc_non_sign_gene:
                asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            continue

        # write significant allele-specific spliced genes to output file
        set_as_gene(gene_to_as_info, gene_id, gene_name, cluster, 'True')
        asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))

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
                gene_to_as_info[gene_id].add_non_as_transcript_cell_type(cluster, trans)
                if inc_non_sign_trans:
                    asts_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
                continue
            gene_to_as_info[gene_id].add_as_transcript_cell_type(cluster, trans, as_hap)
            # write significant allele-specific spliced transcripts to output file
            asts_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))

# phased count filtering criteria for allele-specific expression: (gene only expressed in one haplotype)
# 1. For gene: at least 1 haplotype-resolved transcripts for each gene
# 2. For gene: overall phased ratio is 100%, i.e. all reads are phased
# 3. For gene: total phased reads is at least 20
def is_ae_gene(trans_to_hap_cnt):
    min_ae_max_hap_ratio= 0.8
    min_ae_gene_phased_ratio= 0.8
    min_ae_gene_phased_cnt= 20
    hap1_cnt, hap2_cnt, total_cnt = 0, 0, 0
    for trans, hap_to_cnt in trans_to_hap_cnt.items():
        cnt1, cnt2 = hap_to_cnt[ALL_HAPS[0]], hap_to_cnt[ALL_HAPS[1]]
        hap1_cnt += cnt1
        hap2_cnt += cnt2
        total_cnt += sum(hap_to_cnt.values())
    if total_cnt <= 0:
        return False
    if hap1_cnt + hap2_cnt >= min_ae_gene_phased_cnt and \
        max(hap1_cnt, hap2_cnt)/(hap1_cnt + hap2_cnt) >= min_ae_max_hap_ratio and \
        (hap1_cnt+hap2_cnt)/total_cnt >= min_ae_gene_phased_ratio:
        return True
    return False

def get_gene_to_ae_p(trans_to_hap_cnt, ase_gene_total_cnt, ase_gene_phased_ratio):
    hap1_cnt, hap2_cnt, total_cnt = 0, 0, 0
    for trans, hap_to_cnt in trans_to_hap_cnt.items():
        cnt1, cnt2 = hap_to_cnt[ALL_HAPS[0]], hap_to_cnt[ALL_HAPS[1]]
        hap1_cnt += cnt1
        hap2_cnt += cnt2
        total_cnt += sum(hap_to_cnt.values())

    # only consider genes with phased ratio ≥ 0.8
    if total_cnt >= ase_gene_total_cnt and (hap1_cnt+hap2_cnt)/total_cnt >= ase_gene_phased_ratio:
        return get_binomial_test_p_value(max(hap1_cnt, hap2_cnt), hap1_cnt+hap2_cnt), max(hap1_cnt, hap2_cnt)/(hap1_cnt+hap2_cnt)
    else:
        return None, None

def ae_gene_fdr_corr(ae_gene_list, ae_p_list, ae_gene_to_max_ratio,
                     gene_to_trans_reads, gene_id_to_name, aseg_basic_fp):
    gene_fdrs = get_fdr(ae_p_list)
    for cluster_gene, fdr, p in zip(ae_gene_list, gene_fdrs, ae_p_list):
        cluster, gene_id = cluster_gene[0], cluster_gene[1]
        max_ratio = ae_gene_to_max_ratio[cluster_gene]
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
            aseg_basic_fp.write('\t'.join([cluster, gene_id, gene_name, ','.join(trans_to_hap_cnt.keys()),
                                ','.join(list(map(str, hap1_cnts))),
                                ','.join(list(map(str, hap2_cnts))),
                                ','.join(list(map(str, non_hap_cnts))),
                                ','.join(list(map(str, hap1_ratios))),
                                ','.join(list(map(str, hap2_ratios))),
                                ','.join(list(map(str, non_hap_ratios))),
                                ','.join(list(map(str, total_cnts))),
                                ','.join(list(map(str, total_ratios)))]) + '\n')
        
    
def write_ae_gene(cell_type, gene, trans_to_hap_cnt, gene_id_to_name, aseg_basic_fp):
    hap1_cnts = [trans_to_hap_cnt[trans][ALL_HAPS[0]] for trans in trans_to_hap_cnt]
    hap2_cnts = [trans_to_hap_cnt[trans][ALL_HAPS[1]] for trans in trans_to_hap_cnt]
    non_hap_cnts = [trans_to_hap_cnt[trans]['none'] for trans in trans_to_hap_cnt] 
    hap1_ratios = [cnt / sum(hap1_cnts) if sum(hap1_cnts) > 0 else 0 for cnt in hap1_cnts]
    hap2_ratios = [cnt / sum(hap2_cnts) if sum(hap2_cnts) > 0 else 0 for cnt in hap2_cnts]
    non_hap_ratios = [cnt / sum(non_hap_cnts) if sum(non_hap_cnts) > 0 else 0 for cnt in non_hap_cnts]
    total_cnts = [sum(cnts) for cnts in [hap1_cnts, hap2_cnts, non_hap_cnts]]
    total_ratios = [cnt / sum(total_cnts) if sum(total_cnts) > 0 else 0 for cnt in total_cnts]
    aseg_basic_fp.write('\t'.join([cell_type, gene, gene_id_to_name[gene], ','.join(trans_to_hap_cnt.keys()),
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
def cluster_wise_allele_specific_analysis(cluster_wise_gene_to_trans_reads, gene_id_to_name,
                                          gene_detailed_out, trans_detailed_out, aseg_detailed_out,
                                          inc_non_sign_gene, inc_non_sign_trans,
                                          min_trans_phased_ratio=0.8, 
                                          min_read_cnt=10,
                                          min_gene_phased_ratio=0.8,
                                          ase_gene_total_cnt=_ase_gene_total_cnt, 
                                          ase_gene_phased_ratio=_ase_gene_phased_ratio,
                                          min_trans_cnt=2):
    gene_to_as_info = dd(lambda: None)
    with open(gene_detailed_out, 'w') as gene_detailed_fp, open(trans_detailed_out, 'w') as trans_detailed_fp, \
         open(aseg_detailed_out, 'w') as aseg_detailed_fp:
        gene_detailed_fp.write('\t'.join(['cell_type', 'gene_id', 'gene_name', 'gene_fdr', 'transcripts', 
                                          'hap1_counts', 'hap2_counts', 'hap1_ratios', 'hap2_ratios']) + '\n')
        trans_detailed_fp.write('\t'.join(['cell_type', 'transcript_id', 'transcript_pval', 'gene_id', 'gene_name', 'gene_fdr', 
                                           'hap1_count', 'hap2_count', 'hap1_ratio', 'hap2_ratio']) + '\n')
        aseg_detailed_fp.write('\t'.join(['cell_type', 'gene_id', 'gene_name', 'transcripts',
                                       'hap1_counts', 'hap2_counts', 'non_hap_counts',
                                       'hap1_ratios', 'hap2_ratios', 'non_hap_ratios',
                                       'total_cnts', 'total_ratios']) + '\n')
        # chi-square test on genes for allele-specific splicing
        as_gene_list, as_p_list, as_gene_to_max_ratio = [], [], dict()
        # binomial test on genes for allele-specific expression
        ae_gene_list, ae_p_list, ae_gene_to_max_ratio = [], [], dict()
        # collect valid trans/gene for each cluster
        for cluster, gene_to_trans_reads in cluster_wise_gene_to_trans_reads.items():
            # only use isos with read count >= min_read_cnt in at least one allele
            # only use genes with trans >= min_trans_cnt
            genes = list(gene_to_trans_reads.keys())
            for gene in genes:
                if gene == 'ENSG00000089127' and cluster == 'Mono':
                    print('ok')
                trans_to_hap_cnt = gene_to_trans_reads[gene]
                ae_p, ae_max_ratio = get_gene_to_ae_p(trans_to_hap_cnt, ase_gene_total_cnt, ase_gene_phased_ratio)
                if ae_p != None:
                    ae_gene_list.append((cluster, gene))
                    ae_p_list.append(ae_p)
                    ae_gene_to_max_ratio.append(ae_max_ratio)
                # if is_ae_gene(trans_to_hap_cnt):
                    # write_ae_gene(cluster, gene, trans_to_hap_cnt, gene_id_to_name, aseg_basic_fp)
                transs = list(trans_to_hap_cnt.keys())
                if len(transs) < min_trans_cnt: # gene with <2 transcripts
                    del gene_to_trans_reads[gene]
                    continue
                gene_total_cnt = 0 # gene with gene_total_cnt < 20: None (low expression)
                gene_total_phased_cnt = 0 # total phased cnt of all trans that have >= 10 phased reads in at least one haplotype
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
                if len(trans_to_hap_cnt) < min_trans_cnt or gene_total_phased_cnt / gene_total_cnt < min_gene_phased_ratio: # gene-wise phase ratio < 0.8
                    del gene_to_trans_reads[gene]
                    # set_as_gene(gene_to_as_info, gene, gene_id_to_name[gene], cluster, "Unphased")
                # elif len(trans_to_hap_cnt) < min_trans_cnt: # gene with <2 phased transcripts
                    # del gene_to_trans_reads[gene]
                    # set_as_gene(gene_to_as_info, gene, gene_id_to_name[gene], cluster, "NoSplice")
            for gene, trans_to_hap_cnt in gene_to_trans_reads.items():
                gene_ = (cluster, gene)
                if gene == 'ENSG00000089127' and cluster == 'Mono':
                    print('ok')
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
            # ae_gene_fdr_corr(ae_gene_list, ae_p_list, ae_gene_to_max_ratio, 
            #                  gene_to_trans_reads, gene_id_to_name, aseg_detailed_fp)
            # as_gene_fdr_corr_and_asts(as_gene_list, as_p_list, as_gene_to_max_ratio, min_trans_cnt, 
            #                           cluster_wise_gene_to_trans_reads, gene_id_to_name, 
            #                           gene_detailed_fp, trans_detailed_fp,
            #                           inc_non_sign_gene, inc_non_sign_trans, gene_to_as_info)
            # gene_list, p_list, gene_to_max_ratio = [], [], dict()
        # comment following lines if within cluster FDR is applied
        ae_gene_fdr_corr(ae_gene_list, ae_p_list, ae_gene_to_max_ratio, 
                         gene_to_trans_reads, gene_id_to_name, aseg_detailed_fp)
        as_gene_fdr_corr_and_asts(as_gene_list, as_p_list, as_gene_to_max_ratio, min_trans_cnt, 
                                  cluster_wise_gene_to_trans_reads, gene_id_to_name,
                                  gene_detailed_fp, trans_detailed_fp, 
                                  inc_non_sign_gene, inc_non_sign_trans, 
                                  gene_to_as_info)
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
    phase_par.add_argument('--trans-phased-ratio', type=float, default=_trans_phased_ratio, 
                           help='Min. phase ratio of transcript')
    phase_par.add_argument('--trans-phased-cnt-per-allele', type=int, default=_trans_phased_cnt_per_allele, 
                           help='Min. phased UMI count of transcript in at least one allele')
    phase_par.add_argument('--gene-phased-ratio', type=float, default=_gene_phased_ratio, 
                           help='Min. total phase ratio of gene')
    phase_par.add_argument('--phased-trans-cnt-per-gene', type=int, default=_phased_trans_cnt_per_gene, 
                           help='Min. number of phased transcripts of gene')
    ase_gene_par = parser.add_argument_group('allele-specific gene expression options')
    ase_gene_par.add_argument('--ase-gene-total-cnt', type=int, default=_ase_gene_total_cnt,
                              help='Min. total UMI count of gene to consider for allele-specific expression')
    ase_gene_par.add_argument('--ase-gene-phased-ratio', type=float, default=_ase_gene_phased_ratio,
                              help='Min. phased ratio of gene to consider for allele-specific expression')
    # output 
    out_par = parser.add_argument_group('output options')
    out_par.add_argument('-s', '--stat-out', action='store_true', default=False, help='Output additional allele-specific splicing result for all genes/transcripts')
    out_par.add_argument('-N', '--inc-non-sign-gene', action='store_true', default=False, help='Output both significant and non-significant allele-specific spliced genes')
    out_par.add_argument('-n', '--inc-non-sign-trans', action='store_true', default=False, help='Output both significant and non-significant allele-specific spliced transcripts for each significant gene')
    
    return parser.parse_args()

def main():
    args = parser_argv()

    scrl_bu_fn, bc_to_clu_tsv, hap_list, out_pre = args.scrl_bu, args.bc_cluster, args.hap_list, args.out_pre
    trans_phased_ratio, trans_phased_cnt_per_allele = args.trans_phased_ratio, args.trans_phased_cnt_per_allele
    gene_phased_ratio, phased_trans_cnt_per_gene = args.gene_phased_ratio, args.phased_trans_cnt_per_gene
    ase_gene_total_cnt, ase_gene_phased_ratio = args.ase_gene_total_cnt, args.ase_gene_phased_ratio
    stat_out, inc_non_sign_gene, inc_non_sign_trans = args.stat_out, args.inc_non_sign_gene, args.inc_non_sign_trans

    aseg_detailed_out = out_pre + '_ASE_gene.tsv'
    assg_detailed_out = out_pre + '_ASS_genes.tsv'
    asst_detailed_out = out_pre + '_ASS_transcripts.tsv'

    bc_to_cluster = parse_bc_cell_list(bc_to_clu_tsv)
    all_cell_types = set(bc_to_cluster.values())
    reads_to_hap = get_reads_to_hap(hap_list)
    cluster_wise_gene_to_trans_reads, gene_id_to_name = get_cluster_wise_gene_to_trans_reads(scrl_bu_fn, bc_to_cluster, reads_to_hap)
    # allele-specific spliced gene/transcript
    gene_to_allele_specific = cluster_wise_allele_specific_analysis(cluster_wise_gene_to_trans_reads, gene_id_to_name, 
                                                                    assg_detailed_out, asst_detailed_out, aseg_detailed_out,
                                                                    inc_non_sign_gene, inc_non_sign_trans,
                                                                    trans_phased_ratio, trans_phased_cnt_per_allele, gene_phased_ratio, phased_trans_cnt_per_gene,
                                                                    ase_gene_total_cnt, ase_gene_phased_ratio)
    # write stat of all genes/transcripts to file
    if stat_out:
        assg_stat_out = out_pre + '_all_genes_stat.tsv'
        asst_stat_out = out_pre + '_all_transcripts_stat.tsv'
        output_allele_specific_splicing(all_cell_types, gene_to_allele_specific, assg_stat_out, asst_stat_out)

if __name__ == '__main__':
    main()