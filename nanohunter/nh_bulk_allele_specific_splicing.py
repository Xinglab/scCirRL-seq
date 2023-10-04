import sys
import os
import argparse
from collections import defaultdict as dd
# R stats
from rpy2 import robjects as rb
from rpy2.robjects import r
from rpy2.robjects.packages import importr
stats = importr('stats')
import math

fdr_thres = 0.05
p_thres = 0.05
delta_ratio = 0.05

ALL_HAPS = ['H1', 'H2']  # 'none'

# gene wise allele specific info: 
#   as_gene -> cell types
#   as_transcript -> cell types
class gene_wise_allele_specific_info:
    def __init__(self, gene, gene_name):
        # 5 categories:
        # LowExp: low expression, skipped in allele-specific splicing analysis
        # Unphased: overall phased ratio < 0.8, skipped in allele-specific splicing analysis
        # NoSplice: only one haplotype-resolved transcript, skipped in allele-specific splicing analysis
        # True: allele-specific splicing
        # False: non-allele-specific splicing
        self.as_gene_cate= dd(lambda: "LowExp") # {sample}, LowExp: low expression
        self.as_trans_cate = dd(lambda: dd(lambda: "None")) # {transcript: {sample: hap}}
        self.gene = gene
        self.gene_name = gene_name
        self.chrom = ''

    def add_as_for_sample(self, sample):
        self.as_gene_cate[sample] = 'True'
    def add_non_as_for_sample(self, sample):
        self.as_gene_cate[sample] = 'False'
    def add_unphased_for_sample(self, sample):
        self.as_gene_cate[sample] = 'Unphased'
    def add_nosplice_for_sample(self, sample):
        self.as_gene_cate[sample] = 'NoSplice'
        
    def add_as_transcript_for_sample(self, sample, trans, hap):
        self.as_trans_cate[trans][sample] = hap
    def add_non_as_transcript_for_sample(self, sample, trans):
        self.as_trans_cate[trans][sample] = 'False'
    def set_chrom(self, chrom):
        self.chrom = chrom

    def get_as_gene_to_cate(self):
        return self.as_gene_cate
    def get_as_trans_to_cate(self):
        return self.as_trans_cate
    def get_all_as_trans_hap(self):
        trans_hap = set()
        for trans in self.as_trans_cate:
            for cell_type, hap in self.as_trans_cate[trans].items():
                trans_hap.add((trans, hap))
        return trans_hap

def get_chisq_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.chisq_test(df)
    return list(out[2])[0]

def get_fisher_test_p_value(cnt=[]):
    cnt1, cnt2 = cnt[0], cnt[1]
    r1 = rb.IntVector(cnt1)
    r2 = rb.IntVector(cnt2)
    df = r.rbind(r1, r2)
    out = stats.fisher_test(df)
    return list(out[0])[0]

def get_fdr(ps=[]):
    rp = rb.FloatVector(ps)
    return list(stats.p_adjust(rp, method='fdr'))

def get_reads_to_hap(hap_list):
    read_to_hap = dd(lambda: '')
    for fn in hap_list.rsplit(','):
        with open(fn) as fp:
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

def get_sample_wise_gene_to_trans_reads(esp_cmpt_tsvs, trans_to_gene_id, read_to_hap):
    sample_wise_gene_to_trans_reads = dd(lambda: dd(lambda: dd(lambda: dd(lambda: 0))))  # {sample: {gene: {trans: {hap: cnt}}}}
    for fn in esp_cmpt_tsvs.rsplit(','):
        with open(fn) as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                read_name, sample, cate, trans_ids = line.rstrip().rsplit('\t')
                trans_id = trans_ids.rsplit(',')
                if '' in trans_id:
                    trans_id.remove('')
                if len(trans_id) > 1 or trans_id[0] == 'NA':  # only keep single-match reads
                    continue
                trans_id = trans_id[0]
                if trans_id not in trans_to_gene_id:
                    continue
                gene_id = trans_to_gene_id[trans_id]
                if read_name not in read_to_hap:
                    hap = 'none'
                else:
                    hap = read_to_hap[read_name]
                sample_wise_gene_to_trans_reads[sample][gene_id][trans_id][hap] += 1
    return sample_wise_gene_to_trans_reads

# True/False: gene_total_count >= 10, n_trans >= 2, phased_count_per_trans >= 10, phased_ratio_per_trans >= 0.8
# Unable to phase (NA):  gene_total_count >= 20, but not pass the above criteria
# Low expression (None): gene_total_count < 20
def set_as_gene(gene_to_as_info, gene, gene_name, sample, as_type):
    if gene not in gene_to_as_info:
        gene_to_as_info[gene] = gene_wise_allele_specific_info(gene, gene_name)
    if as_type == 'True':
        gene_to_as_info[gene].add_as_for_sample(sample)
    elif as_type == 'False':
        gene_to_as_info[gene].add_non_as_for_sample(sample)
    elif as_type == 'Unphased':
        gene_to_as_info[gene].add_unphased_for_sample(sample)
    elif as_type == 'NoSplice':
        gene_to_as_info[gene].add_nosplice_for_sample(sample)

def gene_fdr_corr_and_asts(gene_list, p_list, 
                           gene_to_max_ratio, min_trans_cnt, 
                           sample_wise_gene_to_trans_reads, gene_id_to_name,
                           asg_detailed_fp, asts_detailed_fp, gene_to_as_info):
    gene_fdrs = get_fdr(p_list)
    for gene_, fdr, p in zip(gene_list, gene_fdrs, p_list):
        sample, gene = gene_.rsplit('+')
        gene_name = gene_id_to_name[gene]
        gene_to_trans_reads = sample_wise_gene_to_trans_reads[sample]
        transs = gene_to_trans_reads[gene].keys()
        hap1_cnts = [gene_to_trans_reads[gene][trans][ALL_HAPS[0]] for trans in transs]
        hap2_cnts = [gene_to_trans_reads[gene][trans][ALL_HAPS[1]] for trans in transs]
        hap1_ratios = [cnt / sum(hap1_cnts) if sum(hap1_cnts) > 0 else 0 for cnt in hap1_cnts]
        hap2_ratios = [cnt / sum(hap2_cnts) if sum(hap2_cnts) > 0 else 0 for cnt in hap2_cnts]
        out_str = [sample, gene, gene_id_to_name[gene], fdr, 
                   ','.join(transs), 
                   ','.join(list(map(str, hap1_cnts))), 
                   ','.join(list(map(str, hap2_cnts))), 
                   ','.join(list(map(str, hap1_ratios))), 
                   ','.join(list(map(str, hap2_ratios)))]

        if math.isnan(fdr):
            # write nan fdr, where gene is allele-specific-expressed
            asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            continue
        if len(gene_to_trans_reads[gene]) < min_trans_cnt:   # genes with < 2 haplotype-resolved transcripts
            set_as_gene(gene_to_as_info, gene, gene_name, sample, 'NoSplice')
            continue
        elif fdr > fdr_thres or gene_to_max_ratio[gene_] < delta_ratio:
            set_as_gene(gene_to_as_info, gene, gene_name, sample, 'False')
            continue

        asg_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
        set_as_gene(gene_to_as_info, gene, gene_name, sample, 'True')

        # fisher exact test
        tot_cnt = dd(lambda: 0.0)  # {hap: tot_cnt of all trans}
        for hap in ALL_HAPS:
            for trans in gene_to_trans_reads[gene]:
                tot_cnt[hap] += gene_to_trans_reads[gene][trans][hap]
        for trans in gene_to_trans_reads[gene]:
            cnts, ratios = [], []
            for hap in ALL_HAPS:
                cnt = [gene_to_trans_reads[gene][trans][hap], tot_cnt[hap] - gene_to_trans_reads[gene][trans][hap]]
                cnts.append(cnt)
                if tot_cnt[hap] == 0:
                    ratios.append(0.0)
                else:
                    ratios.append(gene_to_trans_reads[gene][trans][hap] / tot_cnt[hap])
            p = get_fisher_test_p_value(cnts)
            if abs(ratios[0] - ratios[1]) < delta_ratio or math.isnan(p) or p > p_thres:
                gene_to_as_info[gene].add_non_as_transcript_for_sample(sample, trans)
                continue

            out_str = [sample, trans, p, gene, gene_id_to_name[gene], fdr, 
                       gene_to_trans_reads[gene][trans][ALL_HAPS[0]], 
                       gene_to_trans_reads[gene][trans][ALL_HAPS[1]], ratios[0], ratios[1]]
            as_hap = ALL_HAPS[0] if ratios[0] > ratios[1] else ALL_HAPS[1]
            asts_detailed_fp.write('{}\n'.format('\t'.join(list(map(str, out_str)))))
            gene_to_as_info[gene].add_as_transcript_for_sample(sample, trans, as_hap)

def is_ae_gene(trans_to_hap_cnt, min_ae_max_hap_ratio, min_ae_gene_phased_ratio, min_ae_gene_phased_cnt):
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

def write_ae_gene(sample, gene, trans_to_hap_cnt, gene_id_to_name, aseg_basic_fp):
    hap1_cnts = [trans_to_hap_cnt[trans][ALL_HAPS[0]] for trans in trans_to_hap_cnt]
    hap2_cnts = [trans_to_hap_cnt[trans][ALL_HAPS[1]] for trans in trans_to_hap_cnt]
    non_hap_cnts = [trans_to_hap_cnt[trans]['none'] for trans in trans_to_hap_cnt] 
    hap1_ratios = [cnt / sum(hap1_cnts) if sum(hap1_cnts) > 0 else 0 for cnt in hap1_cnts]
    hap2_ratios = [cnt / sum(hap2_cnts) if sum(hap2_cnts) > 0 else 0 for cnt in hap2_cnts]
    non_hap_ratios = [cnt / sum(non_hap_cnts) if sum(non_hap_cnts) > 0 else 0 for cnt in non_hap_cnts]
    total_cnts = [sum(cnts) for cnts in [hap1_cnts, hap2_cnts, non_hap_cnts]]
    total_ratios = [cnt / sum(total_cnts) if sum(total_cnts) > 0 else 0 for cnt in total_cnts]
    aseg_basic_fp.write('\t'.join([sample, gene, gene_id_to_name[gene], ','.join(trans_to_hap_cnt.keys()),
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
# phased count filtering criteria for allele-specific expression: (gene only expressed in one haplotype)
# 1. For gene: at least 1 haplotype-resolved transcripts for each gene
# 2. For gene: overall phased ratio is 100%, i.e. all reads are phased
# 3. For gene: total phased reads is at least 20
def sample_wise_allele_specific_analysis(sample_wise_gene_to_trans_reads, gene_id_to_name, 
                                          gene_detailed_out, trans_detailed_out, aseg_basic_out,
                                          min_trans_phased_ratio=0.8, min_gene_phased_ratio=0.8,
                                          min_ae_max_hap_ratio=0.8, min_ae_gene_phased_ratio=0.8, min_ae_gene_phased_cnt=20,
                                          min_trans_cnt=2, min_read_cnt=10, low_exprestion_thres=20):
    gene_to_as_info = dd(lambda: None)
    with open(gene_detailed_out, 'w') as gene_detailed_fp, open(trans_detailed_out, 'w') as trans_detailed_fp, \
         open(aseg_basic_out, 'w') as aseg_basic_fp:
        gene_detailed_fp.write('\t'.join(['sample', 'gene_id', 'gene_name', 'gene_fdr', 'transcripts', 
                                          'hap1_counts', 'hap2_counts', 'hap1_ratios', 'hap2_ratios']) + '\n')
        trans_detailed_fp.write('\t'.join(['sample', 'transcript_id', 'transcript_p', 'gene_id', 'gene_name', 'gene_fdr', 
                                           'hap1_count', 'hap2_count', 'hap1_ratio', 'hap2_ratio']) + '\n')
        aseg_basic_fp.write('\t'.join(['sample', 'gene_id', 'gene_name', 'transcripts',
                                       'hap1_counts', 'hap2_counts', 'non_hap_counts',
                                       'hap1_ratios', 'hap2_ratios', 'non_hap_ratios',
                                       'total_cnts', 'total_ratios']) + '\n')
        # chi-square test on genes
        gene_list, p_list, gene_to_max_ratio = [], [], dict()
        # collect valid trans/gene for each sample
        for sample, gene_to_trans_reads in sample_wise_gene_to_trans_reads.items():
            # only use isos with read count >= min_read_cnt in at least one allele
            # only use genes with trans >= min_trans_cnt
            genes = list(gene_to_trans_reads.keys())
            for gene in genes:
                trans_to_hap_cnt = gene_to_trans_reads[gene]
                if is_ae_gene(trans_to_hap_cnt, min_ae_max_hap_ratio, min_ae_gene_phased_ratio, min_ae_gene_phased_cnt):
                    write_ae_gene(sample, gene, trans_to_hap_cnt, gene_id_to_name, aseg_basic_fp)
                transs = list(trans_to_hap_cnt.keys())
                if len(transs) < min_trans_cnt:
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
                    if filtered_trans or phased_cnt / sum(hap_to_cnt.values()) < min_trans_phased_ratio:
                        del trans_to_hap_cnt[trans]
                    #else:
                        #gene_total_phased_cnt += phased_cnt
                print(gene, gene_id_to_name[gene], trans_to_hap_cnt, gene_total_cnt, gene_total_phased_cnt)
                if gene_total_cnt >= low_exprestion_thres:
                    if gene_total_cnt == 0 or gene_total_phased_cnt / gene_total_cnt < min_gene_phased_ratio:
                        set_as_gene(gene_to_as_info, gene, gene_id_to_name[gene], sample, "Unphased")
                        del gene_to_trans_reads[gene]
                    elif len(trans_to_hap_cnt) < min_trans_cnt:
                        set_as_gene(gene_to_as_info, gene, gene_id_to_name[gene], sample, "NoSplice")
                        del gene_to_trans_reads[gene]
                else:
                    del gene_to_trans_reads[gene]
            for gene, trans_to_hap_cnt in gene_to_trans_reads.items():
                gene_ = sample + '+' + gene
                cnts = []
                for hap in ALL_HAPS:
                    cnt = []
                    for trans, hap_to_cnt in trans_to_hap_cnt.items():
                        cnt.append(hap_to_cnt[hap])
                    cnts.append(cnt)
                # ratios of each iso
                rat1 = [x / sum(cnts[0]) if sum(cnts[0]) > 0 else 0 for x in cnts[0]]
                rat2 = [x / sum(cnts[1]) if sum(cnts[1]) > 0 else 0 for x in cnts[1]]
                gene_to_max_ratio[gene_] = max([abs(x-y) for x, y in zip(rat1, rat2)])
                p = get_chisq_test_p_value(cnts)
                gene_list.append(gene_)
                p_list.append(p)
            # within sample FDR correction
            # comment following 2 lines if across sample FDR is applied
            # gene_fdr_corr_and_asts(gene_list, p_list, gene_to_max_ratio, min_trans_cnt, 
            #                        sample_wise_gene_to_trans_reads, gene_id_to_name, 
            #                        gene_detailed_fp, trans_detailed_fp, gene_to_as_info)
            # gene_list, p_list, gene_to_max_ratio = [], [], dict()
        # comment following 1 line if within sample FDR is applied
        gene_fdr_corr_and_asts(gene_list, p_list, gene_to_max_ratio, min_trans_cnt, 
                               sample_wise_gene_to_trans_reads, gene_id_to_name,
                               gene_detailed_fp, trans_detailed_fp, gene_to_as_info)
    return gene_to_as_info

def output_allele_specific_splicing(all_samples, gene_to_allele_specific, assg_basic_out, asst_basic_out):
    gene_header = ['gene_name', 'gene_id']
    gene_header.extend(all_samples)
    trans_header = ['transcript_id', 'gene_name', 'gene_id']
    trans_header.extend(all_samples)
    with open(assg_basic_out, 'w') as gene_fp, open(asst_basic_out, 'w') as trans_fp:
        gene_fp.write('\t'.join(gene_header) + '\n')
        trans_fp.write('\t'.join(trans_header) + '\n')
        for gene in gene_to_allele_specific:
            # gene
            gene_as_info = gene_to_allele_specific[gene]
            gene_name = gene_as_info.gene_name
            gene_fp.write('\t'.join([gene_name, gene, '\t'.join([gene_as_info.as_gene_cate[sample] for sample in all_samples])]) + '\n')
            # transcript
            for trans, trans_as_info in gene_as_info.as_trans_cate.items():
                if ALL_HAPS[0] in trans_as_info.values() or ALL_HAPS[1] in trans_as_info.values():
                    trans_fp.write('\t'.join([trans, gene_name, gene, '\t'.join([trans_as_info[sample] for sample in all_samples])]) + '\n')


def get_trans_to_gene_id(gtf, ref_gtf):
    trans_to_gene_id = dict()
    gene_id_to_name = dict()
    with open(gtf) as fp, open(ref_gtf) as ref_fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rstrip().rsplit('\t')
            # search for transcript_id and gene_id
            trans_id, gene_id, gene_name = '', '', ''
            if ele[2] == 'transcript':
                if 'transcript_id' in line:
                    trnas_ids = line[15+line.index('transcript_id'):]
                    trans_id = trnas_ids[:trnas_ids.index('"')]
                if 'gene_id' in line:
                    gene_ids = line[9+line.index('gene_id'):]
                    gene_id = gene_ids[:gene_ids.index('"')]
                if 'gene_name' in line:
                    gene_name = line[11+line.index('gene_name'):]
                    gene_name = gene_name[:gene_name.index('"')]
                if trans_id == '' or gene_id == '' or trans_id == 'NA' or gene_id == 'NA' or ',' in trans_id or ',' in gene_id:
                    continue
                trans_to_gene_id[trans_id] = gene_id
                if gene_name == '' or gene_name == 'NA' or ',' in gene_name:
                    continue
                gene_id_to_name[gene_id] = gene_name
        for line in ref_fp:
            if line.startswith('#'):
                continue
            ele = line.rstrip().rsplit('\t')
            # search for transcript_id and gene_id
            gene_id, gene_name = '', ''
            if ele[2] == 'transcript':
                if 'gene_id' in line:
                    gene_ids = line[9+line.index('gene_id'):]
                    gene_id = gene_ids[:gene_ids.index('"')]
                if 'gene_name' in line:
                    gene_name = line[11+line.index('gene_name'):]
                    gene_name = gene_name[:gene_name.index('"')]
                if gene_id in gene_id_to_name or gene_id == '' or gene_name == '':
                    continue
                gene_id_to_name[gene_id] = gene_name
    return trans_to_gene_id, gene_id_to_name

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: allele-specific splicing analysis".format(os.path.basename(__file__)))
    parser.add_argument('esp_cmpt', metavar='esp_cmpt.tsv', help='Read-isoform compatible file, generated by ESPRESSO. Use \',\' to separate mutliple input files', type=str)
    parser.add_argument('gtf', metavar='updated.gtf', type=str, help='Updated gtf file, generated by ESPRESSO')
    parser.add_argument('ref_gtf', metavar='ref.gtf', type=str, help='Reference annotation gtf file')
    parser.add_argument('hap_list', metavar='hap_list.tsv', type=str, help='Haplotype assignment of each read, generated by WhatsHap. Use \',\' to separate mutliple input files')
    parser.add_argument('out_pre', metavar='out_pre', type=str, help='Output prefix')
    
    return parser.parse_args()

def main():
    args = parser_argv()

    esp_cmpt_tsv, gtf, ref_gtf, hap_list, out_pre = args.esp_cmpt, args.gtf, args.ref_gtf, args.hap_list, args.out_pre

    assg_detailed_out = out_pre + '_allele_spliced_gene.tsv'
    asst_detailed_out = out_pre + '_allele_spliced_transcript.tsv'
    assg_basic_out = out_pre + '_allele_spliced_genes_basic.tsv'
    asst_basic_out = out_pre + '_allele_spliced_transcripts_basic.tsv'
    aseg_basic_out = out_pre + '_allele_expressed_gene.tsv'

    reads_to_hap = get_reads_to_hap(hap_list)
    trans_to_gene_id, gene_id_to_name = get_trans_to_gene_id(gtf, ref_gtf)
    sample_wise_gene_to_trans_reads = get_sample_wise_gene_to_trans_reads(esp_cmpt_tsv, trans_to_gene_id, reads_to_hap)
    # allele-specific spliced gene/transcript
    gene_to_allele_specific = sample_wise_allele_specific_analysis(sample_wise_gene_to_trans_reads, 
                                                                   gene_id_to_name,
                                                                   assg_detailed_out, 
                                                                   asst_detailed_out,
                                                                   aseg_basic_out)
    # write to file
    output_allele_specific_splicing(list(sample_wise_gene_to_trans_reads.keys()),  gene_to_allele_specific, assg_basic_out, asst_basic_out)

if __name__ == '__main__':
    main()