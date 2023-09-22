import sys
import argparse
from collections import defaultdict as dd
import gzip
import os
import itertools
import math
from rpy2 import robjects as rb
from rpy2.robjects import r
from rpy2.robjects.packages import importr

stats = importr('stats')

fdr_thres=0.05
p_thres = 0.05
delta_ratio = 0.05


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

def get_bc_from_barcode_tsv(barcode_tsv):
    bcs = dict()
    if barcode_tsv.endswith('gz'):
        fp = gzip.open(barcode_tsv, 'rt')
    else:
        fp = open(barcode_tsv)
    bc_i = 1
    for line in fp:
        bc = line.rstrip().rsplit('-')[0]
        bcs[bc_i] = bc
        bc_i += 1
    return bcs

def get_feature_from_feature_tsv(feature_tsv):
    features = dict()
    if feature_tsv.endswith('gz'):
        fp = gzip.open(feature_tsv, 'rt')
    else:
        fp = open(feature_tsv)
    feature_i = 1
    for line in fp:
        feature = line.rstrip().rsplit()[1]
        features[feature_i] = feature
        feature_i += 1
    return features

def get_feature_to_bc_cnt_from_sc_matrix(in_folder):
    sys.stderr.write("Loading feature-barcode count matrix from {} ...\n".format(in_folder))
    bcs, features = dict(), dict()
    feature_to_bc_cnt = dd(lambda: dd(lambda: 0.0))
    if os.path.exists(in_folder + '/barcodes.tsv.gz'):
        bc_fn = in_folder + '/barcodes.tsv.gz'
    elif os.path.exists(in_folder + '/barcodes.tsv'):
        bc_fn = in_folder + '/barcodes.tsv'
    else:
        sys.stderr.write('Cannot find \'barcodes.tsv.gz\' under \'{}\'\n'.format(in_folder))
        sys.exit(1)

    if os.path.exists(in_folder + '/features.tsv.gz'):
        feature_fn = in_folder + '/features.tsv.gz'
    elif os.path.exists(in_folder + '/features.tsv'):
        feature_fn = in_folder + '/features.tsv'
    else:
        sys.stderr.write('Cannot find \'features.tsv.gz\' under \'{}\'\n'.format(in_folder))
        sys.exit(1)

    if os.path.exists(in_folder + '/matrix.mtx.gz'):
        mtx_fn = in_folder + '/matrix.mtx.gz'
    elif os.path.exists(in_folder + '/matrix.mtx'):
        mtx_fn = in_folder + '/matrix.mtx'
    else:
        sys.stderr.write('Cannot find \'matrix.mtx.gz\' under \'{}\'\n'.format(in_folder))
        sys.exit(1)

    bcs = get_bc_from_barcode_tsv(bc_fn)
    features = get_feature_from_feature_tsv(feature_fn)

    if mtx_fn.endswith('gz'):
        fp = gzip.open(mtx_fn, 'rt')
    else:
        fp = open(mtx_fn)
    n_bc, n_features = len(bcs), len(features)
    header = True
    non_zeros, n_non_zeros = 0, 0
    for line in fp:
        if line.startswith('%'):
            continue
        ele = line.rsplit()
        if header:
            header = False
            _n_features, _n_bc, n_non_zeros = int(ele[0]), int(ele[1]), int(ele[2])
            if _n_features != n_features or _n_bc != n_bc:
                sys.stderr.write('{} ({}, {}) file does not match {}({})/{}({})\n'.format(mtx_fn, _n_features, _n_bc, feature_fn, n_features, bc_fn, n_bc))
                sys.exit(1)
            continue
        non_zeros += 1
        feature_i, bc_i, cnt = int(ele[0]), int(ele[1]), float(ele[2])
        if bc_i > n_bc or bc_i < 0 or feature_i > n_features or feature_i < 0:
            sys.stderr.write('Unexpected matrix numbers: {} {} {}\n'.format(feature_i, bc_i, cnt))
            sys.exit(1)
        bc, feature = bcs[bc_i], features[feature_i]
        feature_to_bc_cnt[feature][bc] += cnt
    if non_zeros != n_non_zeros:
        sys.stderr.write('{} ({}) file does not match itself({})\n'.format(mtx_fn, n_non_zeros, non_zeros))
        sys.exit(1)
    sys.stderr.write("Loading done\n")
    return feature_to_bc_cnt, bcs

def get_trans_to_gene_from_sc_matrix(in_folder):
    sys.stderr.write("Loading transcript-gene mapping from {} ...\n".format(in_folder))
    if os.path.exists(in_folder + '/gene_to_trans.tsv.gz'):
        in_fn = in_folder + '/gene_to_trans.tsv.gz'
    elif os.path.exists(in_folder + '/gene_to_trans.tsv'):
        in_fn = in_folder + '/gene_to_trans.tsv'
    else:
        sys.stderr.write('Cannot find \'gene_to_trans.tsv.gz\' under \'{}\'\n'.format(in_folder))
        sys.exit(1)
    if in_fn.endswith('gz'):
        fp = gzip.open(in_fn, 'rt')
    else:
        fp = open(in_fn)
    trans_to_gene = dict()
    header = True
    for line in fp:
        if header:
            header = False
            continue
        ele = line.rsplit()
        gene_id, gene_name, transcript = ele[0], ele[1], ele[2]
        if transcript in trans_to_gene:
            # transcript with multiple gene_id/gene_name
            continue
        trans_to_gene[transcript] = gene_id+','+gene_name
    sys.stderr.write("Loading done\n")
    return trans_to_gene

    
def get_bc_to_clu(in_fn):
    bc_to_clu = dict()
    with open(in_fn) as fp:
        for line in fp:
            ele = line.rstrip().rsplit('\t')
            if len(ele) == 2:
                bc_to_clu[ele[0].rsplit('-')[0]] = ele[1]
    return bc_to_clu

# perform single-cell-wise and cluster-wise filtering:
#   min_trans_cell_frac_per_type: 
#     minimum fraction of cells in a cluster that express a transcript with at least min_trans_cnt to keep it for all clusters
#     This is to filter out transcripts that are only expressed in a few cells in a cluster, which are likely to be noise
# . min_trans_total_cnt_per_clu:
# .   minimum total counts of a transcript in a cluster to keep it for all clusters
# . min_trans_cnt:
# .   minimum number of transcripts to keep a gene
# NOT filtering cells
def old_get_gene_trans_clu_cnt(bc_to_clu, in_trans_mtx_dir, 
                           min_trans_cell_frac_per_clu=0.00001, min_trans_cnt_per_cell=1,
                           min_trans_total_cnt_per_clu=2, min_trans_cnt=2):
    trans_to_bc_cnt, all_bcs = get_feature_to_bc_cnt_from_sc_matrix(in_trans_mtx_dir)
    trans_to_gene = get_trans_to_gene_from_sc_matrix(in_trans_mtx_dir)
    print('before filtering: {} genes, {} transcripts'.format(len(set(trans_to_gene.values())), len(trans_to_gene.keys())))
    gene_trans_clu_cnt = dd(lambda: dd(lambda: dd(lambda: 0.0)))
    clu_to_cell_cnt = dd(lambda: 0)
    for bc in all_bcs.values():
        if bc not in bc_to_clu:
            continue
        clu = bc_to_clu[bc]
        clu_to_cell_cnt[clu] += 1
    all_clus = list(set(clu_to_cell_cnt.keys()))
    for trans, bc_cnt in trans_to_bc_cnt.items():
        if trans not in trans_to_gene:
            continue
        trans_to_clu_n_cells = dd(lambda: 0)
        trans_to_clu_cnts = dd(lambda: 0.0)
        for bc, cnt in bc_cnt.items():
            if bc not in bc_to_clu:
                continue
            clu = bc_to_clu[bc]
            trans_to_clu_cnts[clu] += cnt
            if cnt >= min_trans_cnt_per_cell:
                trans_to_clu_n_cells[clu] += 1
        single_cell_keep, clu_keep = False, False
        for clu, clu_cell_cnt in clu_to_cell_cnt.items():
            if trans_to_clu_n_cells[clu] >= int(clu_cell_cnt * min_trans_cell_frac_per_clu):
                single_cell_keep = True
                break
        for clu, clu_total_cnt in trans_to_clu_cnts.items():
            if clu_total_cnt >= min_trans_total_cnt_per_clu:
                clu_keep = True
                break
        if single_cell_keep and clu_keep:
            gene = trans_to_gene[trans]
            for clu, trans_cnt in trans_to_clu_cnts.items():
                gene_trans_clu_cnt[gene][trans][clu] += trans_cnt
    kept_trans_n = 0
    with open('cluster_filtered_trans.tsv', 'w') as fp:
        for gene in list(gene_trans_clu_cnt.keys()):
            if len(gene_trans_clu_cnt[gene]) < min_trans_cnt:
                del gene_trans_clu_cnt[gene]
            else:
                kept_trans_n += len(gene_trans_clu_cnt[gene])
                for trans, cnt in gene_trans_clu_cnt[gene].items():
                    fp.write('{}\t{}\t{}\t{}\n'.format(gene, trans, cnt[all_clus[0]], cnt[all_clus[1]]))
    print('after filtering: {} genes {} transcripts'.format(len(gene_trans_clu_cnt), kept_trans_n))
    return gene_trans_clu_cnt

def get_gene_trans_clu_cnt(bc_to_clu, in_trans_mtx_dir):
    trans_to_bc_cnt, all_bcs = get_feature_to_bc_cnt_from_sc_matrix(in_trans_mtx_dir)
    trans_to_gene = get_trans_to_gene_from_sc_matrix(in_trans_mtx_dir)
    gene_trans_clu_cnt = dd(lambda: dd(lambda: dd(lambda: 0.0)))
    for trans, bc_cnt in trans_to_bc_cnt.items():
        if trans not in trans_to_gene:
            continue
        gene = trans_to_gene[trans]
        for bc, cnt in bc_cnt.items():
            if bc not in bc_to_clu:
                continue
            clu = bc_to_clu[bc]
            gene_trans_clu_cnt[gene][trans][clu] += cnt
    return gene_trans_clu_cnt

def get_gene_clu_cnt(bc_to_clu, in_gene_mtx_dir):
    gene_to_bc_cnt, all_bcs = get_feature_to_bc_cnt_from_sc_matrix(in_gene_mtx_dir)
    gene_clu_cnt = dd(lambda: dd(lambda: 0.0))
    for gene in gene_to_bc_cnt:
        for bc in gene_to_bc_cnt[gene]:
            if bc in bc_to_clu:
                clu = bc_to_clu[bc]
                gene_clu_cnt[gene][clu] += gene_to_bc_cnt[gene][bc]
    return gene_clu_cnt

def get_gene_list(gene_list_fn):
    gene_list = dict()
    with open(gene_list_fn, 'r') as fp:
        for line in fp:
            gene_list[line.strip()] = 1
    return gene_list

def cluster_wise_differential_splice1(gene_trans_clu_cnt, ct1, ct2, gene_fp, trans_fp,
                                      gene_list, trans_min_cnt=2, min_trans_per_gene=2):
    used_genes = []
    used_transs = []
    used_transs_cnt1s = []
    used_transs_cnt2s = [] 
    used_gene_chi_square_ps = []
    used_gene_to_max_trans_ratios = []
    used_gene_fdrs = []

    for gene, trans_clu_cnts in gene_trans_clu_cnt.items():
        used_trans = []
        used_trans_cnt1s = []
        used_trans_cnt2s = []
        for trans, clu_cnts in trans_clu_cnts.items():
            cnt1, cnt2 = clu_cnts[ct1], clu_cnts[ct2]
            if cnt1 >= trans_min_cnt or cnt2 >= trans_min_cnt:
                used_trans.append(trans)
                used_trans_cnt1s.append(cnt1)
                used_trans_cnt2s.append(cnt2)
        if len(used_trans) < min_trans_per_gene:
            continue
        
        # chi-square test for gene
        used_trans_ratio1s = [x / sum(used_trans_cnt1s) if sum(used_trans_cnt1s) > 0 else 0 for x in used_trans_cnt1s]
        used_trans_ratio2s = [x / sum(used_trans_cnt2s) if sum(used_trans_cnt2s) > 0 else 0 for x in used_trans_cnt2s]
        
        used_genes.append(gene)
        used_transs.append(used_trans)
        used_transs_cnt1s.append(used_trans_cnt1s)
        used_transs_cnt2s.append(used_trans_cnt2s)
        used_gene_chi_square_ps.append(get_chisq_test_p_value([used_trans_cnt1s, used_trans_cnt2s]))
        used_gene_to_max_trans_ratios.append(max([abs(x-y) for x,y in zip(used_trans_ratio1s, used_trans_ratio2s)]))

    used_gene_fdrs = get_fdr(used_gene_chi_square_ps)
    for gene, transs, trans_cnt1s, trans_cnt2s, gene_fdr, max_ratio in zip(used_genes, used_transs, used_transs_cnt1s, used_transs_cnt2s, used_gene_fdrs, used_gene_to_max_trans_ratios):
        gene_id, gene_name = gene.rsplit(',')
        if gene_list and gene_id not in gene_list:
            continue
        if math.isnan(gene_fdr) or gene_fdr > fdr_thres or max_ratio < delta_ratio:
            continue
        # gene_fp.write('gene_id\tgene_name\tgene_fdr\tmax_delta_ratio\ttranscripts\tcell1_counts\tcell2_counts\tcell1_ratios\tcell2_ratios\t\tcell_type1\tcell_type2\n')
        ratios1, ratios2 = [cnt / sum(trans_cnt1s) for cnt in trans_cnt1s], [cnt / sum(trans_cnt2s) for cnt in trans_cnt2s]
        gene_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(gene.rsplit(',')), gene_fdr, max_ratio, 
                                                            ','.join(transs), 
                                                            ','.join(list(map(str, trans_cnt1s))), ','.join(list(map(str, trans_cnt2s))),
                                                            ','.join(list(map(str, ratios1))), ','.join(list(map(str, ratios2))),
                                                            ct1, ct2))
        for trans, cnt1, cnt2 in zip(transs, trans_cnt1s, trans_cnt2s):
            # fisher exact test for trans
            fisher_trans_cnt1s = [cnt1, sum(trans_cnt1s) - cnt1]
            fisher_trans_cnt2s = [cnt2, sum(trans_cnt2s) - cnt2]
            ratio1 = cnt1 / sum(trans_cnt1s) if sum(trans_cnt1s) > 0 else 0
            ratio2 = cnt2 / sum(trans_cnt2s) if sum(trans_cnt2s) > 0 else 0
            trans_p = get_fisher_test_p_value([fisher_trans_cnt1s, fisher_trans_cnt2s])
            if math.isnan(trans_p) or trans_p > p_thres or abs(ratio1 - ratio2) < delta_ratio:
                continue
            trans_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(gene.rsplit(',')), trans, gene_fdr, trans_p, abs(ratio1-ratio2), cnt1, cnt2, ratio1, ratio2, ct1, ct2))

def cluster_wise_differential_splice(gene_trans_clu_cnt, all_cell_types, 
                                     out_prefix, gene_list, trans_min_cnt, min_trans_per_gene):
    pairwise_cell_types = list(itertools.combinations(all_cell_types, 2))
    with open(f'{out_prefix}_diff_splice_genes.tsv', 'w') as gene_fp, \
         open(f'{out_prefix}_diff_splice_transcripts.tsv', 'w') as trans_fp:
        # gene_trans_clu_cnt[gene][trans][clu] += cnt
        for gene, trans_clu_cnt in gene_trans_clu_cnt.items():
            clu_to_cnts = dd(lambda: 0.0) 
            for trans, clu_cnt in trans_clu_cnt.items():
                for clu, cnt in clu_cnt.items():
                    clu_to_cnts[clu] += cnt
        gene_fp.write('gene_id\tgene_name\tgene_fdr\tmax_delta_ratio\ttranscript_ids\tcell1_counts\tcell2_counts\tcell1_ratios\tcell2_ratios\tcell_type1\tcell_type2\n')
        trans_fp.write('gene_id\tgene_name\ttranscript_id\tgene_fdr\ttranscript_p\tdelta_ratio\tcount1\tcount2\tratio1\tratio2\tcell_type1\tcell_type2\n')
        for (ct1, ct2) in pairwise_cell_types:
            if ct2 < ct1:
                ct1, ct2 = ct2, ct1
            cluster_wise_differential_splice1(gene_trans_clu_cnt, ct1, ct2, gene_fp, trans_fp, gene_list, trans_min_cnt, min_trans_per_gene)

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="{}: identify cell type-wise differentially spliced genes and transcripts from single-cell long-read data".format(os.path.basename(__file__)))
    parser.add_argument('trans_mtx_dir', metavar='trans_mtx_dir', type=str, help='Transcript count matrix directory, generated by NanoHunter')
    parser.add_argument('bc_to_cluster', metavar='bc_to_cluster.tsv', type=str, help='List of cell barcode and corresponding cluster, generated by Seurat/Azimuth')
    parser.add_argument('out_prefix', metavar='out_prefix', type=str, help='Prefix of output files (differentially spliced gene/transcript)')
    parser.add_argument('-c', '--trans-min-count', type=int, default=10, help='Minimum number of read count of a kept transcript in at least one cell cluster')
    parser.add_argument('-t', '--min-trans-per-gene', type=int, default=2, help='Minimum number of transcripts of a kept gene')
    parser.add_argument('-l', '--gene-list', type=str, default=None, help='List of genes of interest. Only DS genes from the list will be outputted. If not specified, all genes are considered')
    return parser.parse_args()

def main():
    # input: 
    #   nanohunter transcript matrix directory, barcode to cluster tsv, output prefix
    # output:
    #   cell type specific spliced gene tsv
    #   gene_details.out:  gene/gene-FDR/transcript isoforms/cell type1/cell type2/trans count1/trans count2/ratio1/ratio2
    #   trans_details.out: trans/trans-pvalue/gene/gene-pvalue/gene-FDR/cell type1/cell type2/count1/count2/ratio1/ratio2
    args = parser_argv()
   
    in_trans_mtx_dir, bc_to_clu_tsv = args.trans_mtx_dir, args.bc_to_cluster, 
    out_prefix, trans_min_cnt, min_trans_per_gene = args.out_prefix, float(args.trans_min_count), float(args.min_trans_per_gene)
    gene_list_fn = args.gene_list
    bc_to_clu = get_bc_to_clu(bc_to_clu_tsv)
    gene_list = {} if gene_list_fn is None else get_gene_list(gene_list_fn)
    gene_trans_clu_cnt = get_gene_trans_clu_cnt(bc_to_clu, in_trans_mtx_dir)
    cluster_wise_differential_splice(gene_trans_clu_cnt, list(set(bc_to_clu.values())), \
                                     out_prefix, gene_list, trans_min_cnt, min_trans_per_gene)
    
if __name__ == '__main__':
    main()