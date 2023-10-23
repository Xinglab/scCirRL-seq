import sys
import os
import gzip
from collections import defaultdict as dd
from nh_cell_type_specific_splicing import *

def write_gene_trans_to_clu_cnt(gene_trans_clu_cnt, all_clus, out_prefix):
    gene_out, trans_out = out_prefix + '.clu_wise_gene_cnt.tsv', out_prefix + '.clu_wise_trans_cnt.tsv'
    with open(gene_out, 'w') as gene_fp, open(trans_out, 'w') as trans_fp:
        gene_fp.write('gene_id\tgene_name\t{}\n'.format('\t'.join(all_clus)))
        trans_fp.write('gene_id\tgene_name\ttranscript_id\t{}\n'.format('\t'.join(all_clus)))
        for gene, trans_to_clu_cnt in gene_trans_clu_cnt.items():
            gene_id, gene_name = gene.split(',')
            gene_to_clu_cnt = dd(lambda: 0)
            gene_fp.write(f'{gene_id}\t{gene_name}')
            for trans, clu_cnt in trans_to_clu_cnt.items():
                trans_fp.write(f'{gene_id}\t{gene_name}\t{trans}')
                for clu in all_clus:
                    trans_fp.write('\t{}'.format(clu_cnt[clu]))
                    gene_to_clu_cnt[clu] += clu_cnt[clu]
                trans_fp.write('\n')
            for clu in all_clus:
                gene_fp.write('\t{}'.format(gene_to_clu_cnt[clu]))
            gene_fp.write('\n')

def write_gene_to_clu_cnt(gene_clu_cnt, all_clus, out_prefix):
    gene_out = out_prefix + '.clu_wise_gene_cnt.tsv'
    with open(gene_out, 'w') as gene_fp:
        gene_fp.write('gene_id\t{}\n'.format('\t'.join(all_clus)))
        for gene, clu_cnt in gene_clu_cnt.items():
            gene_fp.write(f'{gene}')
            for clu in all_clus:
                gene_fp.write('\t{}'.format(clu_cnt[clu]))
            gene_fp.write('\n')

def main():
    if len(sys.argv) != 4:
        print('{} gene/trans_mtx_dir bc_to_clu.tsv out_prefix'.format(sys.argv[0]))
        print('\tgenerate cluster-wise gene/transcript count/CPM matrix')
        sys.exit(1)
    feature_mtx_dir = sys.argv[1]
    bc_to_clu_tsv = sys.argv[2]
    out_prefix = sys.argv[3]
    bc_to_clu = get_bc_to_clu(bc_to_clu_tsv)
    all_clus = set(list(bc_to_clu.values()))

    if os.path.exists(feature_mtx_dir + '/gene_to_trans.tsv.gz'):
        gene_trans_clu_cnt = get_gene_trans_clu_cnt(bc_to_clu, feature_mtx_dir)
        write_gene_trans_to_clu_cnt(gene_trans_clu_cnt, all_clus, out_prefix)
    else:
        gene_clu_cnt = get_gene_clu_cnt(bc_to_clu, feature_mtx_dir)
        write_gene_to_clu_cnt(gene_clu_cnt, all_clus, out_prefix)

if __name__ == '__main__':
    main()