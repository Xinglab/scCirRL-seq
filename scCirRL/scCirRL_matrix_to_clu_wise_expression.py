import sys
import os
import gzip
from collections import defaultdict as dd
from .scCirRL_cell_type_specific_splicing import *

def write_gene_to_clu_cnt(gene_clu_cnt, gene_clu_pct, all_clus, out_prefix):
    gene_cnt_out = out_prefix + '.clu_wise_gene_cnt.tsv'
    gene_cpm_out = out_prefix + '.clu_wise_gene_cpm.tsv'
    gene_pct_out = out_prefix + '.clu_wise_gene_pct.tsv'

    clu_to_total_cnt = dd(lambda: 0.0)
    with open(gene_cnt_out, 'w') as gene_cnt_fp, open(gene_pct_out, 'w') as gene_pct_fp:
        gene_cnt_fp.write('gene_name\t{}\n'.format('\t'.join(all_clus)))
        gene_pct_fp.write('gene_name\t{}\n'.format('\t'.join(all_clus)))
        for gene, clu_cnt in gene_clu_cnt.items():
            gene_cnt_fp.write(f'{gene}')
            gene_pct_fp.write(f'{gene}')
            for clu in all_clus:
                gene_cnt_fp.write('\t{}'.format(clu_cnt[clu]))
                clu_to_total_cnt[clu] += clu_cnt[clu]
                gene_pct_fp.write('\t{}'.format(gene_clu_pct[gene][clu]))
            gene_cnt_fp.write('\n')
            gene_pct_fp.write('\n')

    # calculate cluster-wise CPM
    with open(gene_cpm_out, 'w') as gene_cpm_fp:
        gene_cpm_fp.write('gene_name\t{}\n'.format('\t'.join(all_clus)))
        for gene, clu_cnt in gene_clu_cnt.items():
            gene_cpm_fp.write(f'{gene}')
            for clu in all_clus:
                cnt = clu_cnt[clu]
                gene_cpm_fp.write('\t{}'.format(0.0 if clu_to_total_cnt[clu] == 0.0 else 1000000 * cnt / clu_to_total_cnt[clu]))
            gene_cpm_fp.write('\n')

def write_gene_trans_to_clu_cnt(gene_trans_clu_cnt, gene_clu_pct, all_clus, out_prefix):
    gene_cnt_out, trans_cnt_out = out_prefix + '.clu_wise_gene_cnt.tsv', out_prefix + '.clu_wise_trans_cnt.tsv'
    gene_cpm_out, trans_cpm_out = out_prefix + '.clu_wise_gene_cpm.tsv', out_prefix + '.clu_wise_trans_cpm.tsv'
    gene_pct_out = out_prefix + '.clu_wise_gene_pct.tsv'
    
    gene_to_clu_cnt = dd(lambda: dd(lambda: 0.0))
    clu_to_total_cnt = dd(lambda: 0.0)
    with open(gene_cnt_out, 'w') as gene_cnt_fp, open(gene_pct_out, 'w') as gene_pct_fp, open(trans_cnt_out, 'w') as trans_cnt_fp:
        gene_cnt_fp.write('gene_id\tgene_name\t{}\n'.format('\t'.join(all_clus)))
        gene_pct_fp.write('gene_id\tgene_name\t{}\n'.format('\t'.join(all_clus)))
        trans_cnt_fp.write('gene_id\tgene_name\ttranscript_id\t{}\n'.format('\t'.join(all_clus)))
        for gene, trans_to_clu_cnt in gene_trans_clu_cnt.items():
            gene_id, gene_name = gene.split(',')
            if gene_name == 'PTPRC':
                print('ok')
            gene_cnt_fp.write(f'{gene_id}\t{gene_name}')
            gene_pct_fp.write(f'{gene_id}\t{gene_name}')
            for trans, clu_cnt in trans_to_clu_cnt.items():
                trans_cnt_fp.write(f'{gene_id}\t{gene_name}\t{trans}')
                for clu in all_clus:
                    trans_cnt_fp.write('\t{}'.format(clu_cnt[clu]))
                    gene_to_clu_cnt[gene][clu] += clu_cnt[clu]
                    clu_to_total_cnt[clu] += clu_cnt[clu]
                trans_cnt_fp.write('\n')
            for clu in all_clus:
                gene_cnt_fp.write('\t{}'.format(gene_to_clu_cnt[gene][clu]))
                gene_pct_fp.write('\t{}'.format(gene_clu_pct[gene][clu]))
            gene_cnt_fp.write('\n')
            gene_pct_fp.write('\n')
    with open(gene_cpm_out, 'w') as gene_cpm_fp, open(trans_cpm_out, 'w') as trans_cpm_fp:
        gene_cpm_fp.write('gene_id\tgene_name\t{}\n'.format('\t'.join(all_clus)))
        trans_cpm_fp.write('gene_id\tgene_name\ttranscript_id\t{}\n'.format('\t'.join(all_clus)))
        for gene, trans_to_clu_cnt in gene_trans_clu_cnt.items():
            gene_id, gene_name = gene.split(',')
            gene_cpm_fp.write(f'{gene_id}\t{gene_name}') 
            for clu in all_clus:
                cnt = gene_to_clu_cnt[gene][clu]
                gene_cpm_fp.write('\t{}'.format(0.0 if clu_to_total_cnt[clu] == 0.0 else 1000000 * cnt / clu_to_total_cnt[clu]))
            gene_cpm_fp.write('\n')
            for trans, clu_cnt in trans_to_clu_cnt.items():
                trans_cpm_fp.write(f'{gene_id}\t{gene_name}\t{trans}') 
                for clu in all_clus:
                    cnt = clu_cnt[clu]
                    trans_cpm_fp.write('\t{}'.format(0.0 if clu_to_total_cnt[clu] == 0.0 else 1000000 * cnt / clu_to_total_cnt[clu]))
                trans_cpm_fp.write('\n')

def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: generate cluster-wise gene/transcript count/percent.expressed/CPM matrix".format(os.path.basename(__file__)))
    parser.add_argument('mtx_dir', metavar='gene/trans_mtx_dir', type=str, help='Input gene/transcript count matrix directory')
    parser.add_argument('bc_to_clu', metavar='bc_to_clu.tsv', type=str, help='Input barcode to cluster mapping file')
    parser.add_argument('out_prefix', metavar='out_prefix', type=str, help='Output file prefix')
    return parser.parse_args()

def main():
    args = parser_argv()
    feature_mtx_dir = args.mtx_dir
    bc_to_clu_tsv = args.bc_to_clu
    out_prefix = args.out_prefix
    bc_to_clu = get_bc_to_clu(bc_to_clu_tsv)
    all_clus = set(list(bc_to_clu.values()))

    if os.path.exists(feature_mtx_dir + '/gene_to_trans.tsv.gz'):
        gene_trans_clu_cnt, gene_clu_pct = get_gene_trans_clu_cnt(bc_to_clu, feature_mtx_dir)
        write_gene_trans_to_clu_cnt(gene_trans_clu_cnt, gene_clu_pct, all_clus, out_prefix)
    else:
        gene_clu_cnt, gene_clu_pct = get_gene_clu_cnt(bc_to_clu, feature_mtx_dir)
        write_gene_to_clu_cnt(gene_clu_cnt, gene_clu_pct, all_clus, out_prefix)

if __name__ == '__main__':
    main()