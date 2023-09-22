import argparse
import sys
import os
import gzip
from scipy.sparse import csr_matrix  # for sparse matrix
from collections import defaultdict as dd

from .__init__ import __version__
from .__init__ import __program__

skipped_cates = []
delta_thres = 1e-6
theta_base = 1e6
max_iter = 10000

default_cluster = 'cluster_0'

# return: all_bc, bu_to_trans
def get_bu_cmpt(in_bu_cmpt, bc_to_cluster):
    cand_bc = list(bc_to_cluster.keys())
    bu_to_trans = dd(lambda: dd(lambda: ''))
    gene_to_trans = dd(lambda: set())
    all_bc = dict()
    with open(in_bu_cmpt) as fp:
        i = 0
        for line in fp:
            ele = line.rsplit()
            bc, umi, read_cnt, trans_str, gene_id_str, gene_name_str, reads = ele
            if (cand_bc != [] and bc not in cand_bc) or trans_str == 'NA' or gene_id_str == 'NA' or ',' in gene_id_str:
                continue
            umi = umi + str(i)
            trans = trans_str.rsplit(',')
            gene = gene_id_str + ',' + gene_name_str
            bu_to_trans[gene][(bc, umi)] = trans
            gene_to_trans[gene] = gene_to_trans[gene].union(set(trans))
            all_bc[bc] = 1
            i += 1
    return list(all_bc.keys()), bu_to_trans, gene_to_trans


# read:bc:cellType
def parse_cell_cluster(in_fn):
    if not in_fn:
        return dd(lambda: default_cluster)

    bc_to_cluter = dd(lambda: None)  # barcodes may not be in any clusters
    with open(in_fn) as fp:
        for line in fp:
            ele = line.rstrip().rsplit('\t')
            if len(ele) == 1:
                bc = ele[0].rsplit('-')[0]
                cluster = default_cluster
            elif len(ele) >= 2:
                bc, cluster = ele[0].rsplit('-')[0], ele[1]
            else:
                continue
            bc_to_cluter[bc] = cluster
    return bc_to_cluter


def get_gene_name(in_gtf):
    anno_gene_id_to_name = dict()
    with open(in_gtf) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            if ele[2] != 'gene':
                continue
            gene_id, gene_name = '', ''
            if 'gene_id' in line:
                gene_ids = line[9+line.index('gene_id'):]
                gene_id = gene_ids[:gene_ids.index('"')]

            if 'gene_name' in line:
                gene_names = line[11+line.index('gene_name'):]
                gene_name = gene_names[:gene_names.index('"')]

            if gene_id and gene_name:
                anno_gene_id_to_name[gene_id] = gene_name
    return anno_gene_id_to_name


def get_trans_to_gene(updated_gtf, anno_gtf):
    gene_id_to_name = get_gene_name(anno_gtf)
    trans_to_gene = dict()
    gene_to_trans = dd(lambda: [])
    with open(updated_gtf) as fp:
        for line in fp:
            ele = line.rsplit()
            if ele[2] != 'transcript':
                continue
            trans_id, gene_id, gene_name = '', '', ''
            if 'transcript_id' in line:
                trans_ids = line[15+line.index('transcript_id'):]
                trans_id = trans_ids[:trans_ids.index('"')]
            if 'gene_id' in line:
                gene_ids = line[9+line.index('gene_id'):]
                gene_id = gene_ids[:gene_ids.index('"')]
            if gene_id == 'NA' or ',E' in gene_id:
                continue
            if gene_id in gene_id_to_name:
                gene_name = gene_id_to_name[gene_id]
            else:
                gene_name = gene_id
            if trans_id and gene_name and gene_id:
                trans_to_gene[trans_id] = gene_name
                gene_to_trans[gene_id + ',' + gene_name].append(trans_id)
    return trans_to_gene, gene_to_trans


def parse_cmpt(in_cmpt, read_to_bu, bc_to_cluster, trans_to_gene):
    trans_to_reads = dd(lambda: dd(lambda: []))
    read_to_cluster = dict()
    with open(in_cmpt) as fp:
        for line in fp:
            ele = line.rsplit()
            qname, sample, cate, transs = ele
            if qname not in read_to_bu:
                continue
                # sys.stderr.write('Unexpectecd query name: {}\n'.format(qname))
                # sys.exit(1)
            bc = read_to_bu[qname]['bc']
            cluster = bc_to_cluster[bc]
            if transs == 'NA' or cate in skipped_cates:
                continue
            transs = transs.rsplit(',')
            if '' in transs:
                transs.remove('')
            genes = []
            for trans in transs:
                if trans not in trans_to_gene:
                    continue
                gene = trans_to_gene[trans]
                if gene not in genes:
                    genes.append(gene)
            if len(genes) > 1:
                # sys.stderr.write('{}\t{}\t{}\n'.format(qname, ','.join(transs), ','.join(genes)))
                continue
            for trans in transs:
                trans_to_reads[cluster][trans].append(qname)
            read_to_cluster[qname] = cluster
    return trans_to_reads, read_to_cluster


def em_abun_estimation1(gene, all_reads, all_trans, trans_to_reads, read_to_group):
    n_trans = len(all_trans)
    read_to_trans = dd(lambda: [])  # XXX record read_to_trans each-read-wise,
                                    # instead of simply using the count sum
                                    # this is for furture use of P(read|iso)
    # initial theta
    cur_thetas = dict()
    for trans in all_trans:
        cur_thetas[trans] = theta_base / n_trans
        for read in trans_to_reads[trans]:
            read_to_trans[read].append(trans)

    # EM loop
    n_iter = 0
    cur_trans_abun = dd(lambda: 0.0)
    group_trans_abun = dd(lambda: dd(lambda: 0.0))
    if len(read_to_trans) == 0:
        return cur_trans_abun, group_trans_abun
    while True:
        n_iter += 1
        # E-step:
        # calculate fractional counts based on current thetas
        cur_trans_abun = dd(lambda: 0.0)
        group_trans_abun = dd(lambda: dd(lambda: 0.0))
        for read, cmpt_trans in read_to_trans.items():
            # if 'ENST00000619225.1' in cmpt_trans:
                # print('ok')
            total = sum([cur_thetas[trans] for trans in cmpt_trans])
            for trans in cmpt_trans:
                cur_trans_abun[trans] += theta_base * (cur_thetas[trans] / total)
                group_trans_abun[read_to_group[read]][trans] += theta_base * (cur_thetas[trans] / total)

        # M-step:
        # update the probability of each isoform
        total = sum(cur_trans_abun.values())
        new_thetas = dict()
        for trans in all_trans:
            new_thetas[trans] = theta_base * cur_trans_abun[trans] / total

        # iteration stop check
        delta = 0
        for trans in all_trans:
            delta += abs(cur_thetas[trans] - new_thetas[trans])
            cur_thetas[trans] = new_thetas[trans]
        if delta < delta_thres * theta_base:
            return cur_trans_abun, group_trans_abun
            break

        if n_iter > max_iter:
            # sys.stderr.write('{}:\t{}\n'.format(gene, '\t'.join(all_trans)))
            break
    return cur_trans_abun, group_trans_abun


def em_abun_estimation(read_to_bu, group_dict, all_groups, trans_to_reads, gene_to_trans, read_to_group):
    for gene, all_trans in gene_to_trans.items():
        group_em_res = dd(lambda: ())
        for group in all_groups:  # EM within each group and each gene
            all_reads = dict()
            for trans in all_trans:
                for read in trans_to_reads[group][trans]:
                    all_reads[read] = 1
            abun, sample_abun = em_abun_estimation1(gene, all_reads, all_trans, trans_to_reads[group], read_to_group)
            group_em_res[group] = (abun, sample_abun)
        for trans in all_trans:
            sys.stdout.write('{}\t{}'.format(trans, gene))
            for group in all_groups:
                sys.stdout.write('\t{}'.format(group_em_res[group][0][trans] / theta_base))
                for sample in group_dict:
                    if group_dict[sample] == group:
                        sys.stdout.write('\t{}'.format(group_em_res[group][1][sample][trans] / theta_base))
            sys.stdout.write('\n')


# umi-based count
def em_frac_abun1(gene, all_trans, bu_to_trans, cluster, bc_to_cluster):  # trans_to_reads, read_to_bu):
    n_trans = len(all_trans)
    # XXX record read_to_trans each-read-wise,
    # instead of simply using the count sum
    # this is for furture use of P(read|iso)

    # initial theta
    cur_thetas = dict()
    for trans in all_trans:
        cur_thetas[trans] = theta_base / n_trans
        # for read in trans_to_reads[trans]:
        #     if trans not in read_to_trans:
        #         read_to_trans[read].append(trans)

    # EM loop
    n_iter = 0
    cur_trans_abun = dd(lambda: 0.0)
    bc_trans_abun = dd(lambda: dd(lambda: 0.0))
    if len(bu_to_trans) == 0:
        return bc_trans_abun
    while True:
        n_iter += 1
        # E-step:
        # calculate fractional counts based on current thetas
        cur_trans_abun = dd(lambda: 0.0)
        bc_trans_abun = dd(lambda: dd(lambda: 0.0))
        for (bc, umi), cmpt_trans in bu_to_trans.items():
            if bc_to_cluster[bc] != cluster:
                continue
            total = sum([cur_thetas[trans] for trans in cmpt_trans])
            for trans in cmpt_trans:
                cur_trans_abun[trans] += theta_base * (cur_thetas[trans] / total)
                bc_trans_abun[bc][trans] += theta_base * (cur_thetas[trans] / total)

        # M-step:
        # update the probability of each isoform
        total = sum(cur_trans_abun.values())
        if total == 0.0:
            return bc_trans_abun
        new_thetas = dict()
        for trans in all_trans:
            new_thetas[trans] = theta_base * cur_trans_abun[trans] / total

        # iteration stop check
        delta = 0
        for trans in all_trans:
            delta += abs(cur_thetas[trans] - new_thetas[trans])
            cur_thetas[trans] = new_thetas[trans]
        if delta < delta_thres * theta_base:
            return bc_trans_abun
            break

        if n_iter > max_iter:
            # sys.stderr.write('{}:\t{}\n'.format(gene, '\t'.join(all_trans)))
            break
    return bc_trans_abun


def get_full_length_cnt(all_trans, trans_to_reads):
    read_to_trans = dd(lambda: [])
    trans_to_full_cnt = dd(lambda: 0.0)
    for trans in all_trans:
        for group in trans_to_reads:
            for read in trans_to_reads[group][trans]:
                if trans not in read_to_trans:
                    read_to_trans[read].append(trans)
    n_full = 0
    for read in read_to_trans:
        if len(read_to_trans[read]) == 1:
            n_full += 1
            trans = read_to_trans[read][0]
            trans_to_full_cnt[trans] += 1
    return n_full, len(read_to_trans), trans_to_full_cnt

def write_gene_to_trans_tsv(gene_to_trans, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_fn = out_dir + '/gene_to_trans.tsv.gz'
    with gzip.open(out_fn, 'wt') as fp:
        fp.write('gene_id\tgene_name\ttranscript\n')
        for gene, all_trans in gene_to_trans.items():
            gene_id, gene_name = gene.rstrip().split(',')
            for trans in all_trans:
                fp.write('{}\t{}\t{}\n'.format(gene_id, gene_name, trans))

def write_10X_sparse_matrix(matrix_row, matrix_col, matrix_data, names, all_bc, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_mtx_fn, out_feature_fn, out_bc_fn = out_dir + '/matrix.mtx.gz', out_dir + '/features.tsv.gz', out_dir + '/barcodes.tsv.gz'

    # matrix
    with gzip.open(out_mtx_fn, 'wt') as fp:
        fp.write('%%MatrixMarket matrix coordinate integer general\n')
        fp.write('%metadata_json: {{"software_version": "{}-{}", "format_version": 2}}\n'.format(__program__, __version__))
        fp.write('{} {} {}\n'.format(len(names), len(all_bc), len(matrix_row)))
        for r, c, d in zip(matrix_row, matrix_col, matrix_data):
            fp.write('{} {} {}\n'.format(r, c, d))

    # feature
    with gzip.open(out_feature_fn, 'wt') as fp:
        for name in names:
            fp.write('{}\n'.format(name))

    # barcode
    with gzip.open(out_bc_fn, 'wt') as fp:
        for bc in all_bc:
            fp.write('{}\n'.format(bc))


def em_frac_abun(all_bc, bc_to_cluster, bu_to_trans, gene_to_trans, out_dir):
    # with open(gene_out_fn, 'w') as gene_out_fp, open(iso_out_fn, 'w') as iso_out_fp:
    # gene_out_fp.write('Gene\t{}\n'.format('\t'.join(all_bc)))
    # iso_out_fp.write('Isoform\t{}\n'.format('\t'.join(all_bc)))

    gene_matrix_row, gene_matrix_col, gene_matrix_data, gene_names = [], [], [], []
    iso_matrix_row, iso_matrix_col, iso_matrix_data, iso_ratio_matrix_data, iso_names = [], [], [], [], []
    all_clusters = list(set(bc_to_cluster.values()))
    if all_clusters == []:
        all_clusters = [default_cluster]
    
    gene_i, trans_i = 1, 1  # matrix_row
    for gene, all_trans in gene_to_trans.items():
        all_trans = list(all_trans)
        gene_id, gene_name = gene.rsplit(',')
        bc_to_gene_abun = dd(lambda: 0)
        bc_to_trans_abun = dd(lambda: dd(lambda: 0))
        # gene_full_cnt, gene_total_cnt, trans_full_cnt = get_full_length_cnt(all_trans, trans_to_reads)
        trans_total_cnt = dd(lambda: 0.0)
        for cluster in all_clusters:  # EM within each group and each gene
            bc_trans_abun = em_frac_abun1(gene_name, all_trans, bu_to_trans[gene], cluster, bc_to_cluster)
            for bc in bc_trans_abun:
                for trans in bc_trans_abun[bc]:
                    bc_to_gene_abun[bc] += bc_trans_abun[bc][trans]
                    bc_to_trans_abun[bc][trans] += bc_trans_abun[bc][trans]
                    trans_total_cnt[trans] += bc_trans_abun[bc][trans]/theta_base
        # for bc_i, bc in enumerate(all_bc):
        for bc, gene_abun in bc_to_gene_abun.items():
            if gene_abun == 0.0:
                continue
            gene_matrix_row.append(gene_i)
            gene_matrix_col.append(all_bc.index(bc)+1)
            gene_matrix_data.append(gene_abun/theta_base)
        gene_names.append(gene_id + '\t' + gene_name)
        gene_i += 1

        for bc in bc_to_trans_abun:
            gene_abun = bc_to_gene_abun[bc]
            if gene_abun == 0.0:
                continue
            for trans in bc_to_trans_abun[bc]:
                trans_abun = bc_to_trans_abun[bc][trans]
                if trans_abun == 0.0:
                    continue
                iso_matrix_row.append(trans_i + all_trans.index(trans))
                iso_matrix_col.append(all_bc.index(bc)+1)
                iso_matrix_data.append(trans_abun/theta_base)
                iso_ratio_matrix_data.append(trans_abun/gene_abun)
        for trans in all_trans:
            # iso_names.append(gene_id + ',' + trans + '\t' + gene_name + ',' + trans)
            iso_names.append(trans + '\t' + trans)
        trans_i += len(all_trans)
    write_10X_sparse_matrix(gene_matrix_row, gene_matrix_col, gene_matrix_data, gene_names, all_bc, out_dir+'/gene')
    write_10X_sparse_matrix(iso_matrix_row, iso_matrix_col, iso_matrix_data, iso_names, all_bc, out_dir+'/transcript')
    write_10X_sparse_matrix(iso_matrix_row, iso_matrix_col, iso_ratio_matrix_data, iso_names, all_bc, out_dir+'/transcript_ratio')
        # gene_out_fp.write('{}'.format(gene))
        # for bc in all_bc:
        #     gene_out_fp.write('\t{}'.format(bc_to_gene_abun[bc]/theta_base))
        # gene_out_fp.write('\n')
        # for trans in all_trans:
        #     iso_out_fp.write('{},{}'.format(gene, trans))
        #     for bc in all_bc:
        #         iso_out_fp.write('\t{}'.format(bc_to_trans_abun[bc][trans]/theta_base))
        #     iso_out_fp.write('\n')


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="{}: make sparsed count matrix at both gene and transcript levels".format(os.path.basename(__file__)))
    parser.add_argument('bu_out', metavar='barcode_umi.out', type=str, help='Output file of NanoHunter')
    parser.add_argument('out_dir', metavar='matrix_output_dir', type=str, help='Folder to output single-cell gene/transcript matrix')
    parser.add_argument('-c', '--bc-cluster', type=str, default='', help='List of cell barcode and corresponding cluster (can be used for cell cluster-wise transcript EM abundance estimation)')

    return parser.parse_args()


def make_10X_matrix(bu_fn, cell_to_cluster_tsv, out_dir):
    bc_to_cluster = parse_cell_cluster(cell_to_cluster_tsv)
    all_bc, bu_to_trans, gene_to_trans = get_bu_cmpt(bu_fn, bc_to_cluster)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    write_gene_to_trans_tsv(gene_to_trans, out_dir+'/transcript')
    write_gene_to_trans_tsv(gene_to_trans, out_dir+'/transcript_ratio')
    em_frac_abun(all_bc, bc_to_cluster, bu_to_trans, gene_to_trans, out_dir)



def main():
    args = parser_argv()
    in_bu_fn, out_dir = args.bu_out, args.out_dir
    bc_to_cluster_tsv = args.bc_cluster
    make_10X_matrix(in_bu_fn, bc_to_cluster_tsv, out_dir)

if __name__ == '__main__':
    main()
