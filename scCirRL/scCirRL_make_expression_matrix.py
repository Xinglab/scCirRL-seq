import argparse
import sys
import os
import gzip
from scipy.sparse import csr_matrix  # for sparse matrix
from collections import defaultdict as dd
from .scCirRL_allele_specific_splicing import get_reads_to_hap
from .parameter import scrl_para_class
from .utils import err_log_format_time

from .__init__ import __version__
from .__init__ import __program__

skipped_cates = []
delta_thres = 1e-6
theta_base = 1e6
max_iter = 10000

default_cluster = 'cluster_0'

class scrl_matrix_data:
    def __init__(self, row, col, data, bc, names):
        self.row = row
        self.col = col
        self.data = data
        self.bc = bc
        self.names = names
    def write_matrix_data(self, out_dir):
         write_10X_sparse_matrix(self.row, self.col, self.data, self.names, self.bc, out_dir)

# return: all_bc, bu_to_trans
def get_bu_cmpt(in_bu_cmpt, trans_cate, only_splice, bc_to_cluster):
    cand_bc = list(bc_to_cluster.keys())
    bu_to_trans = dd(lambda: dd(lambda: ''))
    gene_to_trans = dd(lambda: set())
    all_bc = dict()
    with open(in_bu_cmpt) as fp:
        i = 0
        header = True
        for line in fp:
            if header:
                header = False
                continue
            ele = line.rsplit()
            bc, umi, read_cnt, trans_str, gene_id_str, gene_name_str, _reads, _read_cates, _read_splice_tags = ele
            reads, read_cates, read_splice_tags, = filter_by_read_cate(_reads, _read_cates, _read_splice_tags, trans_cate, only_splice)
            if len(reads) == 0:
                continue
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

def get_hap(qnames, read_to_hap):
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

def filter_by_read_cate(_reads, _read_cates, _read_splice_tags, read_cate, only_splice):
    if read_cate == [] and only_splice == False:
        return _reads, _read_cates, _read_splice_tags
    reads, read_cates, tags = [], [], []
    for qname, cate, tag in zip(_reads.rsplit(','), _read_cates.rsplit(','), _read_splice_tags.rsplit(',')):
        # print(qname, cate, tag, read_cate, only_splice)
        if read_cate != [] and cate not in read_cate:
            continue
        if only_splice and tag == 'N':
            continue
        reads.append(qname)
        read_cates.append(cate)
        tags.append(tag)
    return ','.join(reads), ','.join(read_cates), ','.join(tags)
    
def get_bu_hap_cmpt(in_bu_cmpt, trans_cate, only_splice, bc_to_cluster, read_to_hap):
    cand_bc = list(bc_to_cluster.keys())
    bu_to_hap_trans = dd(lambda: dd(lambda: ''))
    gene_to_trans = dd(lambda: set())
    all_bc = dict()
    with open(in_bu_cmpt) as fp:
        i = 0
        header = True
        for line in fp:
            if header:
                header = False
                continue
            ele = line.rsplit()
            bc, umi, read_cnt, trans_str, gene_id_str, gene_name_str, _reads, _read_cates, _read_splice_tags = ele
            reads, read_cates, read_splice_tags, = filter_by_read_cate(_reads, _read_cates, _read_splice_tags, trans_cate, only_splice)
            if len(reads) == 0:
                continue
            if (cand_bc != [] and bc not in cand_bc) or trans_str == 'NA' or gene_id_str == 'NA' or ',' in gene_id_str:
                continue
            hap = get_hap(reads, read_to_hap)
            umi = umi + str(i)
            trans = trans_str.rsplit(',')
            gene = gene_id_str + ',' + gene_name_str
            bu_to_hap_trans[gene][(bc, umi, hap)] = trans
            gene_to_trans[gene] = gene_to_trans[gene].union(set(trans))
            all_bc[bc] = 1
            i += 1
    return list(all_bc.keys()), bu_to_hap_trans, gene_to_trans

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

def hap_em_frac_abun1(gene, all_trans, bu_to_hap_trans, cluster, bc_to_cluster):  # trans_to_reads, read_to_bu):
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
    bc_hap_trans_abun = dd(lambda: dd(lambda: dd(lambda: 0.0)))
    if len(bu_to_hap_trans) == 0:
        return bc_hap_trans_abun
    while True:
        n_iter += 1
        # E-step:
        # calculate fractional counts based on current thetas
        cur_trans_abun = dd(lambda: 0.0)
        bc_hap_trans_abun = dd(lambda: dd(lambda: dd(lambda: 0.0)))
        for (bc, umi, hap), cmpt_trans in bu_to_hap_trans.items():
            if bc_to_cluster[bc] != cluster:
                continue
            total = sum([cur_thetas[trans] for trans in cmpt_trans])
            for trans in cmpt_trans:
                cur_trans_abun[trans] += theta_base * (cur_thetas[trans] / total)
                bc_hap_trans_abun[bc][trans][hap] += theta_base * (cur_thetas[trans] / total)

        # M-step:
        # update the probability of each isoform
        total = sum(cur_trans_abun.values())
        if total == 0.0:
            return bc_hap_trans_abun
        new_thetas = dict()
        for trans in all_trans:
            new_thetas[trans] = theta_base * cur_trans_abun[trans] / total

        # iteration stop check
        delta = 0
        for trans in all_trans:
            delta += abs(cur_thetas[trans] - new_thetas[trans])
            cur_thetas[trans] = new_thetas[trans]
        if delta < delta_thres * theta_base:
            return bc_hap_trans_abun
            break

        if n_iter > max_iter:
            # sys.stderr.write('{}:\t{}\n'.format(gene, '\t'.join(all_trans)))
            break
    return bc_hap_trans_abun

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

def make_matrix_data_with_hap_em(all_bc, bc_to_cluster, bu_to_hap_trans, gene_to_trans):
    all_clusters = list(set(bc_to_cluster.values()))
    hap_to_gene_mtx_data = dd(lambda: dd(lambda: scrl_matrix_data)) # {hap: {gene: mxt_data}}
    hap_to_trans_mtx_data = dd(lambda: dd(lambda: scrl_matrix_data)) # {hap: {gene: mxt_data}}
    if all_clusters == []:
        all_clusters = [default_cluster]
    
    gene_i, trans_i = 1, 1  # matrix_row

    hap_gene_matrix_row, hap_gene_matrix_col, hap_gene_matrix_data, gene_names = dd(lambda: []), dd(lambda: []), dd(lambda: []), []
    hap_iso_matrix_row, hap_iso_matrix_col, hap_iso_matrix_data, iso_names = dd(lambda: []), dd(lambda: []), dd(lambda: []), []
    for gene, all_trans in gene_to_trans.items():
        all_trans = list(all_trans)
        gene_id, gene_name = gene.rsplit(',')
        bc_to_hap_gene_abun = dd(lambda: dd(lambda: 0))
        bc_to_hap_trans_abun = dd(lambda: dd(lambda: dd(lambda: 0)))
        hap_trans_total_cnt = dd(lambda: dd(lambda: 0.0))
        for cluster in all_clusters:  # EM within each group and each gene
            bc_hap_trans_abun = hap_em_frac_abun1(gene_name, all_trans, bu_to_hap_trans[gene], cluster, bc_to_cluster)
            for bc in bc_hap_trans_abun:
                for trans in bc_hap_trans_abun[bc]:
                    for hap in ['H1', 'H2', 'none']:
                        bc_to_hap_gene_abun[bc][hap] += bc_hap_trans_abun[bc][trans][hap]
                        bc_to_hap_trans_abun[bc][trans][hap] += bc_hap_trans_abun[bc][trans][hap]
                        hap_trans_total_cnt[trans][hap] += bc_hap_trans_abun[bc][trans][hap]/theta_base
        for bc in bc_to_hap_gene_abun:
            for hap in ['H1', 'H2', 'none']:
                gene_abun = bc_to_hap_gene_abun[bc][hap]
                if gene_abun == 0.0:
                    continue
                hap_gene_matrix_row[hap].append(gene_i)
                hap_gene_matrix_col[hap].append(all_bc.index(bc)+1)
                hap_gene_matrix_data[hap].append(gene_abun/theta_base)
        gene_names.append(gene_id + '\t' + gene_name)
        gene_i += 1

        for bc in bc_to_hap_trans_abun:
            for hap in ['H1', 'H2', 'none']:
                gene_abun = bc_to_hap_gene_abun[bc][hap]
                if gene_abun == 0.0:
                    continue
                for trans in bc_to_hap_trans_abun[bc]:
                    trans_abun = bc_to_hap_trans_abun[bc][trans][hap]
                    if trans_abun == 0.0:
                        continue
                    hap_iso_matrix_row[hap].append(trans_i + all_trans.index(trans))
                    hap_iso_matrix_col[hap].append(all_bc.index(bc)+1)
                    hap_iso_matrix_data[hap].append(trans_abun/theta_base)
        for trans in all_trans:
            iso_names.append(trans + '\t' + trans)
        trans_i += len(all_trans)
    for hap in ['H1', 'H2', 'none']:
        hap_to_gene_mtx_data[hap] = scrl_matrix_data(hap_gene_matrix_row[hap], hap_gene_matrix_col[hap], hap_gene_matrix_data[hap], all_bc, gene_names)
        hap_to_trans_mtx_data[hap] = scrl_matrix_data(hap_iso_matrix_row[hap], hap_iso_matrix_col[hap], hap_iso_matrix_data[hap], all_bc, iso_names)
    return hap_to_gene_mtx_data, hap_to_trans_mtx_data

def make_matrix_data_with_em(all_bc, bc_to_cluster, bu_to_trans, gene_to_trans):
    gene_matrix_row, gene_matrix_col, gene_matrix_data, gene_names = [], [], [], []
    iso_matrix_row, iso_matrix_col, iso_matrix_data, iso_ratio_matrix_data, iso_names = [], [], [], [], []
    # print('bc_to_cluster:', bc_to_cluster)
    all_clusters = list(set(bc_to_cluster.values()))
    if all_clusters == []:
        all_clusters = [default_cluster]
    
    gene_i, trans_i = 1, 1  # matrix_row
    for gene, all_trans in gene_to_trans.items():
        all_trans = list(all_trans)
        gene_id, gene_name = gene.rsplit(',')
        bc_to_gene_abun = dd(lambda: 0)
        bc_to_trans_abun = dd(lambda: dd(lambda: 0))
        trans_total_cnt = dd(lambda: 0.0)
        for cluster in all_clusters:  # EM within each group and each gene
            bc_trans_abun = em_frac_abun1(gene_name, all_trans, bu_to_trans[gene], cluster, bc_to_cluster)
            for bc in bc_trans_abun:
                for trans in bc_trans_abun[bc]:
                    bc_to_gene_abun[bc] += bc_trans_abun[bc][trans]
                    bc_to_trans_abun[bc][trans] += bc_trans_abun[bc][trans]
                    trans_total_cnt[trans] += bc_trans_abun[bc][trans]/theta_base
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
            iso_names.append(trans + '\t' + trans)
        trans_i += len(all_trans)
    gene_mtx_data = scrl_matrix_data(gene_matrix_row, gene_matrix_col, gene_matrix_data, all_bc, gene_names)
    trans_mtx_data = scrl_matrix_data(iso_matrix_row, iso_matrix_col, iso_matrix_data, all_bc, iso_names)
    return gene_mtx_data, trans_mtx_data
    # write_10X_sparse_matrix(gene_matrix_row, gene_matrix_col, gene_matrix_data, gene_names, all_bc, out_dir+'/gene')
    # write_10X_sparse_matrix(iso_matrix_row, iso_matrix_col, iso_matrix_data, iso_names, all_bc, out_dir+'/transcript')
    # write_10X_sparse_matrix(iso_matrix_row, iso_matrix_col, iso_ratio_matrix_data, iso_names, all_bc, out_dir+'/transcript_ratio')


def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="{}: make sparsed count matrix at both gene and transcript levels".format(os.path.basename(__file__)))
    parser.add_argument('bu_tsv', metavar='bc_umi.tsv', type=str, help='Output file of NanoHunter')
    parser.add_argument('out_dir', metavar='matrix_output_dir', type=str, help='Folder to output single-cell gene/transcript matrix')
    parser.add_argument('-c', '--bc-cluster', type=str, default='', help='List file of cell barcode and corresponding cluster (can be used for cell cluster-wise transcript EM abundance estimation)')
    parser.add_argument('-a', '--hap-list', type=str, default='', help='List file of assigned haplotype for each read, generated by whatshap (can be used for haplotype-wise gene/transcript matrix generation)')
    parser.add_argument('-t', '--read-cate', type=str, default='', help='Only keep reads with specific transcript category, i.e., FSM. Use comma to separete multiple types. Default: keep all reads')
    parser.add_argument('-s', '--splice', default=False, action='store_true', help='Only keep reads with splice tag (ts from minimap2). Default: keep all reads')

    return parser.parse_args()


def make_10X_matrix_main(cell_to_cluster_tsv, hap_list_tsv, bu_fn, trans_cate, only_splice, out_mtx_dir, log_fn):
    err_log_format_time(log_fn, str='Writing cell-gene/transcript expression matrix ... ')
    bu_fn, out_mtx_dir = bu_fn, out_mtx_dir
    # add log info printout
    bc_to_cluster = parse_cell_cluster(cell_to_cluster_tsv)
    if not os.path.exists(out_mtx_dir):
        os.mkdir(out_mtx_dir)
    # write_gene_to_trans_tsv(gene_to_trans, out_dir+'/transcript_ratio')
    if hap_list_tsv:
        read_to_hap = get_reads_to_hap(hap_list_tsv)
        all_bc, bu_to_hap_trans, gene_to_trans = get_bu_hap_cmpt(bu_fn, trans_cate, only_splice, bc_to_cluster, read_to_hap)
    else:
        all_bc, bu_to_trans, gene_to_trans = get_bu_cmpt(bu_fn, trans_cate, only_splice, bc_to_cluster)
    err_log_format_time(log_fn, str='Reading barcode-UMI-transcript/gene count done!')

    if hap_list_tsv:
        hap_to_gene_mtx_data, hap_to_trans_mtx_data = make_matrix_data_with_hap_em(all_bc, bc_to_cluster, bu_to_hap_trans, gene_to_trans)
        err_log_format_time(log_fn, str="Collecting haplotype-wise gene/transcript expression matrix done!")
        for hap in ['H1', 'H2', 'none']:
            hap_to_gene_mtx_data[hap].write_matrix_data(out_mtx_dir+'/gene_haplotype_'+hap)
            hap_to_trans_mtx_data[hap].write_matrix_data(out_mtx_dir+'/transcript_haplotype_'+hap)
            write_gene_to_trans_tsv(gene_to_trans, out_mtx_dir+'/transcript_haplotype_'+hap)
    else:
        gene_mtx_data, trans_mtx_data = make_matrix_data_with_em(all_bc, bc_to_cluster, bu_to_trans, gene_to_trans) 
        err_log_format_time(log_fn, str='Collecting cell-gene/transcript expression matrix done!')
        gene_mtx_data.write_matrix_data(out_mtx_dir+'/gene')
        trans_mtx_data.write_matrix_data(out_mtx_dir+'/transcript')
        write_gene_to_trans_tsv(gene_to_trans, out_mtx_dir+'/transcript')
    err_log_format_time(log_fn, str='Writing cell-gene/transcript expression matrix done!')

def make_10X_matrix(cell_to_cluster_tsv, hap_list_tsv, scrl_para):
    bu_fn, out_mtx_dir = scrl_para.out_bu_fn, scrl_para.out_mtx_dir
    # XXX read type parameter in scrl_para
    make_10X_matrix_main(cell_to_cluster_tsv, hap_list_tsv, bu_fn, [], False, out_mtx_dir, scrl_para.log_fn)
    

def main():
    args = parser_argv()
    bu_fn, out_dir = args.bu_tsv, args.out_dir
    read_cate = args.read_cate.split(',')
    if '' in read_cate:
        read_cate.remove('')
    only_splice = args.splice
    bc_to_cluster_tsv = args.bc_cluster
    hap_list_tsv = args.hap_list
    make_10X_matrix_main(bc_to_cluster_tsv, hap_list_tsv, bu_fn, read_cate, only_splice, out_dir, None)

if __name__ == '__main__':
    main()
