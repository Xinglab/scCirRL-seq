import sys
import os
from collections import defaultdict as dd
from .parameter import nanohunter_para
from .utils import err_log_format_time

filtered_cates = []  # single-exon/NCD/NIC/NNC/ISM/FSM

# return: {trans: (chrom, [[exon1], [exon2] ... [exonN]])}
def get_trans_to_coor(gtf_fn):
    trans_to_coor = dict()
    with open(gtf_fn) as fp:
        trans_id, chrom, coors = '', '', []
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            if ele[2] == 'transcript':
                # end of last transcript
                if trans_id != '' and chrom != '' and coors != []:
                    # sort coors
                    coors.sort(key=lambda a: a[0])
                    trans_to_coor[trans_id] = (chrom, coors)
                trans_id, chrom, coors = '', '', []
                if 'transcript_id' in line:
                    trnas_ids = line[15+line.index('transcript_id'):]
                    trans_id = trnas_ids[:trnas_ids.index('"')]
            elif ele[2] == 'exon':
                chrom, start, end = ele[0], ele[3], ele[4]
                coors.append([int(start), int(end)])
            else:
                continue
    return trans_to_coor


def get_gene_name(in_gtf):
    if in_gtf == '':
        return dict()
    anno_gene_id_to_name = dict()
    with open(in_gtf) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            # if ele[2] != 'gene':
                # continue
            gene_id, gene_name = '', ''
            if 'gene_id ' in line:
                gene_ids = line[9+line.index('gene_id '):]
                gene_id = gene_ids[:gene_ids.index('"')]
            if gene_id in anno_gene_id_to_name:
                continue
            if 'gene_name ' in line:
                gene_names = line[11+line.index('gene_name '):]
                gene_name = gene_names[:gene_names.index('"')]
            if gene_id and gene_name:
                anno_gene_id_to_name[gene_id] = gene_name
    return anno_gene_id_to_name


def collect_trans_to_gene(nh_para=nanohunter_para()):
    anno_gtf = nh_para.anno_gtf
    updated_gtf = nh_para.updated_gtf

    
    anno_gene_id_to_name = dict()
    trans_to_gene_id_name = dict()
    
    if os.path.exists(anno_gtf):
        err_log_format_time(nh_para.log_fn, __name__, 'Collecting gene information from {}'.format(anno_gtf))
        anno_gene_id_to_name = get_gene_name(anno_gtf)
    else:
        err_log_format_time(nh_para.log_fn, 'Warning', 'No annotation GTF file found.')
    if os.path.exists(updated_gtf):
        err_log_format_time(nh_para.log_fn, __name__, 'Collecting transcript information from {}'.format(updated_gtf))
        with open(updated_gtf) as fp:
            for line in fp:
                if line.startswith('#'):
                    continue
                ele = line.rsplit()
                if len(ele) < 3:
                    continue
                if ele[2] != 'transcript':
                    continue
                trans_id, gene_id, gene_name = '', '', ''
                if 'transcript_id ' in line:
                    trnas_ids = line[15+line.index('transcript_id '):]
                    trans_id = trnas_ids[:trnas_ids.index('"')]
                if 'gene_id ' in line:
                    gene_ids = line[9+line.index('gene_id '):]
                    gene_id = gene_ids[:gene_ids.index('"')]
                    if ',' in gene_id: # filter out transcript with multiple gene_id
                        continue
                if 'gene_name ' in line:
                    gene_names = line[11+line.index('gene_name '):]
                    gene_name = gene_names[:gene_names.index('"')]
                elif gene_id in anno_gene_id_to_name:
                    gene_name = anno_gene_id_to_name[gene_id]
                else:
                    gene_name = gene_id  # no gene name available
                if trans_id and gene_id and gene_name:
                    trans_to_gene_id_name[trans_id] = {'id': gene_id, 'name': gene_name}
    else:
        err_log_format_time(nh_para.log_fn, 'Warning', 'No updated GTF file found, no gene/transcript information will be output.')
    return trans_to_gene_id_name


# for ISOQUANT, read may show up in multiple lines with different isoforms
# for ESPRESSO, read only show up in one line, may follow by multiple isoforms
def collect_read_to_trans(trans_to_gene_id_name, nh_para=nanohunter_para()):
    cmpt_iso_fn = nh_para.cmpt_tsv
    is_isoquant = nh_para.isoquant
    # NA_idx = 0
    if not os.path.exists(cmpt_iso_fn):
        err_log_format_time(nh_para.log_fn, 'Warning', 'No read-isoform compatible file found, no gene/transcript quantification will be output.')
        return None
    read_to_trans = dd(lambda: [])
    err_log_format_time(nh_para.log_fn, __name__, 'Collecting compatible transcripts from {}'.format(cmpt_iso_fn))
    with open(cmpt_iso_fn) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            if len(ele) < 2:
                continue
            if not is_isoquant:  # ESPRESSO file
                if len(ele) != 4:
                    continue
                qname, cate, cmpt_trans = ele[0], ele[2], ele[3]
                if cate in filtered_cates:
                    continue
                if cmpt_trans == 'NA':
                    continue
                    # cmpt_trans += str(NA_idx)
                    # NA_idx += 1
                cmpt_trans = cmpt_trans.rsplit(',')
                if '' in cmpt_trans:
                    cmpt_trans.remove('')
                # remove trans not in trans_to_gene_id_name
                cmpt_trans_set = []
                for trans in cmpt_trans:
                    if trans in trans_to_gene_id_name:
                        cmpt_trans_set.append(trans)
                if cmpt_trans_set:
                    read_to_trans[qname] = cmpt_trans_set
            else:  # ISOQUANT
                if len(ele) != 2:
                    continue
                qname, trans = ele[0], ele[1]
                if trans == '*' or trans not in trans_to_gene_id_name: # remove trans not in trans_to_gene_id_name
                    continue
                read_to_trans[qname].append(trans)
    return read_to_trans
