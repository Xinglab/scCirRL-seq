import sys
import os
from collections import defaultdict as dd
from .parameter import scrl_para_class
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


def collect_trans_to_gene(scrl_para=scrl_para_class()):
    anno_gtf = scrl_para.anno_gtf
    updated_gtf = scrl_para.updated_gtf

    
    anno_gene_id_to_name = dict()
    trans_to_gene_id_name = dict()
    
    if os.path.exists(anno_gtf):
        err_log_format_time(scrl_para.log_fn, str='Collecting gene information from {}'.format(anno_gtf))
        anno_gene_id_to_name = get_gene_name(anno_gtf)
    else:
        err_log_format_time(scrl_para.log_fn, 'Warning', 'No annotation GTF file found.')
    if os.path.exists(updated_gtf):
        err_log_format_time(scrl_para.log_fn, str='Collecting transcript information from {}'.format(updated_gtf))
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
        err_log_format_time(scrl_para.log_fn, 'Warning', 'No updated GTF file found, no gene/transcript information will be output.')
    return trans_to_gene_id_name


def collect_trans_to_chrom_span(updated_gtf):
    """Return {trans_id: (chrom, start, end)} from transcript lines in *updated_gtf*.

    Used for region-aware filtering in collect_read_to_trans: only reads whose
    compatible transcripts overlap the requested genomic region are loaded.
    Coordinates are 1-based closed intervals as written in the GTF.
    """
    trans_to_span = {}
    if not updated_gtf or not os.path.exists(updated_gtf):
        return trans_to_span
    with open(updated_gtf) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.split('\t')
            if len(ele) < 5 or ele[2] != 'transcript':
                continue
            chrom, start, end = ele[0], int(ele[3]), int(ele[4])
            trans_id = ''
            if 'transcript_id ' in line:
                tmp = line[15 + line.index('transcript_id '):]
                trans_id = tmp[:tmp.index('"')]
            if trans_id:
                trans_to_span[trans_id] = (chrom, start, end)
    return trans_to_span


def _trans_in_region(trans_ids, trans_to_span, reg_chrom, reg_start, reg_end):
    """True if any transcript in *trans_ids* overlaps the region.

    *reg_start* / *reg_end* are 0-based half-open (BAM/pysam convention).
    GTF coordinates are 1-based closed, so overlap holds when
    gtf_start <= reg_end and gtf_end >= reg_start+1.
    When *reg_start* is None the check is chrom-only.
    """
    for tid in trans_ids:
        span = trans_to_span.get(tid)
        if span is None:
            continue
        chrom, t_start, t_end = span
        if chrom != reg_chrom:
            continue
        if reg_start is None:
            return True
        if t_start <= reg_end and t_end >= reg_start + 1:
            return True
    return False


# for Bambu or ESPRESSO, read only show up in one line, may follow by multiple isoforms
def collect_read_to_trans(trans_to_gene_id_name, scrl_para=scrl_para_class(),
                          region=None, trans_to_span=None):
    """Load read→transcript compatibility from *scrl_para.cmpt_tsv*.

    Parameters
    ----------
    region : str or None
        Genomic region string, e.g. ``"chr1"`` or ``"chr1:10000-50000"``.
        When provided, only reads whose compatible transcripts overlap the
        region are loaded.  Requires *trans_to_span*.
    trans_to_span : dict or None
        Mapping ``{trans_id: (chrom, start, end)}`` from the updated GTF,
        as returned by :func:`collect_trans_to_chrom_span`.  Only used when
        *region* is set.
    """
    cmpt_iso_fn = scrl_para.cmpt_tsv
    # cates = scrl_para.cate
    # NA_idx = 0
    read_to_trans = dd(lambda: [])
    read_to_cate = dd(lambda: 'NA')
    if not os.path.exists(cmpt_iso_fn):
        err_log_format_time(scrl_para.log_fn, 'Warning', 'No read-isoform compatible file found, no gene/transcript quantification will be output.')
        return read_to_trans, read_to_cate

    # parse region once
    reg_chrom, reg_start, reg_end = None, None, None
    if region and trans_to_span:
        if ':' in region:
            reg_chrom, span = region.split(':', 1)
            reg_start, reg_end = map(int, span.split('-'))
        else:
            reg_chrom = region

    err_log_format_time(scrl_para.log_fn, str='Collecting compatible transcripts from {}'.format(cmpt_iso_fn))
    with open(cmpt_iso_fn) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rstrip('\n').rsplit('\t')
            cmpt_trans_set = []
            cate = 'NA'
            # ESPRESSO: 4: readID, sample/bam, category, transcriptID(s)
            # Bambu: 3: readID, equalMatchesTrans, compatibleMatchesTrans
            if len(ele) == 3: # Bambu
                qname, equal_trans, cmpt_trans = ele[0], ele[1], ele[2]
                if qname == 'readId' or equal_trans == 'equalMatches' or cmpt_trans == 'compatibleMatches':
                    continue
                if equal_trans != 'NA':
                    cate = 'FullLen'
                    cmpt_trans_set.extend([equal_trans.rsplit(',')[0]]) # only keep one trans if multiple equal
                elif cmpt_trans != 'NA':
                    if cate == 'NA':
                        cate = 'Partial'
                    cmpt_trans_set.extend(cmpt_trans.rsplit(','))
                else:
                    continue
                if reg_chrom and not _trans_in_region(cmpt_trans_set, trans_to_span, reg_chrom, reg_start, reg_end):
                    continue
                read_to_trans[qname] = cmpt_trans_set
                read_to_cate[qname] = cate
            elif len(ele) == 4: # ESPRESSO
                qname, cate, cmpt_trans = ele[0], ele[2], ele[3]
                if cmpt_trans == 'NA':
                    continue
                # cmpt_trans += str(NA_idx)
                # NA_idx += 1
                cmpt_trans = cmpt_trans.rsplit(',')
                if '' in cmpt_trans:
                    cmpt_trans.remove('')
                # remove trans not in trans_to_gene_id_name
                for trans in cmpt_trans:
                    if trans in trans_to_gene_id_name:
                        cmpt_trans_set.append(trans)
                if reg_chrom and not _trans_in_region(cmpt_trans_set, trans_to_span, reg_chrom, reg_start, reg_end):
                    continue
                if cmpt_trans_set:
                    read_to_trans[qname] = cmpt_trans_set
                read_to_cate[qname] = cate
            else:
                err_log_format_time(scrl_para.log_fn, 'Warning', 'Unrecognized read-isoform compatible file format, no gene/transcript quantification will be output.')
                break
                # continue
    return read_to_trans, read_to_cate
