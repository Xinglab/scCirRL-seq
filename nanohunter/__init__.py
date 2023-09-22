import sys

# from .utils import err_format_time, err_fatal_format_time
# from .seq_utils import _bc_max_ed, _bc_len, _umi_len, get_mp_fetch_set
# from .collect_candidate_bc import collect_cand_ref_bc
# from .assign_bc import assign_ref_bc
# from .parse_quant_file import get_read_to_trans, get_trans_to_gene

__program__ = 'nanohunter'
__module__ = 'nanohunter'
__scripts__ = ['nh_cell_type_specific_splicing',
               'nh_allele_specific_splicing',
               'nh_bulk_allele_specific_splicing',
               'nh_gene_with_gwas_disease_snp',
               'nh_make_expression_matrix',
               'nh_split_bam_by_cluster']
__version__ = '1.0.0'
__author__  = 'Yan Gao'
__email__   = 'gaoy1@chop.edu'
__description__ = 'long-read single-cell RNA-seq analysis'
__url__ = 'https://github.com/Xinglab/NanoHunter'
__cmd__     = ' '.join(sys.argv)

