import sys

__program__ = 'scCirRL'
__module__ = 'scCirRL'
__scripts__ = [
               'scCirRL_allele_specific_expression',
               'scCirRL_allele_specific_splicing',
            #    'scCirRL_bulk_allele_specific_splicing',
               'scCirRL_cell_type_specific_splicing',
               'scCirRL_gene_with_gwas_disease',
               'scCirRL_matrix_to_clu_wise_expression',
               'scCirRL_split_bam_by_cluster',
               'scCirRL_split_chimeric_read',
               'scCirRL_make_expression_matrix',
               # split-parallel-merge barcode calling workflow
               'scCirRL_extract_reads',
               'scCirRL_collect_ref_bc',
               'scCirRL_assign_bc_from_tsv',
               'scCirRL_merge_bc_umi',
               ]
__version__ = '0.0.1'
__author__  = 'Yan Gao'
__email__   = 'gaoy1@chop.edu'
__description__ = 'Haplotype-resolved full-length transcriptome analysis in single cells by scCirRL-seq'
__url__ = 'https://github.com/Xinglab/scCirRL-seq'
__cmd__     = ' '.join(sys.argv)

