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
               ]
__version__ = '1.0.0.beta3'
__author__  = 'Yan Gao'
__email__   = 'gaoy1@chop.edu'
__description__ = 'Haplotype-resolved full-length transcriptome analysis in single cells by scCirRL-seq'
__url__ = 'https://github.com/Xinglab/scCirRL-seq'
__cmd__     = ' '.join(sys.argv)

