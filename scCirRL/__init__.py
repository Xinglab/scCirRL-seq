import sys

__program__ = 'scCirRL'
__module__ = 'scCirRL'
__scripts__ = ['cell_type_specific_splicing',
               'allele_specific_splicing',
               'bulk_allele_specific_splicing',
               'gene_with_gwas_disease',
               'make_expression_matrix',
               'split_bam_by_cluster']
__version__ = '1.0.0.beta3'
__author__  = 'Yan Gao'
__email__   = 'gaoy1@chop.edu'
__description__ = 'scCirRL: single-cell Characterization of isoform regulation using RNA-seq Long reads'
__url__ = 'https://github.com/Xinglab/scRMATS-long'
__cmd__     = ' '.join(sys.argv)

