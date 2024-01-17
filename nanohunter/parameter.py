import sys

_bc_len = 16
_bc_max_ed = 2
_umi_len = 12
_umi_max_ed = 1
_max_clip_len = 200
_five_ada = 'CTACACGACGCTCTTCCGATCT'

class nanohunter_para:
    def __init__(self):
        # input
        self.bc_len = 16
        self.umi_len = 10
        self.bc_max_ed = 2
        self.umi_max_ed = 1
        self.five_ada = 'CTACACGACGCTCTTCCGATCT'
        self.five_max_ed = int(len(self.five_ada) * 0.3)
        self.max_clip_len = 200
        self.all_alignments = False
        self.long_bam = ""
        self.updated_gtf = ""
        self.anno_gtf = ""
        self.cmpt_tsv = ""
        self.isoquant = False
        self.bc_list = ""

        # output
        self.out_dir = ""
        self.high_qual_bu_fn = ""
        self.out_ref_bc_fn = ""
        self.out_bu_fn = ""
        self.out_bu_bam = ""
        self.out_mtx_dir = ""
        self.log_fn = ""