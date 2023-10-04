import sys
import os
import argparse
from collections import defaultdict as dd
import vcf # pyvcf
# for GWAS/LD database
import duckdb as db

filtered_out_gene_prefix = 'HLA-'
ld_populations = ['EUR']
# ld_populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
ld_table_prefix = 'tkg_p3v5a_ld' # tkg_p3v5a_ld_chr1_AFR
ld_threshold = 0.8
promoter_upstream_len = 500
# flanking_len = 0

ALL_HAPS = ['H1', 'H2']  # 'none'

# order of GWAS disease/traits super-term in the database
gwas_disease_list = [
'Cancer',
'Immune system disorder',
'Neurological disorder',
'Cardiovascular disease',
'Digestive system disorder',
'Metabolic disorder',
'Response to drug',
'Other disease',
# 'Biological process',
# 'Body measurement',
# 'Cardiovascular measurement',
# 'Hematological measurement',
# 'Inflammatory measurement',
# 'Lipid or lipoprotein measurement',
# 'Liver enzyme measurement',
# 'Other measurement',
# 'Other trait',
# 'NR', # not reported
# 'NA'
]
n_gwas_disease_traits = len(gwas_disease_list)
gwas_traits_order = dd(lambda: len(gwas_disease_list))
for i, trait in enumerate(gwas_disease_list):
    gwas_traits_order[trait] = i

# gwas_trait_info : gwas_snp/ld_gwas_snp -> trait_info
# each SNP may map to multiple traits
class gwas_trait_info:
    def __init__(self, snp, trait='', parent='', pval=None):
        # snp id
        self.snp = snp
        # leading parent term
        self.leading_parent_term_order = None
        self.leading_trait_term = None
        self.leading_parent_term = None
        self.leading_pval = None
        # all traits
        self.n_traits = 0
        self.n_parent_terms = 0
        self.traits = dd(lambda: dd(lambda: None)) # trait -> {parent: parent, pval: pval}
        self.parents = dd(lambda: None) # parent -> parent_term_order
        
        if trait or parent or pval:
            if parent not in gwas_traits_order:
                sys.stderr.write(f'Warning: parent=\"{parent}\" not in gwas_traits_order\n')
            self.leading_parent_term_order = gwas_traits_order[parent]
            self.leading_trait_term = trait
            self.leading_parent_term = parent
            self.leading_pval = pval
            
            self.n_traits = 1
            self.traits[trait]['parent'] = parent
            self.traits[trait]['pval'] = pval
            self.n_parent_terms = 1
            self.parents[parent] = gwas_traits_order[parent]
        
    def __str__(self):
        return(f'{self.leading_parent_term}({self.leading_trait_term}), P={self.leading_pval}')
    def add_trait(self, trait='', parent='', pval=None):
        if not trait or not parent or not pval:
            sys.stderr.write(f'Error: trait={trait}, parent={parent}, pval={pval}\n')
            return
        if parent not in gwas_traits_order:
            sys.stderr.write(f'Warning: parent=\"{parent}\" not in gwas_traits_order\n')
        
        if trait not in self.traits:
            self.traits[trait]['parent'] = parent
            self.traits[trait]['pval'] = pval
            self.n_traits += 1
            
            # update parent term
            if parent not in self.parents:
                self.parents[parent] = gwas_traits_order[parent]
                self.n_parent_terms += 1
                # update leading parent term
                if not self.leading_parent_term_order or \
                    gwas_traits_order[parent] < self.leading_parent_term_order:
                    self.leading_parent_term_order = gwas_traits_order[parent]
                    self.leading_trait_term = trait
                    self.leading_parent_term = parent
                    self.leading_pval = pval

        elif pval < self.traits[trait]['pval']:
            if parent != self.traits[trait]['parent']:
                sys.stderr.write(f'''Warning: \"new parent={parent}\" != \"old parent={self.traits[trait]['parent']}\"\n''')
                sys.stderr.write(f'''         old parent term will be kept\n''')
            self.traits[trait]['pval'] = pval
            # update leading parent term pval
            if self.traits[trait]['parent'] == self.leading_parent_term and pval < self.leading_pval:
                self.leading_pval = pval

# ld_snp_info     : ld_snp -> ld_info
class ld_snp_info:
    def __init__(self, snp):
        self.snp = snp
        self.ld_snps = dd(lambda: dd(lambda: 0)) # {ld_snp: {population: r2}}
    def add_ld_snp(self, ld_snp, ld_population='', r2=None):
        if ld_population not in self.ld_snps[ld_snp]:
            self.ld_snps[ld_snp][ld_population] = r2
        return

# snp_info    : all hap snps -> GWAS/LD_snp
class snp_info:
    def __init__(self, snp, gene='', chrom='', pos=-1, hap=None):
        self.snp = snp
        self.gene = gene
        self.chrom = chrom
        self.pos = pos
        self.hap = hap
        self._is_gwas_snp = None # None if not searched yet
        self.has_ld_snps = False # None if not searched yet
        self.ld_info = ld_snp_info(snp)
    def is_gwas_snp(self, region_gwas_snps):
        if self._is_gwas_snp is None:
            if self.snp in region_gwas_snps:
                self._is_gwas_snp = True
            else:
                self._is_gwas_snp = False
        return(self._is_gwas_snp)
    def add_ld_snp(self, ld_snp, population, r2):
        if ld_snp and population and r2:
            self.has_ld_snps = True
            self.ld_info.add_ld_snp(ld_snp, population, r2)
    def get_snp_id(self):
        return(self.snp)
    def get_chrom(self):
        return(self.chrom)
    def get_pos(self):
        return(self.pos)
    def get_hap(self):
        return(self.hap)
    # trait/trait_parent_term: combine all traits from the snp and its LD snps, keep the one with highest order
    # gwas_trait: combine all traits from the snp, if it is GWAS
    # ld_snps/ld_traits/ld_scores: all LD snps with LD score > 0.8, for each LD snps, keep the trait with highest order
    def get_leading_trait(self, region_gwas_snps):
        gwas_trait, gwas_trait_order = self.get_gwas_trait(region_gwas_snps)
        ld_trait, ld_trait_order = self.get_ld_trait(region_gwas_snps)
        return gwas_trait if gwas_trait_order <= ld_trait_order else ld_trait
    
    def get_leading_trait_parent(self, region_gwas_snps):
        gwas_trait_parent, gwas_trait_order = self.get_gwas_trait_parent(region_gwas_snps)
        ld_trait_parent, ld_trait_order = self.get_ld_trait_parent(region_gwas_snps)
        return gwas_trait_parent if gwas_trait_order <= ld_trait_order else ld_trait_parent
    def get_gwas_trait(self, region_gwas_snps):
        if self.is_gwas_snp(region_gwas_snps):
            return(region_gwas_snps[self.snp].leading_trait_term, 
                   region_gwas_snps[self.snp].leading_parent_term_order)
        else:
            return "NA", len(gwas_traits_order)
    def get_all_lds(self, region_gwas_snps):
        if self.has_ld_snps:
            ld_info = self.ld_info.ld_snps
            ld_snps, ld_traits, ld_scores, ld_parent_terms = [], [], [], []
            for ld_snp in ld_info:
                ld_snps.append(ld_snp)
                ld_scores.append(max(ld_info[ld_snp].values()))
                ld_traits.append(region_gwas_snps[ld_snp].leading_trait_term) \
                    if ld_snp in region_gwas_snps \
                    else ld_traits.append("NA")
                ld_parent_terms.append(region_gwas_snps[ld_snp].leading_parent_term) \
                    if ld_snp in region_gwas_snps \
                    else ld_parent_terms.append("NA")
            return ld_snps, ld_traits, ld_scores, ld_parent_terms
        else:
            return ["NA"], ["NA"], ["NA"], ["NA"]
    def get_ld_trait(self, region_gwas_snps):
        min_order = len(gwas_traits_order)
        min_order_trait = 'NA'
        ld_info = self.ld_info.ld_snps
        for ld_snp in ld_info:
            if ld_snp in region_gwas_snps:
                order = region_gwas_snps[ld_snp].leading_parent_term_order
                trait = region_gwas_snps[ld_snp].leading_trait_term
                if order < min_order:
                    min_order = order
                    min_order_trait = trait
        return min_order_trait, min_order
    def get_ld_trait_parent(self, region_gwas_snps):
        min_order = len(gwas_traits_order)
        min_order_trait_parent = 'NA'
        ld_info = self.ld_info.ld_snps
        for ld_snp in ld_info:
            if ld_snp in region_gwas_snps:
                order = region_gwas_snps[ld_snp].leading_parent_term_order
                trait = region_gwas_snps[ld_snp].leading_parent_term
                if order < min_order:
                    min_order = order
                    min_order_trait_parent = trait
        return min_order_trait_parent, min_order
    def get_gwas_trait_parent(self, region_gwas_snps):
        if self.is_gwas_snp(region_gwas_snps):
            return(region_gwas_snps[self.snp].leading_parent_term, 
                   region_gwas_snps[self.snp].leading_parent_term_order)
        else:
            return 'NA', len(gwas_traits_order)
    def get_gwas_pval(self, region_gwas_snps):
        if self.is_gwas_snp(region_gwas_snps):
            return(region_gwas_snps[self.snp].leading_pval)
        else:
            return(None)
    # 1. gwas trait of itself
    # 2. leading gwas trait of each of its LD snps
    def get_all_traits(self, region_gwas_snps):
        all_traits = set()
        if self.is_gwas_snp(region_gwas_snps):
            all_traits.add(region_gwas_snps[self.snp].leading_trait_term)
        ld_info = self.ld_info.ld_snps
        for ld_snp in ld_info:
            if ld_snp in region_gwas_snps:
                all_traits.add(region_gwas_snps[ld_snp].leading_trait_term)
        return all_traits
    
    def get_all_parents(self, region_gwas_snps):
        all_parents = set()
        if self.is_gwas_snp(region_gwas_snps):
            all_parents.add(region_gwas_snps[self.snp].leading_parent_term)
        ld_info = self.ld_info.ld_snps
        for ld_snp in ld_info:
            if ld_snp in region_gwas_snps:
                all_parents.add(region_gwas_snps[ld_snp].leading_parent_term)
        return all_parents

def open_vcf_file(vcf_fn):
    if not os.path.exists(vcf_fn + '.tbi'):
        cmd = 'tabix -p vcf {}'.format(vcf_fn)
        os.system(cmd)
    vcf_fp = open(vcf_fn, 'rb') if vcf_fn.endswith('.gz') else open(vcf_fn, 'r')
    return vcf.Reader(vcf_fp)

def set_gwas_trait(gwas_query_res, all_gwas_snp_to_trait, region_prefix):
    for region_ids, snps, trait, parent, pvalue in gwas_query_res:
        if not region_ids or not snps or not trait or not parent or not pvalue:
            continue
        if parent not in gwas_disease_list: # only keep disease related traits
            continue
        region_ids = region_ids.replace('-', ',')
        region_ids = region_ids.replace('x', ',')
        snps = snps.replace('-', ',')
        snps = snps.replace('x', ',')
        for region_id in ''.join(region_ids.split()).split(','):
            for snp in ''.join(snps.split()).split(','):
                region_id = region_prefix + region_id
                if snp not in all_gwas_snp_to_trait[region_id]:
                    all_gwas_snp_to_trait[region_id][snp] = gwas_trait_info(snp, trait, parent, pvalue)
                else:
                    all_gwas_snp_to_trait[region_id][snp].add_trait(trait, parent, pvalue)

def parse_gwas_file(gwas_fn, genes=[]):
    gwas_tb_fn = 'gwas_tb'
    db.sql(f"CREATE OR REPLACE TABLE {gwas_tb_fn} AS SELECT * FROM read_csv_auto('{gwas_fn}', nullstr='NA')")

    all_gwas_snp_to_trait = dd(lambda: dd(lambda: None)) # chrom -> snp -> gwas_trait_info
    
    if not genes: # collec all genes
        gwas_cols = ['CHR_ID', 'SNPS', '"EFO term"', '"Parent term"', '"P-VALUE"'] # MAPPED_GENE: gene_name
        region_prefix = 'chr'
        
        gwas_query_cmd = f'''SELECT {','.join(gwas_cols)} FROM gwas_tb '''
        res = db.sql(gwas_query_cmd).fetchall()
        set_gwas_trait(res, all_gwas_snp_to_trait, region_prefix)
    else:
        gwas_cols = ['"MAPPED_GENE"', 'SNPS', '"EFO term"', '"Parent term"', '"P-VALUE"']
        region_prefix = ''
        for gene in genes:
            gwas_query_cmd = f'''SELECT {','.join(gwas_cols)} FROM gwas_tb WHERE MAPPED_GENE LIKE '%{gene}%' '''
            res = db.sql(gwas_query_cmd).fetchall()
            set_gwas_trait(res, all_gwas_snp_to_trait, region_prefix)
    return all_gwas_snp_to_trait

def get_gene_to_coor(gtf_fn, all_genes=[]):
    gene_id_to_name = dd(lambda: '')
    gene_to_coor = dd(lambda: ('', sys.maxsize, -sys.maxsize), '')
    with open(gtf_fn) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.rsplit()
            # 'transcript' for gene
            # 'exon' for transcript
            gene_id, gene_name, trans_id = '', '', ''
            if ele[2] != 'transcript':
                continue
            if 'gene_id' in line:
                gene_ids = line[9+line.index('gene_id'):]
                gene_id = gene_ids[:gene_ids.index('"')]
            if 'gene_name' in line:
                gene_names = line[11+line.index('gene_name'):]
                gene_name = gene_names[:gene_names.index('"')]
            if all_genes != [] and (gene_id not in all_genes and gene_name not in all_genes):
                continue
            if 'transcript_id' in line:
                trans_ids = line[15+line.index('transcript_id'):]
                trans_id = trans_ids[:trans_ids.index('"')]
            if gene_id == 'NA' or not gene_id or not trans_id or not gene_name:
                continue
            chrom, start, end, strand = ele[0], int(ele[3]), int(ele[4]), ele[6]
            gene_id_to_name[gene_id] = gene_name
            gene_to_coor[gene_id] = (chrom, 
                                     min(gene_to_coor[gene_id][1], start), 
                                     max(gene_to_coor[gene_id][2], end), 
                                     strand)
    return gene_id_to_name, gene_to_coor

snp_id_pos_table = 'tkg_p3v5a_hg38'
def get_snp_id_from_ld_table(ld_con, chrom, pos):
    try:
        res = ld_con.execute(f"SELECT rsid FROM {snp_id_pos_table} WHERE chr = '{chrom}' AND pos = {pos}").fetchone()
    except:
        res = None
    if res:
        return res[0]
    else:
        return None

def get_ld_snps(ld_cand_snps, gwas_snps, ld_con, ld_table_prefix, ld_populations, chrom, ld_threshold):
    gwas_snps_list = list(gwas_snps.keys())
    gwas_snp_str = ','.join(["'"+snp_id+"'" for snp_id in gwas_snps_list])
    ld_cand_snp_str = ','.join(["'"+snp_id+"'" for snp_id in ld_cand_snps])
    
    hap_snp_to_ld_gwas_snps = set()
    # for i in range(0, len(gwas_snps_list), max_gwas_snp_chunk):
        # gwas_snp_str = ','.join(["'"+snp_id+"'" for snp_id in gwas_snps_list[i:i+max_gwas_snp_chunk]])
    for population in ld_populations:
        if len(ld_cand_snps) == 1:
            ld_query_cmd = f"SELECT SNP_A, SNP_B, R2 FROM {ld_table_prefix}_{chrom}_{population} WHERE \
                        R2 > {ld_threshold} AND \
                        ( (SNP_A = \'{ld_cand_snps[0]}\' AND SNP_B in ({gwas_snp_str})) OR \
                          (SNP_B = \'{ld_cand_snps[0]}\' AND SNP_A in ({gwas_snp_str})) )"
        else:
            ld_query_cmd = f"SELECT SNP_A, SNP_B, R2 FROM {ld_table_prefix}_{chrom}_{population} WHERE \
                        R2 > {ld_threshold} AND \
                        ( (SNP_A in ({ld_cand_snp_str}) AND SNP_B in ({gwas_snp_str})) OR \
                          (SNP_B in ({ld_cand_snp_str}) AND SNP_A in ({gwas_snp_str})) )"
        ld_res = ld_con.execute(ld_query_cmd).fetchall()
        for snp_a, snp_b, r2 in ld_res:
            if snp_a in ld_cand_snps:
                hap_snp_to_ld_gwas_snps.add((snp_a, snp_b, population, r2))
            if snp_b in ld_cand_snps:
                hap_snp_to_ld_gwas_snps.add((snp_b, snp_a, population, r2))
    return hap_snp_to_ld_gwas_snps

def get_ld_scores(snps1, snps2, ld_con, ld_table_prefix, ld_populations, chrom, ld_threshold):
    if len(snps1) == 0 or len(snps2) == 0:
        return None
    ld_scores = dd(lambda: dd(lambda: 0))
    snp_str1 = ','.join(["'"+snp_id+"'" for snp_id in snps1])
    snp_str2 = ','.join(["'"+snp_id+"'" for snp_id in snps2])
    for population in ld_populations:
        if len(snps1) == 1 and len(snps2) == 1:
            ld_query_cmd = f"SELECT SNP_A, SNP_B, R2 FROM {ld_table_prefix}_{chrom}_{population} WHERE \
                             R2 > {ld_threshold} AND \
                             ( (SNP_A = \'{snps1[0]}\' AND SNP_B = \'{snps2[0]}\') OR \
                               (SNP_B = \'{snps1[0]}\' AND SNP_A = \'{snps2[0]}\') )"
        elif len(snps1) == 1 and len(snps2) > 1:
            ld_query_cmd = f"SELECT SNP_A, SNP_B, R2 FROM {ld_table_prefix}_{chrom}_{population} WHERE \
                             R2 > {ld_threshold} AND \
                             ( (SNP_A = \'{snps1[0]}\' AND SNP_B in ({snp_str2})) OR \
                               (SNP_B = \'{snps1[0]}\' AND SNP_A in ({snp_str2})) )"
        elif len(snps1) > 1 and len(snps2) == 1:
            ld_query_cmd = f"SELECT SNP_A, SNP_B, R2 FROM {ld_table_prefix}_{chrom}_{population} WHERE \
                             R2 > {ld_threshold} AND \
                             ( (SNP_A in ({snp_str1}) AND SNP_B = \'{snps2[0]}\') OR \
                               (SNP_B in ({snp_str1}) AND SNP_A = \'{snps2[0]}\') )"
        else:
            ld_query_cmd = f"SELECT SNP_A, SNP_B, R2 FROM {ld_table_prefix}_{chrom}_{population} WHERE \
                             R2 > {ld_threshold} AND \
                             ( (SNP_A in ({snp_str1}) AND SNP_B in ({snp_str2})) OR \
                               (SNP_B in ({snp_str1}) AND SNP_A in ({snp_str2})) )"
        ld_res = ld_con.execute(ld_query_cmd).fetchall()
        for snp_a, snp_b, r2 in ld_res:
            if r2 > ld_scores[snp_a][snp_b]:
                ld_scores[snp_a][snp_b] = r2
                ld_scores[snp_b][snp_a] = r2
                ld_scores[snp_a][snp_a] = 1
                ld_scores[snp_b][snp_b] = 1
    return ld_scores

def set_gwas_ld(snp_to_snp_info, gwas_snps, chrom, ld_populations, ld_table_prefix, ld_threshold, ld_con):
    ld_cand_snps = []
    for snp, snp_info in snp_to_snp_info.items():
        ld_cand_snps.append(snp)
    # check if SNP is in GWAS
    for snp in snp_to_snp_info:
        snp_info = snp_to_snp_info[snp]
        if snp in gwas_snps:
            snp_info._is_gwas_snp = True
        else:
            snp_info._is_gwas_snp = False
    # check if SNP is in LD with GWAS
    # print(f"LD searching for chromosome {chrom}")
    if len(ld_cand_snps) > 0 and len(gwas_snps) > 0:
        ld_gwas_snps = get_ld_snps(ld_cand_snps, gwas_snps, ld_con, ld_table_prefix, ld_populations, chrom, ld_threshold)
        for cand_snp, ld_snp, population, r2 in ld_gwas_snps:
            snp_to_snp_info[cand_snp].add_ld_snp(ld_snp, population, r2)

#   collect all SNPs within the gene for each haplotype
def collect_gene_snps(all_gene_id_to_name, gene_to_coor, vcf_fn, ld_duckdb_fn, sample_name):
    all_gene_to_snps = dd(lambda: dd(lambda: dd(lambda: None))) # chr -> gene -> snp_id -> snp_info
    # 1. open vcf file
    vcf_reader = open_vcf_file(vcf_fn)
    if not sample_name:
        if len(vcf_reader.samples) < 1:
            sys.stderr.write('Error: no sample found in the VCF file.\n')
        sample_name = vcf_reader.samples[0]
        sys.stderr.write('No sample name is provided. The first sample ({}) in VCF file will be used.\n'.format(sample_name))
    else:
        if sample_name not in vcf_reader.samples:
            sys.stderr.write('Error: sample name {} is not in VCF file.\n'.format(sample_name))
            sys.exit(1)
    with db.connect(ld_duckdb_fn, read_only=True) as ld_con:
        for gene_id, gene_name in all_gene_id_to_name.items():
            chrom, start, end, strand = gene_to_coor[gene_id]
            if strand == '+':
                start -= promoter_upstream_len
            elif strand == '-':
                end += promoter_upstream_len
            # 2. collect all SNPs within the gene
            try:
                snp_records = vcf_reader.fetch(chrom, max(0, start), end)
            except:
                snp_records = []
            for snp_record in snp_records:
                if not snp_record.ID:
                    chrom_without_chr = chrom.replace('chr', '') if 'chr' in chrom else chrom
                    snp_id = get_snp_id_from_ld_table(ld_con, chrom_without_chr, snp_record.POS)
                    if not snp_id:
                        continue
                else:
                    snp_id = snp_record.ID
                pos = snp_record.POS
                genotype = snp_record.genotype(sample_name)['GT']
                if genotype == '1|0':
                    hap = ALL_HAPS[0]
                elif genotype == '0|1':
                    hap = ALL_HAPS[1]
                else:
                    continue
                all_gene_to_snps[chrom][gene_name][snp_id] = snp_info(snp_id, gene_name, chrom, pos, hap)
    return all_gene_to_snps

# for each perform GWAS analysis
# . 1. check if its SNP is in GWAS
# . 2. if not, check if as-trans-SNP is in LD with GWAS
def collect_gene_to_gwas_traits(all_gene_to_snps, chr_to_all_gwas_snp_trait, ld_duckdb_fn):
    gene_to_gwas_traits = dd(lambda: dd(lambda: False)) # gene_name ->  [traits] : T/F
    gene_to_gwas_parent_traits = dd(lambda: dd(lambda: False)) # gene_name -> [traits]: T/F
    with db.connect(ld_duckdb_fn, read_only=True) as ld_con:
        for chrom, gene_to_snps in all_gene_to_snps.items():
            # for each gene, collect associated disease traits
            gwas_snp_to_trait = chr_to_all_gwas_snp_trait[chrom]
            for gene, all_snps in gene_to_snps.items():
                set_gwas_ld(all_snps, gwas_snp_to_trait, chrom, ld_populations, ld_table_prefix, ld_threshold, ld_con)
                for snp in all_snps:
                    snp_info = all_snps[snp]
                    # collect GWAS traits & parent terms
                    traits = snp_info.get_all_traits(gwas_snp_to_trait)
                    parents = snp_info.get_all_parents(gwas_snp_to_trait)
                    for trait in traits:
                        gene_to_gwas_traits[gene][trait] = True
                    for parent in parents:
                        gene_to_gwas_parent_traits[gene][parent] = True
    return gene_to_gwas_traits, gene_to_gwas_parent_traits

def output_gwas(all_gene_id_to_name, gene_to_gwas_traits, gene_to_gwas_parent_traits, as_gwas_basic_out, as_gwas_parent_out):
    all_traits, all_parents = set(), set()
    all_gene_name_to_id = {v: k for k, v in all_gene_id_to_name.items()}
    for gene in gene_to_gwas_traits:
        all_traits.update(gene_to_gwas_traits[gene].keys())
        all_parents.update(gene_to_gwas_parent_traits[gene].keys())
    with open(as_gwas_basic_out, 'w') as gwas_fp, open(as_gwas_parent_out, 'w') as gwas_parent_fp:
        gwas_fp.write('gene_name\tgene_id\t'+ '\t'.join(all_traits) + '\n')
        gwas_parent_fp.write('gene_name\tgene_id\t' + '\t'.join(all_parents) + '\n')
        for gene in gene_to_gwas_traits:
            gwas_fp.write('\t'.join([gene, all_gene_name_to_id[gene]])+'\t')
            gwas_fp.write('\t'.join([str(gene_to_gwas_traits[gene][trait]) for trait in all_traits]) + '\n')

            gwas_parent_fp.write('\t'.join([gene, all_gene_name_to_id[gene]])+'\t')
            gwas_parent_fp.write('\t'.join([str(gene_to_gwas_parent_traits[gene][parent]) for parent in all_parents]) + '\n')

def output_ld(all_gene_id_to_name, all_gene_to_snps, chr_to_all_gwas_snp_trait, as_gwas_ld_out_dir, ld_duckdb_fn):
    if not os.path.exists(as_gwas_ld_out_dir):
        os.makedirs(as_gwas_ld_out_dir)
    ld_header = ['haplotype', 'chrom', 'pos', 'type', 'gwas_trait', 'gwas_trait_efo_term', 'gwas_pvalue']
    gwas_detail_header = ['snp_id', 'gene_name', 'gene_id', 'haplotype', 'gwas_trait', 'gwas_trait_efo_term', 'ld_trait', 'ld_efo_term', 'ld_snps', 'ld_scores']
    all_gene_name_to_id = {v: k for k, v in all_gene_id_to_name.items()}
    with db.connect(ld_duckdb_fn, read_only=True) as ld_con: #, open(as_gwas_detail_out, 'w') as gwas_detail_fp:
        # gwas_detail_fp.write('{}\n'.format('\t'.join(gwas_detail_header)))
        for chrom, gene_to_snps in all_gene_to_snps.items():
            gwas_snp_to_trait = chr_to_all_gwas_snp_trait[chrom]
            for gene_name, all_snps in gene_to_snps.items():
                gene_id = all_gene_name_to_id[gene_name]
                snps = list(all_snps.keys())
                ld_scores = get_ld_scores(snps, snps, ld_con, ld_table_prefix, ld_populations, chrom, ld_threshold=0.1)
                # write ld scores
                with open(os.path.join(as_gwas_ld_out_dir, '{}.ld'.format(gene_name)), 'w') as ld_fp:
                    ld_fp.write('snp_id\t{}\t{}\n'.format('\t'.join(snps), '\t'.join(ld_header)))
                    for snp1 in snps:
                        ld_fp.write(f'{snp1}\t')
                        ld_fp.write('{}\t'.format('\t'.join([str(ld_scores[snp1][snp2]) for snp2 in snps])))
                        snp_info = all_snps[snp1]
                        hap = snp_info.get_hap()
                        pos = snp_info.get_pos()
                        if snp_info.is_gwas_snp(gwas_snp_to_trait):
                            snp_type = 'gwas'
                        elif snp_info.has_ld_snps:
                            snp_type = 'ld_with_gwas'
                        else:
                            snp_type = 'None'
                        gwas_trait = snp_info.get_gwas_trait(gwas_snp_to_trait)[0]
                        gwas_trait_parent = snp_info.get_gwas_trait_parent(gwas_snp_to_trait)[0]
                        gwas_pvalue = snp_info.get_gwas_pval(gwas_snp_to_trait)
                        ld_fp.write('\t'.join([hap, chrom, str(pos), snp_type, gwas_trait, gwas_trait_parent, str(gwas_pvalue)]) + '\n')
                # write gwas details for each snp
                # for snp in snps:
                #     snp_info = all_snps[snp]
                #     hap = snp_info.get_hap()
                #     gwas_trait = snp_info.get_gwas_trait(gwas_snp_to_trait)[0]
                #     gwas_trait_parent = snp_info.get_gwas_trait_parent(gwas_snp_to_trait)[0]
                #     ld_snps, ld_traits, ld_scores, ld_parents = snp_info.get_all_lds(gwas_snp_to_trait)
                #     gwas_detail_fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                #         snp, gene_name, gene_id, hap, gwas_trait, gwas_trait_parent, \
                #         ','.join(ld_traits), ','.join(ld_parents), ','.join(ld_snps), ','.join([str(score) for score in ld_scores])))


def get_gene_list(gene_list_file):
    all_genes = dict()
    header_idx = dict()
    with open(gene_list_file) as fp:
        header = True
        for line in fp:
            ele = line.rsplit('\t')
            if header:
                header = False
                header_idx = {ele[i].upper(): i for i in range(len(ele))}
                continue
            for h1 in ['GENE', 'GENE_ID', 'GENEID', 'GENE_NAME', 'GENENAME']:
                if h1 in header_idx:
                    all_genes[ele[header_idx[h1]]] = True
                    break
    return list(all_genes.keys())

# if no gene_list_fn provided, collect all genes from GTF file
# 1. collect gene_id to gene_name
# 2. collect gene_id to gene_coordinate
def get_gene_infor(gene_list_fn, gtf_fn):
    all_genes = []
    if gene_list_fn:
        all_genes = get_gene_list(gene_list_fn)
    all_gene_id_to_name, all_gene_id_to_coor = get_gene_to_coor(gtf_fn, all_genes)
    return all_gene_id_to_name, all_gene_id_to_coor
    
def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description="{}: identify disease-associated GWAS SNPs".format(os.path.basename(__file__)))
    # parser.add_argument('gene_list', metavar='gene_list.tsv', type=str, help='gene list, use comma to separate multiple gene list files')
    parser.add_argument('gtf', metavar='gtf', type=str, help='GTF annotation')
    parser.add_argument('vcf', metavar='snp.vcf', type=str, help='VCF file containing SNPs of the sample')
    parser.add_argument('gwas', metavar='gwas.tsv', type=str, help='GWAS file with gene/SNP/trait information')
    parser.add_argument('ld_duckdb', metavar='ld.duckdb', type=str, help='SNPs LD file in duckdb format')
    parser.add_argument('out_pre', metavar='out_pre', type=str, help='Output prefix')
    
    parser.add_argument('-g', '--gene-list', metavar='gene_list', type=str, default='',
                             help='File of gene list (ID or name). Only SNPs within the genes will be considered.\n'+'All genes/SNPs will be considered if no gene list is provided')
    parser.add_argument('-n', '--sample-name', metavar='sample_name', type=str,
                             help='Sample name. If no sample name is provided, the first sample in VCF file will be used')
    return parser.parse_args()


# 4. collect all SNNPs, GWAS, LD scores
# 5. for each gene from the gene list
#   1. collect all SNPs within the gene (including promoter region)
#   2. check if SNP is in GWAS (only disease traits)
#   3. collect other SNPs in LD with the GWAS SNP
# 6. output:
# . 1. gene-wise GWAS disease traits
# . 2. gene-wise LD scores (for LD plot)

def main():
    args = parser_argv()

    gene_list_fn, gtf_fn, vcf_fn, gwas_fn, ld_duckdb_fn, sample_name = \
        args.gene_list, args.gtf, args.vcf, args.gwas, args.ld_duckdb, args.sample_name
    out_pre = args.out_pre

    as_gwas_basic_out = out_pre + '_allele_spliced_gene_gwas.tsv'
    as_gwas_basic_efo_out = out_pre + '_allele_spliced_gene_gwas_efo_term.tsv'
    # as_gwas_detail_out = out_pre + '_allele_spliced_gene_gwas_detailed.tsv'
    as_gwas_ld_out_dir = out_pre + '_allele_spliced_gene_gwas_ld' # directory, each gene has a separate file for all SNPs

    all_gene_id_to_name, all_gene_id_to_coor = get_gene_infor(gene_list_fn, gtf_fn)
    # print(all_gene_id_to_coor)
    chr_to_all_gwas_snp_trait = parse_gwas_file(gwas_fn) #, all_gene_names) # chr -> snp -> trait_info
    all_gene_to_snps = collect_gene_snps(all_gene_id_to_name, all_gene_id_to_coor, vcf_fn, ld_duckdb_fn, sample_name) # chr -> gene_name -> snp_id -> snp_info
    
    # gene_name -> [traits] : T/F
    gene_to_gwas_traits, gene_to_gwas_parent_traits = collect_gene_to_gwas_traits(all_gene_to_snps, chr_to_all_gwas_snp_trait, ld_duckdb_fn)

    # write to file
    output_gwas(all_gene_id_to_name, gene_to_gwas_traits, gene_to_gwas_parent_traits, as_gwas_basic_out, as_gwas_basic_efo_out)
    # output_ld(all_gene_id_to_name, all_gene_to_snps, chr_to_all_gwas_snp_trait, as_gwas_detail_out, as_gwas_ld_out_dir, ld_duckdb_fn)
    output_ld(all_gene_id_to_name, all_gene_to_snps, chr_to_all_gwas_snp_trait, as_gwas_ld_out_dir, ld_duckdb_fn)


if __name__ == '__main__':
    main()