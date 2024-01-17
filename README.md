# NanoHunter: single-cell long-read RNA-seq analysis package

<!-- [![Latest Release](https://img.shields.io/github/release/xinglab/NanoHunter.svg?label=Release)](https://github.com/xinglab/NanoHunter/releases/latest) -->
<!-- [![Github All Releases](https://img.shields.io/github/downloads/xinglab/NanoHunter/total.svg?label=Download)](https://github.com/xinglab/NanoHunter/releases) -->
<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/nanohunter.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/nanohunter) -->
<!-- [![PyPI](https://img.shields.io/pypi/dm/nanohunter.svg?label=pip%20install)](https://pypi.python.org/pypi/nanohunter) -->
<!-- [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963) -->
<!-- [![GitHub Issues](https://img.shields.io/github/issues/xinglab/NanoHunter.svg?label=Issues)](https://github.com/xinglab/NanoHunter/issues) -->
<!-- [![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/xinglab/NanoHunter/blob/main/LICENSE) -->
## Updates (v1.0.0)

- First version

## What is NanoHunter
NanoHunter is an analysis pipeline designed for long-read single-cell RNA-seq data.
It mainly consists of two parts:
1. Barcode/UMI calling and gene/transcript quantification
2. Cell type-specific & allele-specific splicing analysis

<img src="figs/nanohunter-github-workflow.png" width="70%">

NanoHunter is very fast, the barcode/UMI calling step usually takes ~1 hour for 10 million long reads.

## Table of Contents
- [NanoHunter: single-cell long-read RNA-seq analysis package](#nanohunter-single-cell-long-read-rna-seq-analysis-package)
  - [Updates (v1.0.0)](#updates-v100)
  - [What is NanoHunter](#what-is-nanohunter)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Operating system and python version](#operating-system-and-python-version)
    - [Via conda (install locally for now)](#via-conda-install-locally-for-now)
    - [Via pip (not work yet)](#via-pip-not-work-yet)
    - [From source files](#from-source-files)
  - [(Optional) Consensus calling for RCA single-cell long-read data with TideHunter(≥1.5.4)](#optional-consensus-calling-for-rca-single-cell-long-read-data-with-tidehunter154)
    - [Input](#input)
    - [Command](#command)
  - [0. Pre-process: long-read mapping and transcript identification/quantification](#0-pre-process-long-read-mapping-and-transcript-identificationquantification)
    - [0.1 Mapping](#01-mapping)
    - [0.2 Transcript identification and quantification](#02-transcript-identification-and-quantification)
  - [1. Barcode \& UMI calling](#1-barcode--umi-calling)
    - [1.1 Input](#11-input)
    - [1.2 Command](#12-command)
    - [1.3 Output](#13-output)
  - [2. Cell clustering and annotating](#2-cell-clustering-and-annotating)
  - [3. Cell type-specific splicing analysis](#3-cell-type-specific-splicing-analysis)
    - [3.1 Input](#31-input)
    - [3.2 Command](#32-command)
    - [3.3 Output](#33-output)
  - [4. Allele-specific splicing analysis](#4-allele-specific-splicing-analysis)
    - [4.1 Phase long reads with `whatshap`](#41-phase-long-reads-with-whatshap)
    - [4.2 Identifiy allele-specific splicing](#42-identifiy-allele-specific-splicing)
    - [4.3 Identify disease-associated GWAS SNPs from allele-specific spliced genes](#43-identify-disease-associated-gwas-snps-from-allele-specific-spliced-genes)
  - [5. Visualization](#5-visualization)
    - [1. UMAP plot of gene VIM and EPCAM](#1-umap-plot-of-gene-vim-and-epcam)
    - [2. UMAP plot of CD44 transcripts](#2-umap-plot-of-cd44-transcripts)
    - [3. UMAP plot of both CD44 gene and transcripts](#3-umap-plot-of-both-cd44-gene-and-transcripts)

## Installation
### Operating system and python version
NanoHunter was written in python3 and tested on Linux/Unix systems, so it may not work well with python2 and/or other systems(Windows/Mac).

<!-- 
### <a name="dep"></a>~~Dependencies~~ installed via conda
* [minimap2](https://github.com/lh3/minimap2) >= 2.24
* [ESPRESSO](https://github.com/Xinglab/espresso) >= 1.3.1
* [python3](https://www.python.org/) >= 3.8
* [pysam](https://pysam.readthedocs.io/en/latest/) >= 0.21.0
* [edlib](https://pypi.org/project/edlib/) >= 1.3.9
* [kneed](https://pypi.org/project/kneed/) >= 0.8.3
* [scipy](https://scipy.org/) >= 1.10.1
* [rpy2](https://github.com/rpy2/rpy2) >= 3.5.11
* [PyVCF3](https://github.com/dridk/PyVCF3) >= 1.0.3
* [pyinterval](https://pyinterval.readthedocs.io/en/latest/) >= 1.2.0
* [duckdb](https://duckdb.org/docs/api/python/overview.html) >= 0.7.1
* [tabix](https://anaconda.org/bioconda/tabix)
* [whatshap](https://whatshap.readthedocs.io/en/latest/index.html) >= 1.7

(Optional) For rolling circle amplification (RCA) single-cell long-read data:
* [TideHunter](https://github.com/yangao07/TideHunter) >= 1.5.4 
-->

### Via conda (install locally for now)
```
git clone git@github.com:Xinglab/NanoHunter.git
cd NanoHunter
conda create --prefix ./conda_env
conda activate ./conda_env
conda install -c conda-forge -c bioconda python=3.8 --file conda_requirements.txt
python setup.py install
```

<!-- ```
conda create -n nanohunter python=3.8 nanohunter
conda activate nanohunter
``` -->

### Via pip (not work yet)
```
pip install nanohunter
```

### From source files
```
git clone git@github.com:Xinglab/NanoHunter.git
cd NanoHunter && pip install .
```

## (Optional) Consensus calling for RCA single-cell long-read data with [TideHunter](https://github.com/yangao07/TideHunter)(≥1.5.4)
### Input
* `rca_long_read.fq`: RCA long-read fasta/fastq file
* `five_prime.fa`, `three_prime.fa`: 5' and 3' sequence of RCA library. If splint sequence was used, please split splint sequence into two halves and use them as 5' and 3' sequence. Note that both 5' and 3' sequence need to have orientation of 5'->3', i.e., 3' sequence need to be reverse-complement

### Command
```
TideHunter rca_long_read.fq \
           -5 five_prime.fa \
           -3 three_prime.fa \
           -lF -t 16 \
           -o consensus.fa
```


## 0. Pre-process: long-read mapping and transcript identification/quantification
NanoHunter relies on the mapping and transcript identification/quantification result of long reads.
### 0.1 Mapping
For mapping, we recommend using [minimap2](https://github.com/lh3/minimap2) in RNA splice mode.
Any other long-read RNA-seq alignment tools can also be used here.

* Input
  * `long_read.fq/fa`: 1D long reads or RCA consnseus sequence in fastq/fasta format
  * `ref.fa`: reference genome
  * (optional) `anno.gtf`: gene annotation file in GTF format
  * (optional) `n_threads`: number of threads to use

  
* Command
```
# 1. convert GTF to BED12 using paftools.js from minimap2
paftools.js gff2bed anno.gtf > anno_junc.bed12

# 2. mapping with minimap2
minimap2 ref.fa long_read.fq/fa         \
         --junc-bed anno_junc.bed12     \
         -ax splice -ub -k14 -w4        \
         --sam-hit-only --secondary=no  \
         -t n_threads -o long_read.sam

# 3. convert sam to sorted bam
samtools view long_read.sam -b long_read.bam
samtools sort long_read.bam -@ n_threads -o long_read.sorted.bam
```
### 0.2 Transcript identification and quantification
For transcript identification and quantification, we recommend using [ESPRESSO(≥1.3.1)](https://github.com/Xinglab/espresso), other tools like [IsoQuant](https://github.com/ablab/IsoQuant) can also be used.
* Input
  * `long_read.sorted.bam`: sorted long-read alignment file in BAM format
  * `ref.fa`: reference genome file in FASTA format
  * `anno.gtf`: gene annotation file in GTF format 
  * `esp_output_dir`: output directory
* Output
  * `esp_output_dir/esp_cmpt.tsv`: compatible isoforms for all reads
  * `esp_output_dir/for_esp_input_N2_R0_updated.gtf`: updated GTF annotation
* Command
```
# 1. create tab-separated input file for ESPRESSO
mkdir esp_output_dir 2> /dev/null
in_tsv=esp_output_dir/for_esp_input.tsv
abs_path=$(realpath long_read.sorted.bam)
base_name=$(basename long_read.sorted.bam)
echo -e "$abs_path\t$base_name" > $in_tsv

# 2. ESPRESSO S step
perl ESPRESSO_S.pl -L esp_output_dir/for_esp_input.tsv  \
                   -F ref.fa -A anno.gtf            \
                   -M notFilterOutchrM -T n_threads \
                   -O esp_output_dir

# 3. ESPRESSO C step
perl ESPRESSO_C.pl -I esp_output_dir -F ref.fa \
                   -X 0 -T n_threads

# 4. ESPRESSO Q step
perl ESPRESSO_Q.pl -L esp_output_dir/for_esp_input.updated \
                   -A anno.gtf -T n_threads            \
                   -V esp_output_dir/esp_cmpt.tsv
```
* Note that `-V esp_cmpt.tsv` is optional in `ESPRESSO Q` step, but it is required if you want NanoHunter to output gene/transcript quantification inforamtion.

## 1. Barcode & UMI calling
NanoHunter identifies barcode and UMI from sorted alignment BAM file of single-cell long reads ***without*** using reference barcode from short-read data. 

For RCA long reads, before barcode/UMI calling, consensus sequences need to be generated (see [above](#optional-consensus-calling-for-rca-single-cell-long-read-data-with-tidehunter≥154)) and then mapped to reference genome.
### 1.1 Input 
* Required:
  * `long_read.sorted.bam`: sorted long-read BAM (recommend using [minimap2](https://github.com/lh3/minimap2))
* Optional:
  * `read_isoform_compatible.tsv`: tabular file of compatible isoforms for all reads, generated by [ESPRESSO(≥1.3.1)](https://github.com/Xinglab/espresso) or IsoQuant
  * `updated.gtf`: updated GTF annotation, generated by [ESPRESSO(≥1.3.1)](https://github.com/Xinglab/espresso) or IsoQuant
  * `annotation.gtf`: reference annotation GTF file, to retrieve gene names if gene names are not provided in `udpated.gtf`
  * `cell_barcode.tsv`: reference cell barcode list. If provided, nanohunter will directly use it to guide the barcode calling, only long reads with cell barcodes in the provided list will be kept
  * barcode sequence length (default: 16)
  * max. allowed edit distance between barcode and reference barcode (default: 2)
  * UMI sequence length (default: 12)
  * max. allowed edit distance between PCR-duplicated UMIs (default: 1)

### 1.2 Command
```
nanohunter long_read.sorted.bam \
           output_dir           \
           -p updated.gtf          \
           -m read_isoform_compatible.tsv \
           -g annotation.gtf  \
```
### 1.3 Output
* `bc_umi.bam`: BAM file with barcode/UMI information for each long read, BAM tags: `CB` and `UB`, for barcode and UMI.  Note that only long read with barcode/UMI called are kept
* `bc_umi.tsv`: tabular file with barcode/UMI/gene/transcript information for all barcode-called long reads
* `expression_matrix/`: folder containing 10X format sparse expression matrix at both gene and transcrpt level, can be directly parsed using standard single-cell analysis tools, like [Seurat](https://satijalab.org/seurat/)/[Azimuth](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html)
* `ref_bc.tsv`: reference cell barcode list identified from long reads. Note that if cell barcode list was provided to NanoHunter, `ref_bc.tsv` file will not be generated
* `high_qual_bc_umi.rank.tsv`: cell barcode list ranked by unique UMI count based on all high-quality long reads. Note that if cell barcode list was provided to NanoHunter, `high_qual_bc_umi.rank.tsv` file will not be generated
  
Example of `bc_umi.tsv`:

| cell barcode | UMI | # reads | compatible transcript id | gene id | gene name | read names (separeted by \`,\`)|
|-|-|-|-|-|-|-|
| ATCACGACACTTTAGG | ATCACATCCATG | 3 | ENST00000407249, ENST00000341832 | ENSG00000248333 | CDK11B | 741aa2c2-5840-4a29-bd90-3bdcb71604ba, 05f8432a-7f7e-446d-a6d5-8ab9e4eb5102, 6bb13ee1-7397-435c-8840-aeb2a90cf4ab |


## 2. Cell clustering and annotating
All the downstream single-cell long-read analysis rely on the cell type clustering result, which could be acomplished by running [Seurat](https://satijalab.org/seurat/) on the the gene expression matrix and annotating the cell clusters manually or based on known reference annotation like [Azimuth](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html).

For example, for human peripheral blood mononuclear cells (PBMC), the clustering and annotating result can be obtained by mapping to Azimuth human PBMC reference dataset. 
<!-- (https://github.com/Xinglab/NanoHunter/tree/main/scripts_for_paper#run_azimuth.R) -->
We provide the [R script](scripts_for_paper/README.md#run_azimuth_pbmcr) for human PBMC data. For data from other species/tissues, this needs to be done manually by users.

## 3. Cell type-specific splicing analysis
NanoHunter performs ***pairwise comparison*** to identify differntially spliced genes/transcripts between each two cell types/clusters.
### 3.1 Input
* Required
  * `expression_matrix/transcript`: transcript count matrix directory, generated by NanoHunter
  * `bc_to_cell_type.tsv`: list of barcode and corresponding cell type/cluster, generated by [Seurat](https://satijalab.org/seurat/)/[Azimuth](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html)
  * `out_prefix`: prefix of output files
* Optional
  * min. number of reads to keep a transcript for each cell type (default: 10)
  * min. number of transcripts to keep a gene for each cell type (default: 2)
  * list of genes of interest. If provided, all genes will still be processed, but only differentially spliced genes that are in this list will be output
### 3.2 Command
```
nh_cell_type_specific_splicing expression_matrix/transcript \
                               bc_to_cell_type.tsv            \
                               out_prefix
```

### 3.3 Output
* `*_cell_spliced_genes.tsv`: differentially spliced genes between 2 cell types. Gene FDR ≤0.05, max. delta ratio ≥0.05
* `*_cell_spliced_transcripts.tsv`: differentially spliced transcripts between 2 cell types. Gene FDR ≤0.05, transcript p value ≤0.05, detal ratio ≥0.05

Example of `*_cell_spliced_genes.tsv`:

| cell_type1| cell_type2|gene_id| gene_name| gene_fdr| max_delta_ratio| transcript_ids| cell1_counts| cell2_counts| cell1_ratios| cell2_ratios|
|-|-|-|-|-|-|-|-|-|-|-|
|Epithelial| Mesenchymal|ENSG00000026508| CD44| 7.96e-65| 0.68| ENST00000263398, ENST00000433892, ENST00000527326, ESPRESSO:chr11:443:2, ESPRESSO:chr11:443:3| 61.43,122.97, 20.0,23.41,13.63| 322.74,4.24,15.0,1.33-05,1.33e-05| 0.25,0.50,0.08,0.09,0.05| 0.94,0.01,0.04,3.89e-08,3.89e-08 |

Example of `*_cell_spliced_transcripts.tsv`:

| cell_type1| cell_type2 | gene_id| gene_name| transcript_id| gene_fdr| transcript_p| delta_ratio| count1| count2| ratio1| ratio2|
|-|-|-|-|-|-|-|-|-|-|-|-|
| Epithelial| Mesenchymal| ENSG00000026508| CD44| ENST00000263398| 7.96e-65| 8.25e-73| 0.68| 61.43| 322.74| 0.25| 0.94|


## 4. Allele-specific splicing analysis
With additional whole-genome sequencing data, NanoHunter can identify allele-specific splicing within each cell type.
The allele-specific splicing analysis mainly consists of three steps:
1. Phase long reads with whatshap
2. Identify allele-specific spliced genes/transcripts within each cell type
3. Search for disease-associated GWAS SNPs from identified allele-specific spliced genes

### 4.1 Phase long reads with `whatshap`
* Input
  * `wgs_phased.vcf`: WGS phased VCF file, generated from WGS short-read data using [GATK](https://gatk.broadinstitute.org/hc/en-us) and [SHAPEIT](https://odelaneau.github.io/shapeit5/)
  * `ref.fa`: reference genome fasta file
  * `bc_umi.bam`: BAM file with barcode/UMI information, generated by NanoHunter
* Command
```
# sort and index `bc_umi.bam`
samtools sort bc_umi.bam -o bc_umi_sorted.bam
samtools index bc_umi_sorted.bam

# identify haplotype
whatshap haplotag wgs_phased.vcf bc_umi_sorted.bam \
                  --reference ref.fa               \
                  --ignore-read-groups             \
                  --skip-missing-contigs           \
                  -o bc_umi_hap.bam                \
                  --output-haplotag-list hap_list.tsv
```
* Output
  * `bc_umi_hap.bam`: BAM file with barcode/UMI and additional haplotype information
  * `hap_list.tsv`: tabular file containing haplotype information for barcode-called long reads
  
### 4.2 Identifiy allele-specific splicing
* Input
  * `bc_umi.tsv`: tabular file with barcode/UMI/gene/transcript information, generated by NanoHunter
  * `hap_list.tsv`: tabular file containing haplotype information, generated by [whatshap](https://whatshap.readthedocs.io/en/latest/index.html), see [above](#41-phase-long-reads-with-whatshap)
  * `bc_to_cell_type.tsv`: list of barcode and corresponding cell type/cluster, generated by [Seurat](https://satijalab.org/seurat/)/[Azimuth](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html) 

* Command
```
nh_allele_specific_splicing bc_umi.tsv \
                            bc_to_cell_type.tsv \
                            hap_list.tsv \
                            output_prefix

```
* Output
  * `*_allele_spliced_genes.tsv`: allele-specific spliced genes. Gene FDR ≤0.05, max. delta ratio ≥0.05
  * `*_allele_spliced_transcripts.tsv`: allele-specific spliced transcripts. Gene FDR ≤0.05, transcript p value ≤0.05, delta ratio ≥0.05
  
Example of `*_allele_spliced_genes.tsv`:
|cell_type| gene_id| gene_name| gene_fdr| transcripts| hap1_counts| hap2_counts| hap1_ratios| hap2_ratios|
|-|-|-|-|-|-|-|-|-|
|Mono| ENSG00000105383| CD33| 4.42e-07| ENST00000262262, ENST00000421133| 58,9| 24,41| 0.86,0.13| 0.36,0.63|

Example of `*_allele_spliced_transcripts.tsv`:
|cell_type| transcript_id| transcript_p| gene_id| gene_name| gene_fdr| hap1_count| hap2_count| hap1_ratio| hap2_ratio|
|-|-|-|-|-|-|-|-|-|-|
|Mono| ENST00000262262| 3.58e-09| ENSG00000105383| CD33| 4.42e-07| 58| 24| 0.86| 0.36|

### 4.3 Identify disease-associated GWAS SNPs from allele-specific spliced genes
* Input
  * `*_allele_spliced_genes.tsv`: allele-specific spliced genes, generated by [nh_allele_specific_splicing](#42-identifiy-allele-specific-splicing)
  * `wgs_phased.vcf`: WGS phased VCF file, generated from WGS short-read data using [GATK](https://gatk.broadinstitute.org/hc/en-us) and [SHAPEIT](https://odelaneau.github.io/shapeit5/)
  * `annotation.gtf`: annotation GTF file
  * `gwas_catalog_efo.tsv`: GWAS catalog database with EFO term
  * `ld.duckdb`: linkage disequilibrium (LD) score database in duckdb format

* Command
```
nh_gene_with_gwas_disease -g allele_spliced_genes.tsv \
                          annotation.gtf       \
                          wgs_phased.vcf       \
                          gwas_catalog_efo.tsv \
                          ld.duckdb            \
                          output_prefix
```

* Output
  <!-- * `*_allele_spliced_gene_gwas.tsv`: GWAS disease-associated genes -->
  * `*_allele_spliced_gene_gwas_efo_term.tsv`: GWAS diesease-associated genes, with traits described in [EFO](https://www.ebi.ac.uk/efo/) term
  <!-- * `*_allele_spliced_gene_gwas_detailed.tsv`: GWAS and LD information of all SNPs -->
  * `*_allele_spliced_gene_gwas_ld`: directory of LD scores and GWAS traits for all SNPs, each gene has one separate `.ld` file
  
<!-- Format of `*_allele_splice_gwas_detailed.tsv`:
|snp_id| gene_name| gene_id| haplotype| gwas_trait| gwas_efo_term| ld_trait| ld_efo_term| ld_snps| ld_scores|
|-|-|-|-|-|-|-|-|-|-|
|rs3865444| CD33| ENSG00000105383| H2| Alzheimer disease| Neurological disorder| NA, Alzheimer disease, late-onset Alzheimers disease, NA, NA| NA, Neurological disorder,Neurological disorder, NA, NA| rs12459419, rs7245846, rs1354106, rs34813869, rs33978622| 1.0,0.89,0.88,0.88,0.96| -->

Example of `*_allele_splice_gwas_efo_term.tsv`:
|gene_name| gene_id| Cancer| Immune system disorder| Neurological disorder| Cardiovascular disease| Digestive system disorder| Metabolic disorder| Response to drug| Other disease|
|-|-|-|-|-|-|-|-|-|-|
|CD33| ENSG00000105383| False| False | True | False| False| False| False| False|

Example of `*_allele_spliced_gene_gwas_ld/gene.ld`:
|snp_id| rs3865444| rs2459141| rs12459419| rs2455069| rs7245846| rs33978622| rs34813869| rs1354106| haplotype| chrom| pos| type| gwas_trait| gwas_trait_efo_term| gwas_pvalue|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|rs3865444| 1.0| 0.326844| 1.0| 0.327703| 0.891224| 0.967195| 0.887773| 0.881092| H2| chr19| 51224706| gwas| Alzheimer disease| Neurological disorder| 2e-09

## 5. Visualization
After obtaining a list of genes/transcripts of interest (cell-type-specific splicing/allele-specific spciling), NanoHunter also offers visualization of your results in the single-cell UMAP space, simlar to the plots from Seurat R package.

```
# gene/transcript expression folder generated by NanoHunter
gene_mtx <- "output_dir/expression_matrix/gene"
transcript_mtx <- "output_dir/expression_matrix/transcript"
# script for UMAP plot
source("nanohunter_visualization.R)
```

### 1. UMAP plot of gene VIM and EPCAM
```
nanohunter_umap_plot(gene_mtx_dir = gene_mtx,
                     feature_mtx_dir = gene_mtx,
                     feature_list = c("VIM", "EPCAM"))
```

<img src="figs/PGmix_VIM_EPCAM_gene_umap.png" width="40%">


### 2. UMAP plot of CD44 transcripts
```
nanohunter_umap_plot(gene_mtx_dir = gene_mtx,
                     # `feature_mtx_dir` needs to be only one folder or same size as `feature_list`
                     feature_mtx_dir = transcript_mtx,
                     feature_list = c("ENST00000433892", "ENST00000263398"))
```

<img src="figs/PGmix_CD44_transcripts_umap.png" width="40%">

### 3. UMAP plot of both CD44 gene and transcripts
```
nanohunter_umap_plot(gene_mtx_dir = gene_mtx,
                     # `feature_mtx_dir` needs to be only one folder or same size as `feature_list`
                     feature_mtx_dir = c(gene_mtx, transcript_mtx, transcript_mtx),
                     feature_list = c("CD44", "ENST00000433892", "ENST00000263398"))
```

<img src="figs/PGmix_CD44_gene_transcripts_umap.png" width="60%">

Note that, for `gene_mtx_dir` and `feature_mtx_dir`, you can also provide with Seurat object instead of maxtrix folder.

```
gene_obj = make_seurat_obj(gene_mtx)
transcript_obj = make_seurat_obj(transcript_mtx)

nanohunter_umap_plot(gene_mtx_dir = gene_obj,
                     # `feature_mtx_dir` needs to be only one folder or same size as `feature_list`
                     feature_mtx_dir = c(gene_obj, transcript_obj, transcript_obj),
                     feature_list = c("CD44", "ENST00000433892", "ENST00000263398"))
```