# Comparison to GWAS
The scripts that were used for all analyses involving Genome Wide Association Studies.

## 1. LDSC
In the section below, _'downloaded from LDSC'_ means they were downloaded from https://askesgroup.broadinstitute.org/LDSCORE
### 1.1 Performing LDSC
**00_munge_sumstats.job:** makes sure the summary statistics are in the correct format. The summary statistics used are described in the 'data files' section of the manuscript.

**01_ldscore_univariate.job:** creates the weight files needed for script 04. The plink bed files used in this script were _downloaded from LDSC_ (downloaded files: 1000G_EUR_Phase3_plink).\
The original weights files downloaded from LDSC only contained HapMap3 SNPs and no HLA regions. We remade the weights files so they contain the SNPs present in the EUR plink files downloaded from LDSC: ~10 million SNPs. These ~10 million SNPs were saved in a file called regression.snplist.

**02_annotationfiles.job:** create the annotation files needed fro script 04. \
This script runs the following script:
- **annotfile_prep_binary_hnsc.R:** adds information on a SNP being an hNSC raQTL to the annotation files _downloaded from LDSC_ (downloaded file: 1000G_Phase3_baselineLD_v2.2_ldscores.tgz). It uses the bim EUR plink files to ensure the SNPs in the annotation file are in the same order as the bim files. 

**03_ldscore_snpfilter.job:** runs LDSC partitioned LD score calculation for all the annotation files generated with 02. Only prints SuRE SNPs (saved as sure.SNPs.txt), so we are not including any SNPs in the LDSC analysis that were not tested by SuRE. The plink bed files used in this script were _downloaded from LDSC_ (downloaded files: 1000G_EUR_Phase3_plink).

**04-binary1_partitioned_heritability:** runs LDSC partitioned heritability using the LD and Annotation files created with 02 and 03; the weight file from 01; and the summary statistics from 00. The frequency files and plink bed files used in this script were _downloaded from LDSC_ (downloaded files: 1000G_Phase3_frq and 1000G_EUR_Phase3_plink).

### 1.2 Generating the plots
**ldsc_barplot_hnsc-raqtls.R:** creates a bargraph with significance showing the heritability for a subset of annotations for one phenotype (Figure 4B). Takes directly the ouput from script 04.

**ldsc_hnsc_summary_plots.R:** creates a bargraph with significance for the hNSC raQTL category for all phenotypes (Figure 4A). You need a file containing the hNSC raQTL LDSC results for all phenotypes.

## 2. SuRE SNPs and raQTLs in GWAS loci
### 2.1 Determine the number of SuRE SNPs and raQTLs within a locus
This is the information that was used to create Figures 4C, 4D, 4E and S3C in PRISM


### 2.2 Highlighted locus
These are the scripts use to create Figures 5A and 5B.

## 3. Distance to the lead SNP
Figure 4F
