# Comparison to GWAS
The scripts that were used for all analyses involving Genome Wide Association Studies.

## 1. LDSC
### 1.1 Performing LDSC
**00_munge_sumstats.job:** makes sure the summary statistics are in the correct format

**01_ldscore_univariate.job:** creates the weight files
- Filters on a SNP list??

**02_annotationfiles.job:** create the annotation files. This script runs the following script:
- **annotfile_prep_binary_hnsc.R:** adds information on a SNP being an hNSC raQTL to the annotation files downloaded from the LDSC website

**03_ldscore_snpfilter.job:** runs LDSC partitioned LD score calculation. Only prints SuRE SNPs, so we are not including any SNPs in the LDSC analysis that were not tested by SuRE.

**04-binary1_partitioned_heritability:** actually runs LDSC using the LD and Annotation files created with 02 and 03; the weight file from 01; and the summary statistics from 00.

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
