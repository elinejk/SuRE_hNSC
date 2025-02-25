# Comparison to GWAS
The scripts in these directories were used for all analyses involving Genome Wide Association Studies.

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
These scripts were used to generate the information that was used to create Figures 4C, 4D, 4E and S3C in PRISM. \
It requires summary statistics and FUMA output files (files containing the loci and candidate SNPs).

**run_01_sure_in_gwas_loci.job:** run the following script:
- **01new_sure_in_gwas_loci.R:** determines for each locus in the FUMA file: (1) which SuRE SNPs are within the boundaries, if they were part of the sumstats, if they are a candidate SNP, and if they are raQTLs; (2) how many sumstats SNPs are within the boundaries; (3) creates a summary file for each locus with the number of sumstat SNPs, number of candidate SNPs, and how many of those were tested by SuRE and were raQTLs.

**02_raqtl_in_locus_listfile.R:** takes the output from part 1 of script 01, and saves only those that are raQTLs. The information in this output file was later used to determine which locus to highlight.

### 2.2 Highlighted locus
**surexgwas_LOCUS-plots_17052024.R:** for a selected locus of a phenotype it will plot the GWAS SNPs (colored on SuRE-tested), raQTLs within the locus, and genes within the locus. It adds a highlight for the region you specify (here this is the region you want to plot in the next script). This script was used to create Figure 5A.

**surexgwas_ZOOM-plots_17052024.R:** plots the same things as the script above, and also adds DNAse I peaks, and the locations of promoters and enhancers. This script was used to create Figure 5B.

## 3. Distance to the lead SNP
**01_final_distanceto-lead.job:** runs the script below and requires a file with SuRE SNPs (this can be raQTLs, control SNPs or all SuRE SNPs), a file containing the paths to sumstats, FUMA lead SNP files, and FUMA locus files (= infofile). You need to run this script three times, once for each SuRE SNP set.
- **01_distance_to_lead_SNP.R:** for each row (phenotype) in the infofile, it determines which SuRE SNPs are in the locus, if they are part of the sumstats, determines their distance to the lead SNP, and filters on genome-wide significance. This last part was used later in script 03 to generate Figure 4F. _Not used in the manuscript:_ it also determines for each lead SNP the distance to genome-wide significant SuRE SNPs within 100kb, and only retains the SuRE SNP with the lowest p-value. 

**02_combine_distance_sets.R:** combines the output files from script 01 for the different phenotypes in a file with 'phenotype' allphenos. You need to run this script three times, once for each SuRE SNP set.

**03-alt_distance_to_lead_snp_plots.R:** used to create Figure 4F. It takes the output files from script 02 (raQTLs, control SNPs, and all SuRE SNPs), and creates a bargraph showing the distribution in kb for each set.



