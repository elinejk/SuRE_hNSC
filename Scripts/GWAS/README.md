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
- **annotfile_prep_binary_hnsc.R:** adds information on a SNP being an hNSC emVar or hNSC control SNP to the annotation files _downloaded from LDSC_ (downloaded file: 1000G_Phase3_baselineLD_v2.2_ldscores.tgz). It uses the bim EUR plink files to ensure the SNPs in the annotation file are in the same order as the bim files. 

**03_ldscore_snpfilter.job:** runs LDSC partitioned LD score calculation for all the annotation files generated with 02. Only prints SuRE SNPs (saved as sure.SNPs.txt), so we are not including any SNPs in the LDSC analysis that were not tested by SuRE. The plink bed files used in this script were _downloaded from LDSC_ (downloaded files: 1000G_EUR_Phase3_plink).

**04-binary_partitioned_heritability:** runs LDSC partitioned heritability using the LD and Annotation files created with 02 and 03; the weight file from 01; and the summary statistics from 00. The frequency files and plink bed files used in this script were _downloaded from LDSC_ (downloaded files: 1000G_Phase3_frq and 1000G_EUR_Phase3_plink).

### 1.2 Generating the plots
**ldsc_phenotype_barplot.R:** creates a bargraph with significance showing the heritability for a subset of annotations for one phenotype. Takes directly the ouput from script 04.

**ldsc_emvar_control_summary_plot.R:** creates a bargraph with significance for the hNSC emVar and hNSC control category for all phenotypes. You need a file containing the hNSC emVar LDSC results for all phenotypes.

## 2. SuRE SNPs and emVars in GWAS loci
### 2.1 Determine the number of SuRE SNPs and emVars within a locus
These scripts require GWAS summary statistics and FUMA output files (files containing the loci and candidate SNPs).

**run_01_sure_in_gwas_loci.job:** run the following script:
- **01_sure_in_gwas_loci.R:** determines for each locus in the FUMA file: (1) which SuRE SNPs are within the boundaries, if they were part of the sumstats, if they are a candidate SNP, and if they are emVars; (2) how many sumstats SNPs are within the boundaries; (3) creates a summary file for each locus with the number of sumstat SNPs, number of candidate SNPs, and how many of those were tested by SuRE and were emVars. The information from 'part 3' was used as input to generate Figures 4C, 4D, 4E and S4C in PRISM.

**02_emVar_in_locus_listfile.R:** takes the output from part 1 of script 01, and saves only those that are emVars. The information in this output file was later used to determine which locus to highlight.

**emvar_among_multisignif_locus.R:** to determine the number of emVars in loci with more than one significant variant.


### 2.2 Highlighted locus
**surexgwas_LOCUS-plots.R:** for a selected locus of a phenotype it will plot the GWAS SNPs (colored on SuRE-tested), emVars within the locus, and genes within the locus. It adds a highlight for the region you specify (here this is the region you want to plot in the next script). 

**surexgwas_ZOOM-plots.R:** plots the same things as the script above, and also adds DNAse I peaks, and the locations of promoters and enhancers.


## 3. PAINTOR fine-mapping
### 3.1 Fine-mapping and credible sets
**01_create-locus-files.job:** was used to run the **locus-file_paintor.R** script, which splits the summary statistics into separate files for each locus

**02_create-ld-files.job:** runs the PAINTOR CalcLD_1KG_VCF_flip.py script to create LD files for each locus file

**03_create-annotation-files.job:** runs the PAINTOR AnnotateLocus.py script to create annotations files for each locus file. The annotation_paths.txt is used to get the paths to the annotations files that you want to use. 

**04_run_paintor.job:** runs PAINTOR with the annotations you list.

**calculate_credible_set.R:** uses the output from script 04 to define the credible sets

### 3.2 Expression plot for example gene
**HPA_RNA_expression_overview_plot.R:** this script was used to recreate the RNA tissue expression plots from the human protein atlas.

## 4. Distance to the lead SNP
**distance_to_lead_fileprep.R:** for each row (phenotype) in the infofile, it determines which SuRE SNPs are in the locus, if they are part of the GWAS summary statistics, determines their distance to the lead SNP, and filters on genome-wide significance. You need to run this script three times: for emVars, for control SNPs, and for all SuRE SNPs

**distance_to_lead_plots.R:** creates a histogram to show the distance to the lead SNP for emVars, control SNPs and all SuRE SNPs, for all phenotypes combined. It also calcualtes p-values


