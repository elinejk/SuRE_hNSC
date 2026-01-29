# SuRE in human neural stem cells (hNSCs)
Github repository containing the scripts related to the manuscript "Koornstra, E., van Zanten, E.S., Klaassen, N.H.M. et al.: Large-scale identification of the regulatory potential of SNPs in human neural stem cells"

## Scripts
This directory contains the scripts used to generate the data and create several of the plots from the manuscript. \
Each subdirectory has it's own README explaining the scripts and in which order to run them

### 1. Processing
Contains the scripts used to 
- process the raw SuRE data;
- call expression modulating variants (emVars, variants that significantly alter SuRE transcriptional activity); 
- generate plots showing SuRE expression

### 2. Basic Analyses
Contains the scripts used to
- determine the replicate correlation;
- create qq-plots for expected and observed p-values;
- minor allele frequency annotation;
- calculate the storey's pi1 between emVar sets;
- compare emVars and control SNPs based on their expression, fragment count and MAFs.

### 3. DHS
Contains the scripts used to 
- process raw WGS data and filtering for heterozygous SNPs (downloading, sorting, annotation);
- assess allelic imbalance of heterozygous SNPs .

### 3. Genomic Annotations
Contains the scripts used to 
- annotate the SuRE data with enhancers, promoters, DHS, and Gencode genetic information;
- annotate emVars with putative targets;
- generate plots showing the fold enrichment or percentage of overlap.

### 4. DHS sensitivity
Contains the scripts to 
- compare SuRE allelic imbalance to DHS sensitivity;
- generate plots showing the results.

### 5. TFBS
Contains the scripts used to 
- annotate the SuRE data with TFBS from SNP2TFBS;
- generate plots showing the overlap and concordance with the first TFBS;
- determine the enrichment of TFBS for emVars vs controls and between cell types;
- investigate the expression of SNPs within TFBS;
- look into nucleotide conservation;
- determine neighboring TFs of enriched TFBS.

### 6. eQTLs
Contains the scripts used to 
- overlap emVars and control SNPs with p > 0.5 (null variants) with GTEX V7 eQTL data for brain tissues, fibroblasts, whole blood, and liver;
- generate plots showing the overlap and concordance with eVariants;
- generate the lolliplot plots;
- generate a 'heatmap' showing the eGenes for a single eVariant with concordance with emVars; 
- generate Figure S4B which shows the transcripts for an eGene with eVariant locations, chromHMM data and DNAse I peaks.

### 7. GWAS
Contains all the scripts used for the various analyses involving overlap with GWAS data. 

**7.1 LDSC:** contains all the scripts used to 
- run LDSC;
- generate the plots showing the Enrichment. 

**7.2 emVars in loci:** contains all the scripts used to
- determine which GWAS SNPs, including candidate SNPs, were tested by SuRE and were emVars;
- determine the corresponding percentages and counts; 
- generate GWAS locus plots.

**7.3 PAINTOR fine-mapping:** contains the scripts used to 
- perform fine-mapping using PAINTOR;
- create a plot to show tissue expression of a specific gene.

**7.4 Distance to the lead SNP:** contains all the scripts used to
- determine the distance to the lead SNP for emVars, control SNPs and all SuRE SNPs that reached genome-wide significance in the GWAS and are located within 100kb of a lead SNP;
- generate a plot showing the distribution of the distance to the lead SNP.

### 8. Statistics
 Contains the script used to perform the Fisher exact test. It was used for both the genomic annotation comparisons, the TFBS comparisons, eGene tissue distribution comparisons, and distance to the lead SNP comparisons.
