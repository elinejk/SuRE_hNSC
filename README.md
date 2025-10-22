# SuRE in human neural stem cells (hNSCs)
Github repository containing the scripts related to the manuscript "Koornstra, E., van Zanten, E.S., Klaassen, N.H.M. et al.: Large-scale identification of the regulatory potential of SNPs in human neural stem cells"

## Scripts
This directory contains the scripts used to generate the data and create several of the plots from the manuscript. \
Each subdirectory has it's own README explaining the scripts and in which order to run them

### 1. Processing
Contains the scripts used to 
- process the raw SuRE data;
- call expression modulating variants (emVars, variants that significantly alter SuRE transcriptional activity); 
- generate plots showing SuRE expression (Figures 1A, 1B, 1C, 5D)

### 2. Basic Analyses
Contains the scripts used to
- determine the replicate correlation;
- create qq-plots for expected and observed p-values (Figure S1A);
- minor allele frequency annotation;
- calculate the storey's pi1 between emVar sets;
- compare emVars and control SNPs based on their expression, fragment count and MAFs (Figures S1C-F)

### 3. Genomic Annotations
Contains the scripts used to 
- annotate the SuRE data with enhancers, promoters, DHS, and Gencode genetic information;
- generate plots showing the fold enrichment (Figures 1E, S2B, S2C) or percentage of overlap (Figure S2A)

### 4. TFBS
Contains the scripts used to 
- annotate the SuRE data with TFBS from SNP2TFBS;
- generate plots showing the overlap and concordance with the first TFBS (Figures 1F, 1G, S2D, S2E);
- determine the enrichment of TFBS for emVars vs controls (Figures 2A, 2B, S3A, S3B) and between cell types (Figures 2E, 2F, S3D, S3E);
- investigate the expression of SNPs within TFBS (Figure 2C);
- look into nucleotide conservation (Figure 2D, S3C);
- determine neighboring TFs of enriched TFBS (Table S3)

### 5. eQTLs
Contains the scripts used to 
- overlap emVars and control SNPs with p > 0.5 (null variants) with GTEX V7 eQTL data for brain tissues, fibroblasts, whole blood, and liver;
- generate plots showing the overlap and concordance with eVariants (Figure 3A);
- generate the lolliplot plots (Figures 3D, 3E, S4A, S4D);
- generate a 'heatmap' showing the eGenes for a single eVariant with concordance with emVars (Figure 5E); 
- generate Figure S4B which shows the transcripts for an eGene with eVariant locations, chromHMM data and DNAse I peaks.

### 6. GWAS
Contains all the scripts used for the various analyses involving overlap with GWAS data. 

**5.1 LDSC:** contains all the scripts used to 
- run LDSC;
- generate the plots showing the Enrichment (Figures 4A, 4B). 

**5.2 emVars in loci:** contains all the scripts used to
- determine which GWAS SNPs, including candidate SNPs, were tested by SuRE and were emVars;
- determine the percentages and counts used to generate Figures 4C, 4D, 4E, S4C; 
- generate GWAS locus plots (Figures 5A, 5B, S5A, S5B);
- create a plot to show tissue expression of a specific gene (Figure S5C).

**5.3 PAINTOR fine-mapping:** contains the scripts used to perform fine-mapping using PAINTOR.

**5.4 Distance to the lead SNP:** contains all the scripts used to
- determine the distance to the lead SNP for emVars, control SNPs and all SuRE SNPs that reached genome-wide significance in the GWAS and are located within 100kb of a lead SNP 
- generate a plot showing the distribution of the distance to the lead SNP (Figure 4F).

### 6. Statistics
 Contains the script used to perform the Fisher exact test. It was used for both the genomic annotation comparisons, the TFBS comparisons, eGene tissue distribution comparisons, and distance to the lead SNP comparisons.
