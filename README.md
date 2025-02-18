# SuRE_hNSC
Github containing the files related to the manuscript "Identifying the regulatory potential of 7 million SNPs in human neural stem cells"

## Scripts
This directory contains the scripts used to generate the data and create several of the plots from the manuscript. \
Each subdirectory has it's own README explaining the scripts and in which order to run them

### 1. Processing
Contains the scripts used to 
- process the raw SuRE data;
- call raQTLs; and
- generate plots showing SuRE expression (Figures 1A, 1B, 1C, and 5D)

### 2. Genomic Annotations
Contains the scripts used to 
- annotate the SuRE data with enhancers, promoters, DHS, and Gencode genetic information;
- calculate significant differences between raQTLs, control SNPs and all SuRE SNPs (was also used for the TFBS section below); and
- generate plots showing the fold enrichment (Figure 1E) or percentage of overlap (Figures S1B, S1C, and S1D)

### 3. TFBS
Contains the scripts used to 
- annotate the SuRE data with TFBS from SNP2TFBS;
- generate plots showing the overlap and concordance (Figures 1F, 1G, S1E, and S1F). \

Also contains the scripts used for the TFBS enrichment analyses: ..... (Figures 2, and S2)

### 4. eQTLs
Contains the scripts used to 
- overlap raQTLs and control SNPs with p > 0.5 with GTEX V7 eQTL data for brain tissues, fibroblasts, whole blood, and liver;
- generate plots showing the overlap and concordance with eQTLs (Figure 3A);
- generate the lolliplot plots (Figures 3D, 3E, S3A, and S3D);
- generate _a 'heatmap' showing the eGenes for a single eQTL with concordance with raQTLs (Figure 5E);_ 
-  _generate Figure S4 which shows the transcripts for an eGene with eQTL locations, chromHMM data and DNAse I peaks. !!!!FIND THESE SCRIPTS!!!!!_ \
-  overlap raQTL-eGene pairs and 1000 fibroblast eGenes with RNA tissue distribution data from the Human Protein Atlas; and
-  generate the plots showing the RNA tissue distribution (Figure S3B).

### 5. GWAS
Contains all the scripts used for the various analyses involving overlap with GWAS data.

**5.1 LDSC:** contains all the scripts used to 
- run LDSC;
- generate the plots showing the Enrichment (Figures 4A, and 4B). \

**5.2 raQTLs in loci:** contains all the scripts used to
- determine which GWAS SNPs were tested by SuRE and were raQTLs;
- determine the percentages and counts used to generate Figures 4C, 4D, 4E, and S3C; and
- generate GWAS locus plots (Figures 5A, and 5B).

**5.3 Distance to the lead SNP:** contains all the scripts used to
- determine the distance to the lead SNP for raQTLs, control SNPs and all SuRE SNPs within 100kb of a lead SNP that reached genome-wide significance in the GWAS
- generate a plot showing the distribution of the distance to the lead SNP (Figure 4F).
