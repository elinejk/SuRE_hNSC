# TFBS analyses
The scripts that were used for all TFBS analyses.

## 1. TFBS overlap
The following two scripts only include information on the first TFBS in the SNP2TFBS file.

**concordance_SNP2TFBS_20210902.R:** this script matches SuRE rsIDs to the rsIDs listed in SNP2TFBS to obtain a list of SNPs overlapping with TFBS. For each SuRE SNP it indicates if it overlaps with a TFBS or not, and if it overlaps with a TFBS if the effect direction is the same (e.g. REF higher SuRE expression and REF stronger TF binding = concordant). 

**tfbs_plots_final_140224.R:** creates bargraphs for TFBS overlap and TFBS concordance with significance indications. Due to the way the plots are generated, you will have to modify the output of the Fisher exact script which is used to get the p-values. You therefore need two files: one file containing all the percentages, and one with the statistics (for the format see the README in the Genomic_Annotations section). \
The script creates two types of bargraphs (for each there will be one plot with TFBS overlap and one with TFBS concordance):
- Per _cell type_ comparing raQTLs to control SNPs and all SuRE SNPs (Figures 1F, 1G);
- Plots showing the results for hNSC, HepG2 and K562 in a single plot (Figures S1E, S1F).


## 2. TFBS enrichment 
The following scripts include information on all TFBS overlapping raQTLs and control SNPs.

### 2.1 Comparison between raQTLs and controls

*TFBS enrichment*  
**040325_SNP2TFBS_enrichment_conc.R** is used to calculate the enrichment of TFBS amongst raQTLs and control SNPs. It uses the SNP2TFBS file to determine which SNPs bind a TFBS, then creates a df indicating what TFBS is bound by which SNPs (note: a SNP can bind multiple TFBS), and then uses this df to determine, for each unique TFBS, the relative proportion of raQTLs over control SNPs. A proportion test is used to calculate a p-value for statistically significant different proportions, indicating the preference of a TFBS to bind either raQTLs or control SNPs. 

*TFBS enrichment plot*  
**100325_SuRE_figures.R - TFBS enrichment plot section** visualizes the enrichment of each TFBS and labels those that are highly enriched in raQTLs (corresponds to figure 2A, S2C and S2D in manuscript). It plots the proportions between raQTLs and control SNPs per TFBS against the FDR-corrected p-value calculated using a proportion test for each TFBS (see 040325_SNP2TFBS_enrichment_conc.R script for specifics on the enrichment calculations).

*Concordance between SuRE signal and TF binding strength*    
**040325_SNP2TFBS_enrichment_conc.R** is used to calculate the concordance between the SuRE signal and TF binding, as indicated by the score difference between ALT and REF alleles for each SNP. ). If the signs of the SuRE signal and the score difference are the same, the TFBS is denoted as "concordant". For each TFBS, the percentage of concordant raQTL SNPs is then calculated.

*Plot of concordance between SuRE signal and TF binding strength*  
**100325_SuRE_figures.R - TFBS concordance plot section** shows the concordance between SuRE signal and the TFBS binding strength, as indicated by the score difference between REF allele binding and ALT allele binding (corresponds to figure 2B in manuscript; concordance calculated in the 040325_SNP2TFBS_enrichment_conc.R script). 

### 2.2 Comparison between cell types

*Cell-type specific TFBS enrichment*  
**260225_SuRE_barplots_enrichment.R:** is used to calcultate the relative TFBS enrichment of each cell type over the others (hNSC over HepG2 / hNSC over K562). It takes the enrichment results of each cell type as input, then calculates the relative enrichment using the proportion of raQTLs within each TFBS in each cell type (Figures S2A, S2B). 

## 3. Regulatory properties of TFBS binding raQTLs 
The following scripts include information on all TFBS overlapping raQTLs and control SNPs.

*Expression of SNPs within TFBS*   
**100325_SuRE_figures.R - Highest expressing alleles plot section** was used to determine which TFBS confers high expression of the SNPs it regulates. For each TFBS, the mean expression of the highest-expressing alleles for all SNPs it binds is shown (corresponds to figure 2C in the manuscript). 

*Correlation between nucleotide conservation and raQTL abundance within TFBS*   
**100325_SuRE_figures.R - Correlation plot section** shows the mean correlation between raQTL/control fraction and nucleotide conservation for each position within each TFBS (corresponds to figure 2D and S2E in the manuscript).

*Relative enrichment of neighboring TFBS amongst raQTLs vs control SNPs*  
**LOLA_JASPAR_twoSets.sh:** this script is ran using the JASPAR enrichment tool from LOLA, see [here](https://bitbucket.org/CBGR/jaspar_enrichment/src/master/) for more details. 
The script takes as input:
- A file containing raQTLs with 100bp flanking regions.
- A file containing control SNPs with 100bp flanking regions.
We selected raQTLs and control SNPs binding one of the five most enriched TFBS: ZBTB33, YY1, GABPA, TP63, and TP53. The script calculates the enrichment of neighboring TFBS in regions flanking raQTLs versus control SNPs for each of these TFBS separately. For each TFBS, the output file "allEnrichments.tsv" contains enrichment results, including odds ratios (OR) for neighboring TFBS.

**main_jaspar.sh:** this script allows for running the above script in parallel for each of the five TFBS. 

**JASPAR_enrichment_neighboring.R:** This script calculates the mean OR of each TF class, providing an overall measure of whether certain TF classes are more frequently associated with raQTLs than control SNPs. A mean OR > 1 suggests a TF class is more enriched in raQTL regions, while a mean OR < 1 suggests depletion


