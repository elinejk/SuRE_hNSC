# TFBS analyses
The scripts that were used for all TFBS analyses.

## 1. TFBS overlap
### 1.1 Overlap with TFBS and comparisons between raQTLs, control SNPs and all SuRE SNPs
**concordance_SNP2TFBS_20210902.R:** this script matches SuRE rsIDs to the rsIDs listed in SNP2TFBs to obtain a list of SNPs overlapping with TFBS. For each SuRE SNP it indicates if it overlaps with a TFBS or not, and if it overlaps with a TFBS if the effect direction is the same (e.g. REF higher SuRE expression and REF stronger TF binding = concordant). 

**tfbs_plots_final_140224.R:** creates bargraphs for TFBS overlap and TFBS concordance with significance indications. Due to the way the plots are generated, you will have to modify the output of the Fisher exact script which is used to get the p-values. You therefore need two files: one file containing all the percentages, and one with the statistics (_for the format see the README in the Genomic_Annotations section)._ \
The script creates two types of bargraphs (for each there with be a TFBS overlap and a TFBS concordance plot):
- Per cell type comparing raQTLs to control SNPs and all SuRE SNPs (Figures 1F, 1G);
- Plots showing the results for hNSC, HepG2 and K562 (Figures S1E, S1F).

### 1.2 Comparison between cell types
Figures S2A, S2B

## 2. TFBS enrichment
Figures 2A, 2B, S2C, S2D

## 3. ADDITIONAL ANALYSES BY EVA
Figures 2C, 2D, S2E 
Table S3
