# TFBS analyses
The scripts that were used for all TFBS analyses.

## 1. TFBS overlap
### 1.1 Overlap with TFBS and comparisons between raQTLs, control SNPs and all SuRE SNPs
**concordance_SNP2TFBSW_20210902.R:** this script matches SuRE rsIDs to the rsIDs listed in SNP2TFBs to obtain a list of SNPs overlapping with TFBS. For each SuRE SNP it indicates if it overlaps with a TFBS or not, and if it overlaps with a TFBS if the effect direction is the same (e.g. REF higher SuRE expression and REF stronger TF binding = concordant). \

**tfbs_plots_final_140224.R:** this script uses summarized TFBS data obtained from the output above to create two types of plots. 
- Two plots per cell type, one with TFBS overlap and one with concordance, Figures 1F, 1G, S1E, S1F

### 1.2 Comparison between cell types
Figures S2A, S2B

## 2. TFBS enrichment
Figures 2A, 2B, S2C, S2D

## 3. ADDITIONAL ANALYSES BY EVA
Figures 2C, 2D, S2E 
Table S3
