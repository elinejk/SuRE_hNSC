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


## 2. TFBS enrichment
Figures 2A, 2B, S2C, S2D

## 3. Regulatory properties of TFBS binding raQTLs 
Figures 2C, 2D


**260225_SuRE_barplots_enrichment.R:** is used to calcultate the relative TFBS enrichment of each cell type over the others (hNSC over HepG2 / hNSC over K562). It takes the enrichment results of each cell type as input, then calculates the relative enrichment using the proportion of raQTLs within each TFBS in each cell type. 

**LOLA_JASPAR_twoSets.sh:** this script is ran using the JASPAR enrichment tool from LOLA see [here](https://bitbucket.org/CBGR/jaspar_enrichment/src/master/) for more details. 
The script takes as input:
- A file containing raQTLs with 100bp flanking regions.
- A file containing control SNPs with 100bp flanking regions.
We selected raQTLs and control SNPs binding one of the five most enriched TFBS: ZBTB33, YY1, GABPA, TP63, and TP53. The script calculates the enrichment of neighboring TFBS in regions flanking raQTLs versus control SNPs for each of these TFBS separately. For each TFBS, the output file "allEnrichments.tsv" contains enrichment results, including odds ratios (OR) for neighboring TFBS.

**main_jaspar.sh:** this script allows for running the above script in parallel for each of the five TFBS. 

**JASPAR_enrichment_neighboring.R:** This script calculates the mean OR of each TF class, providing an overall measure of whether certain TF classes are more frequently associated with raQTLs than control SNPs. A mean OR > 1 suggests a TF class is more enriched in raQTL regions, while a mean OR < 1 suggests depletion
