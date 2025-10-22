# Annotation with genomic regions
The scripts that were used to overlap emVars, control SNPs and all SuRE SNPs with promoters, enhancers, DHS, and gencode annotations; and to generate plots showing the overlaps.

## 1. Overlap with regions
**annotations.job:** bash script where you can specify your input file and output file, and run the following script:
- **add_all_genomic_features_28022024.R:** adds to your input file: 'general' promoters and enhancers from UCSC and Hoffman et al. (_not used in the manuscript_), cell type 'specific' promoters and enhancers from chromHMM predictions, gencode locations (exon, intron, intergenic), closest TSS, and cell type 'specific' DHS.


## 2. Generate plots
It is important that you have already calculated the p-values for each comparison using the Fisher exact test (in the 'statistics' directory). For the Fisher exact test, create a file that contains the total number of emVars, control SNPs and all SuRE SNPs, and the number of SNPs within each annotation category for each. The script will then compare emVars to control SNPs, and emVars to all SuRE SNPs, and outputs the FE, OR, and p-values.
### 2.1 Overlap percentage plots
**location_plots_final_14022024.R:** creates a bargraph with significance indications for the _percentage_ overlap with genomic features. Due to the way the plots are generated, you will have to modify the output of the Fisher exact script. You need two files:
- One file containing all the percentages
![image](https://github.com/user-attachments/assets/adefa383-bc34-41c4-8a89-c4837ffa23ec)
- One file with the statistics in the format below. The y.position can take some tweaking. Make sure the group names and Annotations match between the files.
![image](https://github.com/user-attachments/assets/009c18e4-3dea-4755-8a4e-f3c7fdc9be08)

### 2.2 Fold enrichment plots
**location_fe-or_plots_08042024.R:** creates a bargraph with significant indications for the _fold enrichment_ (Figures 1E, S1C, and S1D). For this script the Fisher exact output file was also slightly modified. The minimal columns required: Fisher p-value, Annotation category, Fold enrichment, Odds ratio, Comparison (e.g. hNSC emVars vs all SuRE SNPs). This information was taken from the following columns in the Fisher exact output file: Annotation, FE_emvarvscontrol, Fisher_OR_emvarvscontrol, Fisher_p_emVarvscontrol, FE_emvarvsall, Fisher_OR_emvarvsall, Fisher_p_emvarvsall.
