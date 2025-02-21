# Annotation with genomic regions
The scripts that were used to overlap raQTLs, control SNPs and all SuRE SNPs with promoters, enhancers, DHS, and gencode annotations; and to generate plots showing the overlaps.

## 1. Overlap with regions
**annotations.job:** bash script where you can specify your input file and output file, and run the following script:
- **add_all_genomic_features_28022024.R:** adds to your input file: 'general' promoters and enhancers from UCSC and Hoffman et al. (_not used in the manuscript_), cell type 'specific' promoters and enhancers from chromHMM predictions, gencode locations (exon, intron, intergenic), closest TSS, and cell type 'specific' DHS.


## 2. Generate plots
It is important that you have already calculated the p-values for each comparison using the Fisher exact test (in the 'statistics' directory).
### 2.1 Overlap percentage plots

### 2.2 Fold enrichment plots
