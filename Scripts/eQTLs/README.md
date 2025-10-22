# eQTL analyses
Scripts used for all analyses involving eQTLs.

## 1. Overlap with eQTLs
**eqtl_overlap_brain.R** and **eqtl_overlap_non-brain.R:** for each emVar it determines if it was tested by GTEx (.all.pairs.txt.gz files) and if it was an eQTL in a tissue. Next it determines the concordance and prints summary information. These scripts generated the data that was used, among other things, to create Figure 3B and 3C in PRISM.

**eqtl_functional_annotations.R:** takes the results from the emVar output file generated with the script above and annotates it with genomic regions (for example DHS) and gencode information.

## 2. Generation of plots
### 2.1 Overview of all eQTL overlaps
**eqtl_overlap_summary-plot.R:** takes a file summarizing the overlap between emVars and eVariants (determined in step 1) and uses this to create a heatmap (Figure 3A).

### 2.2 Lollipop plot for a single gene
**eqtl_lolliplot.R:** for a specific eGene and a specific tissue, it creates a lollipop plot. It plots a specific transcript and the position of the eVariants (colored based on SuRE concordance) relative to the transcript position. As seen in Figures 3D, 3E, S4A, S4D

### 2.3 Concordance heatmap for a single emVar
**eqtl-tissue-conc_heatmap.R:** takes a single SuRE SNP and determines for each tissue if the SNP is present as eVariant in the tissue and if the effect directions are the same. It then combines all the information into a single 'heatmap' (Figure 5E). The locations of exons need to be manually extracted from gencode and entered into the script.

### 2.4 Gene transcript plot with chromHMM annotations and eQTL locations
**eqtl_gene_chromhmm_plot.R:** creates a plot containing all protein-coding transcripts of a specific gene, and adds pannels containing chromHMM predictions for NPCs and fibroblasts, DNAse I peaks, and the locations of eQTLs (as seen in Figure S4b).
