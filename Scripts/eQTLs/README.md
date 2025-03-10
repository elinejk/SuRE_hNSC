# eQTL analyses
Scripts used for all analyses involving eQTLs.

## 1. Overlap with eQTLs
**eqtl_overlap_alternative-pipeline_31012024.R** and **eqtl_overlap_alternative-pipeline_non-brain_31012024.R:** for each raQTL it determines if it was tested by GTEx (.all.pairs.txt.gz files) and if it was an eQTL in a tissue. Next it determines the concordance and prints summary information. These scripts generated the data that was used, among other things, to create Figure 3B and 3C in PRISM.

**eqtl_functional_annotations_02072024.R:** takes the results from the raQTL output file generated with the script above and annotates it with genomic regions (for example DHS) and gencode information.

## 2. Generation of plots
### 2.1 Overview of all eQTL overlaps
**eqtl_overlap_summary-plot_31012024.R:** takes a file summarizing the overlap between raQTLs and eQTLs (determined in step 1) and uses this to create a heatmap (Figure 3A).

### 2.2 Lollipop plot for a single gene
**eqtl_lolliplot_19062024.R:** for a specific gene and a specific eQTL tissue, it creates a lollipop plot. Here it shows a specific transcript and the position of the eQTLs (colored based on SuRE concordance) relative to the transcript position.
As seen in Fiogures 3D, 3E, S3A, and S3D

### 2.3 Concordance heatmap for a single raQTL
**eqtl-tissue-conc_heatmap_17052024.R:** takes a single SNP and determines for each eQTL tissue if the SNP is present and if the effect directions are the same. It then combines all the information into a single 'heatmap' (Figure 5E). Unfortunately, the locations of exons need to be manually extracted from gencode and entered into the script.

### 2.4 HPA overlap and plots
**hnsc-raqtl_eQTL_comparisons_HPA_05062024.R:** uses the annotated output files generated in step 1 and a set of 1000 random eGenes, and determines in which tissues the protein-coding eGenes are expressed (using the Human Protein Atlas data). This was used to create stacked bargraphs (as seen in Figure S3B).

### 2.5 Gene transcript plot with chromHMM annotations and eQTL locations
**eqtl_gene_chromhmm_plot_27062024.R:** creates a plot containing all protein-coding transcripts of a specific gene, and adds pannels containing chromHMM predictions for NPCs and fibroblasts, DNAse I peaks, and the locations of eQTLs (as seen in Figure S4).
