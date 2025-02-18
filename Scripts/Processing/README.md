# Processing
All files related to the processing of the SuRE data, calling of the raQTLs, and generation of fragment expression plots

## 1. Files to process the raw SuRE data
# 1.1 Process the autosomes and chrX separately
**process.sh**: bash script that runs the R script to process the raw SuRE data for the *autosomes*\
**process_chrx.sh**: bash script that runs the R script to process the raw SuRE data for *chromosome X*\
Both scripts run the script below:
- **process_raw_wohnpcbr3_snp-permutation.R**: script to do the processing. It generates for each SNP a single row with information on number of fragments, mean expression for each allele, and p-values. Includes options on recreating the totals files (necessary for normalization) and downsampling if wanted. Also reprocesses the HepG2 and K562 data

# 1.2 Combine all processed autosome and chrX files into one
**combine_ss.sh**: bash script that runs:
- **combine_ss.R**: script to combine all the processed chromosome files (generated with the scripts above) into one big processed file.

## 2. Files to determine raQTLs
**raqtl_identification_31032022.R**: determines raQTLs using the output file generated with **combine_ss.sh** and **combine_ss.R** \
It can easily be modified to identify raQTLs in HepG2 and K562 by changing the cell type specific column names; and can also be used to determine control SNPs by selecting the SNPs _not_ meeting the p-value cut off.

## 3. Files to generate fragment expression plots
**sure_expression_fragments.R:** requires files that still contain information on the fragment level. These files (reads_sure_CHR.txt.gz) are generated while running the processing scripts. You need to provide the chromosome and the position of the SNP you are interested in. It then creates plots with the fragment expression with the positions of the fragments, and the fragment expression split by allele, as seen in Figure 1B and Figure 1C

**sure_expression_plots.R**: includes the plots generated with **sure_expression_fragments.R** and additionally creates SuRE expression plots, such as seen in Figure 1A, showing the SuRE expression split on the plus and minus strand and split on individuals.
