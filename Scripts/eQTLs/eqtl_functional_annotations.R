library(data.table)
library(dplyr)
library(tidyverse)
library(rtracklayer)

### STEP 1: SET PATHS ###
qtl <- "/home/ekoornstra/raqtls/freeze7/annotations/eqtls/raqtls/all/all_Cells_Transformed_fibroblasts_raqtl-eQTL-candidates-per-gene_id_19062024.txt"
outfile <- "/home/ekoornstra/raqtls/freeze7/annotations/eqtls/raqtls/all/all_Cells_Transformed_fibroblasts_raqtl-eQTL-candidates-per-gene_id_19062024_ANNOTATED.txt"
raqtl_path <- "/home/ekoornstra/raqtls/freeze7/annotations/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024_ANNOTATED.txt"
gc_path <- "/home/ekoornstra/resources/gencode.v19.annotation.gtf.gz"

### STEP 2: SUBSET THE EQTL FILE TO ONLY RETAIN RELEVANT COLUMNS ###
df <- fread(qtl)

cols <- c("variant_id", "SNP_ID", "gene_id", "tss_distance", "maf", "pval_nominal", "slope", "slope_se", "concordant", "total_eQTLs", "total_eQTL-SuRE", "total_eQTL-raQTLs", "total_eQTL-raQTLs-conc")
new_cols <- c("variant_id", "SNP_ID", "egene_id", "egene_tss_distance", "maf", "eqtl_pval_nominal", "eqtl_slope", "eqtl_slope_se", "eqtl_concordance", 
				"total_egene_eqtls", "total_egene_eqtls_sure", "total_egene_eqtls_raqtls", "total_egene_eqtls_conc_raqtls")

df <- select(df, all_of(cols))
names(df) <- new_cols


### STEP 3: MERGE EQTL FILE WITH THE ANNOTATED EMVAR FILE ###
# in the step above you remove the sure expression values, but they will be returned in this step #
raqtl <- fread(raqtl_path)

a <- merge(df, raqtl, by = "SNP_ID", all.x=TRUE)


### STEP 4: ADD GENCODE INFORMATION ###
# the eqtl file does not contain gene information, so we will add another column without the .X, and then merge it with a gencode file to get gene names and gene type #
gencode <- readGFF(gc_path)
gencode <- gencode[gencode$type == "gene", ]

gencode <- select(gencode, all_of(c("gene_id", "gene_name", "gene_type", "strand")))
gencode <- gencode %>% separate(gene_id, c("ensl_gene_id", NA), remove=T)
gencode <- select(gencode, all_of(c("ensl_gene_id", "gene_name", "gene_type", "strand")))
names(gencode) <- c("egene_ensembl", "egene_name", "egene_type", "egene_strand") # as preparation for the merge

a <- a %>% separate(egene_id, c("egene_ensembl", NA), remove=F)

b <- merge(a, gencode, by = "egene_ensembl", all.x=T)


### STEP 5: SAVE FILE ###
fwrite(b, outfile, sep='\t', quote=F, row.names=F)


