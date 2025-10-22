library(data.table)
library(dplyr)

## PATHS ##
pheno_path <- "/gpfs/home6/ekoornstra/surexgwas/freeze7/"
phenos <- c("ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ")
raqtls <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")

raqtl_set <- data.table()

## GET LIST OF EMVARS IN LOCI ##
for (pheno in phenos) {
  df <- fread(paste0(pheno_path, pheno, "_hnsc_sure_in_gwas_snps_11042024.txt"))
  df <- df[df$hNSC_emVar == "1",]
  df$phenotype <- pheno
  
  df_coord <- merge(df, raqtls, by = "SNP_ID")
  
  cols <- c("SNP_ID", "chrom", "pos.hg19", "locus", "chr", "start", "end","GWAS_SNP", "GWAS_P","hNSCemVar_GWAScandidate", "phenotype")
  a <- select(df_coord, all_of(cols))
  raqtl_set <- rbind(raqtl_set, a)  

}

fwrite(raqtl_set, paste0(pheno_path, "raqtls-in-loci_list_all-phenos_12042024.txt"), sep='\t', quote=F, row.names=F)

