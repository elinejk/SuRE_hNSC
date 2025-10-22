library(data.table)
library(dplyr)
library(tools)

# your filtered  plink bim files
bimnames <- paste0("/gpfs/home6/ekoornstra/ldsc/freeze7/1000G_EUR_Phase3_plink/1000G.EUR.QC.", 1:22, ".bim")

# baseline annotation files to which you will add the binary sure annotation categories
filenames <- paste0("/gpfs/home6/ekoornstra/ldsc/freeze7/annotation_ld/1000G_Phase3_baselineLD_v2.2_nohm3filter/baselineLD.", 1:22, ".annot.gz")

# raqtl data
hnpc <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")
cont <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024.txt")

# add sure annotation data to your new annotation files, and save as new annot files
# for each new annotation file...
for(i in seq_along(filenames)){
  dat <- read.table(filenames[i], header=T)
  bim <- read.table(bimnames[i], header=F)
  names(bim) <- c("chr", "rsid","cm","bp","ref","alt")

  # add a column to show that it is a raqtl
  new_dat <- dat %>%
    mutate(hnsc_raqtl = if_else(SNP %in% hnpc$SNP_ID, "1", "0"),
	   hnsc_control = if_else(SNP %in% cont$SNP_ID, "1", "0"))

  # check if the order of the snps is still the same as the bim file
  if (all(new_dat$SNP == bim$rsid)) {
    print("rsids match with bim file, save file")
  } else {
    print("rsids do not match, check script")
  }
 
  # save files
  write.table(new_dat, paste0("/gpfs/home6/ekoornstra/ldsc/freeze7/annotation_ld/1000G_Phase3_baselineLD_v2.2_nohm3filter_hnsc_wcont/baselineLD.v2.2.sure.hnsc.",i,".annot"), quote=F, row.names=F) 
}

