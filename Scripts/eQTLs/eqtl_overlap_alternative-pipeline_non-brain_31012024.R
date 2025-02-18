library(data.table)
library(tidyverse)

# paths and settings             
tissues <- c("Cells_Transformed_fibroblasts","Liver", "Whole_Blood")

all_path <- "/projects/0/AdamsLab/Projects/sure/resources/eQTLs/"
sign_path <- "/projects/0/AdamsLab/Projects/sure/resources/eQTLs/"

raqtl_path <- "/gpfs/home6/ekoornstra/raqtls/freeze7/"

anatype <- "np-raqtls"

# prepare raQTL file
raqtl <- fread(paste0(raqtl_path, "annotations/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024_ANNOTATED.txt")) %>% as.data.table()
#raqtl <- fread(paste0(raqtl_path, "annotations/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024_ANNOTATED.txt")) %>% as.data.table()

#raqtl <- raqtl[raqtl$hNPC.wilcoxon.pvalue > 0.5, ]
raqtl <- raqtl[raqtl$NPC_narrowpeak_DHS == "1", ]

raqtl$chrom <- gsub("chr", "", raqtl$chrom)
raqtl[, variant_id := paste0(chrom, "_", pos.hg19, "_", ref.seq, "_", alt.seq, "_b37")]
raqtl <- raqtl[, c("SNP_ID", "variant_id", "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean")]

print(paste0("number of raqtls in file: ", nrow(raqtl)))

# start the loop to get for each tissue the total tested SNPs, significant pairs, and concordance
for (tissue in tissues){
  # get total tested SNPs
  #tested <- fread(paste0(all_path, tissue, ".allpairs.txt.gz")) %>% as.data.table()
  #ov <- merge(raqtl, tested, by = "variant_id")
  #a <- ov$SNP_ID %>% as.data.table()
  #a <- unique(a)
  #rm(tested)
  #rm(ov)
  
  # then determine the significant eqtls
  eqtl <- fread(paste0(sign_path, tissue, ".v7.signif_variant_gene_pairs.txt.gz")) %>% as.data.table()
  ov_sign <- merge(raqtl, eqtl, by="variant_id")
  b <- ov_sign$SNP_ID %>% as.data.table()
  b <- unique(b)
  rm(eqtl)

  # calculate the concordance
  ov_sign$hnsc_log2_alt_ref <- log2(ov_sign$hNPC.cDNA.alt.mean/ov_sign$hNPC.cDNA.ref.mean)
  ov_sign$raqtl_dir <- ifelse(ov_sign$hNPC.cDNA.alt.mean > ov_sign$hNPC.cDNA.ref.mean, "pos", "neg")
  
  ov_sign$slope_dir <- ifelse(ov_sign$slope > 0, "pos", "neg")
  ov_sign$concordant <- ifelse(ov_sign$raqtl_dir == ov_sign$slope_dir, "conc", ifelse(is.na(ov_sign$slope_dir), "no overlap", "disc"))
  
  sign <- ov_sign[ov_sign$concordant == "conc", ]
  d <- sign$SNP_ID %>% as.data.table()
  d <- unique(d)
  
  # get summary information
  print(paste0("analysis of raqtls: ", anatype))
  print(paste0("printing summary information for ", tissue))
  #print(paste0("hNSC had ", nrow(a), " controls tested as eQTL for ", tissue))
  print(paste0(nrow(b), " raqtls were significant eQTLs in ", tissue))
  print(paste0(nrow(d), " raqtls and significant eQTLs have the same effect direction in ", tissue))
  print("[[[[]]]]")
  
  # save files
  fwrite(sign, paste0(raqtl_path, "annotations/eqtls/raqtls/raqtl_eqtl-", tissue, "_", anatype, "_concordant-hits_18062024.txt"), sep = "\t", quote=F, row.names=F)
  fwrite(ov_sign, paste0(raqtl_path, "annotations/eqtls/raqtls/raqtl_eqtl-", tissue, "_", anatype, "_all-hits_18062024.txt"), sep = "\t", quote=F, row.names=F) 
}


