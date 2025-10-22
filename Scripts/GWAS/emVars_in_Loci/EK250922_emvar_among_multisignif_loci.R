library(data.table)
library(tidyr)

## SETTINGS AND PATHS ##
sumpath <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/GWAS/freeze7/sure-in-gwas/"

phenos <- c("ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ")

combined <- data.table()

## START LOOP ##
for (i in phenos){
  print(paste0("running for phenotype: ", i))
  df <- fread(paste0(sumpath, i, "_hnsc_sure_in_gwas_snps_11042024.txt"))
  
  # filter on candidate and significance
  a <- df[df$SuRE_GWAScandidate == "1" & df$GWAS_sign == "1", ]
  
  # get summary
  res <- a %>%
    group_by(locus) %>%
    summarise(n_total = n(),
              n_emvar = sum(hNSCraQTL_GWAScandidate)) %>%
    as.data.table()
  
  # filter to contain multiple significant variants
  b <- res[res$n_total > 1, ]
  
  # get counts
  print(paste0("Number of loci with more than one sign variant in LD: ", nrow(b)))
  print(paste0("Number of loci where one of them is an emVar: ", nrow(b[b$n_emvar > 0, ])))
  print(paste0("Percentage of loci: ", nrow(b[b$n_emvar > 0, ])/nrow(b) * 100))

  # combine with other phenotypes
  a$phenotype <- i
  combined <- rbind(combined, a)
  
}

## Distance to the lead snp
# load distances
dis <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/GWAS/freeze7/distance_to_lead/20250923_allphenos_raqtls_gwassign_distance-snp-to-lead.txt")
dis <- select(dis, all_of(c("SNP_ID", "phenotype", "distance_lead")))

# merge files
combined <- merge(combined, dis, by = c("SNP_ID", "phenotype"))

# get stats on the distance
combined %>%
  group_by(phenotype) %>%
  summarise(mean_distance = mean(distance_lead)/1000,
            median_distance = median(distance_lead)/1000,
            min_distance = min(distance_lead)/1000,
            max_distance = max(distance_lead)/1000) %>%
  as.data.table()
