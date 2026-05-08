library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)

## FILES ##
emvar_path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/emVars/freeze7/overlap/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024_ANNOTATED.txt"

## READ IN FILES ##
emvar <- fread(emvar_path)
emvar$ID <- paste0("chr", emvar$chrom, ":", emvar$pos.hg19, ":", emvar$ref.seq, ":", emvar$alt.seq)

emvar <- emvar[, .(SNP_ID, ID, chrom, pos.hg19, ref.seq, alt.seq, wilcoxon_pvalue = hNPC.wilcoxon.pvalue, cDNA_ref_mean = hNPC.cDNA.ref.mean, cDNA_alt_mean = hNPC.cDNA.alt.mean)]


### ADDITIONAL ANALYSES ###
# read in the files from the other NPC MPRA datasets
mcafee <- read_excel("V:/ddata/CELB/poot/Eline Koornstra/MPRA_project/Variants/Original_tables/McAfee2023_Supplementary-data_S1.xlsx", skip = 1) %>% as.data.table()
lee <- read_excel("C:/Users/084902/Downloads/1-s2.0-S0092867424014351-mmc2.xlsx", sheet = "Table S1B")  %>% as.data.table()
rummel_s <- read_excel("V:/ddata/CELB/poot/Eline Koornstra/MPRA_project/Variants/Original_tables/Rummel2023_supplementary_TableS3.xlsx", sheet = "A. MPRA_Results_allAnnotated")  %>% as.data.table()
rummel_v <- read_excel("V:/ddata/CELB/poot/Eline Koornstra/MPRA_project/Variants/Original_tables/Rummel2023_supplementary_TableS3.xlsx", sheet = "B. DevLib_Results_processed") %>% as.data.table()
guo_t <- read_excel("V:/ddata/CELB/poot/Eline Koornstra/MPRA_project/Variants/Original_tables/Guo2023_SupplementaryTables_1-11.xlsx", sheet = "table_S2A_DeepSea_scores") %>% as.data.table()
guo_s <- read_excel("V:/ddata/CELB/poot/Eline Koornstra/MPRA_project/Variants/Original_tables/Guo2023_Supplementary-data_S3.xlsx", sheet = "Data_S5_mpra_results_annotated") %>% as.data.table()
guo_rs <- fread("V:/ddata/CELB/poot/Eline Koornstra/MPRA_project/Variants/for-masterfile/guo2023_for-masterfile.txt")

# fix the rummel data.table
rummel_v <- rummel_v[, .(eid, NPC_logFC, NPC_fdr)]
rummel <- merge(rummel_s, rummel_v, by = "eid", all.x = TRUE)
rummel <- rummel[!is.na(rummel$NPC_logFC), ] # remove those without NPC logFC
rummel <- rummel[, .(rsId, NPC_logFC, NPC_fdr)]
rummel$paper <- "rummel"

# fix the guo data.tables
guo_t <- guo_t[, .(chrom, pos, id, ref, alt)]
guo_t$alleles <- paste0(guo_t$ref, "/", guo_t$alt)
guo_s <- guo_s[guo_s$mpra_tissue %like% "NPC", ]
guo_s <- guo_s[, .(Linked_SNP, Chr, Position, `Ref/Alt`, mpra_logfc_mean, mpra_pval_mean)]
guo <- merge(guo_t, guo_s,  by.x = c("chrom", "pos", "alleles"), by.y = c("Chr", "Position", "Ref/Alt"), all.x = TRUE)
guo <- merge(guo, guo_rs, by.x = c("chrom", "pos", "ref", "alt"), by.y = c("chrom", "pos_hg19", "ref", "alt"), all.x = TRUE)
guo <- guo[, .(rsID, mpra_logfc_mean, mpra_pval_mean)]
guo$paper <- "guo"

# prepare the other data.tables for comparison
mcafee <- mcafee[, .(rsID, MPRA_logFC, MPRA_FDR)]
mcafee$paper <- "mcafee"
lee <-lee[, .(RSID, MPRA_logFC, MPRA_FDR)]
lee$paper <- "lee"

# Combine files
mpras <- rbind(mcafee, lee, guo, rummel, use.names=FALSE)

# merge with SuRE variants
sure <- merge(emvar, mpras, by.x = "SNP_ID", by.y = "rsID")

# add concordance
sure$SURE_logFC <- log2(sure$cDNA_alt_mean + 0.001) - log2(sure$cDNA_ref_mean + 0.001)
sure$MPRA_logFC <- as.numeric(sure$MPRA_logFC)
sure$CONCORDANCE  <- ifelse((sure$SURE_logFC > 0 & sure$MPRA_logFC > 0) | (sure$SURE_logFC < 0 & sure$MPRA_logFC < 0), "concordant", "discordant")
sure$SIGNIF <- ifelse(sure$paper == "mcafee" & sure$MPRA_FDR < 0.1 & !is.na(sure$MPRA_FDR), "significant", 
                ifelse(!sure$paper == "mcafee" & sure$MPRA_FDR < 0.05 & !is.na(sure$MPRA_FDR), "significant", "not_significant"))

# get stats
uniqueN(sure$SNP_ID)

# significant variants stats
signif <- sure[sure$SIGNIF == "significant", ]
uniqueN(signif$SNP_ID)
table(signif$paper)

sig_sub <- aggregate(CONCORDANCE ~ SNP_ID, data = signif, FUN = function(x) paste(x, collapse = ", "))
table(sig_sub$CONCORDANCE)


# plots
ggplot(sure, aes(x = MPRA_logFC, y = SURE_logFC, color = SIGNIF, shape = paper)) + 
    geom_point() +
    labs(title = "SuRE hNSC emVars versus NPC MPRA datasets", x = "MPRA logFC", y = "hNSC SuRE logFC") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw()