library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# read in sure file with replicates
a <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/Processing/freeze7/fr7_for_replicates/pvalue.freeze7.snp-permutation.WITHREPLICATES.23052025.txt.gz")

## RUN PART BELOW ONLY IF YOU ARE STARTING FROM THE MAIN REPLICATE FILE ##
# remove unnecessary columns
df <- select(df, all_of(c("SNP_ID", 
                          "hNPC.cDNA.rep1.ref.mean", "hNPC.cDNA.rep1.alt.mean",
                          "hNPC.cDNA.rep2.ref.mean", "hNPC.cDNA.rep2.alt.mean")))

# make it long format
df_long <- gather(df, key=allele_rep1, value=mean_expression_rep1, hNPC.cDNA.rep1.ref.mean:hNPC.cDNA.rep1.alt.mean)
df_long <- gather(df_long, key=allele_rep2, value=mean_expression_rep2, hNPC.cDNA.rep2.ref.mean:hNPC.cDNA.rep2.alt.mean)

# remove columns where the alleles don't match
a <- df_long[(df_long$allele_rep1 == "hNPC.cDNA.rep1.ref.mean" & df_long$allele_rep2 == "hNPC.cDNA.rep2.ref.mean") |
               (df_long$allele_rep1 == "hNPC.cDNA.rep1.alt.mean" & df_long$allele_rep2 == "hNPC.cDNA.rep2.alt.mean"), ]

# add allele column
a$allele <- ifelse(a$allele_rep1 == "hNPC.cDNA.rep1.ref.mean" & a$allele_rep2 == "hNPC.cDNA.rep2.ref.mean", "reference",
                   ifelse(a$allele_rep1 == "hNPC.cDNA.rep1.alt.mean" & a$allele_rep2 == "hNPC.cDNA.rep2.alt.mean", "alternative", "error"))

## FROM HERE IT IS THE SAME FOR THE PRE-SUBSET FILE AND THE MAIN REPLICATE FILE ##
# determine the distribution
hist(a$mean_expression_rep1) # skewed --> log transform
hist(log(a$mean_expression_rep1)) # normal

hist(a$mean_expression_rep2) # skewed --> log transform
hist(log(a$mean_expression_rep2)) # normal

# calculate the correlation for the full set of SNPs
# perform pearson on log transformed values
a$log_mean_expression_rep1 <- log(a$mean_expression_rep1)
a$log_mean_expression_rep2 <- log(a$mean_expression_rep2)
df <- a[is.finite(a$log_mean_expression_rep1) & is.finite(a$log_mean_expression_rep2),]
cor.test(df$log_mean_expression_rep1, df$log_mean_expression_rep2, method = "pearson")

