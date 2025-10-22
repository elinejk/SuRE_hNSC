library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# read in sure file with replicates
a <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/Processing/freeze7/fr7_for_replicates/pvalue.freeze7.snp-permutation.WITHREPLICATES.PLOTCOLSONLY.23052025.txt.gz")

## RUN PART BELOW ONLY IF YOU ARE STARTING FROM THE MAIN REPLICATE FILE ##
# remove unnecessary columns
#df <- select(df, all_of(c("SNP_ID", 
#                          "hNPC.cDNA.rep1.ref.mean", "hNPC.cDNA.rep1.alt.mean",
#                          "hNPC.cDNA.rep2.ref.mean", "hNPC.cDNA.rep2.alt.mean")))

# make it long format
#df_long <- gather(df, key=allele_rep1, value=mean_expression_rep1, hNPC.cDNA.rep1.ref.mean:hNPC.cDNA.rep1.alt.mean)
#df_long <- gather(df_long, key=allele_rep2, value=mean_expression_rep2, hNPC.cDNA.rep2.ref.mean:hNPC.cDNA.rep2.alt.mean)

# remove columns where the alleles don't match
#a <- df_long[(df_long$allele_rep1 == "hNPC.cDNA.rep1.ref.mean" & df_long$allele_rep2 == "hNPC.cDNA.rep2.ref.mean") |
#               (df_long$allele_rep1 == "hNPC.cDNA.rep1.alt.mean" & df_long$allele_rep2 == "hNPC.cDNA.rep2.alt.mean"), ]

# add allele column
#a$allele <- ifelse(a$allele_rep1 == "hNPC.cDNA.rep1.ref.mean" & a$allele_rep2 == "hNPC.cDNA.rep2.ref.mean", "reference",
#                   ifelse(a$allele_rep1 == "hNPC.cDNA.rep1.alt.mean" & a$allele_rep2 == "hNPC.cDNA.rep2.alt.mean", "alternative", "error"))


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


# calculate for 100k random SNPs over 100 iterations
set.seed(123)
runs <- 100
results <- data.table(run = integer(runs), correlation = numeric(runs))

for(i in 1:runs) {
  sampled <- sample_n(df, 100000)
  test <- cor.test(sampled$log_mean_expression_rep1, sampled$log_mean_expression_rep2, method = "pearson")
  results[i, `:=`(run = i, correlation = test$estimate)]
}

sumstats <- results[, .(
  min_correlation = min(correlation),
  max_correlation = max(correlation),
  mean_correlation = mean(correlation),
  median_correlation = median(correlation)
)]

print(sumstats)

# create a random subset with 100k SNPs
set.seed(123)
sset <- sample_n(a, 100000)

# make a correlation plot
# each dot represents the mean expression of all fragments covering either the reference or the alternative allele
cplot <- ggplot(sset, aes(x=mean_expression_rep1, y=mean_expression_rep2)) +
  geom_point(aes(color=allele), alpha=0.5) +
  geom_smooth(method = "lm", se=TRUE, color = "blue") +
  scale_color_manual(values = c("black", "grey")) +
  stat_cor(aes(label = after_stat(r.label))) +
  labs(x ="Mean SuRE expression replicate 1", y = "Mean SuRE expression replicate 2") + 
  theme_bw()

# save
ggsave("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/Processing/freeze7/replicate_correlations_100000SNPs_EK02062025.pdf",
       cplot, width = 5, height = 4)


# correlation plot separated on ref and alt
splot <- ggplot(sset, aes(x=mean_expression_rep1, y=mean_expression_rep2)) +
  geom_point(aes(color=allele), alpha=0.5) +
  geom_smooth(method = "lm", se=TRUE, color = "blue") +
  scale_color_manual(values = c("black", "grey")) +
  stat_cor(aes(label = after_stat(r.label))) +
  facet_grid(~allele) +
  labs(x ="Mean SuRE expression replicate 1", y = "Mean SuRE expression replicate 2") + 
  theme_bw()

# save
ggsave("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/Processing/freeze7/replicate_correlations_separated_100000SNPs_EK02062025.pdf",
       splot, width = 9, height = 4)

# also save the file with snp ids
fwrite(sset, "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/Processing/freeze7/replicate_correlations_100000SNPs_EK02062025.txt", quote=F, row.names=F, sep='\t')
