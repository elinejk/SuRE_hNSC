library(data.table)
library(dplyr)
library(ggplot2)

# read in the sure file
df <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/Processing/freeze7/pvalue.freeze7.snp-permutation.04042024.txt.gz")

# remove unnecessary columns
df <- select(df, all_of(c("SNP_ID", "hNPC.wilcoxon.pvalue", "hNPC.wilcoxon.pvalue.random",
                          "HepG2.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue.random",
                          "K562.wilcoxon.pvalue", "K562.wilcoxon.pvalue.random")))

# remove SNPs without the pvalues
df <- df[!is.na(df$hNPC.wilcoxon.pvalue), ]
df <- df[!is.na(df$hNPC.wilcoxon.pvalue.random),]

# create a random subset with 100k SNPs
a <- sample_n(df, 100000)

# create the datatable for the plot
b <- data.table(expected = sort(-log10(a$hNPC.wilcoxon.pvalue.random)), observed = sort(-log10(a$hNPC.wilcoxon.pvalue)))

# create the qq plot
qplot <- ggplot(b, aes(x=expected, y=observed), color = "black") +
  geom_point() +
  geom_abline(slope=1, intercept=0, color="grey") + 
  labs(x="-log10(random wilcoxon p-values)", y="-log10(wilcoxon p-values)") +
  theme_bw()

