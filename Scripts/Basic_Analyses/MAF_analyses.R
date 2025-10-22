library(data.table)
library(tidyverse)
library(ggpubr)

# load file
maf <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/MAF/310524_hNSC_allSNPs_MAF.txt")
names(maf) <- c("SNP_ID", "MAF")

# annotate file with dataset
mafr <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/MAF/310524_hNSC_raqtls_MAF.txt")
names(mafr) <- c("SNP_ID", "MAF")
mafc <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/MAF/310524_hNSC_controls_MAF.txt")
names(mafc) <- c("SNP_ID", "MAF")

maf$Dataset <- ifelse(maf$SNP_ID %in% mafr$SNP_ID, "emVar", 
                     ifelse(maf$SNP_ID %in% mafc$SNP_ID, "control", "inactive"))
maf$Dataset = factor(maf$Dataset, levels = c("inactive", "control", "emVar"))

# add bins
maf$bins <- cut(maf$MAF, breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
                labels=c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25","0.25-0.3", 
                         "0.3-0.35", "0.35-0.4", "0.4-0.45","0.45-0.5"))

# get some stats
sumstats <- data.table(
  bins = "all",
  count = nrow(maf),
  min = min(maf$MAF),
  max = max(maf$MAF),
  mean = mean(maf$MAF),
  median = median(maf$MAF)
)

ss <- maf %>% 
  group_by(bins) %>%
  summarise(count = n(),
            min = min(MAF),
            max = max(MAF),
            mean = mean(MAF),
            median = median(MAF)
  ) %>% as.data.table()

sumstats <- rbind(sumstats, ss)

rs <- maf %>% 
  group_by(Dataset) %>%
  summarise(count = n(),
            min = min(MAF),
            max = max(MAF),
            mean = mean(MAF),
            median = median(MAF)
  ) %>% as.data.table()

# make plots
density <- ggplot(maf, aes(x = MAF)) +
  geom_density(color = "#009E73", fill = "#009E73", alpha = 0.5) +
  theme_bw()

density_split <- ggplot(maf, aes(x = MAF, group = Dataset)) +
  geom_density(aes(color = Dataset, fill = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("grey", "#56B4E9", "#E69F00")) +
  scale_fill_manual(values = c("grey", "#56B4E9", "#E69F00")) +
  theme_bw()

histo <- ggplot(maf, aes(x = MAF)) +
  geom_histogram(color = "#009E73", fill = "#009E73", alpha = 0.5,
                 binwidth = 0.05, center = 0.025) +
  theme_bw()


## emVar analyses
em <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/emVars/freeze7/overlap/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024_ANNOTATED.txt")

em <- merge(em, mafr, by = "SNP_ID")

em %>% 
  filter(E009_enhancer == 1) %>%
  summarise(count = n(),
            min = min(MAF),
            max = max(MAF),
            mean = mean(MAF),
            median = median(MAF)
  ) %>% as.data.table()
