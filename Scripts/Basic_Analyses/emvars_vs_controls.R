library(data.table)
library(tidyverse)
library(ggpubr)

# settings
raqtls <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/raQTLs/freeze7/overlap/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024_ANNOTATED.txt")
controls <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/raQTLs/freeze7/overlap/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024_ANNOTATED.txt")

mafr <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/MAF/310524_hNSC_raqtls_MAF.txt")
names(mafr) <- c("SNP_ID", "MAF")
mafc <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/MAF/310524_hNSC_controls_MAF.txt")
names(mafc) <- c("SNP_ID", "MAF")

# subset files
raqtls <- select(raqtls, all_of(c("SNP_ID", "ref.element.count", "alt.element.count", 
                                  "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean", "hNPC.wilcoxon.pvalue")))
controls <- select(controls, all_of(c("SNP_ID", "ref.element.count", "alt.element.count", 
                                      "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean", "hNPC.wilcoxon.pvalue")))

# annotate with maf, dataset
raqtls <- merge(raqtls, mafr, by="SNP_ID")
raqtls$Dataset <- "emVar"

controls <- merge(controls, mafc, by="SNP_ID")
controls$Dataset <- "control SNP"

# combine the two sets
df <- rbind(raqtls, controls)

# add min max columns
df$max_fragment <- ifelse(df$`ref.element.count` > df$`alt.element.count`, df$`ref.element.count`, df$`alt.element.count`)
df$min_fragment <- ifelse(df$`ref.element.count` < df$`alt.element.count`, df$`ref.element.count`, df$`alt.element.count`)
df$total_fragment <- df$`ref.element.count` + df$`alt.element.count`

df$max_expr <- ifelse(df$`hNPC.cDNA.ref.mean` > df$`hNPC.cDNA.alt.mean`, df$`hNPC.cDNA.ref.mean`, df$`hNPC.cDNA.alt.mean`)
df$min_expr <- ifelse(df$`hNPC.cDNA.ref.mean` < df$`hNPC.cDNA.alt.mean`, df$`hNPC.cDNA.ref.mean`, df$`hNPC.cDNA.alt.mean`)
df$log2_effect <- log2(df$`hNPC.cDNA.ref.mean` / df$`hNPC.cDNA.alt.mean`)


# fold changes and standard deviations
fc <- aggregate(df[, c(7, 9:13)], list(df$Dataset), mean) %>% 
  column_to_rownames(var="Group.1") %>%
  t()

fc <- as.data.table(fc, keep.rownames = "Category")
fc$foldchange = fc$`emVar`/fc$`control SNP`
names(fc) <- c("Category", "mean control SNP", "mean emVar", "foldchange")

standdev <- aggregate(df[, c(7, 9:13)], list(df$Dataset), sd) %>%
  column_to_rownames(var="Group.1") %>% t()

standdev <- as.data.table(standdev, keep.rownames = "Category")
names(standdev) <- c("Category", "std control SNP", "std emVar")

fc <- merge(fc, standdev, by = "Category")

fc$cutoff_con <- fc$`mean control SNP` + 3 * fc$`std control SNP`
fc$cutoff_em <- fc$`mean emVar` + 3 * fc$`std emVar` 

# generate plots with pvalues
mplot <- ggplot(df, aes(x=Dataset, y=MAF)) +
  geom_violin(trim=FALSE, aes(fill = Dataset)) + geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("#56b4e9","#e69f00")) +
  labs(x = "", y = "Minor allele frequency", title = "MAF") +
  theme_bw()


minfplot <- ggplot(df, aes(x=Dataset, y=min_fragment)) +
  geom_violin(trim=FALSE, aes(fill = Dataset)) + geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("#56b4e9","#e69f00")) +
  labs(x = "", y = "Minimal fragment count", title = "Minimal fragment count") +
  theme_bw()


maxfplot <- ggplot(df, aes(x=Dataset, y=max_fragment)) +
  geom_violin(trim=FALSE, aes(fill = Dataset)) + geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("#56b4e9","#e69f00")) +
  labs(x = "", y = "Maximum fragment count", title = "Maximum fragment count") +
  theme_bw()


tmp <- df[(df$Dataset == "control SNP" & df$total_fragment < 534.9094479) |
            +               (df$Dataset == "emVar" & df$total_fragment < 581.3527187), ]
tfplot <- ggplot(tmp, aes(x=Dataset, y=total_fragment)) +
       geom_violin(trim=FALSE, aes(fill = Dataset)) + geom_boxplot(width=0.1) +
       scale_fill_manual(values = c("#56b4e9","#e69f00")) +
       labs(x = "", y = "Total fragment count", title = "Total fragment count (no outliers: mean + 3SD)") +
       theme_bw()


tmp <- df[(df$Dataset == "control SNP" & df$min_expr < 24.0382236) |
            (df$Dataset == "emVar" & df$min_expr < 15.7240142), ]
mineplot <- ggplot(tmp, aes(x=Dataset, y=min_expr)) +
  geom_violin(trim=FALSE, aes(fill = Dataset)) + geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("#56b4e9","#e69f00")) +
  labs(x = "", y = "Minimal expression", title = "Minimal expression (no outliers: mean + 3SD)") +
  theme_bw() 

tmp <- df[(df$Dataset == "control SNP" & df$max_expr < 46.5616448) |
            (df$Dataset == "emVar" & df$max_expr < 50.2751792), ]
maxeplot <- ggplot(tmp, aes(x=Dataset, y=max_expr)) +
  geom_violin(trim=FALSE, aes(fill = Dataset)) + geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("#56b4e9","#e69f00")) +
  labs(x = "", y = "Maximum expression", title = "Maximum expression (no outliers: mean + 3SD)") +
  theme_bw()


# combine the violin plots
combi <- ggarrange(mplot, tfplot, mineplot, maxeplot, common.legend = TRUE)

ggsave("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/raQTLs/freeze7/figures/raQTL_vs_contol_28052025.pdf",
       combi, width = 8, height = 8.2)


# make a correlation plot
cplot <- ggplot(df, aes(x=max_expr, y=max_fragment)) +
  geom_point(aes(color=df$Dataset), alpha=0.3) +
  scale_color_manual(values = c("#56b4e9","#e69f00")) +
  geom_smooth(method = "lm", se=FALSE, color = "black") +
  stat_cor(label.y = 2,aes(label = after_stat(r.label))) +
  labs(y ="Maximum fragment count", x = "Maximum expression", title="Expression versus Fragment count") + 
  theme_bw() +
  facet_wrap(~Dataset)
