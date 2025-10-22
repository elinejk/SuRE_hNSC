library(data.table)
library(dplyr)
library(ggplot2)

# read file and set parameters
path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/raQTLs/freeze7/figures/"
outpath <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/raQTLs/freeze7/figures/"
df <- fread(paste0(path, "raqtl_FE-OR_for-plot_10042024.txt"))
dag <- "07102024"

tissue <- "E009"
cell <- "hNSC"

if (cell == "hNSC") {
  pal <- c("#56B4E9", "#009E73")
} else if (cell == "HepG2") {
  pal <- c("#C46FA0", "#009E73")
} else if (cell == "K562") {
  pal <- c("#ECE242", "#009E73")
}

global_size = 12

# subset on cell type
a <- df[df$Cell == cell, ]

# set comparison levels
a$Comparison <- factor(a$Comparison, levels = c("raQTLs vs controls", "raQTLs vs all SuRE SNPs"))

# add pvalue asterisks
a <- a %>%
 mutate(label = case_when(
   `Fisher p-value` > 0.05 ~ "", `Fisher p-value` > 0.01 ~ "*",
   `Fisher p-value` > 0.001 ~"**", `Fisher p-value` > 0.0001 ~"***", !is.na(`Fisher p-value`) ~ "****",
   TRUE ~ NA_character_)) 

# Fold-enrichment plots
fe <- ggplot(a, aes(x=Annotation, y=`Fold enrichment`)) +
  geom_bar(aes(fill=Comparison), stat="identity", color = "black", position=position_dodge(), alpha = 0.6, width = 0.7) +
  geom_text(aes(label=label, group = Comparison, y = `Fold enrichment` + 0.2), position = position_dodge(width = .9)) +
  geom_hline(yintercept=1, linetype="dashed") +
  #guides(x = guide_axis(angle = 45)) +
  coord_flip() +
  labs(x = "", y="Fold enrichment") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size=global_size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(paste0(outpath, cell, "_locations_FE_", dag, ".pdf"), fe, units="mm", height=200, width=300)

# Odds ratio plot
fo <- ggplot(a, aes(x=Annotation, y=`Odds ratio`)) +
  geom_bar(aes(fill=Comparison), stat="identity", color = "black", position=position_dodge(), alpha = 0.6) +
  geom_text(aes(label=label, group = Comparison, y = `Odds ratio` + 0.2), position = position_dodge(width = .9)) +
  geom_hline(yintercept=1, linetype="dashed") +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="Odds ratio") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size=global_size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(paste0(outpath, cell, "_locations_OR_", dag, ".pdf"), fo, units="mm", height=200, width=170)


