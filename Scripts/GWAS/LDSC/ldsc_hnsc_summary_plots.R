library(ggplot2)
library(dplyr)
library(data.table)
library(forcats)

path <- "V:/ddata/CELB/poot/Eline Koornstra/LDSC/FINAL_FILES/freeze7/"

results <- fread(paste0(path, "LDSC_hNSC_freeze7_summary_17042024.txt"))

ylimit <- max(results$Enrichment + results$Enrichment_std_error + 1)

phenos <- c("CD", "ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ")

results <- results %>% 
    filter(Phenotype %in% phenos) %>%
    mutate(label= ifelse(Prop._h2 > 0, round(Enrichment_p, digits=2), " "))

results$Phenotype= factor(results$Phenotype, levels= phenos)

global_size = 14

ldscplot <- ggplot(results, aes(x=Phenotype, y=Enrichment, fill = Phenotype)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=pmax(Enrichment-Enrichment_std_error, 0), ymax=pmax(0,Enrichment+Enrichment_std_error)),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=label), nudge_y = (results$Enrichment_std_error + 0.5)) +
  #scale_fill_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73", 
  #                             "#ECE242", "#0072B2","#D55E00","#C46FA0")) +
  scale_fill_manual(values = c("black", "#fb8808", "#d84420", "#c06000", 
                               "#7f32bd", "#5770ff","#87BEFF","#fa98c7")) +
  labs(x="Phenotype", y="Enrichment (Prop. h2g / Prop. SNPs)", title = "hNSC raQTLs")+
  theme_classic() +
  geom_hline(yintercept=1, linetype="dashed") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylimit)) +
  theme(text = element_text(size=global_size))

ldscplot

ggsave(paste0(path, "figures/hNSC-summary-CD_heritability-enr_freeze7_17042024.pdf"), ldscplot, units="mm",width=175, height =200)
