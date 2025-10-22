library(ggplot2)
library(dplyr)
library(data.table)
library(forcats)

path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/LDSC/FINAL_FILES/freeze7/"

results <- fread(paste0(path, "LDSC_hNSC_freeze7_summary_withcontrols_20052025.txt"))

ylimit <- max(results$Enrichment + results$Enrichment_std_error + 1)

phenos <- c("CD", "ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ")

results <- results %>% 
    filter(Phenotype %in% phenos) %>%
    mutate(label= ifelse(Prop._h2 > 0, round(Enrichment_p, digits=2), " "))

results$Phenotype= factor(results$Phenotype, levels= phenos)

global_size = 14

ldscplot <- ggplot(results, aes(x=Phenotype, y=Enrichment, fill = Category)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=pmax(Enrichment-Enrichment_std_error, 0), ymax=pmax(0,Enrichment+Enrichment_std_error)),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=label), nudge_y = (results$Enrichment_std_error + 0.5)) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  #scale_fill_manual(values = c("black", "#fb8808", "#d84420", "#c06000", 
  #                             "#7f32bd", "#5770ff","#87BEFF","#fa98c7")) +
  labs(x="Phenotype", y="Enrichment (Prop. h2g / Prop. SNPs)", title = "hNSC raQTLs and controls")+
  theme_classic() +
  geom_hline(yintercept=1, linetype="dashed") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylimit)) +
  theme(text = element_text(size=global_size))

ldscplot

ggsave(paste0(path, "figures/hNSC-summary-CD_heritability-enr_freeze7_wcontrols_20052025.pdf"), ldscplot, units="mm",width=195, height =200)
