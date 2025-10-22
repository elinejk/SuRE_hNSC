library(data.table)
library(ggplot2)
library(dplyr)

## LOAD THE FILE ##
df <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/figures/hnsc-raqtls_gtex-v7_overlap_10042024.txt")

## MAKE THE PLOT ##
# plot settings
df$Tissue = factor(df$Tissue, levels = c("Amygdala", "Anterior Cingulate Cortex BA24", "Caudate Basal Ganglia",
                                         "Cerebellar Hemisphere", "Cerebellum", "Cortex",
                                         "Frontal Cortex BA9", "Hippocampus", "Hypothalamus", 
                                         "Nucleus Accumbens Basal Ganglia", "Putamen Basal Ganglia",
                                         "Spinal Cord Cervical c1", "Substantia Nigra",
                                         "Fibroblasts", "Liver", "Whole blood"))

global_size = 14

# make the plot
pl <- ggplot(df, aes(x=Tissue,y=Set)) +
  geom_point(aes(size=`Significant eQTLs`, color=Concordance), alpha=I(0.99)) +
  geom_text(data = df, aes(x=Tissue,y=Set,label=round(`Significant eQTLs`, 1)), size = 4, nudge_y = - 0.3) +
  geom_text(data = df, aes(x=Tissue,y=Set,label=round(Concordance, 1), fontface = "italic"), size = 4, nudge_y = - 0.5) +
  theme_minimal() +
  #scale_size_area(max_size=9) +
  scale_color_gradientn(colors=colorRampPalette(RColorBrewer::brewer.pal(9,'Reds'))(255)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1.5, hjust=0),
        text = element_text(size=14)) +
  scale_y_discrete(limits = c("hNSC narrowpeak DHS", "hNSC broadpeak DHS", "all hNSC raQTLs"))+
  scale_x_discrete(position = "top") +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 30), 
                        breaks = c(1, 10, 20, 30)) +
  labs(x ="GTEx V7 Tissue", y = "hNSC raQTL set", color = " % Concordance", size = "% raQTLs that are significant eQTLs") +
  guides(color = guide_colorbar(order = 0), size = guide_legend(order = 1))

ggsave('V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/figures/eQTLs/hnsc-raqtls_gtex-v7_overlap_17062024_wlabels.pdf', pl,
       units = "mm", height = 150, width = 350)

