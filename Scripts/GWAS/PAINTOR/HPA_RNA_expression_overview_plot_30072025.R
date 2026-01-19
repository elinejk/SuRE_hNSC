library(data.table)
library(ggplot2)
library(RColorBrewer)

## PATHS ##
hpa <- fread("V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/HPA/rna_tissue_consensus.tsv")
outpath <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/HPA_expression/"

gene <- "C1orf216"

## TISSUE ORDER AND ORGAN ASSIGNMENTS ##
tissue_order <- c("amygdala", "basal ganglia", "cerebellum", "cerebral cortex", "choroid plexus",
                  "hippocampal formation", "hypothalamus", "medulla oblongata", "midbrain", "pons",
                  "spinal cord", "thalamus", "white matter", # brain
                  "retina", # eye
                  "adrenal gland", "parathyroid gland", "pituitary gland", "thyroid gland", # endocrine
                  "lung", # respiratory
                  "esophagus", "salivary gland", "tongue", # proximal digestive
                  "colon", "duodenum", "rectum", "small intestine", "stomach", # gastrointestinal
                  "gallbladder", "liver", # liver and gallbladder
                  "pancreas", # pancreas
                  "kidney", "urinary bladder", # kidney and bladder
                  "epididymis", "prostate", "seminal vesicle", "testis", # male reproductive
                  "breast", "cervix", "endometrium", "fallopian tube", "ovary", "placenta", "vagina", # female reproductive
                  "heart muscle", "skeletal muscle", "smooth muscle", # muscles
                  "adipose tissue", # adipose tissue
                  "skin", # skin
                  "appendix", "bone marrow", "lymph node", "spleen", "thymus", "tonsil" # bone marrow and lymphoid
                  )

organ_order <- c("brain", "eye", "endocrine", "respiratory", "proximal digestive", "gastrointestinal", 
                 "liver & gallbladder",  "pancreas", "kidney & bladder", "reproductive (M)", 
                 "reproductive (F)", "muscle", "adipose tissue", "skin", "bone marrow & lymphoid")

organ_df <- data.table(Tissue = tissue_order,
                       Organ = c("brain", "brain", "brain", "brain", "brain", "brain", "brain", "brain", "brain", "brain",
                                 "brain", "brain", "brain", "eye", "endocrine", "endocrine", "endocrine", "endocrine", 
                                 "respiratory", "proximal digestive", "proximal digestive", "proximal digestive", 
                                 "gastrointestinal", "gastrointestinal", "gastrointestinal", "gastrointestinal", "gastrointestinal",  
                                 "liver & gallbladder", "liver & gallbladder", "pancreas", "kidney & bladder", "kidney & bladder", 
                                 "reproductive (M)", "reproductive (M)", "reproductive (M)", "reproductive (M)", 
                                 "reproductive (F)", "reproductive (F)", "reproductive (F)", "reproductive (F)", 
                                 "reproductive (F)", "reproductive (F)", "reproductive (F)", "muscle", "muscle", "muscle",
                                 "adipose tissue", "skin", "bone marrow & lymphoid", "bone marrow & lymphoid", 
                                 "bone marrow & lymphoid", "bone marrow & lymphoid", "bone marrow & lymphoid", "bone marrow & lymphoid"))


## PREPARE THE DATATABLE FOR YOUR GENE OF INTEREST ##
hpa_sub <- hpa[hpa$`Gene name` == gene, ] # subset for your gene

hpa_sub <- merge(hpa_sub, organ_df, by = "Tissue", all.x = TRUE) # annotate with organs for plotting


## GENERATE YOUR PLOT ##
custom_palette <- c("#A6CEE3", "#1F78B4","#B2DF8A", "#33A02C","#FB9A99", "#E31A1C","#FDBF6F", "#FF7F00",
                    "#FFFF99", "#CAB2D6", "#6A3D9A","#B15928","#C4A484", "#bebebe", "#252525")

hpa_sub$Tissue <- factor(hpa_sub$Tissue, levels = tissue_order)
hpa_sub$Organ <- factor(hpa_sub$Organ, levels = organ_order)

hpa_plot <- ggplot(hpa_sub, aes(x = Tissue, y = nTPM, fill = Organ)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = custom_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(paste0(outpath, gene, "_RNAconsensus_30072025.pdf"), hpa_plot,
       width = 11, height = 5)
