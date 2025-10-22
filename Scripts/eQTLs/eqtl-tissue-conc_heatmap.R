library(data.table)
library(dplyr)
library(ggplot2)

# paths
path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/"
gc_path <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/GENCODE/gencode.v19.genenames.txt.gz"

tissues <- c("Frontal_Cortex_BA9", "Anterior_cingulate_cortex_BA24", "Amygdala", "Cerebellum", 
             "Nucleus_accumbens_basal_ganglia", "Cortex", "Hypothalamus",
             "Hippocampus", "Putamen_basal_ganglia", "Caudate_basal_ganglia",
             "Spinal_cord_cervical_c-1", "Cerebellar_Hemisphere", "Substantia_nigra", 
             "Cells_Transformed_fibroblasts", "Whole_Blood", "Liver")

snp <- "rs4822045"

# read in each file and put in a single datatable
eqtls <- data.table()

for (i in tissues) {
  tissue <- fread(paste0(path, "overlap/eQTLs/hit_files/raqtl_eqtl-", i, "_all-raqtls_all-hits_18062024.txt"))
  tissue <- tissue[tissue$SNP_ID == snp, ]
  tissue$tissue <- i
  cols <- c("tissue", "gene_id", "concordant")
  sub <- select(tissue, all_of(cols))
  
  eqtls <- rbind(eqtls, sub)
}

# add missing values
vals <- expand.grid(tissue = unique(eqtls$tissue),
                    gene_id = unique(eqtls$gene_id))
df <- merge(vals, eqtls, all = TRUE)
df$tissue = factor(df$tissue, levels = c("Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
                                         "Cerebellar_Hemisphere", "Cerebellum", "Cortex",
                                         "Frontal_Cortex_BA9", "Hippocampus", "Hypothalamus", 
                                         "Nucleus_accumbens_basal_ganglia", "Putamen_basal_ganglia",
                                         "Spinal_cord_cervical_c-1", "Substantia_nigra",
                                         "Cells_Transformed_fibroblasts", "Liver", "Whole_Blood"))

# add gene names
gc <- fread(gc_path)
df <- df %>% separate(gene_id, c("ensembl_id", NA), remove=F)

a <- merge(df, gc, by.x = "ensembl_id", by.y="gene_id", all.x=T)

# now we have a subsetted file with all the eGenes, make a heatmap
q <- ggplot(a, aes(x=gene_name, y=fct_rev(as_factor(tissue)), fill = concordant)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("#6C88C4", "#E77577"), na.value = "white") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))

ggsave(paste0(path, "/figures/eQTLs/",snp, "_eQTL-overview_03072024.pdf"), q, height = 4, width = 5)




