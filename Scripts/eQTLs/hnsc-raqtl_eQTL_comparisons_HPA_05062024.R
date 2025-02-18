library(data.table)
library(tidyverse)

### STEP 1: PREPARE THE INITIAL FILES ###
# load files
path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_Project/raQTLs/freeze7/"
genelist <- "all_Cells_Transformed_fibroblasts_raqtl-eqtl-candidates-per-gene_id_19062024_ANNOTATED"
genes <- fread(paste0(path, "overlap/eQTLs/", genelist, ".txt"))

fibs <- fread(paste0(path, "overlap/eQTLs/fibroblast_1000_random_egenes_02072024_GENCODE_ANNOTATED.txt"))

annotpath <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/HPA/"
hpa <- fread(paste0(annotpath, "proteinatlas.tsv"))


# subset hpa file
cols <- c("Ensembl","Gene", "Protein class", "RNA tissue specificity", "RNA tissue distribution", "RNA tissue specific nTPM",
          "RNA single cell type specificity", "RNA single cell type distribution", "RNA single cell type specific nTPM")
hpa <- select(hpa, all_of(cols))


# subset genes file
genes <- select(genes, all_of(c("SNP_ID", "egene_id", "egene_name", "egene_type", "eqtl_concordance",
                                "NPC_broadpeak_DHS", "NPC_narrowpeak_DHS"))) 

a <- genes[genes$egene_type == "protein_coding", ]
#a <- a[a$NPC_broadpeak_DHS == "1", ]
a <- a[a$NPC_narrowpeak_DHS == "1", ]

a <- a %>% separate(egene_id, into = c("egene_ensembl", NA), sep='\\.', remove=FALSE)

cols <- c("SNP_ID", "egene_ensembl", "egene_name", "eqtl_concordance")
df <- select(a, all_of(cols))

# add random fibroblast genes to file
fibs$SNP_ID <- "X"
fibs$eqtl_concordance <- "random_egenes"

fibs <- select(fibs, all_of(cols))

df <- rbind(df, fibs)

# merge gene and hpa files
annot_df <- merge(df, hpa, by.x = "egene_ensembl", by.y = "Ensembl")



### STEP 2: PROTEIN FIGURES ###
# split column into multiple rows
protein <- annot_df  %>% 
  separate_longer_delim(`Protein class` , delim = ", ") # split column into rows

protein <- select(protein, all_of(c("egene_name","SNP_ID","egene_ensembl","eqtl_concordance","Protein class")))

# get the counts for the plot
p_stat <- table(protein$`Protein class`, protein$eqtl_concordance) %>% as.data.table()
names(p_stat) <- c("protein_class", "eqtl_concordance", "count")

# make the plot
p_pl <- ggplot(p_stat, aes(x = eqtl_concordance, y = count, fill = protein_class)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  theme_bw()

p_pl


### STEP 3: TISSUE FIGURES ###
# split column into multiple rows
tissue <- annot_df %>% 
  separate_longer_delim(`RNA tissue specific nTPM` , delim = ";") %>% # split column into rows
  separate(`RNA tissue specific nTPM`, c("Enhanced Tissue", "RNA tissue specific nTPM"), sep = ": ") # split column into columns

tissue <- select(tissue, all_of(c("egene_name","SNP_ID","egene_ensembl","eqtl_concordance",
                                  "RNA tissue specificity", "RNA tissue distribution","Enhanced Tissue")))


# create a file without the tissues, and only retain unique genes (duplicates only for concordance)
b <- select(annot_df, all_of(c("egene_name","RNA tissue specificity","RNA tissue distribution" ,"eqtl_concordance")))
b <- distinct(b)

# specificity figures
tspec <- table(b$`RNA tissue specificity`, b$eqtl_concordance) %>% as.data.table()
names(tspec) <- c("tissue_specificity", "eqtl_concordance", "count")
tspec$tissue_specificity <- factor(tspec$tissue_specificity, levels = c("Not detected", "Group enriched", "Tissue enriched", 
                                                                        "Tissue enhanced", "Low tissue specificity"))

tspec_pl <- ggplot(tspec, aes(x = eqtl_concordance, y = count, fill = tissue_specificity)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  scale_fill_manual(values=c("#117733", "#DDCC77","#CC6677","#88CCEE", "#332288")) +
  labs(y="Percentage") +
  theme_bw()

ggsave(paste0(path,"figures/eqtls/fibroblasts/fibroblast-eqtls_hnsc.np-dhs.protcod.hpa_tissue-specificity_03072024.pdf"),
       tspec_pl, units="mm",width=250, height =200)

# distribution figures
tdis <- table(b$`RNA tissue distribution`, b$eqtl_concordance) %>% as.data.table()
names(tdis) <- c("tissue_distribution", "eqtl_concordance", "count")
tdis$tissue_distribution <- factor(tdis$tissue_distribution, levels = c("Not detected", "Detected in single", "Detected in some",
                                                                        "Detected in many", "Detected in all"))

tdis_pl <- ggplot(tdis, aes(x = eqtl_concordance, y = count, fill = tissue_distribution)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  scale_fill_manual(values=c("#117733", "#DDCC77","#CC6677","#88CCEE", "#332288")) +
  labs(y="Percentage") +
  theme_bw()

ggsave(paste0(path,"figures/eqtls/fibroblasts/fibroblast-eqtls_hnsc.np-dhs.protcod.hpa_tissue-distribution_03072024.pdf"),
       tdis_pl, units="mm",width=250, height =200)

# tissues figures
tis <- table(tissue$`Enhanced Tissue`, tissue$eqtl_concordance) %>% as.data.table()
names(tis) <- c("Enhanced/Enriched tissues", "eqtl_concordance", "count")

tis_pl <- ggplot(tis, aes(x = eqtl_concordance, y = count, fill = `Enhanced/Enriched tissues`)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  theme_bw()



### STEP 3: CELL TYPE FIGURES ###
# split column into multiple rows
cell <- a %>% 
  separate_longer_delim(`RNA single cell type specific nTPM` , delim = ";") %>% # split column into rows
  separate(`RNA single cell type specific nTPM`, c("Enhanced Cell Type", "RNA cell type specific nTPM"), sep = ": ") # split column into columns

cell <- select(cell, all_of(c("egene_name","SNP_ID","egene_ensembl","eqtl_concordance",
                              "RNA single cell type specificity", "RNA single cell type distribution","Enhanced Cell Type")))

# create a file without the cell types, and only retain unique genes (duplicates only for concordance)
c <- select(a, all_of(c("egene_name","RNA single cell type specificity","RNA single cell type distribution", "eqtl_concordance")))
c <- distinct(c)


# specificity figures
cspec <- table(c$`RNA single cell type specificity`, c$eqtl_concordance) %>% as.data.table()
names(cspec) <- c("celltype_specificity", "eqtl_concordance", "count")
cspec$celltype_specificity <- factor(cspec$celltype_specificity, levels = c("Not detected", "Group enriched", "Cell type enriched", 
                                                                        "Cell type enhanced", "Low cell type specificity"))

cspec_pl <- ggplot(cspec, aes(x = eqtl_concordance, y = count, fill = celltype_specificity)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  labs(y="Percentage") +
  theme_bw()

ggsave(paste0(path,"figures/eqtls/fibroblasts/fibroblast-eqtls_hnsc.np-dhs.protcod.hpa_celltype-specificity_03072024.pdf"),
       cspec_pl, units="mm",width=200, height =200)


# distribution figures
cdis <- table(c$`RNA single cell type distribution`, c$eqtl_concordance) %>% as.data.table()
names(cdis) <- c("celltype_distribution", "eqtl_concordance", "count")
cdis$celltype_distribution <- factor(cdis$celltype_distribution, levels = c("Not detected", "Detected in single", "Detected in some",
                                                                            "Detected in many", "Detected in all"))

cdis_pl <- ggplot(cdis, aes(x = eqtl_concordance, y = count, fill = celltype_distribution)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  labs(y="Percentage") +
  theme_bw()

ggsave(paste0(path,"figures/eqtls/fibroblasts/fibroblast-eqtls_hnsc.np-dhs.protcod.hpa_celltype-distribution_03072024.pdf"),
       cdis_pl, units="mm",width=200, height =200)


# cell types figures
cel <- table(cell$`Enhanced Cell Type`, cell$eqtl_concordance) %>% as.data.table()
names(cel) <- c("Enhanced/Enriched cell types", "eqtl_concordance", "count")


cel_pl <- ggplot(cel, aes(x = eqtl_concordance, y = count, fill = `Enhanced/Enriched cell types`)) +
  geom_bar(position="fill", stat = "identity", color = "black") +
  theme_bw() +
  theme(legend.position = "bottom")



### STEP 4: GET PERCENTAGES
tspec %>% mutate(fraction = ifelse(eqtl_concordance == "conc", count/sum(tspec[tspec$eqtl_concordance == "conc", ]$count), 
                                   count/sum(tspec[tspec$eqtl_concordance == "disc", ]$count)))

tdis %>% mutate(fraction = ifelse(eqtl_concordance == "conc", count/sum(tdis[tdis$eqtl_concordance == "conc", ]$count), 
                           ifelse(eqtl_concordance == "disc", count/sum(tdis[tdis$eqtl_concordance == "disc", ]$count),
                                  count/sum(tdis[tdis$eqtl_concordance == "random_egenes", ]$count))))

cspec %>% mutate(fraction = ifelse(eqtl_concordance == "conc", count/sum(cspec[cspec$eqtl_concordance == "conc", ]$count), 
                                   count/sum(cspec[cspec$eqtl_concordance == "disc", ]$count)))

cdis %>% mutate(fraction = ifelse(eqtl_concordance == "conc", count/sum(cdis[cdis$eqtl_concordance == "conc", ]$count), 
                                   count/sum(cdis[cdis$eqtl_concordance == "disc", ]$count)))


