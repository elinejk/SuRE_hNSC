library(data.table)
library(tidyverse)
library(trackViewer)

### STEP 1: LOAD FILES ###
# paths
path_eqtl <- "/gpfs/home6/ekoornstra/resources/eQTLs/"
path_hits <- "/gpfs/home6/ekoornstra/raqtls/freeze7/annotations/eqtls/"
tissue <- "Cells_Transformed_fibroblasts"

gene <- "ENSG00000142599.13"
gn <- "RERE"


all <- fread(paste0(path_eqtl, tissue, ".v7.signif_variant_gene_pairs.txt.gz"))
all <- all[all$gene_id == gene, ]

### STEP 2: ADD ANNOTATIONS FOR COLORS ###
# load sure file
sure <- fread("/gpfs/home6/ekoornstra/sure_data/freeze7/pvalue.freeze7.snp-permutation.04042024.txt.gz")

sure$chrom <- gsub("chr", "", sure$chrom)
sure[, variant_id := paste0(chrom, "_", pos.hg19, "_", ref.seq, "_", alt.seq, "_b37")]
sure <- sure[, c("SNP_ID", "variant_id", "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean")]

# load raQTLs
raqtls <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/annotations/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024_ANNOTATED.txt")
raqtls[, variant_id := paste0(chrom, "_", pos.hg19, "_", ref.seq, "_", alt.seq, "_b37")]


# add color information
all <- all %>%
	mutate(tested = case_when(all$variant_id %in% sure$variant_id & !all$variant_id %in% raqtls$variant_id ~ "SuRE_tested", 
	all$variant_id %in% raqtls$variant_id ~ "raQTL", 
	!all$variant_id %in% sure$variant_id & !all$variant_id %in% raqtls$variant_id ~ "not SuRE"))
	
			
# add concordance
all <- merge(all, sure, by = "variant_id", all.x=TRUE)
all$hnsc_log2_alt_ref <- log2(all$hNPC.cDNA.alt.mean/all$hNPC.cDNA.ref.mean)
all$raqtl_dir <- ifelse(all$hNPC.cDNA.alt.mean > all$hNPC.cDNA.ref.mean, "pos", "neg")
  
all$slope_dir <- ifelse(all$slope > 0, "pos", "neg")
all <- all %>% mutate(concordant = case_when(all$raqtl_dir == all$slope_dir ~ "concordant", is.na(all$raqtl_dir) ~ "not tested",  
					!all$raqtl_dir == all$slope_dir & !is.na(all$raqtl_dir) ~ "discordant"))



### STEP 3: MAKE THE LOLLIPOP PLOT ###
all$fill <- ifelse(all$concordant == "concordant", "#6C88C4", ifelse(all$concordant == "discordant", "#E77577", "#FFD872"))
all$color <- ifelse(all$tested == "SuRE_tested", "#74737A", ifelse(all$tested == "raQTL", "black", "#cccccc"))			
all <- all %>% separate(variant_id, c("chrom", "position", NA, NA, NA), sep= "_", remove=F)	

# save file inbetween to make the plot on local pc
#fwrite(all, paste0(path_hits, "fibroblasts_eQTLs_RERE_19062024.txt"), sep='\t', row.names=F)
all <- fread("~/Documents/ErasmusMC/work_files/fibroblasts_eQTLs_RERE_19062024.txt")

snps <- GRanges(seqnames=Rle(as.numeric(all$chrom)), ranges=IRanges(start=as.numeric(all$position), width =1))
snps$color <- all$fill
snps$border <- all$color
snps$score <- runif(length(snps))*5
snps$SNPsideID <- ifelse(all$concordant == "concordant", "top", "bottom")

# i manually got the gene information from gencode
features <- GRanges("1", IRanges(c(8412457,8415479,8416160,8418256,8419824,8420172,8421823,8422743,8424116,8424806,8425872,8482787,8525985,8555123,8557465,8568686,8601273,8616534,8617477,8674620,8684369,8716032,8852463,8877219),
                                 width = c(2723,180,146,720,222,1378,113,161,199,92,162,80,98,99,124,48,104,96,105,125,70,468,184,305),
                                 height = 0.1, fill = "red"))


legend <- list(labels = c("Concordant", "Discordant", "Not SuRE tested"),
              fill = c("#6C88C4","#E77577","#FFD872"))

pdf("~/Documents/ErasmusMC/work_files/fibroblasts_RERE_lollipopplot_19062024.pdf")
lolli <- lolliplot(snps, features, legend=legend,
                   cex=0.5)
dev.off()
lolli


pdf(paste0(path_hits, "figures/", tissue, "_RERE_lollipopplot_19062024.pdf"), lolli ,units="mm", width = 200, height = 200)




