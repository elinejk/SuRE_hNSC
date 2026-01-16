library(GenomicAlignments)
library(stringr)
library(data.table)
library(ggplot2) 
library(dplyr)

#Set seed for random number generator for control SNPs 
set.seed(123)

#Load heterozygous SNPs 
df <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/120126_Merkle_H9_heterozygous.csv")

#Give colnames
colnames(df) <- c("chr", "pos", "ref", "alt", "AD", "RD", "QUAL", "PL") 

sure_raqtl <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/raqtls/results/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")
sure_ctrl <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/raqtls/results/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024.txt")

#Change chr names
sure_raqtl$chr <- paste0("chr", sure_raqtl$chrom) 
sure_ctrl$chr <- paste0("chr", sure_ctrl$chrom)

#Overlap WGS heterozygotes w/ SuRE data
raqtl_overlap <- merge(df, sure_raqtl, by.x=c("chr", "pos", "ref", "alt"), by.y=c("chr", "pos.hg19", "ref.seq", "alt.seq"))
control_overlap <- merge(df, sure_ctrl, by.x=c("chr", "pos", "ref", "alt"), by.y=c("chr", "pos.hg19", "ref.seq", "alt.seq"))

#Create files for subsetting w/ samtools
hnsc_sign.bed <- data.frame(cbind(raqtl_overlap$chr, raqtl_overlap$pos-1, raqtl_overlap$pos, rep('*', nrow(raqtl_overlap)),rep('0',nrow(raqtl_overlap)))) 
hnsc_ctrl.bed <- data.frame(cbind(control_overlap$chr, control_overlap$pos-1, control_overlap$pos, rep('*', nrow(control_overlap)),rep('0',nrow(control_overlap))))

#Write the tables
write.table(hnsc_sign.bed,file='hnsc.raqtl.hetero.120126.bed',quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
write.table(hnsc_ctrl.bed,file='hnsc.control.hetero.120126.bed',quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)                                 

#Now we continue with the shell script "SuRE_DHS_bam.sh" 

#Get DHS signal for ref/alt alleles
for(i in 1:nrow(raqtl_overlap)) {
  temp.table <- table(c(
    as.character(stackStringsFromBam("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1_hetero_sign_subset.bam", param=GRanges(seqnames = raqtl_overlap$chr[i], ranges=IRanges(start=raqtl_overlap$pos[i], width=1)))),
    as.character(stackStringsFromBam("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2_hetero_sign_subset.bam", param=GRanges(seqnames = raqtl_overlap$chr[i], ranges=IRanges(start=raqtl_overlap$pos[i], width=1))))
  ))

  # Assign the allelic depths
  raqtl_overlap$DHS.ref[i] <- temp.table[raqtl_overlap$ref[i]]
  raqtl_overlap$DHS.alt[i] <- temp.table[raqtl_overlap$alt[i]]  
}


#Same for control SNPs
for(i in 1:nrow(control_overlap)) {
  temp.table <- table(c(
    as.character(stackStringsFromBam("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1_hetero_ctrl_subset.bam", param=GRanges(seqnames = control_overlap$chr[i], ranges=IRanges(start=control_overlap$pos[i], width=1)))),
    as.character(stackStringsFromBam("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2_hetero_ctrl_subset.bam", param=GRanges(seqnames = control_overlap$chr[i], ranges=IRanges(start=control_overlap$pos[i], width=1))))
  ))

  # Assign the allelic depths
  control_overlap$DHS.ref[i] <- temp.table[control_overlap$ref[i]]
  control_overlap$DHS.alt[i] <- temp.table[control_overlap$alt[i]]  
}


fwrite(raqtl_overlap, "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/140126_raqtl_DHS.txt", row.names=F, quote=F, sep="\t") 
fwrite(control_overlap, "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/140126_ctrl_DHS.txt", row.names=F, quote=F, sep="\t") 


#ANALYSIS 
#First subset control SNPs for 2000 random snps
ctrl_overlap <- sample_n(control_overlap, 2000)


analysis_func <- function(df, snp_type, color_code){
  df$AD <- as.character(df$AD)
  
  # Split AD into REF_AD, ALT_AD, NONREF_AD
  df[, c("REF_AD", "ALT_AD", "NONREF_AD") := tstrsplit(AD, ",", type.convert = TRUE)]

  # Create genomic DNA ratio (REF/ALT)
  df$gdna.ratio.ref.over.alt <- df$REF_AD / df$ALT_AD

  # Plot density for sanity check
  dens_plot <- ggplot(df, aes(x = gdna.ratio.ref.over.alt)) +
    geom_density() + 
    labs(title = paste0("Genomic DNA ratio REF / ALT allele for ", snp_type))

  ggsave(paste0("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/plots/140126_density_H9_ESC_REF_vs_ALT_", snp_type, ".png"), plot = dens_plot, device = "png")

  # Handle missing values for DHS counts
  df$DHS.ref[is.na(df$DHS.ref)] <- 0  # Replace NAs with 0
  df$DHS.alt[is.na(df$DHS.alt)] <- 0  # Replace NAs with 0

  # Create total DHS count column
  df$sum.dhs.counts <- df$DHS.ref + df$DHS.alt

  # Filter loci with at least 10 reads in DNase
  df2 <- df[df$sum.dhs.counts >= 10,]

  # Compute the DHS ref/alt ratio
  df2$dhs.ratio.ref.over.alt <- df2$DHS.ref / df2$DHS.alt

  # Avoid infinity and apply pseudo-counts to handle zero values
  df2$dhs.ratio.ref.over.alt.pseudo <- (df2$DHS.ref + 1) / (df2$DHS.alt + 1)

  # Create the ratio for DHS over genomic DNA
  df2$dhs_over_gdna <- log2(df2$dhs.ratio.ref.over.alt.pseudo / df2$gdna.ratio.ref.over.alt)

  # Create the SuRE ratio
  df2$sure_ratio <- log2(df2$hNPC.cDNA.ref.mean / df2$hNPC.cDNA.alt.mean)

  # Check allelic differences using a DHS threshold of 2 
  df2$dhs_log <- log2(df2$dhs.ratio.ref.over.alt.pseudo)
  df3 <- df2[df2$dhs_log >= 2 | df2$dhs_log <= -2,]

  print(paste0("Number of ", snp_type, " with more than 10 reads and >2 allelic difference:", nrow(df3)))
  # Compute Odds Ratio (OR)
  or <- round((sum(df3$dhs_over_gdna > 0 & df3$sure_ratio > 0) + sum(df3$dhs_over_gdna < 0 & df3$sure_ratio < 0)) / 
              (sum(df3$dhs_over_gdna > 0 & df3$sure_ratio < 0) + sum(df3$dhs_over_gdna < 0 & df3$sure_ratio > 0)), digits = 2)



plot_dhs <- ggplot(df3, aes(x = sure_ratio, y = dhs_over_gdna)) +
  geom_point(color = color_code) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    x = "log2 SuRE signal ratio REF/ALT",
    y = "log2(DHS [REF/ALT] / DNA [REF/ALT])"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("OR =", or),
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black", fontface = "bold"
  ) +
  coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10)) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
  )

  ggsave(paste0("/gpfs/work4/0/AdamsLab/Projects/sure/wgs/plots/160126_", snp_type, "_DHS_10reads_DHS2.pdf"), plot = plot_dhs, device = "pdf")
}

# Apply the function to raqtl and ctrl snps
ctrl_final <- analysis_func(ctrl_overlap, "control SNPs", "#56B4E9")
raqtl_final <- analysis_func(raqtl_overlap, "emVars", "#E69F00")
