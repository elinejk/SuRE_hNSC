library(data.table)

# load files, put them into one data frame, create a column that indicates from which file it came
files <- sort(list.files(path = 'V:/ddata/CELB/poot/Eline Koornstra/PRS/LDSC/1kg_phase3/',
                        pattern= '_wraqtls_1kgphase3.results', full.name=T))
df <- rbindlist(lapply(files, fread))
df$file <- gsub('_partitionedheritability_baselineLDv2.2_wraqtls_1kgphase3.results','',rep(basename(files),each=100))

rows <- c("Coding_UCSCL2_0","Conserved_LindbladTohL2_0","CTCF_HoffmanL2_0", "DGF_ENCODEL2_0","DHS_peaks_TrynkaL2_0",
          "Enhancer_AnderssonL2_0","Enhancer_HoffmanL2_0","FetalDHS_TrynkaL2_0","H3K27ac_HniszL2_0","H3K27ac_PGC2L2_0",
          "H3K4me1_peaks_TrynkaL2_0","H3K4me3_peaks_TrynkaL2_0","H3K9ac_peaks_TrynkaL2_0","Intron_UCSCL2_0",
          "PromoterFlanking_HoffmanL2_0","Promoter_UCSCL2_0","Repressed_HoffmanL2_0","SuperEnhancer_HniszL2_0",
          "TFBS_ENCODEL2_0","Transcr_HoffmanL2_0","TSS_HoffmanL2_0","UTR_3_UCSCL2_0","UTR_5_UCSCL2_0",
          "WeakEnhancer_HoffmanL2_0","hnpc_ziw_raqtlL2_0","hepg2_raqtlL2_0","k562_raqtlL2_0")

df <- df[df$Category %in% rows, ]

# determine z score enrichment
df$Prop_lowbd = pmax(0,df$`Prop._h2` - 1. * df$`Prop._h2_std_error`)

df$Enrichment_z = qnorm(data.matrix(df$Enrichment_p)/2) * (2*(data.matrix(df$Enrichment)<0)-1)
df$Enrichment_z = pmax(0, df$Enrichment_z )

# change category names
df$class = gsub('L2_0','',df$Category)

# make the plot
library(ggplot2)

ggplot(df) +
  geom_point(aes(x=class,y=file, size=Prop_lowbd, color=Enrichment_z), alpha=I(0.99)) +
  theme_minimal() +
  scale_size_area(max_size=9) +
  scale_color_gradientn(colors=colorRampPalette(RColorBrewer::brewer.pal(9,'Reds'))(255)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x ="Annotation category", y = "Trait", color = "enrichment z-score", size = "proportion of explained heritability")



ggsave('V:/ddata/CELB/poot/Eline Koornstra/PRS/LDSC/1kg_phase3/gwas.ldsr.partition.pdf')
