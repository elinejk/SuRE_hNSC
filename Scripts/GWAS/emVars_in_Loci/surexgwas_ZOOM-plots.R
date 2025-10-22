library(rtracklayer)
library(Gviz)
library(data.table)
library(biomaRt)
library(dplyr)

### SET PATHS ###
sure_path <- "V:/ddata/CELB/poot/Eline Koornstra/raQTLsxGWAS/freeze7/sure-data_cp-locus184.txt"
raqtl_path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt"
gwas_path <- "C:/Users/084902/Documents/GWAS_FILES/ea/GWAS_CP_all.txt.gz"

DNase_path <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/ENCODE/DHS/Stamatoyannopoulos_Lab/2018_broadpeak/ENCFF175IJP_DHS_2018_hg19.bigWig"
enh_path <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/ROADMAP/E009_25_imputed12marks_dense_ACTIVEENHANCERS.bed.gz"
prom_path <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/ROADMAP/E009_25_imputed12marks_dense_PROMOTERS.bed.gz"

out_path <- "V:/ddata/CELB/poot/Eline Koornstra/raQTLsxGWAS/freeze7/locus_plots/CP_locus184_zoomlocus_"

### SET PARAMETERS ###
# zoom plot
chr <- "chr22"
st <- 41945000
ed <- 42026000

GR<-makeGRangesFromDataFrame(data.frame(chr=chr, start=st, end=ed, score=1))

# highlight location
hls <- 41985908
hle <- 41986008


### GENERAL TRACKS ###
gtrack <- GenomeAxisTrack(labelPos = "below", col = "black", fontcolor = "black", lwd = 1, exponent= 0)


###############
### GWAS TRACKS ###
gwas <- fread(gwas_path, fill = TRUE)
cols <- c("CHR", "POS", "MarkerName", "Pval")
cols2 <- c("CHR", "BP", "SNP", "P")
gwas <- select(gwas, all_of(cols))
setnames(gwas, cols2)
gwas$CHR <- paste0("chr", gwas$CHR)
locus <- gwas[gwas$CHR == chr & gwas$BP >= st & gwas$BP <= ed, ]

rm(gwas)
gc()

sure <- fread(sure_path)
raqtls <- fread(raqtl_path)

## MAKE MANHATTAN PLOT WITH SIGN NON SIGN
# prepare your data
df <- locus %>%
  # determine chromosome length, and save as chr_len
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  # determine the cumulative position of each chromosome, saved as tot
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  # join it to the original dataset = data
  left_join(locus, ., by=c("CHR"="CHR")) %>%
  # sort and add the cumulative position of each SNP to prevent overlapping sections
  arrange(CHR,BP) %>%
  mutate(BPcum=BP+tot) %>%
  # add annotation information
  mutate(is_sign = if_else(P <= 5*10^-8, "significant", "not significant")) %>%
  mutate(sure_tested = if_else(SNP %in% sure$SNP_ID & SNP %in% raqtls$SNP_ID, "emVar", 
                               if_else(SNP %in% sure$SNP_ID & !SNP %in% raqtls$SNP_ID, "SuRE tested", "not tested")))

# prepare your axis
axisdf = df %>%
  group_by(CHR) %>%
  summarize(center=(max(BPcum) + min(BPcum)) / 2)


# create the plot
df$logP <- -log10(df$P)

df_qtl <- filter(df, sure_tested == "emVar")
df_sign <- filter(df, sure_tested == "SuRE tested")
df_snp <- filter(df, sure_tested == "not tested")
#df_sign <- filter(df, is_sign == "significant")
#df_snp <- filter(df, is_sign == "not significant")

# create granges
manh_sign <- GRanges(seqnames = Rle(df_sign$CHR), ranges = IRanges(start=df_sign$BP, width=1))
manh_snp <- GRanges(seqnames = Rle(df_snp$CHR), ranges = IRanges(start=df_snp$BP, width=1))
manh_qtl <- GRanges(seqnames = Rle(df_qtl$CHR), ranges = IRanges(start=df_qtl$BP, width=1))

manh_sign$p <- df_sign$logP
manh_snp$p <- df_snp$logP
manh_qtl$p <- df_qtl$logP

################
# create the separate manhattan tracks
mtrack1 <- DataTrack(manh_snp, name="-log10(P)", 
                     chromosome = chr, 
                     start = st, 
                     end = ed, 
                     col="#C7C7C7", 
                     type = "p", 
                     ylim=c(min(df$logP), max(df$logP)), 
                     baseline = 7.30103, col.baseline = "black", lty.baseline = "dashed",
                     cex = 0.7,
                     background.title = "white", col.title = "black", col.axis = "black")

if(exists('manh_sign') == TRUE && length(manh_sign) > 0) {
  get('manh_sign')
  mtrack2 <- DataTrack(manh_sign, name="-log10(P)", 
                       chromosome = chr, 
                       start = st, 
                       end = ed, 
                       col="#8D908D", 
                       type = "p", 
                       ylim=c(min(df$logP), max(df$logP)), 
                       cex = 0.7,
                       background.title = "white", col.title = "black", col.axis = "black")
}

mtrack3 <- DataTrack(manh_qtl, name="-log10(P)", 
                     chromosome = chr, 
                     start = st, 
                     end = ed, 
                     col="#E69F00", 
                     type = "p", 
                     ylim=c(min(df$logP), max(df$logP)), 
                     baseline = 7.30103, col.baseline = "black", lty.baseline = "dashed",
                     cex = 0.7,
                     background.title = "white", col.title = "black", col.axis = "black")

list = list(if(exists('mtrack1') == TRUE) {get('mtrack1')} else {NA}, 
            if(exists('mtrack2') == TRUE) {get('mtrack2')} else {NA},
            if(exists('mtrack3') == TRUE) {get('mtrack3')} else {NA})

mtrack <- OverlayTrack(trackList=list[!is.na(list)], name="GWAS -log10(P)",
                       background.title = "white", col.title = "black", col.axis = "black")


### SuRE TRACK ###
#sure$chrom <- paste0("chr", sure$chrom)
sure <- sure[sure$chrom == chr & sure$pos.hg19 >= st & sure$pos.hg19 <= ed, ]
sure <- sure[order(chrom, pos.hg19), ]
sure <- sure[sure$SNP_ID %in% locus$SNP]
sure$status <- ifelse(sure$SNP_ID %in% raqtls$SNP_ID, "emvar", "no emvar")

raqtl <- sure[sure$status == "emvar", ]
rm(raqtls)

sr <- GRanges(seqnames =Rle(raqtl$chrom), 
              ranges=IRanges(start=raqtl$pos.hg19, width=1))
so <- GRanges(seqnames =Rle(sure$chrom), 
              ranges=IRanges(start=sure$pos.hg19, width=1))

rtrack <- AnnotationTrack(sr, name = "emVars", chromosome = chr, start = st, end = ed, 
                          col="#E69F00", 
                          fill = "#E69F00",
                          lwd = 0.1, 
                          background.title = "white", 
                          col.title = "black", 
                          col.axis = "black")

strack <- AnnotationTrack(so, name = "SuRE-tested GWAS SNPs", chromosome = chr, start = st, end = ed, 
                          col="#000000", 
                          fill = "black",
                          lwd = 0.1, 
                          background.title = "white", 
                          col.title = "black", 
                          col.axis = "black")

#strack <- OverlayTrack(list(strack1, strack2), 
#                       name = "SuRE SNPs", 
#                       chromosome = chr, 
#                       start = lst, 
#                       end = led, 
#                       background.title = "white", col.title = "black", col.axis = "black")




###############
### GENE TRACK ###
bm <- useEnsembl(host = "https://grch37.ensembl.org", 
                 biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "hsapiens_gene_ensembl")
biomTrack<-BiomartGeneRegionTrack(genome="hg19",
                                  chromosome=chr,
                                  start=st,
                                  end=ed,
                                  biomart=bm,
                                  filter = list(biotype = "protein_coding"),
                                  collapseTranscripts="meta",
                                  transcriptAnnotation="symbol",
                                  just.group="right",
                                  fontcolor.group = "black",
                                  cex.group=0.5,
                                  shape="smallArrow",
                                  col.line=NULL, col=NULL, fill="darkred",
                                  lwd = 0.1,
                                  #background.panel = "#EEEEEE",
                                  name="Protein-coding genes")





###############
### DNASE TRACK ###
DNase <- import(DNase_path, format = "bigWig", object = "GRanges", selection = BigWigSelection(ranges=GR), which=DNase_path)
dnaseTrack<-DataTrack(range=DNase,
                      genome="hg19",
                      chromosome=chr,
                      start=st,
                      end=ed,
                      type="h",
                      #ylim=c(0,1000),
                      col="darkgreen",
                      fill= "darkgreen",
                      name="DNaseI")

###############
### CHROMHMM TRACK ###
enh <- import(enh_path, format = "bed", which=GR, genome = "hg19")
prom <- import(prom_path, format = "bed", which=GR, genome = "hg19")

etrack <- AnnotationTrack(enh, 
                          chromosome = chr, 
                          start = st, 
                          end = ed, 
                          stacking="dense",
                          col = "#b86c7d",
                          fill = "#b86c7d",
                          col.line=NULL,
                          name = "Enhancer")

ptrack <- AnnotationTrack(prom, 
                          chromosome = chr, 
                          start = st, 
                          end = ed, 
                          stacking="dense",
                          col = "#6c7db8",
                          fill = "#6c7db8",
                          col.line=NULL,
                          name = "Promoter")

###############
### HIGHLIGHT TRACK ###
ht <- HighlightTrack(trackList = list(mtrack, rtrack, biomTrack, ptrack, etrack, dnaseTrack),
                     start=hls, end=hle, chromosome = chr,
                     col="#d8d4e9", fill="#d8d4e9")

###############
### FINAL TRACK ###
plotTracks(list(ht, gtrack),
           chromosome=chr,
           from=st,
           to=ed,
           sizes=c(1, 0.2, 0.8, 0.6, 0.6, 0.5, 1),
           fontcolor.title="black",
           col.border.title="transparent",
           background.title="transparent",
           col.axis="black"
)

pdf(paste0(out_path, chr, "_", st, "_", ed, "_annotationplot.pdf"), width=8, height=4)
plotTracks(list(ht, gtrack),
           chromosome=chr,
           from=st,
           to=ed,
           sizes=c(1, 0.3, 0.6, 0.4, 0.4, 0.4, 0.4),
           fontcolor.title="black",
           col.border.title="transparent",
           background.title="transparent",
           col.axis="black"
)
dev.off()

