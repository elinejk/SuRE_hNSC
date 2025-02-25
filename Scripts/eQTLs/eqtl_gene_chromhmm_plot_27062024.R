library(rtracklayer)
library(Gviz)
library(data.table)
library(biomaRt)
library(dplyr)

### PATHS
enh_path <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/ROADMAP/"
npc <- "E009_25_imputed12marks_dense.bed.gz"
fib_skin <- "E126_25_imputed12marks_dense.bed.gz"
fib_lung <- "E128_25_imputed12marks_dense.bed.gz"

DNase_path <- "V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/ENCODE/DHS/Stamatoyannopoulos_Lab/2018_broadpeak/ENCFF175IJP_DHS_2018_hg19.bigWig"

eqtl_path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/figures/eQTLs/fibroblasts/lolliplot/fibroblasts_eQTLs_RERE_01072024.txt"
out_path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/figures/eQTLs/fibroblasts/"
### SET PARAMETERS
gn <- "RERE"

# coordinates plot
chr <- "chr1"
st <- 8402457
ed <- 8887702

GR<-makeGRangesFromDataFrame(data.frame(chr=chr, start=st, end=ed, score=1))

### eQTL TRACK
eqtl <- fread(eqtl_path)
er <- GRanges(seqnames =Rle(eqtl$chrom), 
              ranges=IRanges(start=eqtl$position, width=1),
              fill = eqtl$fill)

etrack <- AnnotationTrack(er, name = "eQTLs", chromosome = chr, start = st, end = ed, 
                          col=NULL, 
                          fill = er$fill,
                          lwd = 0.1, 
                          background.title = "white", 
                          col.title = "black", 
                          col.axis = "black",
                          stacking = "dense")

### GENE TRACK 
bm <- useEnsembl(host = "https://grch37.ensembl.org", 
                 biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "hsapiens_gene_ensembl")
biomTrack<-BiomartGeneRegionTrack(genome="hg19",
                                  chromosome=chr,
                                  start=st,
                                  end=ed,
                                  biomart=bm,
                                  filter = list(biotype = "protein_coding",
                                                external_gene_name = gn),
                                  #symbol = gn,
                                  transcriptAnnotation="transcript",
                                  just.group="right",
                                  fontcolor.group = "black",
                                  cex.group=0.7,
                                  shape="smallArrow",
                                  col.line="black", col="black", fill="darkred", 
                                  lwd = 0.7,
                                  #background.panel = "#EEEEEE",
                                  name=gn,
                                  background.title = "white", 
                                  col.title = "black", 
                                  col.axis = "black")


### CHROM HMM
# NPC
e009 <- import(paste0(enh_path, npc), format = "bed", which=GR, genome = "hg19")
e126 <- import(paste0(enh_path, fib_skin), format = "bed", which=GR, genome = "hg19")
e128 <- import(paste0(enh_path, fib_lung), format = "bed", which=GR, genome = "hg19")


ntrack <- AnnotationTrack(e009, 
                          chromosome = chr, 
                          start = st, 
                          end = ed, 
                          stacking="dense",
                          fill = e009$itemRgb, col = NULL, col.line = NULL,
                          id = e009$name,
                          name = "E009",
                          background.title = "white", 
                          col.title = "black", 
                          col.axis = "black")

ftrack <- AnnotationTrack(e126, 
                          chromosome = chr, 
                          start = st, 
                          end = ed, 
                          stacking="dense",
                          fill = e126$itemRgb, col = NULL, col.line = NULL,
                          id = e126$name,
                          name = "E126",
                          background.title = "white", 
                          col.title = "black", 
                          col.axis = "black")

ltrack <- AnnotationTrack(e128, 
                          chromosome = chr, 
                          start = st, 
                          end = ed, 
                          stacking="dense",
                          fill = e128$itemRgb, col = NULL, col.line = NULL,
                          id = e128$name,
                          name = "E128",
                          background.title = "white", 
                          col.title = "black", 
                          col.axis = "black")




### DNASE TRACK
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
                      name="DNaseI",
                      background.title = "white", 
                      col.title = "black", 
                      col.axis = "black")

# legend track
gtrack <- GenomeAxisTrack(labelPos = "below", col = "black", fontcolor = "black", lwd = 1, exponent= 0)

### FINAL TRACK
plotTracks(list(etrack, biomTrack, ntrack, ftrack, ltrack, dnaseTrack, gtrack),
           chromosome=chr,
           from=st,
           to=ed,
           #sizes=c(1, 0.2, 0.8, 0.6, 0.6, 0.5, 1),
           fontcolor.title="black",
           col.border.title="transparent",
           background.title="transparent",
           col.axis="black"
)

pdf(paste0(out_path, gn, "_", chr, "_", st, "_", ed, "_annotationplot_03072024.pdf"), width=8, height=8)
plotTracks(list(etrack, biomTrack, ntrack, ftrack, ltrack, dnaseTrack, gtrack),
           chromosome=chr,
           from=st,
           to=ed,
           #sizes=c(1, 0.2, 0.8, 0.6, 0.6, 0.5, 1),
           fontcolor.title="black",
           col.border.title="transparent",
           background.title="transparent",
           col.axis="black"
)
dev.off()
















