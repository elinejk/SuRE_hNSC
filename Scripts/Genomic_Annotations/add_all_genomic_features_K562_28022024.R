library(GenomicRanges)
library(IRanges)
library(data.table)
library(dplyr)
library(rtracklayer)

### STEP 1: LOAD THE FILE AND CONVERT TO GRRANGES OBJECT ###
# where possible, use files that contain the tfbs overlap information
args <- commandArgs(trailingOnly=TRUE)
sure <- args[[1]]
outfile <- args[[2]]
an_path <- "/gpfs/home6/ekoornstra/resources/"

df <- fread(sure)
df_range <- GRanges(seqnames = Rle(df$chrom), 
                    ranges=IRanges(start=df$pos.hg19,width =1)) 

df_range$SNP_ID <- df$SNP_ID
#df_range$eqtl_variant <- df$eqtl_variant

### STEP 2: ADD THE GENERAL GENOMIC FEATURES ###
## 2.1 PROMOTERS ##
print("determining UCSC promoters")
p <- fread(paste0(an_path, "LDSC/Promoter_UCSC.bed"))
names(p) <- c("chr", "start", "end")
p[, chr := sub("chr","",chr)]
pr <- GRanges(seqnames= Rle(p$chr), ranges = IRanges(start = p$start, end = p$end))

closest <- distanceToNearest(df_range, subject = pr)
df_range$UCSC_prom_dist <- NA
df_range$UCSC_prom_dist[queryHits(closest)] <- mcols(closest)$distance

head(df_range)
rm(p)

## 2.2 ENHANCERS ##
print("determining Hoffman enhancers")
e <- fread(paste0(an_path, "LDSC/Enhancer_Hoffman.bed"))
names(e) <- c("chr", "start", "end")
e[, chr := sub("chr","",chr)]
er <- GRanges(seqnames= Rle(e$chr), ranges = IRanges(start = e$start, end = e$end))

closest <- distanceToNearest(df_range, subject = er)
df_range$Hoffman_enh_dist <- NA
df_range$Hoffman_enh_dist[queryHits(closest)] <- mcols(closest)$distance

head(df_range)
rm(e)

## 2.3 GENOME LOCATION ##
print("determining gencode location type")
gc <- rtracklayer::import(paste0(an_path, "gencode.v19.annotation.gtf.gz"))
seqlevels(gc) <- gsub("chr","",seqlevels(gc))

ov <- findOverlaps(df_range, gc)
ov_raqtls <- df_range[queryHits(ov)]
ov_genes <- gc[subjectHits(ov)]
ov_raqtls$GENCODE_location <- ov_genes$type

gc_hit <- as.data.table(ov_raqtls)
gc_hit <- unique(gc_hit)
gc_hit <- select(gc_hit, c("SNP_ID", "GENCODE_location"))
gc_hit <- aggregate(GENCODE_location ~ SNP_ID, gc_hit, FUN = paste, collapse=",")

head(gc_hit)
rm(gc)

## 2.4 CLOSEST PROTEIN CODING GENE AND TSS ##
print("determining closest gene and tss")
gencode <- fread(paste0(an_path, "gencode.v19.annotation.genes.bed.gz")) %>% as.data.table()
gene <- gencode[gencode$gene_type == "protein_coding", ]
gene[, chr := sub("chr","",chr)]
gene$tss <- ifelse(gene$strand == "+", gene$start, gene$end)

# determine the closest gene
genelist <- GRanges(seqnames = Rle(gene$chr), 
                    ranges=IRanges(start=gene$start, end=gene$end))
genelist$gene_id <- gene$gene_id
genelist$gene_strand <- gene$strand

closest <- distanceToNearest(df_range, subject = genelist)
df_range$closest_gene <- NA
df_range$gene_distance <- NA
df_range$gene_strand <- NA
df_range$closest_gene[queryHits(closest)] <- genelist[subjectHits(closest)]$gene_id
df_range$gene_distance[queryHits(closest)] <- mcols(closest)$distance
df_range$gene_strand[queryHits(closest)] <- genelist[subjectHits(closest)]$gene_strand

head(df_range)

# determine the closest tss
tss <- GRanges(seqnames = Rle(gene$chr), 
                     ranges=IRanges(start=gene$tss, width = 1))
tss$gene_id <- gene$gene_id
tss$gene_strand <- gene$strand

closesttss <- distanceToNearest(df_range, subject = tss)
df_range$closest_tss <- NA
df_range$tss_distance <- NA
df_range$tss_strand <- NA
df_range$closest_tss[queryHits(closesttss)] <- tss[subjectHits(closesttss)]$gene_id
df_range$tss_distance[queryHits(closesttss)] <- mcols(closesttss)$distance
df_range$tss_strand[queryHits(closesttss)] <- tss[subjectHits(closesttss)]$gene_strand

head(df_range)
rm(gencode)

### STEP 3: ADD CELL TYPE SPECIFIC ANNOTIONS ##
## 3.1 CHROMHMM MODEL ##
print("determining chromHMM")
tissue <- "E123"
#chromhmm <- fread(paste0(an_path, "chromHMM/E123_15_coreMarks_dense.bed.gz"))
chromhmm <- fread(paste0(an_path, "chromHMM/E123_25_imputed12marks_dense.bed.gz"), skip=1)
names(chromhmm) <- c("chrom", "start", "end", "mark", "score", "strand", "thickStart", "thickEnd", "itemRgb")
chromhmm$chrom <- sub("chr", "", chromhmm$chrom)

hmm_range <- GRanges(seqnames=Rle(chromhmm$chrom),
                     ranges=IRanges(start=chromhmm$start, end = chromhmm$end))
hmm_range$chromHMM_mark <- chromhmm$mark

ov <- findOverlaps(df_range, hmm_range)
ov_raqtls <- df_range[queryHits(ov)]
ov_hmm <- hmm_range[subjectHits(ov)]
ov_raqtls$chromHMM_mark <- ov_hmm$chromHMM_mark

hit <- as.data.table(ov_raqtls)
#hit <- hit %>%
#	mutate(chromHMM_region = ifelse(hit$chromHMM_mark == "1_TssA" | hit$chromHMM_mark == "2_TssAFlnk" | hit$chromHMM_mark == "3_TxFlnk" | hit$chromHMM_mark == "10_TssBiv", "TSS_Promoter", 
#				ifelse(hit$chromHMM_mark == "8_ZNF/Rpts" | hit$chromHMM_mark == "11_BivFlnk" | hit$chromHMM_mark == "13_ReprPC" | hit$chromHMM_mark == "14_ReprPCWk", "Repressor",
#				ifelse(hit$chromHMM_mark == "6_EnhG" | hit$chromHMM_mark == "7_Enh" | hit$chromHMM_mark == "12_EnhBiv", "Enhancer",
#				ifelse(hit$chromHMM_mark == "4_Tx" | hit$chromHMM_mark == "5_TxWk", "Transcribed", "Other")))))

hit <- hit %>%
  mutate(chromHMM_region = ifelse(hit$chromHMM_mark == "1_TssA" | hit$chromHMM_mark == "2_PromU", "TSS_Promoter",
        ifelse(hit$chromHMM_mark == "13_EnhA1" | hit$chromHMM_mark == "14_EnhA2", "Active_Enhancer", "Other")))

cols <- c("SNP_ID", "chromHMM_mark", "chromHMM_region")
hmm <- select(hit, all_of(cols))
hmm <- aggregate(. ~ SNP_ID, hmm, FUN = paste, collapse=",")

head(hmm)
rm(chromhmm)
rm(hit)

## 3.2 DNASE HYPERSENSITIVE SITES ##
print("determining DHS")
# broadpeak
n <- fread(paste0(an_path, "DHS/ENCFF504CYY_DHS_2018_hg19_K562.bed"))
names(n) <- c("chr", "start", "end", "V4", "V5", "V6", "V7", "V8", "V9")
n[, chr := sub("chr", "", chr)]
nr <- GRanges(seqnames=Rle(n$chr), ranges = IRanges(start=n$start, end=n$end))

closest <- distanceToNearest(df_range, subject = nr)
df_range$K562_broadpeak_DHS_dist <- NA
df_range$K562_broadpeak_DHS_dist[queryHits(closest)] <- mcols(closest)$distance

head(df_range)

# narrowpeak
n <- fread(paste0(an_path, "DHS/K562_hg19_2020_DHS.bed"))
names(n) <- c("chr", "start", "end", "V4", "V5", "V6", "V7", "V8", "V9", "V10")
n[, chr := sub("chr", "", chr)]
nr <- GRanges(seqnames=Rle(n$chr), ranges = IRanges(start=n$start, end=n$end))

closest <- distanceToNearest(df_range, subject = nr)
df_range$K562_narrowpeak_DHS_dist <- NA
df_range$K562_narrowpeak_DHS_dist[queryHits(closest)] <- mcols(closest)$distance

head(df_range)

### STEP 4:  CREATE FINAL FILE ###
print("creating final file")
## 4.1 ADD ADDIDIONAL COLUMNS ##
# convert to datatable
dt <- as.data.table(df_range)

# add chromhmm columns
dt <- merge(dt, hmm, by = "SNP_ID", all.x=TRUE) 
dt <- merge(dt, gc_hit, by = "SNP_ID", all.x = TRUE)
dt$GENCODE_location[is.na(dt$GENCODE_location)] <- "intergenic"


# add yes no columns
dt <- dt %>%
  mutate(UCSC_promoter = ifelse(UCSC_prom_dist == 0, "1","0")) %>%
  mutate(Hoffman_enhancer = ifelse(Hoffman_enh_dist == 0, "1","0")) %>%
  mutate(E123_promoter = ifelse(base::grepl("TSS_Promoter",dt$chromHMM_region), "1", "0")) %>%
  mutate(E123_enhancer = ifelse(base::grepl("Active_Enhancer",dt$chromHMM_region), "1", "0")) %>%
  mutate(K562_broadpeak_DHS = ifelse(K562_broadpeak_DHS_dist == 0, "1", "0")) %>%
  mutate(K562_narrowpeak_DHS = ifelse(K562_narrowpeak_DHS_dist == 0, "1", "0"))


head(dt)

## 4.2 FIX COLUMN ORDERS ##
newcols <- c("SNP_ID", "UCSC_prom_dist", "UCSC_promoter", "Hoffman_enh_dist", "Hoffman_enhancer", "GENCODE_location", 
		"closest_gene", "gene_distance", "gene_strand", "closest_tss", "tss_distance", "tss_strand",
		"chromHMM_mark", "chromHMM_region", "E123_promoter", "E123_enhancer", 
		"K562_broadpeak_DHS_dist", "K562_broadpeak_DHS", "K562_narrowpeak_DHS_dist", "K562_narrowpeak_DHS")

dt <- select(dt, all_of(newcols))

## 4.3 PRINT STATISTICS ##
print(paste0("total SNPs used: ", nrow(dt)))

# enhancers
a <- dt[dt$Hoffman_enhancer == "1", ]
b <- dt[dt$Hoffman_enh_dist < 5000, ]
print(paste0("SNPs located in Hoffman enhancers: ", nrow(a)))
print(paste0("SNPs located <5kb to Hoffman enhancers: ", nrow(b)))

# promoters
a <- dt[dt$UCSC_promoter == "1", ]
b <- dt[dt$UCSC_prom_dist < 5000, ]
print(paste0("SNPs located in UCSC promoters: ", nrow(a)))
print(paste0("SNPs located <5kb to Hoffman enhancers: ", nrow(b)))

# gencode
table(dt$GENCODE_location)

# tss
b <- dt[dt$tss_distance < 5000, ]
print(paste0("number of SNPs with a protein coding TSS < 5kb: ", nrow(b)))
print(paste0("mean distance to CLOSEST PROTEIN CODING GENE: ", mean(dt$gene_distance), " | with a range of: ", min(dt$gene_distance), " - ", max(dt$gene_distance), " | and a median of: ", median(dt$gene_distance)))
print(paste0("mean distance to CLOSEST PROTEIN CODING TSS: ", mean(dt$tss_distance), " | with a range of: ", min(dt$tss_distance), " - ", max(dt$tss_distance), " | and a median of: ", median(dt$tss_distance)))

# dhs
b <- dt[dt$K562_broadpeak_DHS == "1", ]
print(paste0("number of SNPs in broadpeak DHS: ", nrow(b)))
b <- dt[dt$K562_narrowpeak_DHS == "1", ]
print(paste0("number of SNPs in narrowpeak DHS: ", nrow(b)))

#chromhmm
table(dt$chromHMM_region)

## 4.4 SAVE FILE ##
final <- merge(df, dt, by = "SNP_ID")
#fwrite(final, outfile, sep="\t", quote=F, row.names=F)
                           
