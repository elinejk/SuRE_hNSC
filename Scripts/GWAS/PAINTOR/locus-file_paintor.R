library(data.table)
library(GenomicRanges)
library(rtracklayer)

## PATHS ##
args <- commandArgs(TRUE)

file_path <- args[1]
loci_file <- args[2]
gwas_file <- args[3]

## PREPARE LOCUS FILE ##
coord <- fread(loci_file)
coord_gr <- GRanges(seqnames=Rle(coord$CHR), ranges=IRanges(start=coord$start, end=coord$end))

## PREPARE GWAS SUMMARY STATISTICS ##
gwas <- fread(gwas_file, fill = TRUE)

#standardize the certain headers
names(gwas)[names(gwas) == 'ID' | names(gwas) == 'MarkerName' | names(gwas) == 'variant_id' | names(gwas) == 'RSID'] <- "SNP"
names(gwas)[names(gwas) == 'chromosome' | names(gwas) == 'Chromosome' | names(gwas) == 'CHROM'] <- "CHR"
names(gwas)[names(gwas) == 'BP' | names(gwas) == 'base_pair_location' | names(gwas) == 'Position'] <- "POS"
names(gwas)[names(gwas) == 'PVAL' | names(gwas) == 'Pval' | names(gwas) == 'p_value' | names(gwas) == 'P-value'] <- "P"

# if necessary add a Z score | check this manually
# z score can be calculated as effect size/SE or beta/SE or log(OR)/SE
head(gwas)
#gwas$Beta <- log(gwas$OR)
gwas$Zscore = gwas$BETA / gwas$SE

# convert to granges | this needs to be changed for each phenotypes because the summary statistics headers are not the same!!
gwas <- gwas[!is.na(gwas$POS), ]

# make sure the chromosome is numbers only, and then create the granges object
gwas[, CHR := sub("chr", "", CHR)]

gwas_gr <- GRanges(seqnames = Rle(gwas$CHR),
  ranges = IRanges(start=gwas$POS, width=1),
  zscore = gwas$Zscore,
  rsid = gwas$SNP,
  A1 = gwas$A1,
  A2 = gwas$A2
  )

## SUBSET THE GWAS FILE TO SEE IF COORDINATES ARE WITHIN LOCI ##
for (i in seq_along(coord_gr)){
  ov_gr <- subsetByOverlaps(gwas_gr, coord_gr[i])
  ov <- as.data.table(ov_gr)
  ov <- ov[, !c("end","width","strand")]
  names(ov) <- c("CHR", "pos","Zscore","rsid","A1","A2")
  ov <- setcolorder(ov, c("CHR", "pos","rsid","A1","A2","Zscore"))
  print(paste0("creating locus", i))
  fwrite(ov, paste0(file_path, "files/gwas-loci/locus", i), sep="\t",row.names=FALSE)
}


