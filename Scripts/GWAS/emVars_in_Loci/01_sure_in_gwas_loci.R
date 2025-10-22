library(data.table)
library(dplyr)
library(GenomicRanges)
library(IRanges)

args <- commandArgs(TRUE)

# paths
path <- "/gpfs/home6/ekoornstra/"
fuma <- args[1]
sumstats <- args[2]
snp_file <- args[3]
outfile <- args[4]
outname <- args[5]
sure <- "sure_data/freeze7/pvalue.freeze7.snp-permutation.04042024.txt.gz"
hnsc <- "raqtls/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt"

### STEP 1: PREPARE FILES ###
## STEP 1.1: prepare locus files
loci <- fread(paste0(path, fuma))
loci$GenomicLocus <- paste0("locus", loci$GenomicLocus)

grlocus <- GRanges(seqnames = Rle(loci$CHR), 
                  ranges=IRanges(start=loci$start, end = loci$end)) 
grlocus$locus <- loci$GenomicLocus
grlocus$nSNPs <- loci$nSNPs
grlocus$LeadSNPs <- loci$LeadSNPs

## STEP 1.2: prepare sure file
cols <- c("SNP_ID", "chrom", "pos.hg19")
df <- fread(paste0(path, sure))
df <- select(df, all_of(cols))

grsure <- GRanges(seqnames = Rle(df$chrom), ranges = IRanges(start = df$pos.hg19, width = 1))
grsure$SNP_ID <- df$SNP_ID

rm(df)

## STEP 1.3: prepare hNSC emVars and GWAS significant snps
raqtl <- fread(paste0(path, hnsc))
raqtl <- raqtl[, "SNP_ID"]

snps <- fread(paste0(path, snp_file))
snps <- snps[, "rsID"]


### STEP 2: START WITH THE OVERLAPS ###
## STEP 2.1: SuRE SNPs within GWAS loci based on the coordinates of the GWAS locus, not the actual GWAS SNPs tested
# 2.1.1: SuRE SNPs in the FUMA locus file
ov <- findOverlaps(grlocus, grsure)
ov_locus <- grlocus[queryHits(ov)]
ov_sure <- grsure[subjectHits(ov)]
ov_locus$SNP_ID <- ov_sure$SNP_ID

inlocus <- as.data.table(ov_locus)
names(inlocus)[names(inlocus) == 'seqnames'] <- "chr"
inlocus <- unique(inlocus)

# 2.1.2: load GWAS sumstats
gwas <- fread(paste0(path, sumstats), fill=TRUE)
names(gwas)[names(gwas) == 'ID' | names(gwas) == 'MarkerName' | names(gwas) == 'variant_id' | names(gwas) == 'RSID'] <- "SNP"
names(gwas)[names(gwas) == 'chromosome' | names(gwas) == 'Chromosome' | names(gwas) == 'CHROM'] <- "CHR"
names(gwas)[names(gwas) == 'BP' | names(gwas) == 'base_pair_location' | names(gwas) == 'Position'] <- "POS"
names(gwas)[names(gwas) == 'PVAL' | names(gwas) == 'Pval' | names(gwas) == 'p_value' | names(gwas) == 'P-value'] <- "P"

gwas <- gwas[, c("SNP", "CHR", "POS", "P")]

gwas <- gwas[!is.na(gwas$POS)]

grgwas <- GRanges(seqnames = Rle(gwas$CHR), ranges = IRanges(start=gwas$POS, width = 1))
grgwas$SNP <- gwas$SNP
grgwas$P <- gwas$P


# 2.1.3: get some statistics
# sure snps in locus that were also tested in the gwas
inlocus <- inlocus %>% mutate(GWAS_SNP = ifelse(inlocus$SNP_ID %in% gwas$SNP, "1", "0"))
inlocus$GWAS_SNP <- as.numeric(inlocus$GWAS_SNP)

# sure snps in locus that are emVars
inlocus <- inlocus %>% mutate(hNSC_emVar = ifelse(inlocus$SNP_ID %in% raqtl$SNP_ID, "1", "0"))
inlocus$hNSC_emVar <- as.numeric(inlocus$hNSC_emVar)

# emVars in the locus that were part of the gwas sumstats (sure in locus vs emVars vs gwas snps)
inlocus <- inlocus %>% mutate(hNSC_emVar_GWAS = ifelse(inlocus$SNP_ID %in% raqtl$SNP_ID & inlocus$SNP_ID %in% gwas$SNP, "1", "0"))
inlocus$hNSC_emVar_GWAS <- as.numeric(inlocus$hNSC_emVar_GWAS)

# sure snps in the locus are "in LD with lead SNPs/indsig" gwas locus SNPs (nSNPs column), named here as candidate
inlocus <- inlocus %>% mutate(SuRE_GWAScandidate = ifelse(inlocus$SNP_ID %in% snps$rsID, "1", "0"))
inlocus$SuRE_GWAScandidate <- as.numeric(inlocus$SuRE_GWAScandidate)

# emVars are also "in LD with lead SNPs/indsig" locus snps
inlocus <- inlocus %>% mutate(hNSCemVar_GWAScandidate = ifelse(inlocus$SNP_ID %in% snps$rsID & inlocus$SNP_ID %in% raqtl$SNP_ID, "1", "0"))
inlocus$hNSCemVar_GWAScandidate <- as.numeric(inlocus$hNSCemVar_GWAScandidate)

# 2.1.4: add GWAS pvalues
a <- gwas[, c("SNP", "P")]
names(a) <- c("SNP_ID", "GWAS_P")
inlocus <- merge(inlocus, a, by = "SNP_ID", all.x=TRUE)
inlocus <- inlocus %>% mutate(GWAS_sign = ifelse(inlocus$GWAS_P <= 5*10^-8, "1", "0"))
inlocus$GWAS_sign[is.na(inlocus$GWAS_sign)] <- 0
inlocus$GWAS_sign <- as.numeric(inlocus$GWAS_sign)


# 2.1.5: save
fwrite(inlocus, paste0(path, outfile), quote=F, row.names=F, sep='\t')


## STEP 2.2: determine number of GWAS SNPs within each locus: irrespective of significance and LD for final stats
# 2.2.1: the overlap
ol <- findOverlaps(grlocus, grgwas)
ol_locus <- grlocus[queryHits(ol)]
ol_gwas <- grgwas[subjectHits(ol)]
ol_locus$SNP <- ol_gwas$SNP
ol_locus$P <- ol_gwas$P

ollocus <- as.data.table(ol_locus)
ollocus <- ollocus %>% count(locus)
names(ollocus) <- c("locus", "nGWAS")

# 2.2.2: merge this information for the final file
final_file <- merge(loci, ollocus, by.x="GenomicLocus", by.y="locus", all.x=TRUE)


### STEP 3: START GENERATION OF THE FINAL SUMMARY FILE ###
## STEP 3.1: first get a summary of the total SuRE SNPs in each locus
a <- inlocus %>% count(locus)
names(a) <- c("locus",  "nSuRE")

loc <- inlocus[ , -c("SNP_ID", "chr", "start", "end", "width", "strand", "nSNPs", "LeadSNPs")]

fullfile <- merge(final_file, a, by.x = "GenomicLocus", by.y = "locus", all.x=TRUE)

## STEP 3.2: then get the statistics for the other overlaps
b <- loc %>% 
  group_by(loc$locus) %>% 
  summarize(nBoth = sum(GWAS_SNP), nhNSC_emVar = sum(hNSC_emVar), nhNSC_emVar_GWAS = sum(hNSC_emVar_GWAS),
            nSuRE_nSNPs = sum(SuRE_GWAScandidate), nhNSC_emVar_nSNPs = sum(hNSCemVar_GWAScandidate), nGWASsign_SuRE = sum(GWAS_sign))
names(b) <- c("GenomicLocus", "nBoth", "nhNSC_emVar", "nhNSC_emVar_GWAS", "nSuRE_nSNPs", "nhNSC_emVar_nSNPs", "nGWASsign_SuRE")

fullfile <- merge(fullfile, b, by = "GenomicLocus", all.x=TRUE)


## STEP 3.3: get a better readable files
subfile <- fullfile[,c("GenomicLocus", "CHR", "start", "end", "LeadSNPs", "nGWAS", "nSNPs", 
                       "nSuRE", "nBoth", "nSuRE_nSNPs", "nGWASsign_SuRE", "nhNSC_emVar", "nhNSC_emVar_GWAS", "nhNSC_emVar_nSNPs")]
subfile <- unique(subfile)


## STEP 3.4: save file
fwrite(subfile, paste0(path, outname), quote=F, row.names=F, sep='\t')



