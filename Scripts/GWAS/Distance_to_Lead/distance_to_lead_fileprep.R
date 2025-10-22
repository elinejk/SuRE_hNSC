library(data.table)
library(dplyr)
library(GenomicRanges)
library(IRanges)


### STEP 1: set paths
args <- commandArgs(trailingOnly=TRUE)
set_path <- args[[1]] # path that leads to the file containing rsIDs and locations
setname <- args[[2]] # indicates which of the sets your are running: emVars, control SNPs, or all SuRE SNPs
info_path <- args[[3]] # this is a file with on each row phenotype specific GWAS paths: phenotype, summary statistics, lead SNPs, and candidate SNPs
out_path <- args[[4]] # where to save the results
dag <- format(Sys.Date(),"%Y%m%d") # date to add to the output file name


### STEP 2: read in the files
# Read in the SuRE file for the set you want to run
dt <- fread(set_path)
cols <- c("SNP_ID", "chrom", "pos.hg19", "ref.seq", "alt.seq", "ref.element.count", "alt.element.count",
	          "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean", "hNPC.wilcoxon.pvalue")
dt <- select(dt, all_of(cols))


# Read in the file with the paths for all GWAS information
infofile <- fread(info_path)



### STEP 3: start loop to determine distance raqtl to closest lead snp
for (row in 1:nrow(infofile)) {
  ## STEP 3.1: set paths for the GWAS files
  set <- dt
  
  pheno <- infofile$pheno[row]
  sumstat_path <- infofile$sumstats[row]
  lead_path <- infofile$leadsnps[row]
  candidate_path <- infofile$candidate[row]
  locus_path <- infofile$locus[row]
  
  
  ## STEP 3.2: is the SNP within a GWAS locus?
  print("testing the overlap with the gwas loci")
  # read the locus file and generate granges
  locus <- fread(locus_path)
  grlocus <- GRanges(seqnames = Rle(locus$CHR), 
                  ranges=IRanges(start=locus$start, end = locus$end)) 
				  
  grsure <- GRanges(seqnames = Rle(set$chrom), ranges = IRanges(start = set$pos.hg19, width = 1))
  grsure$SNP_ID <- set$SNP_ID
  
  # find the overlaps
  ov <- findOverlaps(grlocus, grsure)
  ov_sure <- grsure[subjectHits(ov)]
  inlocus <- as.data.table(ov_sure)
  
  # add to the set file
  set$inlocus <- ifelse(set$SNP_ID %in% inlocus$SNP_ID, 1, 0)
  
  
  ## STEP 3.3: add GWAS information to set file
  print("testing if the snps are part of the gwas summary statistics")
  # read gwas file and change column names
  gwas <- fread(sumstat_path, fill=TRUE)
  names(gwas)[names(gwas) == 'ID' | names(gwas) == 'MarkerName' | names(gwas) == 'variant_id' | names(gwas) == 'RSID' | names(gwas) == 'SNP'] <- "SNP_ID"
  names(gwas)[names(gwas) == 'PVAL' | names(gwas) == 'Pval' | names(gwas) == 'p_value' | names(gwas) == 'P-value' | names(gwas) == 'P'] <- "GWAS_P"
  gwas <- gwas[, c("SNP_ID", "GWAS_P")]
  
  # add information to the set file
  set <- set %>% mutate(ingwas = ifelse(set$SNP_ID %in% gwas$SNP_ID, 1, 0))
  set <- merge(set, gwas, by="SNP_ID", all.x=TRUE)
  
  # add information about candidate SNPs
  can <- fread(candidate_path)
  set <- set %>% mutate(iscandidate = ifelse(set$SNP_ID %in% can$rsID, 1, 0))
  
  
  ## STEP 3.4: determine the distance to the lead snp
  print("determining the distance to the lead snp")
  # read the lead snp file and generate granges 
  leadsnp <- fread(lead_path)
  ld <- GRanges(seqnames = Rle(leadsnp$chr), ranges=IRanges(start=leadsnp$pos, width=1))
  ld$leadsnp <- leadsnp$rsID

  # determine the distance
  dist <- distanceToNearest(grsure, subject = ld)
  grsure$closest_lead <- NA
  grsure$distance_lead <- NA
  grsure$closest_lead[queryHits(dist)] <- ld[subjectHits(dist)]$leadsnp
  grsure$distance_lead[queryHits(dist)] <- mcols(dist)$distance

  # generate a new file where you add the distance to the set file
  distance <- as.data.table(grsure)
  colss <- c("SNP_ID", "closest_lead", "distance_lead")
  distance <- select(distance, all_of(colss))
  snpdist <- merge(set, distance, by = "SNP_ID") 
  
  # save this file
  fwrite(snpdist, paste0(out_path, dag, "_", pheno, "_" , setname, "_distance-snp-to-lead.txt"), quote=F, row.names=F, sep="\t")
  
  # also generate a file with only the raqtls significant in the gwas
  sign <- snpdist[snpdist$ingwas == 1, ]
  sign <- sign[sign$GWAS_P <= 5e-8, ]
  
  fwrite(sign, paste0(out_path, dag, "_", pheno, "_" , setname, "_gwassign_distance-snp-to-lead.txt"), quote=F, row.names=F, sep="\t")  
  

}


### STEP 4: combine the different phenotype files
phenos <- c("ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ")

a <- data.table()

## STEP 4.1: run the loop
for (pheno in phenos){
	# load and prepare the SNP file
	df <- fread(paste0(out_path, dag, "_", pheno, "_" , setname, "_gwassign_distance-snp-to-lead.txt"))
	df$phenotype <- pheno
	
	# combine files
	a <- rbind(a, df)
	fwrite(a, paste0(out_path, dag, "_allphenos_", setname, "_gwassign_distance-snp-to-lead.txt"), quote = F, row.names=F, sep='\t')
	
	# sort on rsID and shortest distance, then retain the first row
	b <- a %>%
	arrange(SNP_ID, distance_lead) %>%
	distinct(SNP_ID, .keep_all = T)
	
	fwrite(b, paste0(out_path, dag, "_allphenos-shortestdist_", setname, "_gwassign_distance-snp-to-lead.txt"), quote = F, row.names=F, sep='\t')

}
