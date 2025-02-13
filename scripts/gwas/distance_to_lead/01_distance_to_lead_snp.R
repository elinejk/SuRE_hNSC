library(data.table)
library(dplyr)
library(GenomicRanges)
library(IRanges)

### STEP 1: set paths
args <- commandArgs(trailingOnly=TRUE)
set_path <- args[[1]] # path that leads to the file containing the rsids and locations: raqtls, controls
setname <- args[[2]]
info_path <- args[[3]] # this is a file with on each row phenotype specific GWAS paths: pheno, sumstats, leadsnps
out_path <- args[[4]]

### STEP 2: read files
infofile <- fread(info_path)

dt <- fread(set_path)
cols <- c("SNP_ID", "chrom", "pos.hg19", "ref.seq", "alt.seq", "ref.element.count", "alt.element.count",
	          "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean", "hNPC.wilcoxon.pvalue")
dt <- select(dt, all_of(cols))

### STEP 3: start loop to determine distance raqtl to closest lead snp
for (row in 1:nrow(infofile)) {
  ## STEP 3.1: set paths
  set <- dt
  
  pheno <- infofile$pheno[row]
  sumstat_path <- infofile$sumstats[row]
  lead_path <- infofile$leadsnps[row]
  locus_path <- infofile$locus[row]
  
  ## STEP 3.2: overlap with locus coordinates
  print("testing the overlap with the gwas loci")
  # read the files and generate granges
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
  
  
  ## STEP 3.3: add gwas information to set file
  print("testing if the snps are part of the gwas summary statistics")
  # read gwas file and change column names
  gwas <- fread(sumstat_path, fill=TRUE)
  names(gwas)[names(gwas) == 'ID' | names(gwas) == 'MarkerName' | names(gwas) == 'variant_id' | names(gwas) == 'RSID' | names(gwas) == 'SNP'] <- "SNP_ID"
  names(gwas)[names(gwas) == 'PVAL' | names(gwas) == 'Pval' | names(gwas) == 'p_value' | names(gwas) == 'P-value' | names(gwas) == 'P'] <- "GWAS_P"
  gwas <- gwas[, c("SNP_ID", "GWAS_P")]
  
  # add information to the set file
  set <- set %>% mutate(ingwas = ifelse(set$SNP_ID %in% gwas$SNP_ID, 1, 0))
  set <- merge(set, gwas, by="SNP_ID", all.x=TRUE)
  
  
  ## STEP 3.4: determine the distance to the lead snp
  print("determining the distance to the lead snp")
  # read the lead snp file and generate granges 
  leadsnp <- fread(lead_path)
  ld <- GRanges(seqnames = Rle(leadsnp$chr), ranges=IRanges(start=leadsnp$pos, width=1))
  ld$leadsnp <- leadsnp$rsID

 
  # gdetermine the distance
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
  fwrite(snpdist, paste0(out_path, pheno, "_" , setname, "_distance-snp-to-lead_16042024.txt"), quote=F, row.names=F, sep="\t")
  
  # also generate a file with only the raqtls significant in the gwas
  sign <- snpdist[snpdist$ingwas == 1, ]
  sign <- sign[sign$GWAS_P <= 5e-8, ]
  
  fwrite(sign, paste0(out_path, pheno, "_" , setname, "_gwassign_distance-snp-to-lead_16042024.txt"), quote=F, row.names=F, sep="\t")  
  
  ## STEP 3.5: now the distance from the lead snp to the snps
  print("determining the distance from the lead snp to the gwas-tested snp with the lowest sure pvalue < 100kb")
  # prepare the lead snp file
  leadsnp$start_window <- leadsnp$pos - 100000
  leadsnp$end_window <- leadsnp$pos + 100000
  leadsnp$chr <- as.character(leadsnp$chr)
  
  # subset the file file to only retain snps tested by the gwas and prepare file for overlap
  subs <- set[set$ingwas == 1, ]
  subs <- set[set$GWAS_P <= 5e-8, ]
  subs$end <- subs$pos.hg19
  
  # determine if the snp falls within the 100kb window
  setkey(leadsnp, chr, start_window, end_window)
  wind <- foverlaps(subs, leadsnp, by.x = c("chrom", "pos.hg19", "end"), type = "any", nomatch=NULL)
  
  # only retain the match with the lowest sure pvalue
  wind <- wind %>%
	 arrange(rsID, hNPC.wilcoxon.pvalue) %>%
	 distinct(rsID, .keep_all = T)
   
   # calculate the distance
   wind$distance <- abs(wind$pos.hg19 - wind$pos)
   
   # save the file
   print(paste0("lead snps with a GWAS-tested SuRE SNP within 100kb: ", nrow(wind)))
   fwrite(wind, paste0(out_path, pheno,  "_" , setname, "_distance-lead-to-snp_16042024.txt"), quote=F, row.names=F, sep="\t")
}




