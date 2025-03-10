library(data.table)
library(dplyr)

### paths
path <- "/gpfs/home6/ekoornstra/surexgwas/freeze7/distance-to-lead/"

phenos <- c("ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ")
set <- "controls"

a <- data.table()

### make a loop for our method
for (pheno in phenos){
	## load an prepare the SNP file
	df <- fread(paste0(path, set, "/", pheno, "_", set, "_gwassign_distance-snp-to-lead_16042024.txt"))
	df$Phenotype <- pheno
	
	# combine files
	a <- rbind(a, df)
	fwrite(a, paste0(path, set, "/allphenos_", set, "_gwassign_distance-snp-to-lead_16042024.txt"), quote = F, row.names=F, sep='\t')
	
	# sort on rsid and shortest distance, then retain the first row
	b <- a %>%
	arrange(SNP_ID, distance_lead) %>%
	distinct(SNP_ID, .keep_all = T)
	
	fwrite(b, paste0(path, set, "/allphenos-shortestdist_", set, "_gwassign_distance-snp-to-lead_16042024.txt"), quote = F, row.names=F, sep='\t')

}
