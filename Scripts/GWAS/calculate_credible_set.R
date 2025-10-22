library(data.table)

# Set paths
path <- "/gpfs/home6/ekoornstra/surexgwas/freeze7/finemapping/SCZ/files/results/"
output_path <- "/gpfs/home6/ekoornstra/surexgwas/freeze7/finemapping/SCZ/files/"
credible_percentage <- 95
annot <- "ldsc"
running_date <- "23072025"

# Get list of files 
pattern <- "^locus[0-9]+\\.results$"
files <- list.files(path = path, pattern = pattern, full.names = TRUE)

# Set necessary information
credible_set <- data.table()
threshold <- credible_percentage / 100


# Loop through each file
for (file in files) {
  # Read the file
  dt <- fread(file)
  
  # Confirm that it countains the column
  if (!"Posterior_Prob" %in% colnames(dt)) {
    warning(paste("Missing 'Posterior_Prob' column in file:", file))
    next
  }

  # Order by decreasing posterior probability
  setorder(dt, -Posterior_Prob)

  # Calculate cumulative sum of posterior probabilities
  dt[, cum_prob := cumsum(Posterior_Prob)]
  total_post_prob <- sum(dt$Posterior_Prob)

  # Select SNPs meeting the threshold
  cs <- dt[cum_prob <= threshold * total_post_prob]

  # Add the name of the locus, aka the filename
  cs[, source_file := basename(file)]

  # Combine all loci into one file
  credible_set <- rbind(credible_set, cs, fill = TRUE)
}

head(credible_set)


# Annotate with emVars
emvars <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/annotations/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_ANNOTATED_TFBS_120325.txt")

credible_set <- merge(credible_set, emvars, by.x = "rsid", by.y = "SNP_ID", all.x=TRUE)

# Save results
fwrite(credible_set, paste0(output_path, "credible_set_", credible_percentage, "percentage_", annot, "_", running_date, ".txt"), quote=F, row.names=F, sep='\t') 
