library(qvalue)
library(ggplot2)
library(reshape2)

# set up the files
cells <- c("hepg2", "hnsc", "k562") # for the file names
cell_type <- c("HepG2", "hNPC", "K562") # for the column names

data_list <- lapply(paste0("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/emVars/freeze7/", 
                    cells, "_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt"), read.delim, header=TRUE)

names(data_list) <- cell_type

active <- read.delim("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/emVars/freeze7/freeze7_cell_line_shared_active_snps.txt")

# function to compute storeys pi1
compute_pi1 <- function(pval) {
  qobj <- qvalue(p = pval)
  pi1 <- 1 - qobj$pi0
  return(pi1)
}

# save results in a matrix
results <- matrix(NA, nrow=3, ncol=3, dimnames=list(cell_type, cell_type))

# calculate the pi1 for all combinations of cell types
for (i in seq_along(cell_type)) {
  for (j in seq_along(cell_type)) {
    if (i != j) {
      cell_A <- cell_type[i]
      cell_B <- cell_type[j]
      
      df_A <- data_list[[cell_A]]
      col_name_B <- paste0(cell_B, ".wilcoxon.pvalue")
      
      #df_A <- df_A[df_A$SNP_ID %in% active$SNP_ID, ] # optional, filter on active snps
      
      if (col_name_B %in% names(df_A)) {
        pvals_B_in_A <- df_A[[col_name_B]]
        pi1 <- compute_pi1(pvals_B_in_A)
        results[cell_A, cell_B] <- round(pi1, 4)
      } else {
        warning(paste("Missing", col_name_B, "in", cell_A, "file"))
      }
    }
  }
}

# results
# interpretation: 0 (no replication), 0.2-0.4 (weak/moderate replication), 
# 0.5-0.7 (strong evidence of shared signal), 0.8-1 (very strong replication = same underlying biology)
# 0.7 means that about 70% of pvalues are true signals
print(results)


# plot the results if wanted
results_long <- melt(results, varnames = c("Signif in", "Evaluated by"), value.name = "pi1")

ggplot(results_long, aes(x = `Signif in`, y = `Evaluated by`, fill = pi1)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", pi1)), na.rm = TRUE) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey80") +
  theme_minimal() +
  labs(title = "Storey's pi1 Heatmap",
       x = "Significant in",
       y = "Evaluated in") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

