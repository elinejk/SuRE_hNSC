#Packages
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(xlsx)



#1. Define TF names
tf_names <- c("ZBTB33", "YY1", "TP63", "TP53", "GABPA") 

#2. Define paths
file_paths <- lapply(tf_names, function(tf){
paste0("/gpfs/home6/evzanten/jaspar_enrichment/020524_twosets_raqtl_enr_100bp_", tf, "/allEnrichments.tsv") 
} 
  ) %>% unlist()
  
#3. Read in files
ra_data <- lapply(file_paths, function(read_paths){
fread(read_paths)
}
  ) 

#4. Add name of TFs to each element in the list
TF <- c("ZBTB33", "YY1", "TP63", "TP53", "GABPA")
ra_data <- map2(ra_data, TF, ~cbind(.x, TF=.y))

#5. Make a df of the listed data
ra_df <- rbindlist(ra_data)

#6. Remove Inf values. 
ra_df <- ra_df[is.finite(ra_df$oddsRatio),]

#7. Calculate mean OR per TF / neighboring TF class combo. 
ra_class <- ra_df %>% group_by(TF, class) %>% summarise(means=mean(oddsRatio)) %>%
as.data.frame()


#8. Convert to wide format to have the TFs as columns
ra_wide_class <- ra_class %>%
           pivot_wider(names_from=TF, values_from = means) %>%
           as.data.frame()

           
#9. Save as Excel file
write.xlsx(ra_wide_class, "/gpfs/home6/evzanten/jaspar_enrichment/tables/101024_TFBS_classes_neighboring.xlsx", row.names=F)




