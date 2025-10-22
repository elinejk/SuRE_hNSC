library(ggplot2)
library(ggrepel) 
library(dplyr)
library(data.table)


###COLOR CODES
##emVars: #E69F00
##Controls: #56B4E9
##All: #009E73

##hNSC: #E69F00
##HepG2: #0072B2
##K562: #D55E00	

# Define cell types and paths
cell_types <- c("K562", "HepG2", "hNPC") 
file_paths <- lapply(cell_types, function(cell_type){
  paste0("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/transcription_factors/TF_enrichment/",cell_type,"_results/", cell_type, "_results_firstTFBS_EK100625.txt") 
}
) %>% unlist() 

# Read in files 
df_enrichment <- lapply(file_paths, function(read_paths){
  fread(read_paths)
}
)

# Add name of each cell type to list elements
names(df_enrichment) <- cell_types

# Create new column for each cell type denoting enrichment (part of old script but not really relevant here)
df2_enrichment <- lapply(df_enrichment, function(df){
  df %>% 
    mutate(
      enrichment_status = factor(
      ifelse(enriched == TRUE & -log10(p_fdr) >= 1.30103, "Enriched in emVars",
             ifelse(enriched == FALSE & -log10(p_fdr) >= 1.30103, "Enriched in control SNPs", "Not enriched")),
      levels = c("Enriched in emVars", "Enriched in control SNPs", "Not enriched")))
})



#Cell line Barplot 

#Define color codes for each cell type (colors below are for emvars)
cell_colors <- c("hNPC" = "#E69F00", "HepG2" = "#0072B2", "K562" = "#D55E00")

#Function to calculate enrichment of each cell type over the others (BARPLOT) 
calculate_relative_enrichment <- function(df_list, cell1, cell2) {
  df1 <- df_list[[cell1]]
  df2 <- df_list[[cell2]]
  
  #Merge datasets by TFBS
  df_merged <- merge(df1, df2, by = "TFBS", suffixes = c(".x", ".y"))

  #Add log2 scores
  df_enriched <- df_merged %>%
    #mutate(log2_bar = log2(MPRA_hits_prop.x / MPRA_hits_prop.y)) #emVar
    mutate(log2_bar = log2(MPRA_nonhits_prop.x / MPRA_nonhits_prop.y)) #control SNPs
    #mutate(log2_bar = log2(MPRA_totalhits_prop.x / MPRA_totalhits_prop.y)) #emvar+control(active)

  
  # remove inf and NaN
  df_enriched <- df_enriched[!is.na(df_enriched$log2_bar) & !df_enriched$log2_bar == "-Inf" & 
                               !df_enriched$log2_bar == "Inf", ]
  
  # get the top highest and top lowest enriched
  n <- 10 
  
  df_enriched <- df_enriched %>% 
   arrange(df_enriched$log2_bar) %>% 
   mutate(rnum = row_number()) %>% 
   mutate(bottom_n = ifelse(rnum %in% seq(1, n), 1, 0)) %>% 
   mutate(top_n = ifelse(rnum %in% seq( n()-n+1, n()), 1, 0)) %>% 
   select(-rnum)
  
  df_set <- df_enriched[df_enriched$bottom_n == "1" | df_enriched$top_n == "1", ]
  
  #Combine top TFBS from both cell types
  df_plot <- df_set %>%
    mutate(color_bar = ifelse(log2_bar > 0, cell_colors[[cell1]], cell_colors[[cell2]]))
  
  return(df_plot)
}

#Function to generate plot 
generate_barplot <- function(df_plot, cell1, cell2) {
  eplot <- ggplot(df_plot, aes(x = reorder(TFBS, log2_bar), y = log2_bar)) +
    geom_bar(stat = "identity", fill = df_plot$color_bar, alpha = 1) +
    scale_color_identity() +
    labs(x = "", y = paste0("Enrichment of motifs (control SNPs)\nlog2(", cell1, "/", cell2, ")")) + # change control to set of interest
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
}

#Apply function for all possible comparisons
comparisons <- list(
  c("K562", "HepG2"),
  c("hNPC", "K562"),
  c("hNPC", "HepG2")
)

#Loop through comparisons, calculate enrichment, and generate plots
plot_list <- lapply(comparisons, function(pair) {
  cell1 <- pair[1]
  cell2 <- pair[2]
  
  # Calculate enrichment
  enrichment_df <- calculate_relative_enrichment(df_enrichment, cell1, cell2)
  
  # Generate and return plot
  generate_barplot(enrichment_df, cell1, cell2)
  ggsave(paste0("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/transcription_factors/paper_plots/barplot_rel_enrichment_celltypes/EK20250610_barplot_firstTFBS_control_", cell1, "_vs_", cell2, ".pdf"), width=5, height=4)
})





