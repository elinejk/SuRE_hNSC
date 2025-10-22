## Packages
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# Load SNP2TFBS data
snp2tfbs <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/resources/SNP2TFBS/snp2tfbs_JASPAR_CORE_2014_vert_RESTRUCTURED.txt.gz")

# Keep only SNPs with rs or ss identifiers
snp2tfbs <- snp2tfbs[grepl("rs|ss", snp2tfbs$SNP_ID),]

# Keep only the first TFBS if wanted
firsttfbs <- TRUE

if(firsttfbs == TRUE) {
  basesave <- "_results_firstTFBS_" # part of name of the results file
  snp2tfbs <- snp2tfbs %>%
    mutate(row_id = row_number(),
           abs_scoreDiff = abs(scoreDiff)) %>%
    group_by(SNP_ID) %>%
    slice_max(order_by = abs_scoreDiff, with_ties = TRUE) %>%
    slice(1) %>%   # among ties, keep the first
    ungroup()
} else {
  basesave <- "_results_"
}

# Define cell types
cell_types <- c("k562", "hepg2", "hnsc")

# Cell type names within sure data columns, make aliases
cell_type_data <- list(
    k562 = "K562",
    hepg2 = "HepG2",
    hnsc = "hNPC"
)

# Function to calculate TFBS enrichment & concordance
calculate_tfbs_enrichment <- function(cell_type) {
    
    # Define file paths
    file_path_raqtl <- paste0("/gpfs/work4/0/AdamsLab/Projects/sure/raqtls/results/freeze7/", 
                              cell_type, "_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")
    
    file_path_control <- paste0("/gpfs/work4/0/AdamsLab/Projects/sure/raqtls/results/freeze7/", 
                                cell_type, "_no_downsampling_snp-permutation_freeze7_controls_04042024.txt")

    # Read data
    raqtls <- fread(file_path_raqtl)
    controls <- fread(file_path_control)

    # Merge with snp2tfbs
    raqtls_merged <- merge(raqtls, snp2tfbs, by = "SNP_ID")
    raqtls_merged$raqtl <- TRUE

    controls_merged <- merge(controls, snp2tfbs, by = "SNP_ID")
    controls_merged$raqtl <- FALSE

    # Combine emVar and control SNPs
    all_snp <- rbind(raqtls_merged, controls_merged)

    # Get unique TFBS names
    u_tf <- unique(na.omit(snp2tfbs$TFBS))


    #1. Make new df indicating whether the SNPs are emVars or not
    dd <- all_snp[,c("SNP_ID", "raqtl")] 
    dd <- as.data.frame(dd)
    
    #2. Now for each SNP, we determine which TF it binds (0 or 1)
    dd[,u_tf] <- lapply(u_tf, function(n){
      as.logical(all_snp[,TFBS] == n, na.rm=T)
      })
      
    #3. Now for the emVars and controls we get a count
    raqtl_n <- table(dd$raqtl)
    
    fam_test <- lapply(u_tf, function(n) {
      x <- aggregate(dd[, n], by = list(dd$raqtl), sum)$x
      expected_value <- sum(x) / sum(raqtl_n) * sum(raqtl_n)
      correct <- if (expected_value < 10) TRUE else FALSE
      prop.test(x, raqtl_n, correct = correct)
    }) 
    
    fam_p <- sapply(fam_test, function(x) x$p.value)
    fam_enriched <- sapply(fam_test, function(x) which.max(x$estimate) == 2)
    
    fam_results <- data.frame(tf_name = u_tf,
                              p = fam_p,
                              p_fdr = p.adjust(fam_p, "BH"),
                              enriched = fam_enriched)
                              
      
    ###MAKE SUMSTATS###
    summary_df <- dd %>%
   # Gather all TF columns into key-value pairs (long format)
   tidyr::pivot_longer(cols = 3:ncol(.), names_to = "TF_name", values_to = "TF_value") %>%
   # Group by TF_name and calculate the counts
   group_by(TF_name) %>%
   summarise(
     MPRA_hits = sum(TF_value & raqtl==TRUE),
     MPRA_hits_prop = MPRA_hits / sum(raqtl==TRUE),
     MPRA_nonhits = sum(TF_value & raqtl==FALSE),
     MPRA_nonhits_prop = MPRA_nonhits / sum(raqtl==FALSE)
       )

    ### **Concordance Analysis**
    add_concordance <- function(data_frame, cell_type) {
        # Use the cell_type_data alias for concordance column names
        cell_type_alias <- cell_type_data[[cell_type]]  # Get the alias
        
        data_frame$concordance <- ifelse(
            data_frame[[paste0(cell_type_alias, ".cDNA.ref.mean")]] < data_frame[[paste0(cell_type_alias, ".cDNA.alt.mean")]] & data_frame$scoreDiff > 0 |
            data_frame[[paste0(cell_type_alias, ".cDNA.ref.mean")]] < data_frame[[paste0(cell_type_alias, ".cDNA.alt.mean")]] & data_frame$PWMref == "." & !data_frame$PWMalt %in% c("-", "+") |
            data_frame[[paste0(cell_type_alias, ".cDNA.ref.mean")]] > data_frame[[paste0(cell_type_alias, ".cDNA.alt.mean")]] & data_frame$scoreDiff < 0 |
            data_frame[[paste0(cell_type_alias, ".cDNA.ref.mean")]] > data_frame[[paste0(cell_type_alias, ".cDNA.alt.mean")]] & data_frame$PWMalt == "." & !data_frame$PWMref %in% c("-", "+"),
            "concordant",
            ifelse(data_frame$PWMref == data_frame$PWMalt, "equal", "discordant")
        )
        return(data_frame)
    }

    # Apply concordance function
    df_concordance <- add_concordance(all_snp, cell_type)
    
    df_new <- df_concordance %>% 
        group_by(TFBS, raqtl) %>% 
        summarise(Percentage_Concordant = sum(concordance == "concordant") / n() * 100, .groups="drop") %>%
        spread(key = raqtl, value = Percentage_Concordant) %>%
        rename_with(~ gsub("TRUE", "percentage_concordant_raqtl", .x)) %>%
        rename_with(~ gsub("FALSE", "percentage_concordant_control", .x))


    # Merge enrichment and statistics
    summary_df <- merge(fam_results, summary_df, by.x = "tf_name", by.y = "TF_name")
    summary_df <- summary_df[order(summary_df$p_fdr),]

    # Merge final statistics
    all_sumstats <- merge(df_new, summary_df, by.x = "TFBS", by.y = "tf_name")

    # Save results
    opath <- paste0("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/transcription_factors/TF_enrichment/", 
                    toupper(cell_type), "_SNP-level", basesave, format(Sys.Date(), "%d%m%y"), ".txt")
    fwrite(df_concordance, opath, sep = "\t")
    
    output_path <- paste0("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/transcription_factors/TF_enrichment/", 
                          toupper(cell_type), basesave, format(Sys.Date(), "%d%m%y"), ".txt")
    fwrite(all_sumstats, output_path, sep = "\t")
}


# Apply Function to Each Cell Type
results_list <- lapply(cell_types, calculate_tfbs_enrichment)


