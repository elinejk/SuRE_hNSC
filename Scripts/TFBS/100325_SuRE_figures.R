library(ggplot2)
library(ggrepel) 
library(dplyr)
library(data.table)


###COLOR CODES
##raQTLs: #E69F00
##Controls: #56B4E9
##All: #009E73

##hNSC: #E69F00
##HepG2: #0072B2
##K562: #D55E00	

# Define cell types and paths
cell_types <- c("K562", "HepG2", "hNPC") 
file_paths <- lapply(cell_types, function(cell_type){
              paste0("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/TF_enrichment/",cell_type,"_results_060325.txt")
              }
              ) %>% unlist() 
              
# Read in files 
df_enrichment <- lapply(file_paths, function(read_paths){
                 fread(read_paths)
                 }
                 )
                 
# Add name of each cell type to list elements
names(df_enrichment) <- cell_types

# Create new column for each cell type denoting enrichment
df2_enrichment <- lapply(df_enrichment, function(df){
                  df %>% 
                  mutate(enrichment_status = factor(
                    ifelse(enriched == TRUE & -log10(p_fdr) >= 1.30103, "Enriched in raQTLs",
                      ifelse(enriched == FALSE & -log10(p_fdr) >= 1.30103, "Enriched in control SNPs", "Not enriched")),
                  levels = c("Enriched in raQTLs", "Enriched in control SNPs", "Not enriched")
    ))
})
                  







####PLOTS####


#TFBS enrichment plot 
volcano_plot <- lapply(cell_types, function(volcano_celltype){ 

df <- df2_enrichment[[volcano_celltype]]

ggplot(data = df , aes(x = log2(MPRA_hits_prop / MPRA_nonhits_prop), y = -log10(p_fdr), color = enrichment_status)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Enriched in raQTLs" = "#E69F00", "Enriched in control SNPs" = "#56B4E9", "Not enriched" = "grey"),
    labels = c("Enriched in raQTLs", "Enriched in control SNPs", "Not enriched"),
    name = "Enrichment Status"
  ) +
  geom_label_repel(
    aes(label = ifelse(enriched == TRUE & -log10(p_fdr) >= 30, TFBS, "")),
    size = 6,
    fontface = 1,
    label.size= NA, 
    colour = "black",
    nudge_x = 0.2, nudge_y = 0.2,
    segment.color = "grey50",
    segment.size = 0.5,
    segment.alpha = 0.5
  ) +
  xlab("log2(raQTL proportion/control proportion)") +
  ylab("-log10(FDR-corrected enrichment P-value)") +
  ggtitle(paste0(volcano_celltype, " transcription factor binding site enrichment")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18),
    legend.text=element_text(size=18),
    legend.title=element_text(size=18), 
    plot.title = element_text(size=22)
     
  ) +
  scale_y_continuous(breaks = seq(0, 300, 100), labels = seq(0, 300, 100), limits = c(0, 300))
  
  ggsave(paste0("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/paper_plots/060325_volcano_enrichment_", volcano_celltype, ".pdf"), device="pdf", height=15, width=15)
  }
    )
    
    
    
#TFBS concordance plot 
concordance_plot <- lapply(cell_types, function(concordance_celltype){

df <- df2_enrichment[[concordance_celltype]]
df_remove <- df[df$MPRA_hits <= 5 | df$MPRA_nonhits <= 5,] 
df <- df[!(df$TFBS %in% df_remove$TFBS),] #Remove TFBS that only bind a few SNPs


ggplot(df, aes(x=log2(MPRA_hits_prop / MPRA_nonhits_prop), y=percentage_concordant_raqtl, color=enrichment_status)) + geom_point(mapping=aes(size=MPRA_hits)) +scale_color_manual(
    values = c("Enriched in raQTLs" = "#E69F00", "Enriched in control SNPs" = "#56B4E9", "Not enriched" = "grey"),
    labels = c("Enriched in raQTLs", "Enriched in control SNPs", "Not enriched"),
    name = "Enrichment Status"
  ) + scale_size_continuous(name="number of raQTLs") +
    geom_text_repel(aes(label=ifelse(enrichment_status=="Enriched in raQTLs" & percentage_concordant_raqtl >= 90, TFBS, "")), colour="black", size=2.82, fontface=1) + labs(x="log2(raQTL proportion/control proportion)", y="TFBS concordance (%)", title=paste0("Concordance between SuRE signal for ",concordance_celltype, " raQTLs and TFBS binding prediction")) + theme_bw() +
  theme(
    text=element_text(size=8)
     )    
  ggsave(paste0("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/paper_plots/060325_volcano_concordance_", concordance_celltype, ".pdf"), device="pdf", height=15, width=15)
  }
    )
  

#Barplot 

#Define color codes for each cell type (see above) 
cell_colors <- c("hNPC" = "#E69F00", "HepG2" = "#0072B2", "K562" = "#D55E00")

#Function to calculate enrichment of each cell type over the others (BARPLOT) 
calculate_relative_enrichment <- function(df_list, cell1, cell2) {
  df1 <- df_list[[cell1]]
  df2 <- df_list[[cell2]]
  
  #Merge datasets by TFBS
  df_merged <- merge(df1, df2, by = "TFBS", suffixes = c(".x", ".y"))
  
  #Filter for significantly enriched TFBS in either dataset
  df_enriched <- df_merged %>%
    filter((enriched.x == TRUE & p_fdr.x <= 0.05) | (enriched.y == TRUE & p_fdr.y <= 0.05))
  
  #Add pseudocount to avoid division by zero
  df_enriched <- df_enriched %>%
    mutate(
      MPRA_hits_prop.x = MPRA_hits_prop.x + 1,
      MPRA_hits_prop.y = MPRA_hits_prop.y + 1,
      log2_bar = log2(MPRA_hits_prop.x / MPRA_hits_prop.y)
    )
  
  
  #Combine top TFBS from both cell types
  df_plot <- df_enriched %>%
    mutate(color_bar = ifelse(log2_bar > 0, cell_colors[[cell1]], cell_colors[[cell2]]))
  
  return(df_plot)
}

  #Function to generate plot 
  generate_barplot <- function(df_plot, cell1, cell2) {
  ggplot(df_plot, aes(x = reorder(TFBS, log2_bar), y = log2_bar)) +
    geom_bar(stat = "identity", fill = df_plot$color_bar, alpha = 1) +
    scale_color_identity() +
    labs(x = "", y = paste0("Enrichment of motif disruptions log2(", cell1, "/", cell2, ")")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size = 30), axis.title = element_text(size = 30))
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
  ggsave(paste0("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/paper_plots/060325_barplot_", cell1, "_vs_", cell2, ".pdf"), device="pdf", width=20, height=15)
})







#The following plots were only generated for hNSCs 

#Highest expressing alleles plot
#Load SuRE data
df <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/project_on_processing/results/results_eline/freeze7/pvalue.freeze7.snp-permutation.04042024.txt.gz")

#Subset only based on fragment count (not expression)
df_sub <- df %>%
filter(ref.element.count >= 10 & ref.element.count <1000 & alt.element.count >= 10 & alt.element.count < 1000) %>%
dplyr::select(SNP_ID, chrom, pos.hg19, ref.seq, alt.seq, hNPC.cDNA.ref.mean, hNPC.cDNA.alt.mean, hNPC.wilcoxon.pvalue)

#Annotate TFBS
#Load SNP2TFBS data for determining what SNP binds our focal TFBS
snp2tfbs <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/resources/SNP2TFBS/snp2tfbs_JASPAR_CORE_2014_vert_RESTRUCTURED.txt.gz")
#Keep only rs or ss variants
tf_snp2tfbs <- snp2tfbs[grepl("rs|ss", snp2tfbs$SNP_ID),]


#Merge w/ df_sub
df_tf <- merge(df_sub, tf_snp2tfbs, by="SNP_ID")


#Now for all TFs we have the SNPs they bind. Now we want to know per TF the mean of the highest expressing SNPs 
df_tf2 <- df_tf %>%
mutate(TFBS=toupper(TFBS)) %>%
group_by(TFBS) %>% 
mutate(highest_expr=pmax(hNPC.cDNA.ref.mean, hNPC.cDNA.alt.mean)) %>%
mutate(mean_highest_expr=mean(highest_expr)) %>%
as.data.frame()

df_tf2_sorted <- df_tf2 %>%
  arrange(desc(mean_highest_expr))

#The following 2 steps are arbitrary choices (how many TFs do you want to plot?). We chose to visualize the TFBS with the highest expressing SNPs, going from high to lower (so stopping at "AR"). 
# Find the row number where TFBS == "AR"
ar_index <- which(df_tf2_sorted$TFBS == "AR")

# Subset the data to include all rows up to and including "AR"
df_tf2 <- df_tf2_sorted[1:ar_index, ]
df_tf2 <- df_tf2 %>%
          filter(!duplicated(TFBS)) %>%
          as.data.frame()

#DOTPLOT 
dotplot_highexp <- ggplot(df_tf2, aes(x=reorder(TFBS, -mean_highest_expr), y=mean_highest_expr, colour=mean_highest_expr)) + geom_point() + labs(title="Mean expression of highest-expressing alleles of SNPs binding a TF", x="TF", y="Mean expression of highest expressing alleles", colour="Expression") + scale_colour_gradient(high="red", low="yellow") + theme_minimal() + theme(axis.text.x=element_text(angle=90, size=10))

ggsave("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/paper_plots/060324_hNPC_dotplot.pdf", plot=dotplot_highexp, width=25)


#Correlation plot (between nucleotide conservation and relative raQTL fraction) 

#Load large SNP2TFBS file containing positions of SNPs binding all SNP2TFBS TFs (n=195) to overlap w/ SuRE data --> the difference between this df and the snp2tfbs df is that this one contains start + end positions of each tfbs so that we can determine the position of each snp within the tfbs  
tfs <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/TF_PWMs/210524_all_TFs.txt.gz")
#For the PWM data, we just need rsid, position (to be sure), start of PWM and end of PWM, REF/ALT alleles, motif itself and the TF name. 
tf_data <- select(tfs, V1, V3, V5, V6, V7, V8, V9, V18)
colnames(tf_data) <- c("SNP_ID", "pos.snp", "ref", "alt", "start_PWM", "end_PWM", "sequence", "TF")
#Remove instances where position of TFBS is unknown to avoid NA's later on
tf_data <- tf_data[!(tf_data$start_PWM=="." | tf_data$end_PWM=="."),]


#Load PWM data for all TFBS --> I made this file myself before by simply subsetting SNP2TFBS PWM's for nucleotide frequencies per position + added the TFBS name, then concatenated the data for all 195 TFBS. So this file just contains the nucleotide frequencies and TFBS names
pwm <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/TF_PWMs/matrices/epd.expasy.org/ftp/snp2tfbs/pwms/210524_PWM_logodds.txt")

#Load raw hNSC data (raQTLs + control SNPs) 
hnsc_raqtl <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/raqtls/results/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")
hnsc_raqtl$is_raqtl <- TRUE
    
hnsc_control <- fread("/gpfs/work4/0/AdamsLab/Projects/sure/raqtls/results/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024.txt")
hnsc_control$is_raqtl <- FALSE
hnsc_all <- rbind(hnsc_raqtl, hnsc_control) 

#Merge with TFBS data
hnsc_df <- merge(tf_data, hnsc_all, by="SNP_ID")
hnsc_df$TF <- toupper(hnsc_df$TF) 

#Add position per nucleotide within each TFBS in the PWM data
pwm2 <- pwm %>%
        group_by(TF_name) %>%
        mutate(pos_within_TFBS = row_number()) %>% #Because the row number corresponds with the position of the nucleotide within the TFBS
        ungroup() %>%
        as.data.frame()
        
pwm2$TF_name <- toupper(pwm2$TF_name)  

      

#Function that determines the position of each SuRE SNP within the PWM, then calculates the proportion of raQTLs vs controls at that position, then calculates the correlation between the highest PWM score for that position and the raQTL/control proportions. 
annot_SNP <- function(df, pwm){
             df$pos.sure.snp <- as.numeric(df$pos.hg19) - as.numeric(df$start_PWM) + 1 #Because position is 0-based we add +1 
             df <- merge(df, pwm, by.x=c("TF", "pos.sure.snp"), by.y=c("TF_name", "pos_within_TFBS"))
             df <- select(df, SNP_ID, chrom, pos.hg19, ref.seq, alt.seq, is_raqtl, pos.sure.snp, TF, start_PWM, end_PWM, A, C, G, T) 
             
             df <- df %>%
                   mutate(conservation_score=pmax(A,C,G,T)) %>%
                   group_by(TF) %>%
                   mutate(total_controls=sum(is_raqtl==FALSE), total_raqtls=sum(is_raqtl==TRUE)) %>%
                   ungroup() %>%
                   group_by(TF, pos.sure.snp) %>%
                   mutate(prop_raqtl=sum(is_raqtl==TRUE)/total_raqtls, prop_control=sum(is_raqtl==FALSE)/total_controls) %>%
                   mutate(prop_raqtl=ifelse(is.nan(prop_raqtl), 0, prop_raqtl), prop_control=ifelse(is.nan(prop_control), 0, prop_control),            
                   fraction_raqtl_vs_control=(prop_raqtl+ 1/(prop_control + 1))) %>%  
                   ungroup() %>%
                   group_by(TF) %>%
                   mutate(cor_TF=cor(fraction_raqtl_vs_control, conservation_score, use="complete.obs")) %>%
                   ungroup() %>%
                   as.data.frame() 
                   return(df)
             }
hnsc_df <- annot_SNP(hnsc_df, pwm2)  



#For the main plot, we visualized all TFBS enriched for raQTLs, and the top 10 most enriched for controls (total n=38 TFBS) 
#Load hNSC enrichment data (we loaded that at the beginning of the script already as "df2_enrichment" 
hnsc_enrichment <- as.data.frame(df2_enrichment[["hNPC"]])
df_top <- hnsc_enrichment %>%
              arrange(p_fdr) 
              
df_ctrl <- df_top %>%
           filter(enriched == FALSE & p_fdr <= 0.05) %>%
           slice_head(n=10) %>%
           mutate(experiment="control") %>%
           as.data.frame()

df_raqtl <- df_top %>%
           filter(enriched == TRUE & p_fdr <= 0.05) %>%
           mutate(experiment="raQTL") %>%
           as.data.frame()  
                        
all <- bind_rows(df_ctrl, df_raqtl)
all$TFBS <- toupper(all$TFBS) 

all_merged <- merge(all, hnsc_df, by.x="TFBS", by.y="TF")

all_main <- all_merged %>%
            filter(!duplicated(TFBS)) %>% #Because we did not summarize, so here the output contains duplicate tfbs even though we only need one cor_TF score per TF 
            select(TFBS, p_fdr, enriched, experiment, conservation_score, fraction_raqtl_vs_control, cor_TF) %>% 
            arrange(-cor_TF) %>%
            as.data.frame() 
            
plot_cor <- ggplot(all_main, aes(x=factor(TFBS, levels=unique(TFBS)), y=cor_TF)) +
                    geom_point(aes(color=cor_TF)) +
                    scale_color_gradient2(low = "purple", mid = "grey", high = "red", midpoint = 0) +  # Continuous color scale
                    theme_minimal(base_size = 16) +  # Base size for better readability
                    theme(
                    axis.title.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.text.x= element_text(angle=90,size=14), 
                    axis.text.y = element_text(size = 14),
                    plot.title = element_text(size = 16, face = "bold"),  # Increase title size and make it bold
                    axis.title.y = element_text(size = 14)  # Increase y-axis title size
                     ) +
                    ylim(-1, 1) + 
                    labs(
                    y = "Correlation between raQTL/control fraction with TFBS nucleotide conservation",
                    color = "Correlation"
                    )
            ggsave("/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/paper_plots/060325_hnsc_cor_main.pdf", plot=plot_cor, device="pdf", width=20, height=15)
                           
