library(ggplot2)
library(dplyr)
library(data.table)
library(forcats)

## PATHS AND SETTINGS ##
pheno <- "cp_lee2018"
phenotype <- "Cognitive Performance"
path <- "/gpfs/home6/ekoornstra/ldsc/freeze7/"

results <- fread(paste0(path, "heritability_results/with_controls/", pheno, "_hnsc_heritability_freeze7.results"))

## PREPARE FILE ##
rows <- c("Coding_UCSCL2_0","Conserved_LindbladTohL2_0","CTCF_HoffmanL2_0", "DGF_ENCODEL2_0","DHS_peaks_TrynkaL2_0",
          "Enhancer_AnderssonL2_0","Enhancer_HoffmanL2_0","FetalDHS_TrynkaL2_0","H3K27ac_HniszL2_0","H3K27ac_PGC2L2_0",
          "H3K4me1_peaks_TrynkaL2_0","H3K4me3_peaks_TrynkaL2_0","H3K9ac_peaks_TrynkaL2_0","Intron_UCSCL2_0",
          "PromoterFlanking_HoffmanL2_0","Promoter_UCSCL2_0","Repressed_HoffmanL2_0","SuperEnhancer_HniszL2_0",
          "TFBS_ENCODEL2_0","Transcr_HoffmanL2_0","TSS_HoffmanL2_0","UTR_3_UCSCL2_0","UTR_5_UCSCL2_0",
          "WeakEnhancer_HoffmanL2_0", "hnsc_raqtlL2_0", "hnsc_controlL2_0")

new_rows <-c("Coding (UCSC)","Conserved (Lindblad)","CTCF (Hoffman)", "DGF (ENCODE)","DHS (Trynka)",
             "Enhancer (Andersson)","Enhancer (Hoffman)","Fetal DHS (Trynka)","H3K27ac (Hnisz)","H3K27ac (PGC2)",
             "H3K4me1 (Trynka)","H3K4me3 (Trynka)","H3K9ac (Trynka)", 
             "hNSC control SNP", "hNSC raQTL",
             "Intron (UCSC)","Promoter (UCSC)", "PromoterFlanking (Hoffman)","Repressed (Hoffman)",
             "SuperEnhancer (Hnisz)","TFBS (ENCODE)","Transcribed (Hoffman)","TSS (Hoffman)",
             "3UTR (UCSC)","5UTR (UCSC)","WeakEnhancer (Hoffman)")

results_sub <- results[results$Category %in% rows, ]
results_sub$new_names <- new_rows

ylimit <- max(results_sub$Enrichment + results_sub$Enrichment_std_error + 0.5)

# add significance, > 0.05 not sign, > 0.01 indicates p-values 0.01 - 0.05, > 0.001 indicates 0.001-0.01 etc
results_sub <- results_sub %>%
  mutate(label= case_when(
    Enrichment_p > 0.05 ~ "", Enrichment_p > 0.01 ~ "*",
    Enrichment_p > 0.001 ~"**", Enrichment_p > 0.0001 ~"***", !is.na(Enrichment_p) ~ "****",
    TRUE ~ NA_character_
  ))

## MAKE THE PLOT ##
ldscplot <- ggplot(results_sub, aes(x=Category, y=Enrichment)) +
  geom_bar(stat="identity", color="black", fill = if_else(results_sub$Category == "hnsc_raqtlL2_0", "#56B4E9",
						if_else(results_sub$Category == "hnsc_controlL2_0", "#56B4E9", "#0072B2")),
           position=position_dodge()) +
  geom_errorbar(aes(ymin=pmax(Enrichment-Enrichment_std_error, 0), ymax=pmax(0,Enrichment+Enrichment_std_error)),
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=label), nudge_y = (results_sub$Enrichment_std_error + 0.5)) +
  guides(x = guide_axis(angle = 90))+
  labs(x="Annotation Category", y="Enrichment (Prop. h2g / Prop. SNPs)", title=phenotype)+
  scale_x_discrete(labels=new_rows) +
  theme_classic() +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0, ylimit)

ldscplot

ggsave(paste0(path, "figures/", pheno, "_hnsc_heritability-enr_freeze7_witcontrol.pdf"), ldscplot, units="mm",width=214, height =175)

