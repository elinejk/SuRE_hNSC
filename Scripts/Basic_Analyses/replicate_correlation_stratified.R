library(data.table)
library(dplyr)
library(ggplot2)

### PREPARE FILES ###
# read in relevant files
a <- fread("/projects/0/AdamsLab/Projects/sure/project_on_processing/results/results_eline/fr7_for_replicates/pvalue.freeze7.snp-permutation.withreplicates.23052025.txt.gz")
emvar <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")
controls <- fread("/gpfs/home6/ekoornstra/raqtls/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024.txt")
active <- rbind(emvar, controls)

# subset the full sure file to only retain the relevant columns
df <- select(a, all_of(c("SNP_ID", "ref.element.count", "alt.element.count", 
                         "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean", 
                         "hNPC.wilcoxon.pvalue", "hNPC.wilcoxon.pvalue.random",
                         "hNPC.cDNA.rep1.ref.mean", "hNPC.cDNA.rep1.alt.mean",
                         "hNPC.cDNA.rep2.ref.mean", "hNPC.cDNA.rep2.alt.mean")))

# Add emvar and active SNP annotation
df <- df %>%
  mutate(IS_EMVAR = ifelse(SNP_ID %in% emvar$SNP_ID, "1", "0"),
         IS_ACTIVE = ifelse(SNP_ID %in% active$SNP_ID, "1", "0"))

# add effect and concordance columns
df$EFFECT <- log2(df$hNPC.cDNA.ref.mean/df$hNPC.cDNA.alt.mean)
df$REP1_EFFECT <- log2(df$hNPC.cDNA.rep1.ref.mean/df$hNPC.cDNA.rep1.alt.mean)
df$REP2_EFFECT <- log2(df$hNPC.cDNA.rep2.ref.mean/df$hNPC.cDNA.rep2.alt.mean)
df$CONCORDANCE <- ifelse(df$REP1_EFFECT > 0 & df$REP2_EFFECT > 0, 1, ifelse(df$REP1_EFFECT < 0 & df$REP2_EFFECT < 0, 1, 0))

# Create separate data.tables for emvars and active snps
e <- df[df$IS_EMVAR == "1", ]
b <- df[df$IS_ACTIVE == "1", ]


### GET CONCORDANCE STATS ###
table(e$CONCORDANCE)
table(b$CONCORDANCE)


### FUNCTIONs TO GET CORRELATION AND CONCORDANCE SPLIT ON BINS ###
# Function to compute the correlation and concordance
compute_metrics <- function(subdt) {
  
  safe_cor <- function(x, y) {
    
    df <- data.frame(x = x, y = y)
    df <- df[is.finite(df$x) & is.finite(df$y), ]
    
    if (nrow(df) < 2) return(NA_real_)
    if (sd(df$x) == 0 || sd(df$y) == 0) return(NA_real_)
    
    cor(df$x, df$y, method = "pearson")
  }
  
  safe_cor2 <- function(x, y) {
    
    df <- data.frame(x = x, y = y)
    df <- df[is.finite(df$x) & is.finite(df$y), ]
    
    if (nrow(df) < 2) return(NA_real_)
    if (sd(df$x) == 0 || sd(df$y) == 0) return(NA_real_)
    
    cor(df$x, df$y, method = "spearman")
  }
  
  rep1_all <- c(subdt$hNPC.cDNA.rep1.ref.mean, subdt$hNPC.cDNA.rep1.alt.mean)
  rep2_all <- c(subdt$hNPC.cDNA.rep2.ref.mean, subdt$hNPC.cDNA.rep2.alt.mean)
  fin <- subdt[is.finite(subdt$REP1_EFFECT) & is.finite(subdt$REP2_EFFECT), ]
  effect_diff <- fin$REP1_EFFECT - fin$REP2_EFFECT
  
  list(
    n = nrow(subdt),
    cor_ref = safe_cor2(log(subdt$hNPC.cDNA.rep1.ref.mean), log(subdt$hNPC.cDNA.rep2.ref.mean)),
    cor_alt = safe_cor2(log(subdt$hNPC.cDNA.rep1.alt.mean), log(subdt$hNPC.cDNA.rep2.alt.mean)),
    cor_combined = safe_cor2(log(rep1_all), log(rep2_all)),
    cor_effect = safe_cor(fin$REP1_EFFECT, fin$REP2_EFFECT),
    concordance = if (nrow(subdt) == 0) NA_real_ else mean(subdt$CONCORDANCE, na.rm = TRUE)
  )
}

# Function to split the bins and calculate the metrics
analyze_bins <- function(dt, n_bins = 10) {
  dt <- copy(dt)
  
  # Remove problematic rows
  dt <- dt[!is.na(hNPC.wilcoxon.pvalue) & !is.na(EFFECT) & is.finite(EFFECT)]
  
  ## Create quantile bins ##
  # p-value (use -log10)
  dt[, pval_bin := cut(
    -log10(hNPC.wilcoxon.pvalue),
    breaks = quantile(-log10(hNPC.wilcoxon.pvalue), probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
  )]
  
  # absolute effect size
  dt[, effect_bin := ntile(abs(EFFECT), 10)]
  
  # allele coverage bins
  dt[, coverage_bin := ntile(ref.element.count + alt.element.count, 10)]
  dt[, ref_coverage_bin := ntile(ref.element.count, 10)]
  dt[, alt_coverage_bin := ntile(alt.element.count, 10)]
  
  ## Apply function to bins ##
  # By p-value bins
  res_pval <- dt[, {
    metrics <- compute_metrics(.SD)
    
    c(metrics,
      list(pval_min = min(hNPC.wilcoxon.pvalue, na.rm = TRUE), pval_max = max(hNPC.wilcoxon.pvalue, na.rm = TRUE),
           logp_min = min(-log10(hNPC.wilcoxon.pvalue), na.rm = TRUE),logp_max = max(-log10(hNPC.wilcoxon.pvalue), na.rm = TRUE)))
    
  }, by = pval_bin]
  
  # By effect bins
  res_effect <- dt[, {
    metrics <- compute_metrics(.SD)
    
    c(metrics,
      list(effect_min = min(EFFECT, na.rm = TRUE), effect_max = max(EFFECT, na.rm = TRUE),
           abs_effect_min = min(abs(EFFECT), na.rm = TRUE), abs_effect_max = max(abs(EFFECT), na.rm = TRUE)))
    
  }, by = effect_bin]
  
  # By coverage bins
  res_coverage <- dt[, {
    metrics <- compute_metrics(.SD)
    
    c(metrics,
      list(cov_min = min(ref.element.count + alt.element.count, na.rm = TRUE), 
           cov_max = max(ref.element.count + alt.element.count, na.rm = TRUE)))
    
  }, by = coverage_bin]
  
  res_ref_coverage <- dt[, {
    metrics <- compute_metrics(.SD)
    
    c(metrics,
      list(ref_cov_min = min(ref.element.count, na.rm = TRUE), 
           ref_cov_max = max(ref.element.count, na.rm = TRUE)))
    
  }, by = ref_coverage_bin]
  
  res_alt_coverage <- dt[, {
    metrics <- compute_metrics(.SD)
    
    c(metrics,
      list(alt_cov_min = min(alt.element.count, na.rm = TRUE), 
           alt_cov_max = max(alt.element.count, na.rm = TRUE)))
    
  }, by = alt_coverage_bin]
  
  # Results in list
  return(list(
    pval_bins = res_pval[order(pval_bin)],
    effect_bins = res_effect[order(effect_bin)],
    coverage_bins = res_coverage[order(coverage_bin)],
    ref_coverage_bins = res_ref_coverage[order(ref_coverage_bin)],
    alt_coverage_bins = res_alt_coverage[order(alt_coverage_bin)]
  ))
}

    
### APPLY FUNCTION TO THE DATASETS ###
res_main <- analyze_bins(df)
res_emvar <- analyze_bins(e)
res_active <- analyze_bins(b)

    
### PLOT RESULTS ###
# pvalue
pval_dt <- as.data.table(res_emvar$pval_bins)
pval_dt_long <- melt(pval_dt,
                     measure.vars = c("cor_combined", "concordance"),
                     variable.name = "metric",
                     value.name = "value")

pplot <- ggplot(pval_dt_long, aes(x = pval_bin, y = value, color = metric)) +
  geom_line(aes(group = metric)) +
  geom_point(size = 2) +
  scale_color_manual(name = "Metric",
                     labels = c("Activity correlation (rho)", "Effect concordance"),
                     values = c("#ec9006","#d24e01")) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1),
                     limits=c(0, 1)) +
  labs(x = "-log10(p-value) bin (bin 1 = least significant)", 
       y = "Metric",
       title = "emVars stratified by statistical significance") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("/projects/0/AdamsLab/Projects/sure/emvars_stratified_by_pval.pdf",pplot)

# absolute effects
effect_dt <- as.data.table(res_emvar$effect_bins)
ef_dt_long <- melt(effect_dt,
                     measure.vars = c("cor_combined", "concordance"),
                     variable.name = "metric",
                     value.name = "value")

eplot <- ggplot(ef_dt_long, aes(x = effect_bin, y = value, color = metric)) +
  geom_line(aes(group = metric)) +
  geom_point(size = 2) +
  scale_color_manual(name = "Metric",
                     labels = c("Activity correlation (rho)", "Effect concordance"),
                     values = c("#ec9006","#d24e01")) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1),
                     limits=c(0, 1)) +
  labs(x = "Absolute effect size bin (bin 1 = smallest effects)",
       y = "Metric",
       title = "emVars stratified by effect size") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("/projects/0/AdamsLab/Projects/sure/emvars_stratified_by_effect.pdf",
       eplot)

# coverage
coverage_dt <- as.data.table(res_emvar$coverage_bins)
cov_dt_long <- melt(coverage_dt,
                    measure.vars = c("cor_combined", "concordance"),
                    variable.name = "metric",
                    value.name = "value")


covplot <- ggplot(cov_dt_long, aes(x = coverage_bin, y = value, color = metric)) +
  geom_line(aes(group = metric)) + 
  geom_point(size = 2) +
  scale_color_manual(name = "Metric",
                     labels = c("Activity correlation (rho)", "Effect concordance"),
                     values = c("#ec9006","#d24e01")) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1),
                     limits=c(0, 1)) +
  labs(x = "Coverage bin (bin 1 = lowest coverage)",
       y = "Metric (correlation or concordance)",
       title = "emVars stratified by total coverage") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("/projects/0/AdamsLab/Projects/sure/emvars_stratified_by_totalcoverage.pdf",
       covplot)

ref_cov_dt <- as.data.table(res_emvar$ref_coverage_bins)

refplot <- ggplot(ref_cov_dt, aes(x = ref_coverage_bin, y = cor_ref, color = "Activity correlation (REF; rho)")) +
  geom_line() +  
  geom_point(size = 2) +
  scale_color_manual(name = "Metric",
                     values = c("Activity correlation (REF; rho)" = "#a5a4a4")) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1),
                     limits=c(0, 1)) +
  labs(x = "Reference allele coverage bin (bin 1 = lowest coverage)",
       y = "Spearman correlation",
       title = "emVars stratified by reference allele coverage") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("/projects/0/AdamsLab/Projects/sure/emvars_stratified_by_refcoverage.pdf",
       refplot)


alt_cov_dt <- as.data.table(res_emvar$alt_coverage_bins)

altplot <- ggplot(alt_cov_dt, aes(x = alt_coverage_bin, y = cor_alt, color = "Activity correlation (ALT; rho)")) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(name = "Metric",
                     values = c("Activity correlation (ALT; rho)" = "#666666")) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_y_continuous(breaks=seq(-0.5, 1, 0.1),
                     limits=c(-0.5, 1)) +
  labs(x = "Alternative allele coverage bin (bin 1 = lowest coverage)",
       y = "Spearman correlation",
       title = "emVars stratified by alternative allele coverage") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("/projects/0/AdamsLab/Projects/sure/emvar_stratified_by_altcoverage.pdf",
       altplot)



ggsave("/projects/0/AdamsLab/Projects/sure/activesnp_minfrag_density.pdf", dminplot)
