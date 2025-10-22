library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

### PATH AND PARAMETERS ###
path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/GWAS/freeze7/distance_to_lead/"

sets <- c("controls", "all_sure") # the sets you want to compare your emVars to


### PREPARE FILES ###
## load emVar file
raqtl <- fread(paste0(path, "20250923_allphenos_raqtls_gwassign_distance-snp-to-lead.txt"))
raqtl <- raqtl[raqtl$ingwas == 1, ]
#raqtl <- raqtl[raqtl$iscandidate == 1, ] # optional: filter on candidate SNPs (for LD)
raqtl <- raqtl[!is.na(raqtl$distance_lead), ]
raqtl$Dataset <- "emVar"

# Stats on distances
raqtl %>%
  group_by(phenotype) %>%
  summarise(mean_distance = mean(distance_lead)/1000,
            median_distance = median(distance_lead)/1000,
            min_distance = min(distance_lead)/1000,
            max_distance = max(distance_lead)/1000) %>%
  as.data.table()


# subset for further analyses
cols <- c("Dataset", "distance_lead")
raqtl <- select(raqtl, all_of(cols))
 
b <- raqtl
	
## start a loop for the different control datasets to merge them to the EMVAR set
for (set in sets) {
  # prepare the control dataset
  cont <- fread(paste0(path, "20250923_allphenos_", set, "_gwassign_distance-snp-to-lead.txt"))
  cont <- cont[cont$ingwas == 1, ]
  #cont <- cont[cont$iscandidate == 1, ] # optional: filter on candidate SNPs (for LD)
  cont <- cont[!is.na(cont$distance_lead), ]
  cont$Dataset <- set
  cont <- select(cont, all_of(cols))
        
  # combine the two datasets, calculate the distance in kb, subset to only have within 100kb
  b <- rbind(b, cont)
}


### FILTER ON DISTANCE ###
b$distance_lead <- as.numeric(b$distance_lead)
b$distance_lead <- as.numeric(b$distance_lead/1000)
b$Dataset <- factor(b$Dataset, levels = c("emVar", "controls", "all_sure"))
kb <- b[b$distance_lead <= 100, ]
#kb <- b # if you don't want to filter on distance

        

### MAKE THE PLOT ##
pl <- ggplot(data = kb, aes(x=distance_lead, y = after_stat(density))) +

  geom_histogram(aes(fill=Dataset), position=position_dodge(), binwidth = 20, 
  color="black", boundary=0) +
  scale_x_continuous(breaks=seq(0, 100, by = 20)) +
  labs(x = "Distance (in kb) to lead SNP", y="Density") +
  theme_classic() +
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73"))

ggsave(paste0(path, "allphenos_gwassign_snp-to-lead_emVar_vs_control_vs_allsure_histo_23092025.pdf"), pl, units="mm",width=200, height = 250)


### CALCULATE P-VALUES ###
## first put them into bins
sub <- b[b$distance_lead < 100, ]

sub <- sub %>%
  mutate(distance_bin = ifelse(sub$distance_lead < 20, "0-20",
                        ifelse(sub$distance_lead >= 20 & sub$distance_lead < 40, "20-40",
                        ifelse(sub$distance_lead >= 40 & sub$distance_lead < 60, "40-60",
                        ifelse(sub$distance_lead >= 60 & sub$distance_lead < 80, "60-80", "80-100")))))

res <- table(sub$Dataset, sub$distance_bin) %>% as.data.table()
names(res) <- c("dataset", "bin", "count")

# add total SNPs per dataset
total_snps <- res[, .(total = sum(count)), by = dataset]
res <- merge(res, total_snps, by = "dataset")

# turn into wide format
count_wide <- dcast(res, bin ~ dataset, value.var = "count", fill = 0)
setnames(count_wide, old = c("emVar", "controls", "all_sure"), new = c("count_emVar", "count_control", "count_all_sure"))

total_wide <- dcast(res, bin ~ dataset, value.var = "total", fill = 0)
setnames(total_wide, old = c("emVar", "controls", "all_sure"), new = c("total_emVar", "total_control", "total_all_sure"))

res_wide <- merge(count_wide, total_wide, by = "bin") %>% as.data.table()

# create function to compute Fisher test between two datasets per bin
fisher_compare <- function(dtw, d1, d2) {
  dtw[, {
    d1_in <- get(paste0("count_", d1))
    d2_in <- get(paste0("count_", d2))
    d1_total <- get(paste0("total_", d1))
    d2_total <- get(paste0("total_", d2))
    
    # 2x2 contingency table
    tbl <- matrix(c(
      d1_in, d1_total - d1_in,
      d2_in, d2_total - d2_in
    ), nrow = 2, byrow = TRUE)
    
    list(p.value = fisher.test(tbl)$p.value)
  }, by = bin]
}

# Run comparisons
result_ec <- fisher_compare(res_wide, "emVar", "control")
result_ec[, comparison := "emVar_vs_control"]

result_ea <- fisher_compare(res_wide, "emVar", "all_sure")
result_ea[, comparison := "emVar_vs_all"]

# Combine results
all_results <- rbind(result_ec, result_ea)

print(all_results)



