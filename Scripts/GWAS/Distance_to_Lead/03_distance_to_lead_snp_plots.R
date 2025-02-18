library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

### paths
path <- "/gpfs/home6/ekoornstra/surexgwas/freeze7/distance-to-lead/"


### parameters
phenos <- c("ADHD", "ASD", "BPD", "CP", "IQ", "MDD", "SCZ", "allphenos", "allphenos-shortestdist", "allphenos-gwassign")
sets <- c("controls", "sampled_matched_controls", "matched_controls")


### make a loop for our method
for (pheno in phenos){
	## load an prepare the raQTL file
	raqtl <- fread(paste0(path, "raqtls/", pheno, "_raqtls_distance-snp-to-lead_16042024.txt"))
	raqtl <- raqtl[raqtl$ingwas == 1, ]
  raqtl <- raqtl[!is.na(raqtl$distance_lead), ]
	raqtl$Dataset <- "raQTL"
	cols <- c("Dataset", "distance_lead")
	raqtl <- select(raqtl, all_of(cols))
	
	## start a loop for the different control datasets
	for (set in sets) {
	# prepare the control dataset
	cont <- fread(paste0(path, set, "/", pheno, "_", set, "_distance-snp-to-lead_16042024.txt"))
	cont <- cont[cont$ingwas == 1, ]
  cont <- cont[!is.na(cont$distance_lead), ]
	cont$Dataset <- set
	cont <- select(cont, all_of(cols))
	
	# combine the two datasets, calculate the distance in kb, subset to only have within 100kb
	b <- rbind(raqtl, cont)
	b$distance_lead <- as.numeric(b$distance_lead)
	b$distance_lead <- as.numeric(b$distance_lead/1000)
	b$Dataset <- factor(b$Dataset, levels = c("raQTL", set))
	kb <- b[b$distance_lead <= 100, ]
	
	# make the plots
	pl <- ggplot(data = kb, aes(x=distance_lead, y = after_stat(density))) +
    # the histogram wants to put 0 in the middle of a bin so the bin would be -5 to 5, but you want 0-10. you can change this by updating the boudaries. without the boundary you get 11 bins, and with the boundary you get 10
	# here using blocks of 20 kb
    geom_histogram(aes(fill=Dataset), position=position_dodge(), binwidth = 20, 
                   color="black", boundary=0) +
    scale_x_continuous(breaks=seq(0, 100, by = 20)) +
    labs(x = "Distance (in kb) to lead SNP", y="Density", title=pheno) +
    theme_classic() +
    scale_fill_manual(values=c("#E69F00","#56B4E9"))

	ggsave(paste0(path, "figures/raqtls_vs_", set, "/", pheno, "_snp-to-lead_raqtl_vs_", set, "_histo_17042024.pdf"), pl, units="mm",width=200, height = 250)
	
	}
}
