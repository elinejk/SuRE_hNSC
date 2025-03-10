library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

### paths
path <- "/gpfs/home6/ekoornstra/surexgwas/freeze7/distance-to-lead/"


### parameters
phenos <- c("allphenos") # the phenotypes you want to plot, here I am only plotting the combined phenotype set
sets <- c("controls", "all_sure") # the sets you want to compare your raQTLs to


### make a loop for our method
for (pheno in phenos){
	## load an prepare the raQTL file
	raqtl <- fread(paste0(path, "raqtls/", pheno, "_raqtls_gwassign_distance-snp-to-lead_16042024.txt"))
	raqtl <- raqtl[raqtl$ingwas == 1, ]
        raqtl <- raqtl[!is.na(raqtl$distance_lead), ]
	raqtl$Dataset <- "raQTL"
	cols <- c("Dataset", "distance_lead")
	raqtl <- select(raqtl, all_of(cols))
 
  	b <- raqtl
	
	## start a loop for the different control datasets to merge them to the raqtl set
	for (set in sets) {
        	# prepare the control dataset
		cont <- fread(paste0(path, set, "/", pheno, "_", set, "_gwassign_distance-snp-to-lead_16042024.txt"))
        	cont <- cont[cont$ingwas == 1, ]
        	cont <- cont[!is.na(cont$distance_lead), ]
        	cont$Dataset <- set
        	cont <- select(cont, all_of(cols))
        
        	# combine the two datasets, calculate the distance in kb, subset to only have within 100kb
        	b <- rbind(b, cont)
	
	}

	### retain those within 100 kb
	b$distance_lead <- as.numeric(b$distance_lead)
	b$distance_lead <- as.numeric(b$distance_lead/1000)
	b$Dataset <- factor(b$Dataset, levels = c("raQTL", "controls", "all_sure"))
	kb <- b[b$distance_lead <= 100, ]
        
  	### make the plots
  	pl <- ggplot(data = kb, aes(x=distance_lead, y = after_stat(density))) +
   	# the histogram wants to put 0 in the middle of a bin so the bin would be -5 to 5, but you want 0-10. you can change this by updating the boudaries. without the boundary you get 11 bins, and with the boundary you get 10
  	# here using blocks of 20 kb
   	geom_histogram(aes(fill=Dataset), position=position_dodge(), binwidth = 20, 
   	color="black", boundary=0) +
   	scale_x_continuous(breaks=seq(0, 100, by = 20)) +
   	labs(x = "Distance (in kb) to lead SNP", y="Density", title=pheno) +
   	theme_classic() +
   	scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73"))

   	ggsave(paste0(path, "figures/", pheno, "_gwassign_snp-to-lead_raqtls_vs_control_vs_allsure_histo_25022025.pdf"), pl, units="mm",width=200, height = 250)
   
}
