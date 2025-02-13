library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

args <- commandArgs(TRUE)

## VARIABLES AND FILES ##
# SNP information
chr <- args[1]
snp <- args[2]

print(paste0("generating plots for SNP: ", snp))

# files
df_path <- "/projects/0/AdamsLab/Projects/sure/project_on_processing/results/results_eline/freeze6/for_plots/freeze6_reads/"
out_path <- "/gpfs/home6/ekoornstra/raqtls/freeze7/expression_plots/"

df <- fread(paste0(df_path, "reads_sure_", chr, ".txt.gz"))

# subset df for just the SNP of interest
a <- df[df$SNP_ID == snp, ]
a <- a %>%
  mutate(lib = ifelse(a$library == "42-1" | a$library == "42-2", "42",
                      ifelse(a$library == "43-1" | a$library == "43-2", "43",
                             ifelse(a$library == "44-1" | a$library == "44-2", "44", 
                                    "45")))) %>%
  mutate(genome = ifelse(a$library == "42-1" | a$library == "42-2", "HG02601, Pakistan",
                         ifelse(a$library == "43-1" | a$library == "43-2", "GM18983, Japan",
                                ifelse(a$library == "44-1" | a$library == "44-2", "HG01241, Puerto Rico", 
                                       "HG03464, Sierra Leone"))))

table(a$SNPbaseInf, a$lib, a$genome)
a %>% group_by(genome) %>% summarise_at(vars(cDNA.hNPC.norm.ipcr), list(name=mean))
a %>% group_by(lib) %>% summarise_at(vars(cDNA.hNPC.norm.ipcr), list(name=mean))

snppos <- a$SNPabspos[2]
snppos <- as.numeric(snppos)

# segment plot
segs <- ggplot(a) + 
	geom_segment(aes(x=a$start, xend = a$end, y=a$cDNA.hNPC.norm.ipcr, yend=a$cDNA.hNPC.norm.ipcr, color = a$SNPbaseInf),
	             position = position_jitter(height=0.05, seed = 1)) + 
	geom_vline(xintercept=snppos, colour = "black") +
	scale_x_continuous(breaks=c(snppos - 400, snppos, snppos + 400)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) + 
  scale_color_manual(values = c("black", "grey")) +
	labs(x=paste0("base pair position, chr", chr), y = "Normalized SuRE expression", title = snp) + 
	theme_classic() + 
  theme(legend.title = element_blank(),
        axis.text=element_text(size=11),strip.text = element_text(size = 12),
        axis.title=element_text(size=13,face="bold"))

ggsave(paste0(out_path, snp, "_fragment_expression.pdf"), segs, width = 200, height = 200,units="mm")


# jitter plot
plot.mean <- function(x) {
  m <- mean(x)
  c(y=m, ymin=m, ymax=m)
}
pt <- ggplot(a, aes(x = a$SNPbaseInf, y = a$cDNA.hNPC.norm.ipcr, fill=a$SNPbaseInf, colour = a$SNPbaseInf)) +
	geom_jitter(binaxis = 'y', stackdir='center', binwidth=0.06,
	             stackgroups=TRUE, binpositions="all", stackratio=0.5, alpha=0.8) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) + 
  scale_color_manual(values = c("black", "grey")) +
  scale_fill_manual(values = c("black", "grey")) +
  stat_summary(fun.data="plot.mean", geom="errorbar", colour="red", width=0.5, linewidth=1) +
	labs(x="", y = "Normalized SuRE expression", title = snp) + 
  theme_classic() + 
  theme(legend.title = element_blank(),
        axis.text=element_text(size=11),strip.text = element_text(size = 12),
        axis.title=element_text(size=13,face="bold"))

ggsave(paste0(out_path, snp, "_allelic_expression_jitterplot.pdf"), pt, units="mm",width = 140, height = 200)





