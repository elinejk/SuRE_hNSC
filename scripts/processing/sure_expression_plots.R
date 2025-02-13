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
df_path <- "/projects/0/AdamsLab/Projects/sure/project_on_processing/results/results_eline/freeze7/freeze7_readsfiles/"
#df_path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze6/raQTLs/Figures/FINAL_FIGURES/Expression_plots/"
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

# start section for track plots
snppos <- a$SNPabspos[2]
snppos <- as.numeric(snppos)

bp <- (snppos-5000):(snppos+5000)
b <- data.table(bp)
names(b) <- "BP"

sure <- data.table()

for (row in 1:nrow(b)) {
  r <- b$BP[row]
  sel <- df[r >= df$start & r <= df$end, ]
  sure <- rbind(sure, sel)
  sure <- sure %>% distinct()
}

rm(df)
rm(sel)
gc()


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

# point plot
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


# generate the mean sure expression at each position
b$mean.hNPC.expr <- NA
b$strand <- NA

m <- b

# all libaries combined
for (row in 1:nrow(b)){
  r <- b$BP[row]
  
  # get the means for the plus strand
  p <- sure[r >= sure$start & r <= sure$end & sure$strand == '+', ]
  pexpr <- mean(p$cDNA.hNPC.norm.ipcr)
  b$mean.hNPC.expr[row] <- pexpr
  b$strand[row] <- "plus"
  
  # now for the minus strand
  d <- sure[r >= sure$start & r <= sure$end & sure$strand == '-', ]
  mexpr <- mean(d$cDNA.hNPC.norm.ipcr)
  m$mean.hNPC.expr[row] <- -(mexpr)
  m$strand[row] <- "minus"
}


final_expr <- rbind(b, m)

final_expr <- final_expr %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))

plexpr <- ggplot(final_expr, aes(x=BP, y=mean.hNPC.expr, colour = strand)) +
  geom_line() +
  geom_area(aes(fill=strand, group=strand), position='identity', alpha = 0.5) +
  scale_x_continuous(breaks=c(snppos - 4000, snppos - 2000, snppos, snppos + 2000, snppos + 4000)) +
  geom_vline(xintercept=snppos, colour = "black") +
  labs(x=paste0("base pair position, chr", chr), y = "Mean normalized SuRE expression", title = snp) + 
  scale_color_manual(values=c("#69B79F","#009E73")) +
  scale_fill_manual(values=c("#69B79F","#009E73")) +
  theme(axis.text=element_text(size=11),strip.text = element_text(size = 12),
        axis.title=element_text(size=13,face="bold")) +
  theme_classic()

ggsave(paste0(out_path, snp, "_combined_mean_sure_minus-plus_expression_5kb.pdf"), plexpr,units="mm", width = 200, height = 200)

# libraries split
sure <- sure %>%
  mutate(lib = ifelse(sure$library == "42-1" | sure$library == "42-2", "42",
                      ifelse(sure$library == "43-1" | sure$library == "43-2", "43",
                             ifelse(sure$library == "44-1" | sure$library == "44-2", "44", 
                                    "45")))) %>%
  mutate(genome = ifelse(sure$library == "42-1" | sure$library == "42-2", "HG02601, Pakistan",
                           ifelse(sure$library == "43-1" | sure$library == "43-2", "GM18983, Japan",
                                  ifelse(sure$library == "44-1" | sure$library == "44-2", "HG01241, Puerto Rico", 
                                         "HG03464, Sierra Leone"))))

genomes <- c("HG02601, Pakistan", "GM18983, Japan", "HG01241, Puerto Rico", "HG03464, Sierra Leone")

final <- data.table()

for (val in genomes) {
  sel <- sure[sure$genome == val, ]
  
  for (row in 1:nrow(b)){
    r <- b$BP[row]
    
    # get the means for the plus strand
    p <- sel[r >= sel$start & r <= sel$end & sel$strand == '+', ]
    pexpr <- mean(p$cDNA.hNPC.norm.ipcr)
    b$mean.hNPC.expr[row] <- pexpr
    b$strand[row] <- "plus"
    
    # now for the minus strand
    d <- sel[r >= sel$start & r <= sel$end & sel$strand == '-', ]
    mexpr <- mean(d$cDNA.hNPC.norm.ipcr)
    m$mean.hNPC.expr[row] <- -(mexpr)
    m$strand[row] <- "minus"
  }
  
  fin_expr <- rbind(b, m)
  fin_expr$genome <- val
  
  final <- rbind(final, fin_expr)
  
}

final <- final %>% 
  mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))


plexpr2 <- ggplot(final, aes(x=BP, y=mean.hNPC.expr, colour = strand)) +
  geom_line() +
  geom_area(aes(fill=strand, group=strand), position='identity', alpha = 0.5) +
  scale_x_continuous(breaks=c(snppos - 4000, snppos - 2000, snppos, snppos + 2000, snppos + 4000)) +
  geom_vline(xintercept=snppos, colour = "black") +
  labs(x=paste0("base pair position, chr", chr), y = "Mean normalized SuRE expression", title = snp) + 
  scale_color_manual(values=c("#69B79F","#009E73")) +
  scale_fill_manual(values=c("#69B79F","#009E73")) +
  facet_wrap(~genome, ncol = 1, strip.position = "right") +
  theme(axis.text=element_text(size=11),strip.text = element_text(size = 12),
        axis.title=element_text(size=13,face="bold")) +
  theme_classic()
  
ggsave(paste0(out_path, snp, "_separate_mean_sure_minus-plus_expression_5kb.pdf"),plexpr2,  units="mm",width = 200, height = 225)









