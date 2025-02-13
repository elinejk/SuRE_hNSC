library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)

# STEP 1: read files
path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/figures/"

df <- fread(paste0(path, "hepg2-locations_percentages_for-plot_24042024.txt"))
stat.test <- fread(paste0(path, "hepg2-locations_statistics_for-plot_24042024.txt"))


# STEP 2: prepare the order of the locations
levs <- c("E118 Active enhancer", "E118 Promoter",
          "HepG2 DHS broadpeak","HepG2 DHS narrowpeak")

df$Annotation= factor(df$Annotation, levels = levs)
df$Group = factor(df$Group, levels = c("HepG2 raQTLs", "HepG2 controls", "all SuRE SNPs"))

# STEP 3: prepare the locations of the statistics
x <- factor(stat.test$Annotation, 
            levels=levels(df$Annotation))
stat.test$Annotation <- (as.numeric(x, 1:7))-1

stat.test1 <- stat.test[stat.test$group2 == "HepG2 controls", ]
stat.test2 <- stat.test[stat.test$group2 == "all SuRE SNPs", ]

# STEP 4: set the font size
global_size = 11

# STEP 5: make the plot
p <- ggplot(df, aes(x=Annotation, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="Overlap (%)") +
  scale_fill_manual(values = c("#0072B2", "#C46FA0", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(df$Percentage + 4))) +
  geom_bracket(
    aes(xmin = 0.7 + Annotation, xmax = 1+ Annotation, label=p.signif), data=stat.test1, y.position=stat.test1$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, label.size=4) +
  geom_bracket(
    aes(xmin = 0.7 + Annotation, xmax = 1.3 + Annotation, label=p.signif), data=stat.test2, y.position = stat.test2$y.position,
    inherit.aes=FALSE, position="dodge", tip.length=0.01, size=0.4,label.size=4)


p

# STEP 6: save the plot
ggsave(paste0(path, "locations/hepg2_locations_freeze7_Routput_24042024.pdf"), p, units="mm", height=120, width=150)





