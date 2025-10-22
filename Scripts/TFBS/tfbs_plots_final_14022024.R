library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)

# STEP 1: read files
path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/figures/"

df <- fread(paste0(path, "tfbs_percentages_for-plot_10042024.txt"))
stat.test <- fread(paste0(path, "tfbs_statistics_for-plot_10042024.txt"))

# STEP 2: prepare files the concordance rounds
# hNSC concordance
df.hnsc <- df[df$Cell == "hNSC" & df$Source == "TFBS concordance", ]
df.hnsc$Group = factor(df.hnsc$Group, levels=c("hNSC emVars", "hNSC controls", "All SuRE SNPs"))
stat.test.hnsc <- stat.test[stat.test$Cell == "hNSC" & stat.test$Source == "TFBS concordance", ]

# HepG2 concordance
df.hepg2 <- df[df$Cell == "HepG2" & df$Source == "TFBS concordance", ]
df.hepg2$Group = factor(df.hepg2$Group, levels=c("HepG2 emVars", "HepG2 controls", "All SuRE SNPs"))
stat.test.hepg2 <- stat.test[stat.test$Cell == "HepG2" & stat.test$Source == "TFBS concordance", ]

# K562 concordance
df.k562 <- df[df$Cell == "K562" & df$Source == "TFBS concordance", ]
df.k562$Group = factor(df.k562$Group, levels=c("K562 emVars", "K562 controls", "All SuRE SNPs"))
stat.test.k562 <- stat.test[stat.test$Cell == "K562" & stat.test$Source == "TFBS concordance", ]



# STEP 3: prepare files for the overlap rounds
# hNSC overlap
df.hnsco <- df[df$Cell == "hNSC" & df$Source == "TFBS overlap", ]
df.hnsco$Group = factor(df.hnsco$Group, levels=c("hNSC emVars", "hNSC controls", "All SuRE SNPs"))
stat.test.hnsco <- stat.test[stat.test$Cell == "hNSC" & stat.test$Source == "TFBS overlap", ]

# HepG2 overlap
df.hepg2o <- df[df$Cell == "HepG2" & df$Source == "TFBS overlap", ]
df.hepg2o$Group = factor(df.hepg2o$Group, levels=c("HepG2 emVars", "HepG2 controls", "All SuRE SNPs"))
stat.test.hepg2o <- stat.test[stat.test$Cell == "HepG2" & stat.test$Source == "TFBS overlap", ]

# K562 overlap
df.k562o <- df[df$Cell == "K562" & df$Source == "TFBS overlap", ]
df.k562o$Group = factor(df.k562o$Group, levels=c("K562 emVars", "K562 controls", "All SuRE SNPs"))
stat.test.k562o <- stat.test[stat.test$Cell == "K562" & stat.test$Source == "TFBS overlap", ]



# STEP 4: prepare files for the comparison rounds
# overlap comparison
df.overlap <- df[df$Source == "TFBS overlap", ]
stat.test.overlap <- stat.test[stat.test$Source == "TFBS overlap", ]

# concordance comparison
df.conc <- df[df$Source == "TFBS concordance", ]
stat.test.conc <- stat.test[stat.test$Source == "TFBS concordance", ]



# STEP 4: create plots for the concordance round
global_size = 12

hnsc <- ggplot(df.hnsc, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS concordance (%)") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.hnsc$Percentage + 10)), breaks=seq(0, max(df.hnsc$Percentage + 10), 10)) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.hnsc, y.position=stat.test.hnsc$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)


hepg2 <- ggplot(df.hepg2, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS concordance (%)") +
  scale_fill_manual(values = c("#0072B2", "#C46FA0", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.hepg2$Percentage + 10)), breaks=seq(0, max(df.hepg2$Percentage + 10), 10)) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.hepg2, y.position=stat.test.hepg2$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)


k562 <- ggplot(df.k562, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS concordance (%)") +
  scale_fill_manual(values = c("#D55E00", "#ECE242", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.k562$Percentage + 10)), breaks=seq(0, max(df.k562$Percentage + 10), 10)) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.k562, y.position=stat.test.k562$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)

# plots for the concordance round
ggsave(paste0(path, "hnsc_TFBSconcordance_Routput_14052024.pdf"), hnsc, units="mm", height=200, width=120)
ggsave(paste0(path, "hepg2_TFBSconcordance_Routput_14052024.pdf"), hepg2, units="mm", height=200, width=120)
ggsave(paste0(path, "k562_TFBSconcordance_Routput_14052024.pdf"), k562, units="mm", height=200, width=120)




# STEP 5: create plots for the overlap round
hnsc <- ggplot(df.hnsco, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS overlap (%)") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.hnsco$Percentage + 10)), breaks=seq(0, max(df.hnsco$Percentage + 10), 10)) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.hnsco, y.position=stat.test.hnsco$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)


hepg2 <- ggplot(df.hepg2o, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS overlap (%)") +
  scale_fill_manual(values = c("#0072B2", "#C46FA0", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.hepg2o$Percentage + 10)), breaks=seq(0, max(df.hepg2o$Percentage + 10), 10)) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.hepg2o, y.position=stat.test.hepg2o$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)


k562 <- ggplot(df.k562o, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS overlap (%)") +
  scale_fill_manual(values = c("#D55E00", "#ECE242", "#009E73")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.k562o$Percentage + 10)), breaks=seq(0, max(df.k562o$Percentage + 10), 10)) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.k562o, y.position=stat.test.k562o$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)

# plots for the overlap round
ggsave(paste0(path, "hnsc_TFBSoverlap_Routput_14052024.pdf"), hnsc, units="mm", height=200, width=120)
ggsave(paste0(path, "hepg2_TFBSoverlap_Routput_14052024.pdf"), hepg2, units="mm", height=200, width=120)
ggsave(paste0(path, "k562_TFBSoverlap_Routput_14052024.pdf"), k562, units="mm", height=200, width=120)




# STEP 6: create plots for the comparison round
#group by cell line
tfbs <- ggplot(df.overlap, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS overlap (%)") +
  scale_fill_manual(values = c('#009E73',"#C46FA0", "#0072B2",
                               "#56B4E9", "#E69F00", 
                               "#ECE242", "#D55E00")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.overlap$Percentage + 10)), breaks=seq(0, max(df.overlap$Percentage + 10), 10)) +
  facet_grid(~ Cell, scales="free_x") +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.overlap, y.position=stat.test.overlap$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)

conc <- ggplot(df.conc, aes(x=Group, y=Percentage)) +
  geom_bar(aes(fill=Group), stat="identity", color = "black", position=position_dodge()) +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "", y="TFBS concordance (%)") +
  scale_fill_manual(values = c('#009E73',"#C46FA0", "#0072B2",
                               "#56B4E9", "#E69F00", 
                               "#ECE242", "#D55E00")) +
  theme_bw() +
  theme(text = element_text(size=global_size), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, max(df.conc$Percentage + 10)), breaks=seq(0, max(df.conc$Percentage + 10), 10)) +
  facet_grid(~ Cell, scales="free_x") +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label=p.signif), data=stat.test.conc, y.position=stat.test.conc$y.position,
    inherit.aes=FALSE, tip.length=0.01, size=0.4, vjust=0.5, label.size=4)

# save plots for the comparison round
ggsave(paste0(path, "TFBS_overlap_comparison_Routput_14052024.pdf"), tfbs, units="mm", height=100, width=200)
ggsave(paste0(path, "TFBS_concordance_comparison_Routput_14052024.pdf"), conc, units="mm", height=100, width=200)




