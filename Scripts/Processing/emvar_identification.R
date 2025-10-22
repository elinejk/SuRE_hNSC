library(data.table)
library(dplyr)

# paths and parameters
sure_path <- "/gpfs/home6/ekoornstra/sure_data/freeze7/pvalue.freeze7.snp-permutation.04042024.txt.gz"
out_path <- "/gpfs/home6/ekoornstra/raqtls/freeze7/"
runtype <- "no_downsampling_snp-permutation"
freeze <- "freeze7"
dag <- "04042024"

# read and subset sure file
sure <- fread(sure_path) %>%
  as.data.table

colss <- c("SNP_ID", "ref.seq", "alt.seq", "chrom", "pos.hg19",
           "ref.element.count", "alt.element.count", 
           "hNPC.cDNA.ref.mean", "hNPC.cDNA.alt.mean", "hNPC.wilcoxon.pvalue", "hNPC.wilcoxon.pvalue.random",
           "K562.cDNA.ref.mean", "K562.cDNA.alt.mean", "K562.wilcoxon.pvalue", "K562.wilcoxon.pvalue.random",
           "HepG2.cDNA.ref.mean", "HepG2.cDNA.alt.mean", "HepG2.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue.random")

sure <- select(sure, all_of(colss))

# subset on element counts
a <- sure[sure$ref.element.count >= 10 & sure$ref.element.count < 1000 & sure$alt.element.count >= 10 & sure$alt.element.count < 1000, ]

## HNSC ##
# with a loop
subs <- a[a$hNPC.cDNA.ref.mean > 4 | a$hNPC.cDNA.alt.mean > 4, ]

p1 <- sort(subs$hNPC.wilcoxon.pvalue[subs$hNPC.wilcoxon.pvalue < 0.05], decreasing = TRUE)
p2 <- sort(subs$hNPC.wilcoxon.pvalue.random[subs$hNPC.wilcoxon.pvalue.random < 0.05], decreasing = TRUE)

for (value in p1) {
  k1 <- sum(p1 < value)
  k2 <- sum(p2 < value)
  fdr <- k2/(k1+k2)
  if (fdr < 0.05) {
    print(paste0("For hNSCs, an FDR of 5% is reached when a p-value of ", value, " is used as cut off"))
    break
  }
}

# emVars
raqtls <- filter(subs, hNPC.wilcoxon.pvalue < value)
controls <- filter(subs, hNPC.wilcoxon.pvalue >= value)
paste0(nrow(raqtls), " hNSC raQTLs were called with a p-value <", value)
paste0(nrow(controls), " hNSC controls")
fwrite(raqtls, paste0(out_path, "hnsc_", runtype, "_", freeze, "_wilc-raqtls_", dag, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
fwrite(controls, paste0(out_path, "hnsc_", runtype, "_", freeze, "_controls_", dag, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")


## HEPG2 ##
# subset on mean and sort
subs <- a[a$HepG2.cDNA.ref.mean > 4 | a$HepG2.cDNA.alt.mean > 4, ]

p1 <- sort(subs$HepG2.wilcoxon.pvalue[subs$HepG2.wilcoxon.pvalue < 0.05], decreasing = TRUE)
p2 <- sort(subs$HepG2.wilcoxon.pvalue.random [subs$HepG2.wilcoxon.pvalue.random < 0.05], decreasing = TRUE)

# calculate cut off
for (value in p1) {
  k1 <- sum(p1 < value)
  k2 <- sum(p2 < value)
  fdr <- k2/(k1+k2)
  if (fdr < 0.05) {
    print(paste0("For HepG2, an FDR of 5% is reached when a p-value of ", value, " is used as cut off"))
    break
  }
}


# emVars
raqtls <- filter(subs, HepG2.wilcoxon.pvalue < value)
controls <- filter(subs, HepG2.wilcoxon.pvalue >= value)
paste0(nrow(raqtls), " HepG2 raQTLs were called with a p-value <", value)
paste0(nrow(controls), " HepG2 controls")
fwrite(raqtls, paste0(out_path, "hepg2_", runtype, "_", freeze, "_wilc-raqtls_", dag, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
fwrite(controls, paste0(out_path, "hepg2_", runtype, "_", freeze, "_controls_", dag, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")


## K562 ##
# subset on mean and sort 
subs <- a[a$K562.cDNA.ref.mean > 4 | a$K562.cDNA.alt.mean > 4, ]

p1 <- sort(subs$K562.wilcoxon.pvalue[subs$K562.wilcoxon.pvalue < 0.05], decreasing = TRUE)
p2 <- sort(subs$K562.wilcoxon.pvalue.random[subs$K562.wilcoxon.pvalue.random < 0.05], decreasing = TRUE)

# calculate cut off
for (value in p1) {
  k1 <- sum(p1 < value)
  k2 <- sum(p2 < value)
  fdr <- k2/(k1+k2)
  if (fdr < 0.05) {
    print(paste0("For K562, an FDR of 5% is reached when a p-value of ", value, " is used as cut off"))
    break
  }
}

# emVars
raqtls <- filter(subs, K562.wilcoxon.pvalue < value)
controls <- filter(subs, K562.wilcoxon.pvalue >= value)
paste0(nrow(raqtls), " K562 raQTLs were called with a p-value <", value) 
paste0(nrow(controls), " K562 controls")
fwrite(raqtls,  paste0(out_path, "k562_", runtype, "_", freeze, "_wilc-raqtls_", dag, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
fwrite(controls, paste0(out_path, "k562_", runtype, "_", freeze, "_controls_", dag, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")

