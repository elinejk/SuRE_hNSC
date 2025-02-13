#libraries
library(data.table)
library(tidyverse)


# load raqtl sets
input.dir <- "/gpfs/home6/ekoornstra/sure_data/freeze7/"
raqtl.hnpc <- fread(paste0(input.dir, "pvalue.freeze7.snp-permutation.04042024.txt.gz"))
#raqtl.hepg2 <- fread(paste0(input.dir,"hepg2_no_downsampling_snp-permutation_freeze7_controls_04042024.txt"))
#raqtl.k562 <- fread(paste0(input.dir,"k562_no_downsampling_snp-permutation_freeze7_controls_04042024.txt"))
output.prefix <- "/gpfs/home6/ekoornstra/raqtls/freeze7/"
condition <- "all SuRE"
raqtl.hepg2 <- raqtl.hnpc
raqtl.k562 <- raqtl.hnpc

# load SNP2TFBS file
# data is a gz-zipped txt file downloaded from ftp://ccg.vital-it.ch/snp2tfbs/mapped_files/ at 28-11-2019. file  = snp2tfbs_JASPAR_CORE_2014_vert.txt.gz
snp2tfbs.df.2 <- read.delim("/gpfs/home6/ekoornstra/resources/TFBS/snp2tfbs_JASPAR_CORE_2014_vert.txt.gz", header = FALSE, as.is = TRUE, )

# alternative rows are removed 
snp2tfbs.df.filtered <- snp2tfbs.df.2[grepl("rs|ss", snp2tfbs.df.2$V1),]
rm(snp2tfbs.df.2)

# reformat the snp2tfbs dataframe
df <- snp2tfbs.df.filtered
df$V10 <- as.numeric(df$V10)
df$ref.stronger.binding <- NA
rm(snp2tfbs.df.filtered)

# look at the scoredifferences column
df[df$V10 > 0,"ref.stronger.binding"] <- 0
df[df$V10 < 0,"ref.stronger.binding"] <- 1

# some of them have dots (no PWM)
#df[df$V8 == ".","ref.stronger.binding"] <- 0
#df[df$V9 == ".","ref.stronger.binding"] <- 1
df[df$V8 == "." & !df$V9 == "-" | df$V8 == "." & !df$V9 == "+","ref.stronger.binding"] <- 0
df[df$V9 == "." & !df$V8 == "-" | df$V9 == "." & !df$V8 == "+","ref.stronger.binding"] <- 1
sum(table(df$ref.stronger.binding))

# some have the same value, remove those
#df <- df[df$V8 != df$V9,]
#table(df$ref.stronger.binding)
df[df$V8 == df$V9, "ref.stronger.binding"] <- 2

# add stronger binding column
raqtl.hnpc$ref.stronger.binding <- df[match(raqtl.hnpc$SNP_ID, df$V1),"ref.stronger.binding"]
raqtl.hepg2$ref.stronger.binding <- df[match(raqtl.hepg2$SNP_ID, df$V1),"ref.stronger.binding"]
raqtl.k562$ref.stronger.binding <- df[match(raqtl.k562$SNP_ID, df$V1),"ref.stronger.binding"]

raqtl.hnpc$ref.stronger.binding[is.na(raqtl.hnpc$ref.stronger.binding)] <- "no match"
raqtl.hepg2$ref.stronger.binding[is.na(raqtl.hepg2$ref.stronger.binding)] <- "no match"
raqtl.k562$ref.stronger.binding[is.na(raqtl.k562$ref.stronger.binding)] <- "no match"

# add columns for max expression (EK)
raqtl.hnpc$max <- if_else(raqtl.hnpc$hNPC.cDNA.ref.mean > raqtl.hnpc$hNPC.cDNA.alt.mean, "ref", "alt")
raqtl.hepg2$max <- if_else(raqtl.hepg2$HepG2.cDNA.ref.mean > raqtl.hepg2$HepG2.cDNA.alt.mean, "ref", "alt")
raqtl.k562$max <- if_else(raqtl.k562$K562.cDNA.ref.mean > raqtl.k562$K562.cDNA.alt.mean, "ref", "alt")

# add tf and conc columns
df$tf <- df$V7
raqtl.hnpc$tf <- df[match(raqtl.hnpc$SNP_ID, df$V1),"tf"]
raqtl.hnpc$conc <- if_else(raqtl.hnpc$max == "ref" & raqtl.hnpc$ref.stronger.binding == 1, "conc",
                     if_else(raqtl.hnpc$max == "alt" & raqtl.hnpc$ref.stronger.binding == 0, "conc",
                             if_else(raqtl.hnpc$max == "ref" & raqtl.hnpc$ref.stronger.binding == 0, "disc",
                                     if_else(raqtl.hnpc$max == "alt" & raqtl.hnpc$ref.stronger.binding == 1, "disc", 
                                             if_else(raqtl.hnpc$ref.stronger.binding == 2, "equal","no match")))))

raqtl.hepg2$tf <- df[match(raqtl.hepg2$SNP_ID, df$V1),"tf"]
raqtl.hepg2$conc <- if_else(raqtl.hepg2$max == "ref" & raqtl.hepg2$ref.stronger.binding == 1, "conc",
                           if_else(raqtl.hepg2$max == "alt" & raqtl.hepg2$ref.stronger.binding == 0, "conc",
                                   if_else(raqtl.hepg2$max == "ref" & raqtl.hepg2$ref.stronger.binding == 0, "disc",
                                           if_else(raqtl.hepg2$max == "alt" & raqtl.hepg2$ref.stronger.binding == 1, "disc", 
                                                   if_else(raqtl.hepg2$ref.stronger.binding == 2, "equal","no match")))))

raqtl.k562$tf <- df[match(raqtl.k562$SNP_ID, df$V1),"tf"]
raqtl.k562$conc <- if_else(raqtl.k562$max == "ref" & raqtl.k562$ref.stronger.binding == 1, "conc",
                           if_else(raqtl.k562$max == "alt" & raqtl.k562$ref.stronger.binding == 0, "conc",
                                   if_else(raqtl.k562$max == "ref" & raqtl.k562$ref.stronger.binding == 0, "disc",
                                           if_else(raqtl.k562$max == "alt" & raqtl.k562$ref.stronger.binding == 1, "disc", 
                                                   if_else(raqtl.k562$ref.stronger.binding == 2, "equal","no match")))))


table(raqtl.hnpc$conc)
table(raqtl.hepg2$conc)
table(raqtl.k562$conc)

# get required values for concordance (EK)
## hnsc ##
over_count <- sum(
  sum(raqtl.hnpc$max == "ref" & raqtl.hnpc$ref.stronger.binding == 1, na.rm = TRUE), #conc
  sum(raqtl.hnpc$max == "alt" & raqtl.hnpc$ref.stronger.binding == 0, na.rm = TRUE), #conc
  sum(raqtl.hnpc$max == "ref" & raqtl.hnpc$ref.stronger.binding == 0, na.rm = TRUE), #disconc
  sum(raqtl.hnpc$max == "alt" & raqtl.hnpc$ref.stronger.binding == 1, na.rm = TRUE), #disconc
  sum(raqtl.hnpc$ref.stronger.binding == 2, na.rm = TRUE) # equal
)

raqtls <- nrow(raqtl.hnpc)

print(paste0("THE OVERLAP FOR hNSC ", condition, " IS:"))
print(paste0("overlap count: ", over_count))
print(paste0("overlap percentage: ", (over_count/raqtls)*100))

conc <- sum(
  sum(raqtl.hnpc$max == "ref" & raqtl.hnpc$ref.stronger.binding == 1, na.rm = TRUE),
  sum(raqtl.hnpc$max == "alt" & raqtl.hnpc$ref.stronger.binding == 0, na.rm = TRUE))

print(paste0("THE CONCORDANCE FOR hNSC ", condition, " IS:"))
print(paste0("concordance count: ", conc))
print(paste0("concordance percentage: ", (conc/over_count)*100))

## hepg2 ##
over_count2 <- sum(
  sum(raqtl.hepg2$max == "ref" & raqtl.hepg2$ref.stronger.binding == 1, na.rm = TRUE), #conc
  sum(raqtl.hepg2$max == "alt" & raqtl.hepg2$ref.stronger.binding == 0, na.rm = TRUE), #conc
  sum(raqtl.hepg2$max == "ref" & raqtl.hepg2$ref.stronger.binding == 0, na.rm = TRUE), #disconc
  sum(raqtl.hepg2$max == "alt" & raqtl.hepg2$ref.stronger.binding == 1, na.rm = TRUE), #disconc
  sum(raqtl.hepg2$ref.stronger.binding == 2, na.rm = TRUE) # equal
)

raqtls2 <- nrow(raqtl.hepg2)

print(paste0("THE OVERLAP FOR HepG2 ", condition, " IS:"))
print(paste0("overlap count: ", over_count2))
print(paste0("overlap percentage: ", (over_count2/raqtls2)*100))

conc2 <- sum(
  sum(raqtl.hepg2$max == "ref" & raqtl.hepg2$ref.stronger.binding == 1, na.rm = TRUE),
  sum(raqtl.hepg2$max == "alt" & raqtl.hepg2$ref.stronger.binding == 0, na.rm = TRUE))

print(paste0("THE CONCORDANCE FOR HepG2 ", condition, " IS:"))
print(paste0("concordance count: ", conc2))
print(paste0("concordance percentage: ", (conc2/over_count2)*100))

## k562 ##
over_count3 <- sum(
  sum(raqtl.k562$max == "ref" & raqtl.k562$ref.stronger.binding == 1, na.rm = TRUE), #conc
  sum(raqtl.k562$max == "alt" & raqtl.k562$ref.stronger.binding == 0, na.rm = TRUE), #conc
  sum(raqtl.k562$max == "ref" & raqtl.k562$ref.stronger.binding == 0, na.rm = TRUE), #disconc
  sum(raqtl.k562$max == "alt" & raqtl.k562$ref.stronger.binding == 1, na.rm = TRUE), #disconc
  sum(raqtl.k562$ref.stronger.binding == 2, na.rm = TRUE) # equal
)

raqtls3 <- nrow(raqtl.k562)

print(paste0("THE OVERLAP FOR K562 ", condition, " IS:"))
print(paste0("overlap count: ", over_count3))
print(paste0("overlap percentage: ", (over_count3/raqtls3)*100))

conc3 <- sum(
  sum(raqtl.k562$max == "ref" & raqtl.k562$ref.stronger.binding == 1, na.rm = TRUE),
  sum(raqtl.k562$max == "alt" & raqtl.k562$ref.stronger.binding == 0, na.rm = TRUE))

print(paste0("THE CONCORDANCE FOR K562 ", condition, " IS:"))
print(paste0("concordance count: ", conc3))
print(paste0("concordance percentage: ", (conc3/over_count3)*100))

# SAVE THE DATA AT THIS POINT
fwrite(raqtl.hnpc, paste0(output.prefix,"pvalue.freeze7.snp-permutation.tfbs.overlap.23072024.txt.gz"), sep='\t')
#fwrite(raqtl.hepg2, paste0(output.prefix,"hepg2.freeze7.", condition, "_wTF_09042024.txt"), sep='\t')
#fwrite(raqtl.k562, paste0(output.prefix,"k562.freeze7.", condition, "_wTF_09042024.txt"), sep='\t')


