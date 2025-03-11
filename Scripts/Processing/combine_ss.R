#### PACKAGES
library(data.table)
library(magrittr)

setDTthreads(2)

#### INIT

# summary statistics
d_res <- "/home/ekoornstra/sure_data/freeze7"
chr <- c(1:22, "X")
f_res <- paste0("out_chr", chr, ".txt.gz")
p_res <- file.path(d_res, f_res)

# output
f_out <- "pvalue.freeze7.snp-permutation.04042024.txt.gz"
p_out <- file.path(d_res, f_out)

#### MAIN

# load
ss <- lapply(p_res, fread) %>%
  rbindlist()

# save
fwrite(ss, p_out, sep = '\t')
