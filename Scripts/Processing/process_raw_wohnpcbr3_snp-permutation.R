#### PACKAGES
library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)

setDTthreads(32)

#### FUNCTIONS
fun_tau <- function(x) {
  x <- tabulate(x)
  sum(x^3 - x)
}

sander_log <- function(...) {
  time <- format(Sys.time(), "%H:%M:%S ")
  do.call("message", list(time, ...))
}

support_z_1 <- function(n.x, n.y) n.x * (n.x + 1) / 2 + n.x * n.y / 2
support_z_2 <- function(n.x, n.y, tau) sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - tau/((n.x + n.y) * (n.x + n.y - 1))))

#### ARGS
args <- commandArgs(trailingOnly=TRUE)
chr <- args[[1]]
downsampling_must_happen <- FALSE
p_res <- args[[2]]

#### INIT
redo_totals <- FALSE
chrs <- c(1:22, "X")

#### PATHS
wd <- "/projects/0/AdamsLab/Projects/sure/project_on_processing"
p_data <- file.path(wd, "data/updated_with_hnpc2")

p_temp <- file.path(wd, "temp_permutation")
if (!dir.exists(p_temp)) dir.create(p_temp)

#p_res <- file.path(wd, "results")
if (!dir.exists(p_res)) dir.create(p_res)

totals_rds <- file.path(p_temp, "totals.rds")

# downsampling info
f_ds <- c("k562.downsampling.vector.SuRE42_1.180803.rds",
          "HepG2.downsampling.vector.SuRE42_1.180803.rds",
          "HepG2.downsampling.vector.SuRE43_1.180803.rds")
p_ds <- file.path(p_data, f_ds)
ds <- lapply(p_ds, readRDS)
dsv <- list(`42_1` = c(ds[[1]], ds[[2]]),
            `43_1` = ds[[3]])
#ds_cols_hnpc <- paste0("cDNA.hNPC.B",1:2)
ds_cols_k562 <- paste0("cDNA.K562.B",1:3)
ds_cols_hepg2 <- paste0("cDNA.HepG2.B",1:2)

# find all files
sure_dirs <- list.files(p_data, pattern = "sure", full.names = TRUE)
#sure_dirs <- sure_dirs[!grepl("sure42-2", sure_dirs, fixed = TRUE)]
sure_files <- lapply(sure_dirs, list.files, full.names = TRUE)
names(sure_files) <- sapply(sure_files, 
                            function(x) sub(".*sure([0-9]{1}.*)$", 
                                            "\\1", 
                                            dirname(x[1])))

### FIRST, GET TOTAL COUNTS FOR NORMALISATION
sum_u <- c("hNPC_B1_T1", "hNPC_B2_T1", "K562_B1_T1", "K562_B2_T1", "K562_B3_T1", "HEPG2_B1_T1", "HEPG2_B2_T1")
new_names <- c("cDNA.hNPC.B1","cDNA.hNPC.B2", "cDNA.K562.B1", "cDNA.K562.B2", "cDNA.K562.B3", "cDNA.HepG2.B1", "cDNA.HepG2.B2")

if (!file.exists(totals_rds) | redo_totals == TRUE) {
  sander_log(": Calculating totals...")

  totals <- lapply(seq_along(sure_files), function(i) {
    i_name <- names(sure_files)[i]
    sander_log(": Calculating totals for ", i_name)
    
    s_name <- sub("-", "_", i_name, fixed = TRUE)
    sum_cols <- c("iPCR", paste0("SuRE", s_name, "_", sum_u))
    rowSums(sapply(seq_along(chrs), function(j) colSums(data.table::fread(sure_files[[i]][j], select = sum_cols))))
  }) %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    set_names(c("count", new_names)) 
    
  saveRDS(totals, totals_rds)
  
  message(format(Sys.time(), "%H:%M:%S"), ": Totals calculated.")
} else {
  totals <- readRDS(totals_rds)
}

##################
sander_log(": Running chromosome ", chr)
system("free")

out_sure <- file.path(p_temp, paste0("reads_sure_", chr, ".txt.gz"))

counts.dt <- lapply(seq_along(sure_files), function(i) {

  ### LOAD + SIMPLE CLEAN ###

  i_name <- names(sure_files)[i]
  sander_log(": Running ", i_name)
  sander_log(": Memory status:")
  system("free")

  # load data
  chr_id <- grep(paste0("chr", chr, ".txt.gz"), sure_files[[i]], fixed = TRUE)
  a <- fread(sure_files[[i]][chr_id], drop = c("SNPrelpos", "SNPvar", "PAT_MAT", "SNPidx"))
  
  # remove reads without a SNP in them
  sander_log(": Rows before removing non-SNPs: ", nrow(a))
  a <- subset(a, SNP_ID != "" & !is.na(SNP_ID))
  sander_log(": Rows after removing non-SNPs: ", nrow(a))
  a <- subset(a, SNPabspos != "" & !is.na(SNPabspos))
  sander_log(": Rows after removing rows without position: ", nrow(a))

  # rename things
  setnames(a, "iPCR", "count")
  
  s_name <- sub("-", "_", i_name, fixed = TRUE)
  old_names <- paste0("SuRE", s_name, "_", sum_u)
  
  setnames(a, old_names, new_names)
  
  # remove `chr` from the chr column
  a[, chr := as.numeric(sub("chr", "", chr, fixed = TRUE))]
  
  ### DOWNSAMPLING ###
  if (downsampling_must_happen == TRUE) {

    if (i_name == "42-1") {
      ds_cols <- c(ds_cols_k562, ds_cols_hepg2)
      nr <- nrow(a)
      sander_log(": downsampling library 42_1")
      
      for (j in seq_along(ds_cols)) {
        ds_temp <- rep(1:nr, a[, ds_cols[j], with = FALSE][[1]])
        nds <- round(length(ds_temp) * as.numeric(dsv[["42_1"]][j]))
        ds_rle <- rle(sort(sample(ds_temp, size = nds)))
        a[ds_rle$values, (ds_cols[j]) := ds_rle$lengths]
      }
      
      totals$cDNA.K562.B1[i] <- totals$cDNA.K562.B1[i] * dsv[["42_1"]][1]
      totals$cDNA.K562.B2[i] <- totals$cDNA.K562.B2[i] * dsv[["42_1"]][2]
      totals$cDNA.K562.B3[i] <- totals$cDNA.K562.B3[i] * dsv[["42_1"]][3]
      totals$cDNA.HepG2.B1[i] <- totals$cDNA.HepG2.B1[i] * dsv[["42_1"]][4]
      totals$cDNA.HepG2.B2[i] <- totals$cDNA.HepG2.B2[i] * dsv[["42_1"]][5]
    }
    
    if (i_name == "43-1") {
      ds_cols <- ds_cols_hepg2
      nr <- nrow(a)
      sander_log(": downsampling library 43_1")
  
      for (j in seq_along(ds_cols)) {
        ds_temp <- rep(1:nr, a[, ds_cols[j], with = FALSE][[1]])
        nds <- round(length(ds_temp) * as.numeric(dsv[["43_1"]][j]))
        ds_rle <- rle(sort(sample(ds_temp, size = nds)))
        a[ds_rle$values, (ds_cols[j]) := ds_rle$lengths]
      }
      
      totals$cDNA.HepG2.B1[i] <- totals$cDNA.HepG2.B1[i] * dsv[["43_1"]][1]
      totals$cDNA.HepG2.B2[i] <- totals$cDNA.HepG2.B2[i] * dsv[["43_1"]][2]
    }
  
  }
    
  ### REMOVE HNPC THIRD REPLICATE ###
  if (i_name == "42-2") {
    a[, "SuRE42_2_hNPC_B3_T1" := NULL]
  }

  
  ### LONG FORMAT ###
  sander_log(": Creating long format for ", i_name)
  
  # given that reads can contain the same snp twice, we have to remove the duplicates
  # first, we find the duplicates
  snp.duplicates <- unlist(lapply(strsplit(a$SNP_ID, split = ","), duplicated))
  
  # then we go to long format
  #sure_snp <- a[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE, type.convert = TRUE)))][!is.na(SNP_ID)]
  sure_snp <- as.data.table(separate_rows(a, SNPbase, SNPbaseInf, SNPabspos, SNP_ID, SNPvarInf, sep = ","))
  sure_snp[, SNPabspos := as.numeric(SNPabspos)]
  sure_snp[, SNPvarInf := as.numeric(SNPvarInf)]
  sure_snp[, library := ..i_name]
  
  # remove the duplicates
  sure_snp <- sure_snp[!snp.duplicates, ]
  
  # remove those where SNPvarInf %in% c(2,3)
  sure_snp <- sure_snp[SNPvarInf %in% 0:1]
  
  sander_log(": Number of rows: ", nrow(sure_snp))
  sander_log(": Memory status:")
  system("free")
  
  ### NORMALIZE ###
  sure_snp[, count := round(count / totals$count[i] * 1e9, digits = 2)]
  sure_snp[, cDNA.hNPC.B1 := round(cDNA.hNPC.B1 / totals$cDNA.hNPC.B1[i] * 1e9, digits = 1)]
  sure_snp[, cDNA.hNPC.B2 := round(cDNA.hNPC.B2 / totals$cDNA.hNPC.B2[i] * 1e9, digits = 1)]
  sure_snp[, cDNA.K562.B1 := round(cDNA.K562.B1 / totals$cDNA.K562.B1[i] * 1e9, digits = 1)]
  sure_snp[, cDNA.K562.B2 := round(cDNA.K562.B2 / totals$cDNA.K562.B2[i] * 1e9, digits = 1)]
  sure_snp[, cDNA.K562.B3 := round(cDNA.K562.B3 / totals$cDNA.K562.B3[i] * 1e9, digits = 1)]
  sure_snp[, cDNA.HepG2.B1 := round(cDNA.HepG2.B1 / totals$cDNA.HepG2.B1[i] * 1e9, digits = 1)]
  sure_snp[, cDNA.HepG2.B2 := round(cDNA.HepG2.B2 / totals$cDNA.HepG2.B2[i] * 1e9, digits = 1)]
  
  # average over replicates and normalize for ipcr counts
  sure_snp[, cDNA.hNPC.norm.ipcr := (cDNA.hNPC.B1 + cDNA.hNPC.B2) / 2 / count]
  sure_snp[, cDNA.K562.norm.ipcr := (cDNA.K562.B1 + cDNA.K562.B2 + cDNA.K562.B3) / 3 / count]
  sure_snp[, cDNA.HepG2.norm.ipcr := (cDNA.HepG2.B1 + cDNA.HepG2.B2) / 2 / count]
  
  
  ### OUT ###
  sure_snp[, c("count", "cDNA.hNPC.B1","cDNA.hNPC.B2", "cDNA.K562.B1", 
               "cDNA.K562.B2", "cDNA.K562.B3", "cDNA.HepG2.B1", 
               "cDNA.HepG2.B2") := NULL]
  
  append_sure <- i != 1
  fwrite(sure_snp, out_sure, sep = "\t", row.names = FALSE, quote = FALSE, append = append_sure)
  NULL
})
  
rm(sure_snp)
gc()

######### CREATE SUMMARY STATS
sander_log("reading back in all values for all snps for all libraries")
counts.dt <- fread(out_sure)

sander_log("rows: ", nrow(counts.dt))

# order by SNP
setorder(counts.dt, SNP_ID)             

# get those that have reads for REF and ALT; additionally, remove all non-REF/ALT
counts.dt <- counts.dt[, if(all(0:1 %in% SNPvarInf)) .SD, SNP_ID]
snp_indel.ids.vector <- unique(counts.dt$SNP_ID)
sander_log("unique snps: ", length(snp_indel.ids.vector))

#### RESULTS
# create empty results matrix
colnames <- c("SNP_ID","ref.seq", "alt.seq", "chrom", "pos.hg19", "ref.element.count", "alt.element.count",
              "ref.hNPC.active.element.count", "alt.hNPC.active.element.count",
              "ref.K562.active.element.count", "alt.K562.active.element.count",
              "ref.HepG2.active.element.count", "alt.HepG2.active.element.count", 
              "hNPC.cDNA.ref.mean","hNPC.cDNA.alt.mean", "hNPC.cDNA.ref.median", "hNPC.cDNA.alt.median", 
              "K562.cDNA.ref.mean","K562.cDNA.alt.mean", "K562.cDNA.ref.median", "K562.cDNA.alt.median", 
              "HepG2.cDNA.ref.mean","HepG2.cDNA.alt.mean","HepG2.cDNA.ref.median", "HepG2.cDNA.alt.median",
              "hNPC.wilcoxon.pvalue", "K562.wilcoxon.pvalue", "HepG2.wilcoxon.pvalue",
              "hNPC.wilcoxon.pvalue.random", "K562.wilcoxon.pvalue.random", "HepG2.wilcoxon.pvalue.random")

results.indel <- data.frame(matrix(nrow = length(snp_indel.ids.vector), ncol = length(colnames)))
colnames(results.indel) <- colnames
sander_log("created `results.indel`")

first_rows <- which(!duplicated(counts.dt$SNP_ID))
results <- data.table(SNP_ID = counts.dt$SNP_ID[first_rows])
sander_log("created `results`")

if (chr == "X") {
  results[, chrom := "X"]
  } else {
  results[, chrom := counts.dt$chr[first_rows]]
  }
results[, pos.hg19 := counts.dt$SNPabspos[first_rows]]
results[, ref.seq := counts.dt[SNPvarInf == 0, SNPbaseInf[1], by = SNP_ID][, -1]]
results[, alt.seq := counts.dt[SNPvarInf == 1, SNPbaseInf[1], by = SNP_ID][, -1]]
sander_log("results 1")

# add the number of elements (gDNA fragments) that are present in the library
# we force the integers to numeric to avoid some integer overflow problems in the wilcox test (as integers have a low min/max than doubles)
results[, ref.element.count := as.numeric(counts.dt[SNPvarInf == 0, .N, by = SNP_ID][[2]])]
results[, alt.element.count := as.numeric(counts.dt[SNPvarInf == 1, .N, by = SNP_ID][[2]])]
sander_log("results 2")

# add the number of elements (gDNA fragments) that have a transcript count of > 0
stat_cols <- c("cDNA.hNPC.norm.ipcr", "cDNA.K562.norm.ipcr", "cDNA.HepG2.norm.ipcr")
cell_types <- c("hNPC", "K562", "HepG2")

results[, ref.hNPC.active.element.count  := counts.dt[SNPvarInf == 0, sum(cDNA.hNPC.norm.ipcr > 0),  by = SNP_ID][, -1]]
results[, alt.hNPC.active.element.count  := counts.dt[SNPvarInf == 1, sum(cDNA.hNPC.norm.ipcr > 0),  by = SNP_ID][, -1]]
results[, ref.K562.active.element.count  := counts.dt[SNPvarInf == 0, sum(cDNA.K562.norm.ipcr > 0),  by = SNP_ID][, -1]]
results[, alt.K562.active.element.count  := counts.dt[SNPvarInf == 1, sum(cDNA.K562.norm.ipcr > 0),  by = SNP_ID][, -1]]
results[, ref.HepG2.active.element.count := counts.dt[SNPvarInf == 0, sum(cDNA.HepG2.norm.ipcr > 0), by = SNP_ID][, -1]]
results[, alt.HepG2.active.element.count := counts.dt[SNPvarInf == 1, sum(cDNA.HepG2.norm.ipcr > 0), by = SNP_ID][, -1]]
sander_log("results 3")

# + mean
results[, hNPC.cDNA.ref.mean  := counts.dt[SNPvarInf == 0, mean(cDNA.hNPC.norm.ipcr),  by = SNP_ID][, -1]]
results[, hNPC.cDNA.alt.mean  := counts.dt[SNPvarInf == 1, mean(cDNA.hNPC.norm.ipcr),  by = SNP_ID][, -1]]
results[, K562.cDNA.ref.mean  := counts.dt[SNPvarInf == 0, mean(cDNA.K562.norm.ipcr),  by = SNP_ID][, -1]]
results[, K562.cDNA.alt.mean  := counts.dt[SNPvarInf == 1, mean(cDNA.K562.norm.ipcr),  by = SNP_ID][, -1]]
results[, HepG2.cDNA.ref.mean := counts.dt[SNPvarInf == 0, mean(cDNA.HepG2.norm.ipcr), by = SNP_ID][, -1]]
results[, HepG2.cDNA.alt.mean := counts.dt[SNPvarInf == 1, mean(cDNA.HepG2.norm.ipcr), by = SNP_ID][, -1]]
sander_log("results 4")

# + median
results[, hNPC.cDNA.ref.median  := counts.dt[SNPvarInf == 0, median(cDNA.hNPC.norm.ipcr),  by = SNP_ID][, -1]]
results[, hNPC.cDNA.alt.median  := counts.dt[SNPvarInf == 1, median(cDNA.hNPC.norm.ipcr),  by = SNP_ID][, -1]]
results[, K562.cDNA.ref.median  := counts.dt[SNPvarInf == 0, median(cDNA.K562.norm.ipcr),  by = SNP_ID][, -1]]
results[, K562.cDNA.alt.median  := counts.dt[SNPvarInf == 1, median(cDNA.K562.norm.ipcr),  by = SNP_ID][, -1]]
results[, HepG2.cDNA.ref.median := counts.dt[SNPvarInf == 0, median(cDNA.HepG2.norm.ipcr), by = SNP_ID][, -1]]
results[, HepG2.cDNA.alt.median := counts.dt[SNPvarInf == 1, median(cDNA.HepG2.norm.ipcr), by = SNP_ID][, -1]]  
sander_log("results 5")

# + mean without zero-expressing reads
results[, hNPC.cDNA.ref.mean_excl0  := counts.dt[SNPvarInf == 0, mean(cDNA.hNPC.norm.ipcr[cDNA.hNPC.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, hNPC.cDNA.alt.mean_excl0  := counts.dt[SNPvarInf == 1, mean(cDNA.hNPC.norm.ipcr[cDNA.hNPC.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, K562.cDNA.ref.mean_excl0  := counts.dt[SNPvarInf == 0, mean(cDNA.K562.norm.ipcr[cDNA.K562.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, K562.cDNA.alt.mean_excl0  := counts.dt[SNPvarInf == 1, mean(cDNA.K562.norm.ipcr[cDNA.K562.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, HepG2.cDNA.ref.mean_excl0 := counts.dt[SNPvarInf == 0, mean(cDNA.HepG2.norm.ipcr[cDNA.HepG2.norm.ipcr > 0]), by = SNP_ID][, -1]]
results[, HepG2.cDNA.alt.mean_excl0 := counts.dt[SNPvarInf == 1, mean(cDNA.HepG2.norm.ipcr[cDNA.HepG2.norm.ipcr > 0]), by = SNP_ID][, -1]]
sander_log("results 6")

# + median without zero-expressing reads
results[, hNPC.cDNA.ref.median_excl0  := counts.dt[SNPvarInf == 0, median(cDNA.hNPC.norm.ipcr[cDNA.hNPC.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, hNPC.cDNA.alt.median_excl0  := counts.dt[SNPvarInf == 1, median(cDNA.hNPC.norm.ipcr[cDNA.hNPC.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, K562.cDNA.ref.median_excl0  := counts.dt[SNPvarInf == 0, median(cDNA.K562.norm.ipcr[cDNA.K562.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, K562.cDNA.alt.median_excl0  := counts.dt[SNPvarInf == 1, median(cDNA.K562.norm.ipcr[cDNA.K562.norm.ipcr > 0]),   by = SNP_ID][, -1]]
results[, HepG2.cDNA.ref.median_excl0 := counts.dt[SNPvarInf == 0, median(cDNA.HepG2.norm.ipcr[cDNA.HepG2.norm.ipcr > 0]), by = SNP_ID][, -1]]
results[, HepG2.cDNA.alt.median_excl0 := counts.dt[SNPvarInf == 1, median(cDNA.HepG2.norm.ipcr[cDNA.HepG2.norm.ipcr > 0]), by = SNP_ID][, -1]]
sander_log("finished `results`")

### WILCOXON

# calculate rank
sander_log("calculating ranks")
counts.dt[, r.hNPC := rank(cDNA.hNPC.norm.ipcr), by = SNP_ID]
counts.dt[, r.K562 := rank(cDNA.K562.norm.ipcr), by = SNP_ID]
counts.dt[, r.HepG2 := rank(cDNA.HepG2.norm.ipcr), by = SNP_ID]

# hnpc
sander_log("processing wilcoxon for hnpc")
results[, z := counts.dt[SNPvarInf == 0, sum(r.hNPC), by = SNP_ID][, -1] - support_z_1(ref.element.count, alt.element.count)]
results[, z := (z - sign(z) * 0.5) / support_z_2(ref.element.count, alt.element.count, as.numeric(counts.dt[, fun_tau(r.hNPC), by = SNP_ID][[2]]))]
results[, hNPC.wilcoxon.pvalue := 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE))]

# k562
sander_log("processing wilcoxon for k562")
results[, z := counts.dt[SNPvarInf == 0, sum(r.K562), by = SNP_ID][, -1] - support_z_1(ref.element.count, alt.element.count)]
results[, z := (z - sign(z) * 0.5) / support_z_2(ref.element.count, alt.element.count, as.numeric(counts.dt[, fun_tau(r.K562), by = SNP_ID][[2]]))]
results[, K562.wilcoxon.pvalue := 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE))]

# hepg2
sander_log("processing wilcoxon for hepg2")
results[, z := counts.dt[SNPvarInf == 0, sum(r.HepG2), by = SNP_ID][, -1] - support_z_1(ref.element.count, alt.element.count)]
results[, z := (z - sign(z) * 0.5) / support_z_2(ref.element.count, alt.element.count, as.numeric(counts.dt[, fun_tau(r.HepG2), by = SNP_ID][[2]]))]
results[, HepG2.wilcoxon.pvalue := 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE))]

results[, z := NULL]

### RANDOM WILCOXON

set.seed(1)
counts.dt[, SNPvarInf.random := sample(SNPvarInf), by = SNP_ID]
#counts.dt[, SNPvarInf.random := sample(SNPvarInf)]

# hnpc
sander_log("processing wilcoxon (random) for hnpc")
results[, rz := counts.dt[, sum(r.hNPC[SNPvarInf.random == 0]), by = SNP_ID][, -1] - support_z_1(as.numeric(counts.dt[, sum(SNPvarInf.random == 0), by = SNP_ID][[2]]), 
                                                                                                 as.numeric(counts.dt[, sum(SNPvarInf.random == 1), by = SNP_ID][[2]]))]            
results[, rz := (rz - sign(rz) * 0.5) / support_z_2(as.numeric(counts.dt[, sum(SNPvarInf.random == 0), by = SNP_ID][[2]]), 
                                                    as.numeric(counts.dt[, sum(SNPvarInf.random == 1), by = SNP_ID][[2]]), 
                                                    as.numeric(counts.dt[, fun_tau(r.hNPC), by = SNP_ID][[2]]))]
results[, hNPC.wilcoxon.pvalue.random  := 2 * pmin(pnorm(rz), pnorm(rz, lower.tail = FALSE))]

# k562
sander_log("processing wilcoxon (random) for k562")
results[, rz := counts.dt[, sum(r.K562[SNPvarInf.random == 0]), by = SNP_ID][, -1] - support_z_1(as.numeric(counts.dt[, sum(SNPvarInf.random == 0), by = SNP_ID][[2]]), 
                                                                                                 as.numeric(counts.dt[, sum(SNPvarInf.random == 1), by = SNP_ID][[2]]))]            
results[, rz := (rz - sign(rz) * 0.5) / support_z_2(as.numeric(counts.dt[, sum(SNPvarInf.random == 0), by = SNP_ID][[2]]), 
                                                    as.numeric(counts.dt[, sum(SNPvarInf.random == 1), by = SNP_ID][[2]]), 
                                                    as.numeric(counts.dt[, fun_tau(r.K562), by = SNP_ID][[2]]))]
results[, K562.wilcoxon.pvalue.random  := 2 * pmin(pnorm(rz), pnorm(rz, lower.tail = FALSE))]

# hepg2
sander_log("processing wilcoxon (random) for hepg2")
results[, rz := counts.dt[, sum(r.HepG2[SNPvarInf.random == 0]), by = SNP_ID][, -1] - support_z_1(as.numeric(counts.dt[, sum(SNPvarInf.random == 0), by = SNP_ID][[2]]), 
                                                                                                  as.numeric(counts.dt[, sum(SNPvarInf.random == 1), by = SNP_ID][[2]]))]            
results[, rz := (rz - sign(rz) * 0.5) / support_z_2(as.numeric(counts.dt[, sum(SNPvarInf.random == 0), by = SNP_ID][[2]]), 
                                                    as.numeric(counts.dt[, sum(SNPvarInf.random == 1), by = SNP_ID][[2]]), 
                                                    as.numeric(counts.dt[, fun_tau(r.HepG2), by = SNP_ID][[2]]))]
results[, HepG2.wilcoxon.pvalue.random  := 2 * pmin(pnorm(rz), pnorm(rz, lower.tail = FALSE))]

results[, rz := NULL]

############ WRITE OUTPUT
sander_log("writing output file")
fwrite(results, file.path(p_res, paste0("out_chr", chr, ".txt.gz")), quote = FALSE, row.names = FALSE, sep = "\t")

sander_log("Done.")
