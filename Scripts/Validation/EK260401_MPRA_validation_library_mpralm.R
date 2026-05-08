library(data.table)
library(mpra)

## PATHS ##
# the files below contain the pDNA and cDNA values normalized to reads per million, each biological replicate is a column, SNP identifiers are in the format barcode1_SNP1_allele1
pdna <- "/projects/0/AdamsLab/Projects/sure/validation/mpra_library_pdna.txt" 
cdna <- "/projects/0/AdamsLab/Projects/sure/validation/mpra_library_cdna.txt" 

## PREPARE FILES ##
dna <- fread(pdna)
rna <- fread(cdna)

dna <- data.frame(dna, row.names = 1) 
rna <- data.frame(rna, row.names = 1) 


## INFORMATION ##
E <- 54 # number of SNPs
s <- 3 # number of biological replicates
nalleles <-2 # number of alleles per SNP

## MAKE AN MPRASET ##
# Aggregate counts for barcodes corresponding to the same element and allele
agg_output <- lapply(seq_len(E), function(elem_id) {
    pattern1 <- paste0(paste0("SNP", elem_id), "_allele1")
    bool_rna_allele1 <- grepl(pattern1, rownames(rna))
    pattern2 <- paste0(paste0("SNP", elem_id), "_allele2")
    bool_rna_allele2 <- grepl(pattern2, rownames(rna))
    agg_rna <- c(
        colSums(rna[bool_rna_allele1,]),
        colSums(rna[bool_rna_allele2,])
    )
    names(agg_rna) <- paste0(rep(c("allele1", "allele2"), each = s), "_", names(agg_rna))
    bool_dna_allele1 <- grepl(pattern1, rownames(dna))
    bool_dna_allele2 <- grepl(pattern2, rownames(dna))
    agg_dna <- c(
        colSums(dna[bool_dna_allele1,]),
        colSums(dna[bool_dna_allele2,])
    )
    names(agg_dna) <- paste0(rep(c("allele1", "allele2"), each = s), "_", names(agg_dna))
    list(rna = agg_rna, dna = agg_dna)
})

agg_rna <- do.call(rbind, lapply(agg_output, "[[", "rna"))
agg_dna <- do.call(rbind, lapply(agg_output, "[[", "dna"))
eid <- paste0("SNP", seq_len(E))
rownames(agg_rna) <- eid
rownames(agg_dna) <- eid

# make the MPRAset
mpraset <- MPRASet(DNA = agg_dna, RNA = agg_rna, eid = eid, barcode = NULL)

## RUN MPRALM ##
# create the design matrix
design <- data.frame(intcpt = 1, alleleB = grepl("allele2", colnames(mpraset)))

# create the block vector, indicates which column is which allele
block_vector <- rep(1:3, 2)

# run mpralm
mpralm_allele <- mpralm(object = mpraset, design = design, aggregate = "none", normalize = FALSE, block = block_vector, 
model_type = "corr_groups", plot = FALSE, endomorphic = TRUE, coef=2)

# get all results in table
rownames(mpralm_allele) <- paste0("SNP_", seq_len(nrow(mpralm_allele)))
rowData(mpralm_allele)

results <- rowData(mpralm_allele)
results <- as.data.table(results)

## SAVE RESULTS ##
fwrite(results, file = "/projects/0/AdamsLab/Projects/sure/validation/mpralm_allele_results.txt", quote = F, row.names = T, sep = '\t')
