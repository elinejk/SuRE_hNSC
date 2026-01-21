library(dplyr)
library(data.table)
library(GenomicRanges)
library(IRanges)

## LOAD FILES
# interactions
pir <- fread("V:/ddata/CELB/poot/Eline Koornstra/Resources/Annotations/Song2019_NG_TableS2_HiC_excitatory_neurons.txt")
pir$link_id <- seq_len(nrow(pir))


# sure data
hnsc_emvar <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/emVars/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_wilc-raqtls_04042024.txt")
hnsc_emvar$is_emvar <- TRUE

hnsc_control <- fread("V:/ddata/CELB/poot/Eline Koornstra/SuRE_hNSC_project/emVars/freeze7/hnsc_no_downsampling_snp-permutation_freeze7_controls_04042024.txt")
hnsc_control$is_emvar <- FALSE

hnsc <- rbind(hnsc_emvar, hnsc_control) 
hnsc$chrom <- paste0("chr", hnsc$chrom)



## GENOMIC RANGES
hnsc_gr <- GRanges(seqnames = hnsc$chrom,
                   ranges = IRanges(start = hnsc$`pos.hg19`, width = 1),
                   SNP_ID = hnsc$SNP_ID,
                   is_emvar = hnsc$is_emvar)

lhs_gr <- GRanges(seqnames = pir$lhs_chr,
                  ranges = IRanges(start = pir$lhs_start, end = pir$lhs_end),
                  link_id = pir$link_id)

rhs_gr <- GRanges(seqnames = pir$rhs_chr,
                  ranges = IRanges(start = pir$rhs_start, end = pir$rhs_end),
                  link_id = pir$link_id)


## FIND OVERLAPS
ov_lhs <- findOverlaps(hnsc_gr, lhs_gr)
ov_rhs <- findOverlaps(hnsc_gr, rhs_gr)


## GET MATCHES
matches_lhs <- data.frame(
  SNP_ID = mcols(hnsc_gr[queryHits(ov_lhs)])$SNP_ID,
  is_emvar = mcols(hnsc_gr[queryHits(ov_lhs)])$is_emvar,
  link_id = mcols(lhs_gr[subjectHits(ov_lhs)])$link_id,
  hit = "lhs"
)

matches_rhs <- data.frame(
  SNP_ID = mcols(hnsc_gr[queryHits(ov_rhs)])$SNP_ID,
  is_emvar = mcols(hnsc_gr[queryHits(ov_rhs)])$is_emvar,
  link_id = mcols(rhs_gr[subjectHits(ov_rhs)])$link_id,
  hit = "rhs"
)


all_matches <- bind_rows(matches_lhs, matches_rhs)


# add unmatched SNPs
matched_hnsc_list <- unique(all_matches$SNP_ID)

unmatched_hnsc <- hnsc %>%
  filter(!(SNP_ID %in% matched_hnsc_list)) %>%
  mutate(hit = "none", link_id = NA) %>%
  select(SNP_ID, is_emvar, link_id, hit)


# Combine matched + unmatched
full_result <- bind_rows(all_matches, unmatched_hnsc)


# Merge with pir matches
final_results <- merge(full_result, pir, by = "link_id", all.x = TRUE)


## MERGE BACK WITH THE ORIGINAL emVAR FILES
full_emvar <- merge(final_results, hnsc, by = c("SNP_ID", "is_emvar"))



## SUMMARY FOR PIR
# Get a summary
hnsc_summary <- full_result %>%
  mutate(hit_binary = ifelse(hit == "none", "non-hit", "hit")) %>%
  distinct(SNP_ID, is_emvar, hit_binary)

overlap_table <- table(hnsc_summary$is_emvar, hnsc_summary$hit_binary)
print(overlap_table)








