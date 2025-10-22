#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.vanzanten.1@erasmusmc.nl

cd /gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/



#hNSC
zcat "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/input/310524_1kg_intermediate.txt.gz" | grep -Fwf "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/input/hnsc_raqtl_rsids.txt" > 310524_hNSC.raqtls_annot.txt 


zcat "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/input/310524_1kg_intermediate.txt.gz" | grep -Fwf "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/input/hnsc_control_rsids.txt" > 310524_hNSC_controls.annot.txt

zcat "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/input/310524_1kg_intermediate.txt.gz" |  grep -Fwf "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/input/hnsc_all_snps.txt" > 310524_hNSC_allSNPs.annot.txt


#Just fetch MAF
grep "MAF=" 310524_hNSC_raqtls_annot.txt  | awk -F '[ ;]' '{for(i=1;i<=NF;i++) if ($i ~ /^MAF=/) print $1, substr($i, 5)}' > "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/output/310524_hNSC_raqtls_MAF.txt"

grep "MAF=" 310524_hNSC_controls.annot.txt  | awk -F '[ ;]' '{for(i=1;i<=NF;i++) if ($i ~ /^MAF=/) print $1, substr($i, 5)}' > "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/output/310524_hNSC_controls_MAF.txt"

grep "MAF=" 310524_hNSC_allSNPs.annot.txt  | awk -F '[ ;]' '{for(i=1;i<=NF;i++) if ($i ~ /^MAF=/) print $1, substr($i, 5)}' > "/gpfs/work4/0/AdamsLab/Projects/sure/transcription_factors/SNP_lists_MAF/output/310524_hNSC_allSNPs_MAF.txt"




