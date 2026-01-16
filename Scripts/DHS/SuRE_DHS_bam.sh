#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.vanzanten.1@erasmusmc.nl

module load 2024
module load SAMtools/1.21-GCC-13.3.0

samtools view -b -h -L "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/hnsc.raqtl.hetero.120126.bed" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1.bam" -o "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1_hetero_sign_subset.bam" 

samtools index "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1_hetero_sign_subset.bam" 

#Rep2 raqtl
samtools view -b -h -L "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/hnsc.raqtl.hetero.120126.bed" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2.bam" -o "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2_hetero_sign_subset.bam" 

samtools index "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2_hetero_sign_subset.bam" 


#And for ctrls 
samtools view -b -h -L "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/hnsc.control.hetero.120126.bed" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1.bam" -o "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1_hetero_ctrl_subset.bam" 

samtools index "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep1_hetero_ctrl_subset.bam" 

#Rep2 ctrl
samtools view -b -h -L "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/hnsc.control.hetero.120126.bed" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2.bam" -o "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2_hetero_ctrl_subset.bam" 

samtools index "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/dhs/120126_hNPC_rep2_hetero_ctrl_subset.bam" 

