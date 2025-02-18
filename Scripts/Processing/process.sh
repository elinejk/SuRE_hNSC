#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22
#SBATCH --mail-type=END
#SBATCH --mail-user=e.koornstra@erasmusmc.nl
#SBATCH --output=/home/ekoornstra/sure_data/freeze7/logs/process_reads-%A_%a.out



Rscript process_raw_wohnpcbr3_snp-permutation.R \
  ${SLURM_ARRAY_TASK_ID} \
  /gpfs/home6/ekoornstra/sure_data/freeze7/
