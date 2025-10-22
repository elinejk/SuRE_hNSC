#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200GB
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.vanzanten.1@erasmusmc.nl

cd /gpfs/home6/evzanten/jaspar_enrichment/

tfbs=$1
setone="/gpfs/home6/evzanten/jaspar_enrichment/input/hg38/020524_${tfbs}_raQTL_100bp_hg38_JASPAR.bed"
settwo="/gpfs/home6/evzanten/jaspar_enrichment/input/hg38/020524_${tfbs}_control_100bp_hg38_JASPAR.bed"

module load 2021
module load R/4.1.0-foss-2021a
module load Python/3.9.5-GCCcore-10.3.0
module load cairo/1.16.0-GCCcore-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0


mkdir output

#emVar over control enrichment
mkdir output/twosets_raqtl_enr_100bp_"$tfbs"
bin/JASPAR_enrich.sh twoSets hg38 $setone $settwo 020524_twosets_raqtl_enr_100bp_"$tfbs" 2022 2 

#control over emVar enrichment
mkdir output/twosets_ctrl_enr_100bp_"$tfbs"
bin/JASPAR_enrich.sh twoSets hg38 $settwo $setone 020524_twosets_ctrl_enr_100bp_"$tfbs" 2022 2 


 
