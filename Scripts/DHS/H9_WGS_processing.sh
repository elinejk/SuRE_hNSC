#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e.vanzanten.1@erasmusmc.nl

module load 2024
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0

cd /gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/
#First downloaded sra-toolkit https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
#Now download all H9 ESC WGS data
prefetch SRR1091088
fasterq-dump SRR1091088

prefetch SRR1091091
fasterq-dump SRR1091091

prefetch SRR1091092
fasterq-dump SRR1091092

#Align fastq files to hg37 
 #First index files
 for file in *.fastq ; do samtools faidx $file ; done 
 
 #Index ref genome fasta
 bwa index "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/tools/bwa/GRCh37.p13.genome.fa"
 
 #Now align the files to the ref genome
bwa mem -t 16 "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/tools/bwa/GRCh37.p13.genome.fa" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091088_1.fastq" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091088_2.fastq" > "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091088_aligned.sam"
 
bwa mem -t 16 "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/tools/bwa/GRCh37.p13.genome.fa" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091091_1.fastq" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091091_2.fastq" > "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091091_aligned.sam"
 
bwa mem -t 16 "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/tools/bwa/GRCh37.p13.genome.fa" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091092_1.fastq" "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091092_2.fastq" > "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/raw/SRR1091092_aligned.sam"

#Convert sam into bam (binary format) 
samtools view -bS SRR1091088_aligned.sam > SRR1091088_aligned.bam
samtools view -bS SRR1091091_aligned.sam > SRR1091091_aligned.bam
samtools view -bS SRR1091092_aligned.sam > SRR1091092_aligned.bam

#Sort the BAM files
samtools sort SRR1091088_aligned.bam -o SRR1091088_aligned_sorted.bam
samtools sort SRR1091091_aligned.bam -o SRR1091091_aligned_sorted.bam
samtools sort SRR1091092_aligned.bam -o SRR1091092_aligned_sorted.bam

#Index the sorted BAM files
samtools index SRR1091088_aligned_sorted.bam
samtools index SRR1091091_aligned_sorted.bam
samtools index SRR1091092_aligned_sorted.bam


#Merge the 3 runs 
samtools merge 070126_H9_ESC_merged.bam SRR1091088_aligned_sorted.bam SRR1091091_aligned_sorted.bam SRR1091092_aligned_sorted.bam
samtools index 070126_H9_ESC_merged.bam

#Create tarball of raw files to save space
tar -czvf H9_ESC_wgs_raw.tar.gz *

#Add dummy sample information
java -jar "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/tools/picard.jar"  AddOrReplaceReadGroups   I=/gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/070126_H9_ESC_merged.bam   O=/gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/070126_H9_ESC_merged_with_RG.bam   RGID=1   RGLB=lib1   RGPL=illumina   RGPU=unit1   RGSM=H9_ESC

#Call variants of processed bam file w/ gatk
gatk --java-options "-Xmx4g" HaplotypeCaller \
  -R /gpfs/work4/0/AdamsLab/Projects/sure/wgs/tools/ref_fasta/GRCh37.p13.genome.fa \
  -I /gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/070126_H9_ESC_merged_with_RG.bam \
  -O /gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/080126_H9_ESC_GATK.vcf.gz \
  -ERC GVCF \
  -G StandardAnnotation \
  -G AS_StandardAnnotation

#Filter for heterozygous SNPs
bcftools view -i 'GT="0/1"' "/gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/080126_H9_ESC_GATK.vcf.gz" -Oz -o /gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/120126_H9_heterozygous.vcf.gz
tabix -p vcf /gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/120126_H9_heterozygous.vcf.gz 

#Convert into csv to work in R 
bcftools query \
         -f '%CHROM,%POS,%ID,%REF,%ALT,%INFO/AF,%INFO/AC\n' \
         /gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/120126_H9_heterozygous.vcf.gz > /gpfs/work4/0/AdamsLab/Projects/sure/wgs/processed/120126_H9_heterozygous.csv
