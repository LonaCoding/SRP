#!/bin/bash
#PBS -N STAR
#PBS -l walltime=08:00:00
#PBS -l vmem=100gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M aagt1@student.le.ac.uk
reference_genome=/scratch/spectre/a/aagt1/GRCh38/star/GRCh38/star
outprefix=/scratch/spectre/a/aagt1/SteeredResearch/star_output/
files=/scratch/spectre/a/aagt1/SteeredResearch/fastq_trimmed
module load star/2.7.9a
mkdir $outprefix/SRR1974635
STAR --runThreadN 16 --genomeDir /scratch/spectre/a/aagt1/GRCh38/star --readFilesIn $files/SRR1974635_R1_paired.fastq $files/SRR1974635_R2_paired.fastq --outReadsUnmapped Fastx --outFileNamePrefix $outprefix/SRR1974635/SRR1974635 --outSAMtype BAM SortedByCoordinate