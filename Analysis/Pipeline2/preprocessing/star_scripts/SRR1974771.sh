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
mkdir $outprefix/SRR1974771
STAR --runThreadN 16 --genomeDir /scratch/spectre/a/aagt1/GRCh38/star --readFilesIn $files/SRR1974771_R1_paired.fastq $files/SRR1974771_R2_paired.fastq --outReadsUnmapped Fastx --outFileNamePrefix $outprefix/SRR1974771/SRR1974771 --outSAMtype BAM SortedByCoordinate