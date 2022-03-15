#!/bin/bash 
#PBS -N genome_index
#PBS -l walltime=16:00:00
#PBS -l vmem=32gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M aagt1@student.le.ac.uk

reference_genome=/scratch/spectre/a/aagt1/GRCh38
gtf=$reference_genome/gencode.v39.primary_assembly.annotation.gtf
genome_fasta=$reference_genome/GRCh38.p13.genome.fa

command=/cm/shared/apps/star/2.7.9a/STAR

module load star/2.7.9a 

$command --runThreadN 16 --runMode genomeGenerate --genomeDir $reference_genome/star --genomeFastaFiles $genome_fasta --sjdbGTFfile $gtf --sjdbOverhang 100 


