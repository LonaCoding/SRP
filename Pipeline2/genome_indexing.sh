#!/bin/bash 

#PBS -N ANGELO_TEST
#PBS -l walltime=8:00:00
#PBS -l vmem=32gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M aagt1@student.le.ac.uk

reference_genome=/scratch/spectre/a/aagt1/GRCh38/refdata-gex-GRCh38-2020-A
gtf=$reference_genome/genes/genes.gtf
genome_fasta=$reference_genome/fasta/genome.fa
command=/cm/shared/apps/star/2.7.9a/STAR

export OMP_NUM_THREADS=$PBS_NUM_PPN

module load star/2.7.9a  

$command --runThreadN 8 --runMode genomeGenerate --genomeDir $reference_genome/star --genomeFastaFiles $genome_fasta --sjdbGTFfile $gtf --sjdbOverhang 100

