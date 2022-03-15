#!/bin/bash
#
#PBS -N STAR_2.4.0h
#PBS -l walltime=00:40:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk

STAR=/scratch/spectre/s/sesm2/SRP/STAR-STAR_2.4.0h/bin/Linux_x86_64_static/STAR

# Execute the job code

$STAR --runMode genomeGenerate --runThreadN 16 --genomeDir /scratch/spectre/s/sesm2/SRP/genome_directory --genomeFastaFiles /scratch/spectre/s/sesm2/SRP/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --sjdbGTFfile /scratch/spectre/s/sesm2/SRP/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 74 --genomeChrBinNbits 12
