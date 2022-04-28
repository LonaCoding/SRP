#!/bin/bash
#
#PBS -N STAR_2.4.0h
#PBS -l walltime=06:00:00
#PBS -l vmem=100gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk

STAR=/scratch/spectre/s/sesm2/SRP/STAR-STAR_2.4.0h/bin/Linux_x86_64_static/STAR

# Execute the job code

mkdir /scratch/spectre/s/sesm2/SRP/STAR_out

for i in /scratch/spectre/s/sesm2/SRP/trimgalore_out/*_prin2_1_val_1.fq
do 
fn=$(echo "$(basename "$i")" | awk -F'[_.]' '{print $1}')
$STAR --runThreadN 4 --genomeDir /scratch/spectre/s/sesm2/SRP/genome_directory --readFilesIn /scratch/spectre/s/sesm2/SRP/trimgalore_out/"$fn"_prin2_1_val_1.fq /scratch/spectre/s/sesm2/SRP/trimgalore_out/"$fn"_prin2_2_val_2.fq --outFileNamePrefix /scratch/spectre/s/sesm2/SRP/STAR_out/"$fn"_STAR_ \-outFilterType BySJout \--outFilterMultimapNmax 20 \-- alignSJoverhangMin 8 \--alignSJDBoverhangMin 1 \--outFilterMismatchNmax 999 \--outFilterMismatchNoverLmax 0.04 \--alignIntronMin 20 \--alignIntronMax 1000000 \--alignMatesGapMax 1000000 \--outSAMstrandField intronMotif
done

# STAR (version 2.4.0h) was used to align each paired set of fastq files to produce a SAM output. All parameters were adopted from the methods section of the original publication.
