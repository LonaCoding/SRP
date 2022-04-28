#!/bin/bash
#
#PBS -N Prin2
#PBS -l walltime=80:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


# Execute the job code

mkdir /scratch/spectre/s/sesm2/SRP/prin_2_output

# Prinseq (prinseq-lite-0.20.4) is used to filter for minimum length (30bp) and remove orphan pairs from paired end files as per the original publication.

for i in /scratch/spectre/s/sesm2/SRP/cutadapt_out/*_cutadapt_1.fastq
do 
fn=$(echo "$(basename "$i")" | awk -F'[_.]' '{print $1}') # File prefixes used from previous output to ensure enitre set of reads is processed and called in their respective pairs.
perl /scratch/spectre/s/sesm2/SRP/prinseq-lite-0.20.4/prinseq-lite.pl -fastq /scratch/spectre/s/sesm2/SRP/cutadapt_out/"$fn"_cutadapt_1.fastq -fastq2 /scratch/spectre/s/sesm2/SRP/cutadapt_out/"$fn"_cutadapt_2.fastq -min_len 30 -out_format 3 -out_good /scratch/spectre/s/sesm2/SRP/prin_2_output/"$fn"_prin2 -out_bad null
done

