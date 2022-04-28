#!/bin/bash
#
#PBS -N fastqc
#PBS -l walltime=20:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


# Execute the job code


# Fastqc (FastQC v0.11.2) report on each pair of files from first prinseq step to be returned to prin_1_output. 

for i in /scratch/spectre/s/sesm2/SRP/prin_1_output/*_prin1_1.fastq
do 
fn=$(echo "$(basename "$i")" | awk -F'[_.]' '{print $1}')
#echo $fn"_prin1_1.fastq" $fn"_prin1_2.fastq"
#Test; to confirm prefix is accurate before using it to run the perl script.
perl /scratch/spectre/s/sesm2/SRP/FastQC/fastqc -extract -t 6 /scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_1.fastq" /scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_2.fastq"
done





