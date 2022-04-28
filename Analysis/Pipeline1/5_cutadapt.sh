#!/bin/bash
#
#PBS -N fastqc
#PBS -l walltime=5:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


mkdir cutadapt_out

for i in /scratch/spectre/s/sesm2/SRP/prin_1_output/*_prin1_1_fastqc
do 
fn=$(echo "$(basename "$i")" | awk -F'[_.]' '{print $1}')

#extract -a seqs from overrepresented seqs in /1 file and write fasta file (python script or_1.py).
python /scratch/spectre/s/sesm2/SRP/or_1.py /scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_1_fastqc"/fastqc_data.txt

aseqs=/scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_1_fastqc"/fastqc_data_aseqs.txt

#extract -a seqs from overrepresented seqs in /2 file and write fasta file (python script or_2.py).
python /scratch/spectre/s/sesm2/SRP/or_2.py /scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_2_fastqc"/fastqc_data.txt

Aseqs=/scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_2_fastqc"/fastqc_data_Aseqs.txt

#echo $aseqs $Aseqs #test

cutadapt -a file:$aseqs -A file:$Aseqs -e 0.15 -m 30 -o /scratch/spectre/s/sesm2/SRP/cutadapt_out/$fn"_cutadapt_1.fastq" -p /scratch/spectre/s/sesm2/SRP/cutadapt_out/$fn"_cutadapt_2.fastq" /scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_1.fastq" /scratch/spectre/s/sesm2/SRP/prin_1_output/$fn"_prin1_2.fastq"
done
# Cutadapt (version 2.0) was used to trim overrepresented sequences (identified by FastQC) specific to each file. Parameters used were as per the original publication.


