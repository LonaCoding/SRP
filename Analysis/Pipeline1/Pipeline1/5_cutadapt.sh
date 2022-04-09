#!/bin/bash
#
#PBS -N fastqc
#PBS -l walltime=80:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


mkdir cutadapt_out

#extract -a seqs from overrepresented seqs in /1 file and write with prefix -a (python script or_1seqs.py).

#extract -a seqs from overrepresented seqs in /2 file and write with prefix -A (python script or_2seqs.py).

for i in /scratch/spectre/s/sesm2/SRP/prin_t1/*_prin1_1.fastq
do 
fn=$(echo "$(basename "$i")" | awk -F'[_.]' '{print $1}')
python /scratch/spectre/s/sesm2/SRP/or_1seqs.py /scratch/spectre/s/sesm2/SRP/prin_t1/$fn"_prin1_1_fastqc"/fastqc_data.txt
python /scratch/spectre/s/sesm2/SRP/or_2seqs.py /scratch/spectre/s/sesm2/SRP/prin_t1/$fn"_prin1_2_fastqc"/fastqc_data.txt
cutadapt $fn"_aseqs.txt" $fn"_Aseqs.txt" -e 0.15 -m 30 -o /scratch/spectre/s/sesm2/SRP/cutadapt_out/$fn"_cutadapt_1.fastq" -p /scratch/spectre/s/sesm2/SRP/cutadapt_out/$fn"_cutadapt_2.fastq" /scratch/spectre/s/sesm2/SRP/prin_t1/$fn"_prin1_1.fastq" /scratch/spectre/s/sesm2/SRP/prin_t1/$fn"_prin1_2.fastq"
done

# Cutadapt (version 2.0) was used to trim overrepresented sequences (identified by FastQC) specific to each file. Parameters used were as per the original publication.
