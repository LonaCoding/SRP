#!/bin/bash
#
#PBS -N Prin1
#PBS -l walltime=80:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


# Execute the job code

mkdir /scratch/spectre/s/sesm2/SRP/prin_1_output


for i in /scratch/spectre/s/sesm2/SRP/SRA_dump/*
do
perl /scratch/spectre/s/sesm2/SRP/prinseq-lite-0.20.4/prinseq-lite.pl -fastq /scratch/spectre/s/sesm2/SRP/corrected_format_SRA/"$(basename "$i")"_1_adp.fastq -fastq2 /scratch/spectre/s/sesm2/SRP/corrected_format_SRA/"$(basename "$i")"_2_adp.fastq -min_len 30 -trim_left 10 -trim_qual_right 25 -lc_method entropy -lc_threshold 65 -out_format 3 -out_good /scratch/spectre/s/sesm2/SRP/prin_1_output/"$(basename "$i")"_prin1 -out_bad null
done

# Use the name of directories pulled from NCBI to provide a prefix for all files to be processed in Prinseq (prinseq-lite-0.20.4).
# The above parameters "-min_len 30 -trim_left 10 -trim_qual_right 25 -lc_method entropy -lc_threshold 65" were all stipulated in the original paper and have been adhered to.
