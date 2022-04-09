#!/bin/bash
#
#PBS -N trimgalore
#PBS -l walltime=16:00:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


# Execute the job code

mkdir trimgalore_out
cd trimgalore_out # To return trimmed files to their own folder depsite the absence of an output function in this version of TrimGalore (version 0.4.1).

#Trim parameteres taken from the methods section of the original publication and report file excluded in the output as not required.

for i in /scratch/spectre/s/sesm2/SRP/prin_2_output/*_prin2_1.fastq
do 
fn=$(echo "$(basename "$i")" | awk -F'[_.]' '{print $1}')
perl /scratch/spectre/s/sesm2/SRP/trim-galore-0.4.1-0/bin/trim_galore --nextera --stringency 1 --no_report_file --paired /scratch/spectre/s/sesm2/SRP/prin_2_output/"$fn"_prin2_1.fastq /scratch/spectre/s/sesm2/SRP/prin_2_output/"$fn"_prin2_2.fastq
done

cd .. #To return to working directory.


