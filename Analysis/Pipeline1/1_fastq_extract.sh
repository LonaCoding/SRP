#!/bin/bash
#
#PBS -N fastq_extract
#PBS -l walltime=1:00:00
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=2
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


# Execute the job code

#set up both directories
mkdir SRA_dump
mkdir SRA_files

#extract SRA files
sratoolkit.3.0.0-centos_linux64/bin/./prefetch --option-file acc_list -O SRA_dump

#fastq-dump and split files into pairs
for i in SRA_dump/*
do
sratoolkit.3.0.0-centos_linux64/bin/./fastq-dump -I --split-files $(basename "$i") -O SRA_files
done

