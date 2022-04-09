#!/bin/bash 
#PBS -N fastqdump 
#PBS -l walltime=08:00:00
#PBS -l vmem=32gb 
#PBS -l nodes=1:ppn=16 
#PBS -m bea 
#PBS -M aagt1@student.le.ac.uk 

export OMP_NUM_THREADS=$PBS_NUM_PPN 

SCRATCH=/scratch/spectre/a/aagt1/SteeredResearch/batch_sra/
outdir=/scratch/spectre/a/aagt1/SteeredResearch/batch_fastq
fastqdump=/cm/shared/apps/sratoolkit/2.11.1/bin/fastq-dump


module load sratoolkit/2.11.1

for f in $SCRATCH/*

do 
	$fastqdump --split-3 $f --outdir $outdir 
done



