#!/bin/bash 
#PBS -N prefetch 
#PBS -l walltime=08:00:00
#PBS -l vmem=32gb 
#PBS -l nodes=1:ppn=16 
#PBS -m bea 
#PBS -M aagt1@student.le.ac.uk 

export OMP_NUM_THREADS=$PBS_NUM_PPN 

scratchdir=/scratch/spectre/a/aagt1/SteeredResearch/test_data/batch_sra
sra_accessions=$scratchdir/SRR_Acc_List.txt
prefetch_command=/cm/shared/apps/sratoolkit/2.11.1/bin/prefetch

module load sratoolkit/2.11.1

$prefetch_command --option-file $sra_accessions


