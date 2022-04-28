#!/bin/bash
#
#PBS -N HTSeq
#PBS -l walltime=00:30:00
#PBS -l vmem=60gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M sesm2@student.le.ac.uk


# Execute the job code

mkdir /scratch/spectre/s/sesm2/SRP/HTSeq_out

for i in /scratch/spectre/s/sesm2/SRP/STAR_out/*_STAR_Aligned.out.sam
do 
python /home/s/sesm2/anaconda3/envs/py2/bin/htseq-count -m intersection-nonempty \-s no -c /scratch/spectre/s/sesm2/SRP/HTSeq_out/"$fn".tsv -n 6 /scratch/spectre/s/sesm2/SRP/STAR_out/"$fn"_STAR_Aligned.out.sam /scratch/spectre/s/sesm2/SRP/Homo_sapiens.GRCh37.75.gtf
done

# HTSeq (version 0.6.1) was used to produce counts. All parameters were adopted from the methods section of the original publication.





