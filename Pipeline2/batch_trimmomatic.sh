#!/bin/bash 

#PBS -N trimmomatic
#PBS -l walltime=16:00:00
#PBS -l vmem=100gb
#PBS -l nodes=1:ppn=16
#PBS -m bea
#PBS -M aagt1@student.le.ac.uk



files=/scratch/spectre/a/aagt1/SteeredResearch/batch_fastq
java_file=/cm/shared/apps/java/11.0.2/bin/java

module load java/11.0.2

for i in "$files"/* 
do
	
	if [[ $i == *_1.fastq ]]
		
	then
		prefix="${i::-8}"
		input_1="$i"
		input_2="${i::-8}_2.fastq"
		output_1=$prefix"_R1_paired.fastq"
		output_2=$prefix"_R1_unpaired.fastq"
		output_3=$prefix"_R2_paired.fastq"
		output_4=$prefix"_R2_unpaired.fastq"

		$java_file -jar /home/a/aagt1/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 -trimlog trimlog.txt $input_1 $input_2 $output_1 $output_2 $output_3 $output_4 ILLUMINACLIP:/home/a/aagt1/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:36:10 LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20
		
	else
		:
	fi
done




		
