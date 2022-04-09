# Replace the .1 .2 extensions with /1 /2 so that the reads are recognised as pairs in Prinseq.
#This is quick so does not need too be submitted as a job for the cluster.

mkdir corrected_format_SRA

for i in SRA_files/* # Loop through list of all 466 accession numbers from the NCBI short read archive (SRA), in project SRP057196, stored to a text file.
do
python3 idadaptor.py SRA_files/$(basename "$i") # Python scipt to substitue the extensions for both the identifier and quality identifier.
mv SRA_files/*adp.fastq /scratch/spectre/s/sesm2/SRP/corrected_format_SRA # Relocate output to be called at next step.
done
