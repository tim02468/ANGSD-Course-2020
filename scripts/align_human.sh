#! /bin/bash

#SBATCH  --partition=angsd_class
#SBATCH  --nodes=1
#SBATCH  --ntasks=1
#SBATCH  --job-name=angsd_project925
#SBATCH  --time=15:00:00   
#SBATCH  --mem=40G          

# Usage: align_Gierlinski.sh <STARindex_dir> <FASTQ_dir> <SAMPLE NAME>
# Make sure that STAR has been loaded and is accessible via `STAR` and that
# the genome/transcriptome index has already been built.

# Check that we have our command line argument(s)
arg_count=$#
if [ $arg_count -lt 3 ]; then
	echo "Not enough command line arguments. Exiting ..."
	echo "Usage: align_Gierlinski.sh <STARindex_dir> <FASTQ_dir> <SAMPLE NAME>"
	exit
fi 

# Read in arguments
STAR_DIR=$1
FASTQ_DIR=$2
SAMPLE=$3

# Load packages
spack load fastqc
spack load star@2.7.0e
spack load samtools@1.9%gcc@6.3.0


# FILES=`'ls' ${FASTQ_DIR}/${SAMPLE}/*.fastq.gz | paste -s -d , -`
echo "Aligning files for ${SAMPLE}, files:"
# echo $FILES 
STAR --genomeDir ${STAR_DIR}/ --readFilesIn ./${FASTQ_DIR}/${SAMPLE}/${SAMPLE}-R_1.fastq.gz.1 ./${FASTQ_DIR}/${SAMPLE}/${SAMPLE}-R_2.fastq.gz.1 --readFilesCommand zcat --outFileNamePrefix ./alignments/${SAMPLE}. --outSAMtype BAM SortedByCoordinate --runThreadN 4 --twopassMode Basic

samtools index ./alignments/${SAMPLE}.Aligned.sortedByCoord.out.bam
