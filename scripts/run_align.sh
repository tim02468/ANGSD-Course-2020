#! /bin/bash

#SBATCH  --partition=angsd_class
#SBATCH  --nodes=1
#SBATCH  --ntasks=1
#SBATCH  --job-name=angsd_project925
#SBATCH  --time=24:00:00   
#SBATCH  --mem=40G 


for sample in `echo MSDP075 MSDP080 SDPC082 SDPC087`; do ./align_human.sh ./ref_gene/hg38_index ./raw_reads ${sample}; done 
