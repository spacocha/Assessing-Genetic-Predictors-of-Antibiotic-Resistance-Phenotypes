#!/bin/bash -l

#SBATCH
#SBATCH --job-name=checkm
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem

module load mamba/0.23.0
module load anaconda3/2024.02-1

cd /home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/checkm2

conda activate checkm2

ASSEMBLY=`awk "NR==$SLURM_ARRAY_TASK_ID" ../assembly.txt`

echo "Starting checkm2"
date
checkm2 predict --threads 1 --input ../assemblies/${ASSEMBLY}/contigs.fasta --output-directory ${ASSEMBLY}

echo "Complete"
date


