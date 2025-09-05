#!/bin/bash -l

#SBATCH
#SBATCH --job-name=trim
#SBATCH --time=12:00:00

module load mamba/0.23.0
module load anaconda3/2024.02-1

cd /home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype

conda activate rgi_env

conda activate rgi
ASSEMBLY=`awk "NR==$SLURM_ARRAY_TASK_ID" assembly.txt`

echo "Starting rgi"
date

mkdir ./rgi_output/${ASSEMBLY}/
rgi main -i ./assemblies/${ASSEMBLY}/contigs.fasta --output_file ./rgi_output/${ASSEMBLY}/${ASSEMBLY}_rgi_output --local --clean -a DIAMOND --low_quality --include_nudge

echo "Complete"
date


