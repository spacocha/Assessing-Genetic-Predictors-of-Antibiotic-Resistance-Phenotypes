#!/bin/bash -l

#SBATCH
#SBATCH --job-name=trim
#SBATCH --time=12:00:00

#be in rgi_output
#/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/rgi_output
module load python/3.9.0
source /home/sprehei1/scr4_sprehei1/Sarah/python_envs/genotype-phenotype/bin/activate

echo "Starting rgi conversion"
date
python /home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/Assessing-Genetic-Predictors-of-Antibiotic-Resistance-Phenotypes/Scripts/07_rgi/convert_rgi_output.py

echo "Complete"
date


