#!/bin/bash -l

#SBATCH
#SBATCH --job-name=trim
#SBATCH --time=12:00:00

# Define some variables
# Input sequences
STRAIN=`awk "NR==$SLURM_ARRAY_TASK_ID" strain.txt`
ACC=`awk "NR==$SLURM_ARRAY_TASK_ID" accession.txt`

echo ${ACC}
echo ${STRAIN}
FWD_FQ="/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/NCBI_prefetch_207_isolates/${ACC}_1.fastq"
REV_FQ="/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/NCBI_prefetch_207_isolates/${ACC}_2.fastq"

echo ${FWD_FQ}
echo ${REV_FQ}

# Output sequences
FORWARD_PAIRED_OUTPUT_FILE_NAME=/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/trim_results/${STRAIN}_1_pe.fastq
FORWARD_UNPAIRED_OUTPUT_FILE_NAME=/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/trim_results/${STRAIN}_1_se.fastq
REVERSE_PAIRED_OUTPUT_FILE_NAME=/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/trim_results/${STRAIN}_2_pe.fastq
REVERSE_UNPAIRED_OUTPUT_FILE_NAME=/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/trim_results/${STRAIN}_2_se.fastq

# Output folders
SPADES_OUTPUT_FOLDER=/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/assemblies/${STRAIN}
QUAST_OUTPUT_FOLDER=/home/sprehei1/scr4_sprehei1/Sarah/Genotype-Phenotype/quast/${STRAIN}

TRIMMOMATIC_ADAPTER=/data/apps/extern/anaconda/envs/trimmomatic/0.39_conda/share/trimmomatic/adapters/TruSeq3-SE.fa

# Load module
ml trimmomatic/0.39_conda

# Trimomatic inputs and outputs, parameters
trimmomatic PE -phred33 ${FWD_FQ} ${REV_FQ} $FORWARD_PAIRED_OUTPUT_FILE_NAME $FORWARD_UNPAIRED_OUTPUT_FILE_NAME $REVERSE_PAIRED_OUTPUT_FILE_NAME $REVERSE_UNPAIRED_OUTPUT_FILE_NAME ILLUMINACLIP:${TRIMMOMATIC_ADAPTER}:2:30:10 SLIDINGWINDOW:4:33 LEADING:10 TRAILING:10 MINLEN:100 HEADCROP:18 TAILCROP:9 AVGQUAL:10
# Run Spades </directory/spades.py> inputs, parameters, output
~/scr4_sprehei1/Xingyou/environment/SPAdes-3.15.5-Linux/bin/spades.py     --pe-1 1 $FORWARD_PAIRED_OUTPUT_FILE_NAME     --pe-2 1 $REVERSE_PAIRED_OUTPUT_FILE_NAME --pe-s 1 $FORWARD_UNPAIRED_OUTPUT_FILE_NAME --pe-s 1 $REVERSE_UNPAIRED_OUTPUT_FILE_NAME    -k 21,25,29     --sc     --cov-cutoff auto     --phred-offset 33     -o $SPADES_OUTPUT_FOLDER
# Run Quast </directory/quast.py> input, output
~/scr4_sprehei1/Xingyou/environment/quast-5.2.0/quast.py $SPADES_OUTPUT_FOLDER/scaffolds.fasta  -o $QUAST_OUTPUT_FOLDER

