#!/bin/sh
#SBATCH --job-name=trimmomatic

# Define some variables
# Input sequences
SEQUENCE_ID_FORWARD_FILE_NAME=/working_dir/target_dir/${1}_R1.fastq
SEQUENCE_ID_REVERSE_FILE_NAME=/working_dir/target_dir/${1}_R2.fastq

# Output sequences
FORWARD_PAIRED_OUTPUT_FILE_NAME=/working_dir/work/trim_results/${1}_FW_PD.fastq
FORWARD_UNPAIRED_OUTPUT_FILE_NAME=/working_dir/work/trim_results/${1}_FW_UPD.fastq
REVERSE_PAIRED_OUTPUT_FILE_NAME=/working_dir/work/trim_results/${1}_RE_PD.fastq
REVERSE_UNPAIRED_OUTPUT_FILE_NAME=/working_dir/work/trim_results/${1}_RE_UPD.fastq

# Output folders
SPADES_OUTPUT_FOLDER=/working_dir/work/spades_results/${1}
QUAST_OUTPUT_FOLDER=/working_dir/work/quast_results/${1}

TRIMMOMATIC_ADAPTER=/data/apps/extern/anaconda/envs/trimmomatic/0.39_conda/share/trimmomatic/adapters/TruSeq3-SE.fa

# Load module
ml trimmomatic/0.39_conda

# Trimomatic inputs and outputs, parameters
trimmomatic PE -phred33 $SEQUENCE_ID_FORWARD_FILE_NAME $SEQUENCE_ID_REVERSE_FILE_NAME $FORWARD_PAIRED_OUTPUT_FILE_NAME $FORWARD_UNPAIRED_OUTPUT_FILE_NAME $REVERSE_PAIRED_OUTPUT_FILE_NAME $REVERSE_UNPAIRED_OUTPUT_FILE_NAME ILLUMINACLIP:$TRIMMOMATIC_ADAPTER:2:30:10 SLIDINGWINDOW:4:33 LEADING:10 TRAILING:10 MINLEN:100 HEADCROP:18 TAILCROP:9 AVGQUAL:10
# Run Spades </directory/spades.py> inputs, parameters, output
/working_dir/environment/SPAdes-3.15.5-Linux/bin/spades.py     --pe1-1 $FORWARD_PAIRED_OUTPUT_FILE_NAME     --pe1-2 $REVERSE_PAIRED_OUTPUT_FILE_NAME     -k 21,25,29     --sc     --cov-cutoff auto     --phred-offset 33     -o $SPADES_OUTPUT_FOLDER
# Run Quast </directory/quast.py> input, output
/working_dir/environment/quast-5.2.0/quast.py $SPADES_OUTPUT_FOLDER/scaffolds.fasta     -o $QUAST_OUTPUT_FOLDER
