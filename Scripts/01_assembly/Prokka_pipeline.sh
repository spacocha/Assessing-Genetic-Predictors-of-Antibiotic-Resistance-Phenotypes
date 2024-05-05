#!/bin/sh
#SBATCH --job-name=prokka
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-198
ml prokka/1.14.5

TARGET_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /working_dir/work/prokka_target.txt)
prokka  --outdir  /working_dir/work/prokka_result/$TARGET_NAME  /working_dir/work/spades_results/$TARGET_NAME/contigs.fasta
