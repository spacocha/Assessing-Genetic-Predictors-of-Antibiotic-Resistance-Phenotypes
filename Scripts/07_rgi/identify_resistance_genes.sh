wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json card.json --local

#!/bin/sh
#SBATCH --job-name=prokka
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-198

ml anaconda3/2022.05
ml mamba/0.23.0
source ~/.bashrc
conda env list
mamba activate rgi

TARGET_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /working_dir/work/prokka_target.txt)
mkdir output/${TARGET_NAME}

rgi main --input_sequence /working_dir/work/spades_results/$TARGET_NAME/contigs.fasta --output_file output/${TARGET_NAME}/${TARGET_NAME}_rgi_output --local --clean -a DIAMOND --low_quality --include_nudge
