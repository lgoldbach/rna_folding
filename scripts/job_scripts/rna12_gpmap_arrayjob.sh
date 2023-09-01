#!/bin/bash

#SBATCH --job-name=RNA12
#SBATCH --nodes=1
#SBATCH --array=1-51
#SBATCH --time=02-00:00:00
#SBATCH --output=/home/lgold/rna_folding/data/RNA12/alt_alphabets/four_letters/graph4_1/output-%A_%a.out
#SBATCH --cpus-per-task=1

SCRIPT_DIR="/home/lgold/rna_folding/scripts"
INPUT_DIR="/home/lgold/rna_folding/data/RNA12/alt_alphabets/four_letters/genotypes"
OUTPUT_DIR="/home/lgold/rna_folding/data/RNA12/alt_alphabets/four_letters/graph4_1"

echo START: $(date '+%d/%m/%Y %H:%M:%S')

srun python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA12_gt${SLURM_ARRAY_TASK_ID}.txt -o $OUTPUT_DIR/RNA12_ph${SLURM_ARRAY_TASK_ID}.txt -m 3 -s 0 -b JKLM -p 1

echo END: $(date '+%d/%m/%Y %H:%M:%S')
