#!/bin/bash

#SBATCH --job-name=RNA12
#SBATCH --nodes=1
#SBATCH --array=1-49
#SBATCH --time=02-00:00:00
#SBATCH --output=/home/lgold/rna_folding/data/RNA12/output/output-%A_%a.out
#SBATCH --cpus-per-task=1

SCRIPT_DIR="/home/lgold/rna_folding/scripts"
INPUT_DIR="/home/lgold/rna_folding/data/RNA12/input"
OUTPUT_DIR="/home/lgold/rna_folding/data/RNA12/output"

echo START: $(date '+%d/%m/%Y %H:%M:%S')

srun python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA12_gt.txt.${SLURM_ARRAY_TASK_ID} -o $OUTPUT_DIR/RNA12_ph.txt.${SLURM_ARRAY_TASK_ID} -m 3 -s 0

echo END: $(date '+%d/%m/%Y %H:%M:%S')
