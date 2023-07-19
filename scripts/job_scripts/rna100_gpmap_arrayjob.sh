#!/bin/bash

#SBATCH --job-name=RNA100
#SBATCH --nodes=1
#SBATCH --array=1-50
#SBATCH --time=02-00:00:00
#SBATCH --output=/home/lgold/rna_folding/data/RNA100/output/output-%A_%a.out
#SBATCH --cpus-per-task=1

SCRIPT_DIR="/home/lgold/rna_folding/scripts"
INPUT_DIR="/home/lgold/rna_folding/data/RNA100/input"
OUTPUT_DIR="/home/lgold/rna_folding/data/RNA100/output"

echo START: $(date '+%d/%m/%Y %H:%M:%S')

srun python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA100_gt_sample.txt.${SLURM_ARRAY_TASK_ID} -o $OUTPUT_DIR/RNA100_ph.txt.${SLURM_ARRAY_TASK_ID} -m 3 -s 0 -z 100

echo END: $(date '+%d/%m/%Y %H:%M:%S')
