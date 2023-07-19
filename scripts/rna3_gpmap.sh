#!/bin/bash

#SBATCH --job-name=RNA12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00-01:00:00

SCRIPT_DIR="/home/lgold/rna_folding/scripts"
INPUT_DIR="/home/lgold/rna_folding/data/RNA3/input"
OUTPUT_DIR="/home/lgold/rna_folding/data/RNA3/output"

# Your script goes here
srun --ntasks=1 python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA3_gt.txt -o $OUTPUT_DIR/RNA3_ph.txt -m 1 -s 0
