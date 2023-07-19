#!/bin/bash

#SBATCH --job-name=RNA12
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00-01:00:00

SCRIPT_DIR="/home/lgold/rna_folding/scripts"
INPUT_DIR="/home/lgold/rna_folding/data/RNA12/input"
OUTPUT_DIR="/home/lgold/rna_folding/data/RNA12/output"

# Your script goes here
srun --ntasks=1 python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA12_gt_100_1.txt -o $OUTPUT_DIR/RNA12_ph_100_1.txt -m 3 -s 0
srun --ntasks=1 python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA12_gt_100_2.txt -o $OUTPUT_DIR/RNA12_ph_100_2.txt -m 3 -s 0
srun --ntasks=1 python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA12_gt_100_3.txt -o $OUTPUT_DIR/RNA12_ph_100_3.txt -m 3 -s 0
srun --ntasks=1 python $SCRIPT_DIR/build_gp_map.py -f $INPUT_DIR/RNA12_gt_100_4.txt -o $OUTPUT_DIR/RNA12_ph_100_4.txt -m 3 -s 0
