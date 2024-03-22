#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

srun python make_amplicon_triplets.py --n 327896 --dir data/hcov_global_2023-11-16_09-28_primers_template_ONT --out data/hcov_ONT_327896_triplets --strategy ONT