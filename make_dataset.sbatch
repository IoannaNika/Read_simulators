#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

srun python make_amplicon_tuples.py --n 50 --dir data/hcov_global_2023-11-16_09-28 --out data/dataset