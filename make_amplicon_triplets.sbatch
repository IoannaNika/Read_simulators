#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

srun python make_amplicon_triplets.py --n 500000 --dir data/hcov_global_2023-11-16_09-28 --out data/triplet_dataset_test_500000
