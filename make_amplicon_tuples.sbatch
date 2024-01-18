#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

# srun python make_amplicon_tuples.py --n 500000 --dir data/hcov_global_2023-11-16_09-28 --out data/dataset_500000
srun python make_amplicon_tuples.py --n 400 --dir data/test_datasets/pox_B.1.17 --out data/test_datasets/pox_B.1.17/dataset
