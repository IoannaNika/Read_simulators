#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=144:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

srun python make_amplicon_triplets.py --n 100000 --dir data/Zika/Gisaid_dataset_2024_02_12_12 --out data/Zika/Zika_2024_02_12_12_100000