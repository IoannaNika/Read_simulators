#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2048
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

srun python make_dataset.sbatch --n 10 --dir data/hcov_global_2023-11-16_09-28 --out data/dataset