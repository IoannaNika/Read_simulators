#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

for sample in 03_50 #09_0 07_98 06_95 05_90 04_75 03_50 02_100 01_100 
do 
    for seed in 100
    do 
        echo "python subsample_and_make_all_pairs_lumc.py --input_dir data/lumc_data/natural_mixtures/$sample --n 1000 --seed $seed --gr_start 25712"
        python subsample_and_make_all_pairs_lumc.py --input_dir data/lumc_data/natural_mixtures/$sample --n 1000 --seed $seed --gr_start 25712
    done
done   