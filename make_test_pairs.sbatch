#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 


python make_test_pairs.py --reads Experiments_data/HCV-1b/cov_100x_mixture_pb_hifi/dataset/all_reads.tsv --output Experiments_data/HCV-1b/cov_100x_mixture_pb_hifi/dataset/all_test_pairs.tsv

