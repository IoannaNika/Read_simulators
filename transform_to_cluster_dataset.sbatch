#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

# python transform_to_cluster_dataset.py --tsv data/2022_val_dataset/samples.tsv --out  data/2022_val_dataset/single_samples.tsv
python transform_to_cluster_dataset.py --tsv data/HCV-1b-NCBI/subselection_95/dataset/samples.tsv --out data/HCV-1b-NCBI/subselection_95/dataset/single_samples.tsv
