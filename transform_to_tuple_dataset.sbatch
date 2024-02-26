#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

python transform_to_tuple_dataset.py --tsv data/triplet_dataset_primers_template_excl_neg_identicals_327896/triplets.tsv --out data/triplet_dataset_primers_template_excl_neg_identicals_327896/samples.tsv