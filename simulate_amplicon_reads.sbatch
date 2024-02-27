#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 


python simulate_amplicon_reads.py --dir data/HIV-1-NCBI/subselection_95/NCBI_processed
