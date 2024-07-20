#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5GB

python src/calculate_number_of_mutations.py --dir data/tuples_pacbio_sars_cov_2_rev_compl/dataset