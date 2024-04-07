#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

echo "srun python make_amplicon_triplets.py --n 327896 --dir data/hcov_global_2023-11-16_09-28_primers_template_rev_compl --out data/triplets_pacbio_sars_cov_2_rev_compl/dataset --strategy  pacbio-hifi"
srun python make_amplicon_triplets.py --n 327896 --dir data/hcov_global_2023-11-16_09-28_primers_template_rev_compl --out data/triplets_pacbio_sars_cov_2_rev_compl/dataset --strategy  pacbio-hifi