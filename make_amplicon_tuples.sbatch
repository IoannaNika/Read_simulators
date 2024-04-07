#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048
#SBATCH --gres=gpu:0 

# srun python make_amplicon_tuples.py --n 500000 --dir data/hcov_global_2023-11-16_09-28 --out data/dataset_500000
# srun python make_amplicon_tuples.py --n 300 --dir data/test_datasets/test_set_BA.2 --out data/test_datasets/test_set_BA.2/dataset
# srun python make_amplicon_tuples.py --n 2000 --dir data/HCV-1b-NCBI/subselection_95/subselection_95 --out data/HCV-1b-NCBI/subselection_95/dataset

# echo "srun python make_amplicon_tuples.py --n 655817 --dir data/hcov_global_2023-11-16_09-28_primers_template_rev_compl --out tuples_pacbio_sars_cov_2_rev_compl/dataset --strategy pacbio-hifi"
# srun python make_amplicon_tuples.py --n 655817 --dir data/hcov_global_2023-11-16_09-28_primers_template_rev_compl --out data/tuples_pacbio_sars_cov_2_rev_compl/dataset --strategy pacbio-hifi


echo "srun python make_amplicon_tuples.py --n 655817 --dir data/hcov_global_2023-11-16_09-28_primers_template_rev_compl --out data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset --strategy pacbio-hifi --append_mode True"
srun python make_amplicon_tuples.py --n 655817 --dir data/hcov_global_2023-11-16_09-28_primers_template_rev_compl --out data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset --strategy pacbio-hifi --append_mode True --dup_pos True

# echo "python make_amplicon_tuples.py --n 655817 --dir data/hcov_global_2023-11-16_09-28_primers_template_ONT_rev_compl --out data/tuples_ONT_sars_cov_2_rev_compl/dataset --strategy ONT"
# srun python make_amplicon_tuples.py --n 655817 --dir data/hcov_global_2023-11-16_09-28_primers_template_ONT_rev_compl --out data/tuples_ONT_sars_cov_2_rev_compl/dataset --strategy ONT
