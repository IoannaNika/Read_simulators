#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

python make_lumc_dataset.py --bed data/lumc_data/Wuhan/aln_wuhan_like_1.sorted.bam.bed --fastq  data/lumc_data/Wuhan/trimmed_01_100.fastq --out data/lumc_data/Wuhan/combined.tsv
python make_lumc_dataset.py --bed data/lumc_data/Omicron/aln_omicron_1.sorted.bam.bed --fastq  data/lumc_data/Omicron/trimmed_08_0.fastq --out data/lumc_data/Omicron/combined.tsv
python make_lumc_tuples.py --tsv1 data/lumc_data/Wuhan/combined.tsv --tsv2 data/lumc_data/Omicron/combined.tsv --out data/lumc_data/lumc_dataset.tsv --n 10000