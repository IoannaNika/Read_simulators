#!/bin/sh
#SBATCH --partition=general,insy
#SBATCH --qos=short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 


python make_lumc_dataset.py --bed data/cerbaresearch/Pacbio-nsp12-fragment/P0835SEQ00003-HDF6.bc1005--bc1035_reads.bed --fastq  data/cerbaresearch/Pacbio-nsp12-fragment/P0835SEQ00003-HDF6.bc1005--bc1035.fastq --out data/cerbaresearch/Pacbio-nsp12-fragment/reads_1005.tsv

# python make_lumc_dataset.py --fastq data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/P0835SEQ00003-HDF6.bc1008--bc1035.fastq --bed data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/P0835SEQ00003-HDF6.bc1008--bc1035_reads.bed --out Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/reads_1008.tsv
# python make_lumc_dataset.py --fastq data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/P0835SEQ00003-HDF6.bc1012--bc1035.fastq --bed data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/P0835SEQ00003-HDF6.bc1012--bc1035_reads.bed --out Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/reads_1012.tsv
