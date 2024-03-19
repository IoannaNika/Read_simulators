#!/bin/sh
#SBATCH --partition=general,insy
#SBATCH --qos=short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 


python make_lumc_dataset.py --bed data/cerbaresearch/Pacbio-nsp12-fragment/P0835SEQ00003-HDF6.bc1005--bc1035_reads.bed --fastq  data/cerbaresearch/Pacbio-nsp12-fragment/P0835SEQ00003-HDF6.bc1005--bc1035.fastq --out data/cerbaresearch/Pacbio-nsp12-fragment/grouped_reads.tsv
