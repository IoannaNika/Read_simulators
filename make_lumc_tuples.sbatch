#!/bin/sh

#SBATCH --partition=general,insy
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=30 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=END
#SBATCH --gres=gpu:0 

# python make_lumc_dataset.py --bed data/lumc_data/Wuhan/aln_wuhan_like_1.sorted.bam.bed --fastq  data/lumc_data/Wuhan/trimmed_01_100.fastq --out data/lumc_data/Wuhan/combined.tsv
# python make_lumc_dataset.py --bed data/lumc_data/Omicron/aln_omicron_1.sorted.bam.bed --fastq  data/lumc_data/Omicron/trimmed_08_0.fastq --out data/lumc_data/Omicron/combined.tsv
# python make_lumc_tuples.py --tsv1 data/lumc_data/Wuhan/combined.tsv --tsv2 data/lumc_data/Omicron/combined.tsv --out data/lumc_data/lumc_dataset_rev_compl.tsv --n 100000

# echo python make_lumc_tuples.py --tsv1 data/lumc_data/Wuhan/combined.tsv --tsv2 data/lumc_data/Omicron/combined.tsv --out Experiments_data/lumc/dataset/all_test_pairs.tsv --all True --with_prefix True
# python make_lumc_tuples.py --tsv1 data/lumc_data/Wuhan/combined.tsv --tsv2 data/lumc_data/Omicron/combined.tsv --out Experiments_data/lumc/dataset/all_test_pairs.tsv --all True --with_prefix True


# cerbaresearch


# echo "python make_lumc_tuples.py --tsv1 data/cerbaresearch/Pacbio-nsp12-fragment/reads_1005.tsv --out Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1005/dataset/all_test_pairs.tsv --all True" 
# python make_lumc_tuples.py --tsv1 data/cerbaresearch/Pacbio-nsp12-fragment/reads_1005.tsv --out Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1005/dataset/all_test_pairs.tsv --all True 

echo "make_lumc_tuples.py --tsv1 data/cerbaresearch/Pacbio-nsp12-fragment/reads_1007.tsv --out Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1007/dataset/all_test_pairs.tsv --all True" 
python make_lumc_tuples.py --tsv1 data/cerbaresearch/Pacbio-nsp12-fragment/reads_1007.tsv --out Experiments_data/cerbaresearch/Pacbio-nsp12-fragment/1007/dataset/all_test_pairs.tsv --all True 

# echo "make_lumc_tuples.py --tsv1 Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1008/reads_1008.tsv --out Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1008/dataset/all_test_pairs.tsv --all True" 
# python make_lumc_tuples.py --tsv1 Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1008/reads_1008.tsv --out Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1008/dataset/all_test_pairs.tsv --all True 

# echo "make_lumc_tuples.py --tsv1 Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1012/reads_1012.tsv --out Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1012/dataset/all_test_pairs.tsv --all True" 
# python make_lumc_tuples.py --tsv1 Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1012/reads_1012.tsv --out Experiments_data/cerbaresearch/Pacbio-Spike-5-fragment-CDC/1012/dataset/all_test_pairs.tsv --all True 
