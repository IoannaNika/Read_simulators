#!/bin/sh

#SBATCH --partition=general,insy
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --ntasks=30 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --mail-type=END
#SBATCH --output=08_0%j.out # Set name of output log. %j is the Slurm jobIds
#SBATCH --gres=gpu:0 


ref_seq="data/SARS-CoV-2-NC_045513.fa"
# 01_100 02_100 03_50 04_75 05_90 06_95 07_98 08_0 09_0
for mixture in 08_0
do
    hifi_reads=trimmed_${mixture}
    outdir=data/lumc_data/natural_mixtures/$mixture
    bash files_for_analysis.sh $ref_seq $hifi_reads $outdir
    python make_lumc_dataset.py --bed  ${outdir}/${hifi_reads}_reads.bed --fastq ${outdir}/${hifi_reads}.fastq --out ${outdir}/reads.tsv
    python make_lumc_tuples.py --tsv1 ${outdir}/reads.tsv --out ${outdir}/all_test_pairs.tsv --all True 
done