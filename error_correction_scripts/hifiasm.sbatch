#!/bin/sh

#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --gres=gpu:0 
#SBATCH -o Error_corrected_data/hifiasm/lumc/03_50/output.%j.out 

bash error_correction_scripts/hifiasm.bash
echo "ERROR CORRECTION DONE!"
