
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO

def simulate_hifi_reads(directory, id, cores):
    os.system("pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-SEQUEL.model  --template {directory}/{id}.template --pass-num 10 --seed 20 --prefix {directory}/{id}")
    os.system("samtools view -bS {directory}/{id}.sam > {directory}/{id}.bam")
    os.system("ccs {directory}/{id}.bam -j {cores} {directory}/{id}.fastq")
    return

def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--cores', dest = 'cores', default=8, required=False, type=int, help="cpu cores for the simulation")
    args = parser.parse_args()

    cores = args.cores

    # get directories in data directory
    data_dir = args.dir

    # get all directories in data directory
    dirs = [x[0] for x in os.walk(data_dir)]
    
    if data_dir in dirs:
        dirs.remove(data_dir)

    # for each directory, get all fasta files
    for directory in dirs:
        # get all fasta files in directory
        files = os.listdir(directory)
        # take only fasta files
        files = [file for file in files if file.endswith(".template")]
        
        for file in files:
            id = file.split(".")[0]
            simulate_hifi_reads(directory, id, cores)


if __name__ == "__main__":
    sys.exit(main())