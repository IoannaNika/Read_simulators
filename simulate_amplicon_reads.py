
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO


def simulate_ONT_HQ_reads(directory, identifier, cores): 
    os.system("pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-ONT-HQ.model --template {}/{}.template --prefix {}/{}".format(directory, identifier, directory, identifier))
    return

def simulate_hifi_reads(directory, identifier, cores):
    os.system("pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-SEQUEL.model  --template {}/{}.template --pass-num 10 --seed 20 --prefix {}/{}".format(directory, identifier, directory, identifier))
    os.system("samtools view -bS {}/{}.sam > {}/{}.bam".format(directory, identifier, directory, identifier))
    os.system("ccs {}/{}.bam -j {} {}/{}.fastq".format(directory, identifier, cores, directory, identifier))
    return

def main():
    parser = argparse.ArgumentParser(description="Simulates pacbio hifi amplicon reads")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--cores', dest = 'cores', default=8, required=False, type=int, help="cpu cores for the simulation")
    parser.add_argument('--strategy', dest = 'strategy', default="pacbio-hifi", required=True, type=str, help="ONT or pacbio-hifi")    
    args = parser.parse_args()

    cores = args.cores
    strategy = args.strategy

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
            identifier = file.split(".template")[0]
            if strategy == "pacbio-hifi":
                simulate_hifi_reads(directory, identifier, cores)
            if strategy == "ONT": 
                simulate_ONT_HQ_reads(directory, identifier, cores)


    print("Done")



if __name__ == "__main__":
    sys.exit(main())