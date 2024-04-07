
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Filter NCBI sequences to only include those with complete genomes and write to new directory.")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory, containing metadata.txt and sequences.fasta files")
    args = parser.parse_args()

    # # metadata
    # metadata = args.dir + "/metadata.txt"
    # sequence path
    seq_path = args.dir + "/sequences.fasta"

    # read txt file of the format: 
    # 1. Zika virus isolate Homo sapiens/Haiti/0728/2014, complete genome
    # 10,807 bp linear RNA 
    # OK571913.1 GI:2128382615

    # list_of_identifiers = []

    # with open(metadata) as f:
    #     lines = f.readlines()
    
    # counter = 0

    # while counter < len(lines): 
    #     if  lines[counter].split(", ")[-1].strip() == "complete genome":
    #         counter+=2
    #         identifier = lines[counter].split(" ")[0].strip()
    #         list_of_identifiers.append(identifier)
    #         counter +=2
    #     else: 
    #         counter +=4

    # # print list of identifiers in a output.txt file where each identifier is written below each other without comma
    # out_file = "output.txt"
    # # make file 
    # with open(out_file, "w") as out_file:
    #     for id in list_of_identifiers:
    #         out_file.write(id + "\n")


    # read sequences.fasta file
    fasta_sequences = SeqIO.parse(open(seq_path),'fasta')
    # create output directory in data directory
    parent_out_dir = args.dir + "/NCBI_processed"
    if not os.path.exists(parent_out_dir):
        os.mkdir(parent_out_dir)

    # total number of sequences considered
    total_seqs = 0

    for fasta in fasta_sequences:
        strain, sequence = fasta.id, str(fasta.seq)
        identifier = strain.split(" ")[0].strip()

        
        # if identifier not in list_of_identifiers:
        #     continue
        
        total_seqs += 1

        # create output file name
        output_file = parent_out_dir + "/" + identifier + ".fasta"
        # if file exists, delete and rewrite
        if os.path.exists(output_file):
            os.remove(output_file)
        # write to file
        with open(output_file, "w") as out_file:
            out_file.write(">" + strain + "\n" + sequence + "\n")

    print("Total number of sequences considered: " + str(total_seqs))

if __name__ == "__main__":
    sys.exit(main())





