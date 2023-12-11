

import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO

def create_template(seq, id):
    template = "data/nCoV-2019.bed"
    template_df = pd.read_csv(template, sep='\t', header=None)
    template_df.columns = ["chr", "start", "end", "name_1", "name_2", "strand"]

    final_template = ""
    # read template file two lines at a time (positive & negative strand info)
    for i in range(0, len(template_df), 2):
        # get positive strand info
        pos_strand = template_df.iloc[i]
        # get negative strand info
        neg_strand = template_df.iloc[i+1]

        seq_start = pos_strand["end"]
        seq_end = neg_strand["start"]

        if seq_start > len(seq):
            break

        if seq_end > len(seq):
            seq_end = len(seq)
        
        amplicon = seq[seq_start:seq_end]

        final_template += ">" + id + ":" + str(seq_start) + "_" + str(seq_end) + "\n"
        final_template += str(amplicon) + "\n"
    
    return final_template


def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    args = parser.parse_args()

    # get directories in data directory
    data_dir = args.dir

    # get all directories in data directory
    dirs = [[x[0] for x in os.walk(data_dir)][1]]
    print(dirs)
    
    # for each directory, get all fasta files
    for directory in dirs:
        # get all fasta files in directory
        files = os.listdir(directory)
        # take only fasta files
        files = [file for file in files if file.endswith(".fasta")]
        
        for file in files:
            file_path = directory + "/" + file
            gisaid_id = file.split(".")[0]
            # parse the fasta file
            record = SeqIO.parse(file_path, "fasta")

            seq = next(record).seq

            template = create_template(seq, gisaid_id)
            # write template to file
            template_file = directory + "/" + gisaid_id + ".template"

            with open(template_file, "w") as out_file:
                out_file.write(template)

    
if __name__ == "__main__":
    sys.exit(main())