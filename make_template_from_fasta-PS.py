

import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
import subprocess
def parse_R_output(output):
    # parse output from R script
    output = output.split("\n")
    start = output[1].split(" ")[1]
    end = output[4].split(" ")[1]
    return start, end

def get_amplicon_positions(Fprob, Rprob, seq_path):
    # call r script to get amplicon positions
    results = subprocess.run(["Rscript match_primers.R " + Fprob + " " + Rprob + " " + seq_path], shell=True, capture_output=True, text=True)
    results = parse_R_output(results.stdout)
    start = int(results[0]) + len(Fprob)
    end = int(results[1]) - len(Rprob)
    return start, end

def create_template(seq_path, id, template):
    template_df = pd.read_csv(template, sep='\t', header=None)
    template_df.columns = ["chr", "start", "end", "name_1", "score", "strand", "primer"]
    final_template = "" 

    # read template file two lines at a time (positive & negative strand info)
    for i in range(0, len(template_df), 2):
        cnt = 0
        # get positive strand info for reference genome
        pos_strand = template_df.iloc[i]
        # get negative strand info for reference genome
        neg_strand = template_df.iloc[i+1]

        seq_start = pos_strand["end"]
        seq_end = neg_strand["start"]

        Fprob = pos_strand["primer"].strip()
        # print(Fprob)
        Rprob = neg_strand["primer"].strip()
        # print(Rprob)
        
        try:
            # get amplicon positions
            start, end = get_amplicon_positions(Fprob, Rprob, seq_path)

        except:
            print("Something went wrong with the amplicon positions.")
            continue

        # parse the fasta file
        record = SeqIO.parse(seq_path, "fasta")
        seq = next(record).seq
        
        amplicon = seq[start:end]

        while cnt < 10:
            final_template += ">" + id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt) + "\n"
            final_template += str(amplicon) + "\n"
            cnt += 1

    return final_template

def main():
    parser = argparse.ArgumentParser(description="Creates amplicon fragments")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--bed_file', dest = 'bed_file', required=True, type=str, help="path to bed file with primer locations")
    args = parser.parse_args()

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
        files = [file for file in files if file.endswith(".fasta")]
        
        for file in files:
            file_path = directory + "/" + file
            gisaid_id = file.split(".fasta")[0]
            
            seq_path = directory + "/" + gisaid_id + ".fasta"

            template = create_template(seq_path, gisaid_id, args.bed_file)
            # write template to file
            template_file = directory + "/" + gisaid_id + ".template"
            # if template file already exists, remove it
            if os.path.exists(template_file):
                os.remove(template_file)
            with open(template_file, "w") as out_file:
                out_file.write(template)


    
if __name__ == "__main__":
    sys.exit(main())