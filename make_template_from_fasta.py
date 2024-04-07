

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

def get_primer_sequences(primer_file, index):

    # read primer file
    primer_df = pd.read_csv(primer_file, sep='\t', header=0)
    primer_df.columns = ["name", "pool", "seq", "length", "%gc", "tm (use 65)"]

    # get primer sequences
    Fprob = primer_df.iloc[index]["seq"]
    Rprob = primer_df.iloc[index+1]["seq"]
    
    return Fprob, Rprob


def create_template(seq_path, s_id, template, primer_file, n_templates):
    template_df = pd.read_csv(template, sep='\t', header=None)
    template_df.columns = ["chr", "start", "end", "name_1", "name_2", "strand"]
    final_template = "" 
    # read template file two lines at a time (positive & negative strand info)
    for i in range(0, len(template_df), 2):
        # get positive strand info for reference genome
        pos_strand = template_df.iloc[i]
        # get negative strand info for reference genome
        neg_strand = template_df.iloc[i+1]

        seq_start = pos_strand["end"]
        seq_end = neg_strand["start"]

        try:
            # get primer sequences
            Fprob, Rprob = get_primer_sequences(primer_file, i)

            # get amplicon positions
            start, end = get_amplicon_positions(Fprob, Rprob, seq_path)

        except:
            print("Something went wrong with the amplicon positions.")
            print("Will get amplicon with MSAs instead")

            parent_dir = "/".join(seq_path.split("/")[: -2])
            amplicon = subprocess.run(['python find_amplicons_with_mfft.py --dir {} --start {} --end {} --record_id {}'.format(parent_dir, seq_start, seq_end, s_id)], shell=True, capture_output=True, text=True)
            print("Amplicon: ", amplicon)
            amplicon = str(amplicon.stdout).strip()
            cnt = 0
            while cnt < int(n_templates/2):
                # write hald with + strand and half with - strand
                final_template += ">" + "+_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt) + "\n"
                final_template += amplicon + "\n"
                cnt += 1
            leftover_templates = n_templates - int(n_templates/2)
            cnt2 = 0
            while cnt2 < leftover_templates:
                final_template += ">" + "-_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt+cnt2) + "\n"
                final_template += str(amplicon[::-1].translate(str.maketrans("ATGC", "TACG"))) + "\n"
                cnt2 += 1
            continue

        # parse the fasta file
        record = SeqIO.parse(seq_path, "fasta")
        seq = next(record).seq
        
        amplicon = seq[start:end]
        cnt = 0
        while cnt < int(n_templates/2):
            # write hald with + strand and half with - strand
            final_template += ">" + "+_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt) + "\n"
            final_template += str(amplicon) + "\n"
            cnt += 1
        leftover_templates = n_templates - int(n_templates/2)
        cnt2 = 0
        while cnt2 < leftover_templates:
            final_template += ">" + "-_" + s_id + ":" + str(seq_start) + "_" + str(seq_end) + ":" + str(cnt+cnt2) + "\n"
            final_template += str(str(amplicon)[::-1].translate(str.maketrans("ATGC", "TACG"))) + "\n"
            cnt2 += 1

    return final_template



def main():
    parser = argparse.ArgumentParser(description="Creates amplicon fragments")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--bed_file', dest = 'bed_file', required=False, default="data/nCoV-2019.bed", type=str, help="path to bed file with primer locations")
    parser.add_argument('--primer_file', dest = 'primer_file', required=False, default="data/nCoV-2019.tsv",type=str, help="path to primer file")
    parser.add_argument('--seq_id', dest = 'seq_id', required=False, default=None, type=str, help="create template for specific sequence")
    parser.add_argument('--cnt_templates', dest = 'cnt_templates',  required=False, default=10, type=int, help="number of templates to create per sequence")
    args = parser.parse_args()

    # get directories in data directory
    data_dir = args.dir

    if args.seq_id == None:
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
                gisaid_id = file.split(".")[0]
                
                seq_path = directory + "/" + gisaid_id + ".fasta"
            

                template = create_template(seq_path, gisaid_id, args.bed_file, args.primer_file)
                # write template to file
                template_file = directory + "/" + gisaid_id + ".template"

                with open(template_file, "w") as out_file:
                    out_file.write(template)
    else: 
            
        seq_path = args.dir + "/" + args.seq_id + ".fasta"
        template = create_template(seq_path, args.seq_id, args.bed_file, args.primer_file, args.cnt_templates)
        # write template to file
        template_file = args.dir + "/" + args.seq_id + ".template"
        # if template file already exists, remove it
        if os.path.exists(template_file):
            os.remove(template_file)

        open(template_file, 'x').close()

        with open(template_file, "w") as out_file:
            out_file.write(template)    



    
if __name__ == "__main__":
    sys.exit(main())