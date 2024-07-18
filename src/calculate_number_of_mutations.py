import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO

def load_read(fasta_file): 

    record = SeqIO.read(fasta_file, "fasta")

    return record.seq

def write_to_temp_file(sequence1, sequence2):
    with open("temp.fasta", "w") as f:
        f.write(">sequence1\n")
        f.write(sequence1 + "\n")
        f.write(">sequence2\n")
        f.write(sequence2 + "\n")

def mafft_msa():
    os.system("mafft --auto --quiet --thread 4 temp.fasta > mafft_temp.fasta")
    return 

def clean_up():
    os.system("rm temp.fasta")
    os.system("rm mafft_temp.fasta")

def calc_n_mutations(read1, read2):

    write_to_temp_file(read1, read2)

    # align
    mafft_msa()

    # load aligned sequences expected two sequences
    aligned_sequences = SeqIO.parse("mafft_temp.fasta", "fasta")

    seq1 = next(aligned_sequences).seq
    seq2 = next(aligned_sequences).seq

    n_mutations = 0

    for i in range(len(seq1)):
        if (seq1[i] != seq2[i]) and (seq1[i] != "-") and (seq2[i] != "-"):
            n_mutations += 1
    
    clean_up()

    return n_mutations

    

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="data directory")
    args = parser.parse_args()

    samples = pd.read_csv(args.dir + "/samples.tsv", sep="\t")
    out_file = open(args.dir + "/samples_mutations.tsv", "w")
    out_file.write("read_1\tread_2\tlabel\tmutations\n")
    out_file.close()

    for index, row in samples.iterrows():
        read_1 = row["read_1"] 
        read_2 = row["read_2"]
        label = row["label"]
        

        read1 = load_read(args.dir + "/" + read_1 + ".fasta")
        read2 = load_read(args.dir + "/" + read_2 + ".fasta")

        mutation_n = calc_n_mutations(read1, read2)
        
        out_file = open(args.dir + "/samples_mutations.tsv", "a")
        out_file.write(read_1 + "\t" + read_2 + "\t" + label + "\t" + str(mutation_n) + "\n")
        out_file.close()


if __name__ == "__main__":
    sys.exit(main())