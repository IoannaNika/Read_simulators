import argparse
import sys
import os
from Bio import SeqIO
import pandas as pd


def mafft_msa(dir): 
    if os.path.exists("{}/msa_alignment.fasta".format(dir)):
        return
    os.system("mafft --auto --quiet --thread 4 {}/sequences.fasta > {}/msa_alignment.fasta".format(dir, dir))
    return 


def find_amplicons(dir ,primers, start, end, record_id):

    msa_alignment = "{}/msa_alignment.fasta".format(dir)
    fasta_sequences = SeqIO.parse(open(msa_alignment),'fasta')

    # primer_locations = parse_primer_file(primers)
    
    for record in fasta_sequences:
        if record.id == record_id:
            if start > len(record.seq) or end > len(record.seq):
                raise ValueError("Start or end is out of range, skipping")
            amplicon = record.seq[start:end]
            return amplicon
    
    raise ValueError("Record not found")
      

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="directory with sequences.fasta which contains all sequences to be used for amplicon detection")
    parser.add_argument('--start', dest = 'start', required=True, type=int, help="start of amplicon")
    parser.add_argument('--end', dest = 'end', required=True, type=int, help="end of amplicon")
    parser.add_argument('--record_id', dest = 'record_id', required=True, type=str, help="record id of sequence to be used for amplicon detection")
    args = parser.parse_args()

    dir = args.dir
    start = args.start
    end = args.end
    record_id = args.record_id


    mafft_msa(dir)

    amplicon = find_amplicons(dir, start, end, record_id)
    # remove any gaps
    amplicon = amplicon.replace("-", "").strip()

    print(amplicon)


if __name__ == "__main__":
    sys.exit(main())