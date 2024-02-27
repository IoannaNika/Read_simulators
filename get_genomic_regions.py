import argparse
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Get genomic regions")
    parser.add_argument('--bed_file', dest = 'bed_file', required=True, type=str, help="path to bed file with primer locations")
    args = parser.parse_args()

    # get amplicon positions from the bed file
    bed_file = args.bed_file
    # read as tsv file
    # OQ551959.1	120	140	scheme_1_LEFT	1	+	AGGACTCGGCTTGCTGAAGT
    # OQ551959.1	1081	1107	scheme_1_RIGHT	1	-	ACCGATCTACATAGTCTCTAAAGGGT

    # read bed file
    bed = pd.read_csv(bed_file, sep="\t", header=None)
    bed.columns = ["chr", "start", "end", "name_1", "score", "strand", "primer"]
    result = []
    # read in pairs of 2 rows
    for i in range(0, len(bed), 2):
        # get positive strand info for reference genome
        pos_strand = bed.iloc[i]
        # get negative strand info for reference genome
        neg_strand = bed.iloc[i+1]

        seq_start = pos_strand["end"]
        seq_end = neg_strand["start"]

        result.append((seq_start, seq_end))
    return result

if __name__ == "__main__":
    sys.exit(main())