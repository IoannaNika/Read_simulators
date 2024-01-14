import pandas as pd
import os
import sys
import argparse
import random
import editdistance
import csv


def main():
    parser = argparse.ArgumentParser(description="Transform triplets to tuples")
    parser.add_argument('--tsv', dest = 'tsv_file', required=True, type=str, help="")
    parser.add_argument('--out', dest = 'out', required=True, type=str, help="")
    args = parser.parse_args()

    tsv_file = args.tsv_file
    outfile = args.out

    # make a new tsv file with the following columns:
    # read_id, genomic region, label
    # Read in the tsv file
    tsv = pd.read_csv(tsv_file, sep="\t")

    # check if the outfile exists
    if os.path.exists(outfile):
        os.remove(outfile)
    # create it again
    open(outfile, 'x').close()
    # write header to the outfile
    with open(outfile, "a") as f:
        f.write("read_1\tread_2\tgenomic_region\tlabel\tlineage_r1\tlineage_r2\n")

    # go through the tsv file to map it to the new tsv file
    for indx, row in tsv.iterrows():
        # read 1
        read_anch = row["read_anch"]
        read_pos = row["read_pos"]
        read_neg = row["read_neg"]
        label_neg = row["neg_label"]
        label_pos = "positive"
        genomic_region = row["genomic_region"]
        lineage_pos = row["lineage_pos"]
        lineage_neg = row["lineage_neg"]


        # write to the new tsv file
        with open(outfile, "a") as f:
            f.write(f"{read_anch}\t{read_pos}\t{genomic_region}\t{label_pos}\t{lineage_pos}\t{lineage_pos}\n")
            f.write(f"{read_anch}\t{read_neg}\t{genomic_region}\t{label_neg}\t{lineage_neg}\t{lineage_neg}\n")

if __name__ == "__main__":
    sys.exit(main())