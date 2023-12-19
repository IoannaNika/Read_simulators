import pandas as pd
import os
import sys
import argparse
import random
import editdistance
import csv


def main():
    parser = argparse.ArgumentParser(description="")
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
        f.write("read_id\tgenomic_region\tlabel\n")

    # go through the tsv file to map it to the new tsv file
    for row in tsv.iterrows():
        # read 1
        read_id_1 = row["read_1"]
        genomic_region_1 = row["genomic_region"]
        label_1 = "_".join(read_id_1.split("_")[:-1])

        # read 2
        read_id_2 = row["read_2"]
        genomic_region_2 = row["genomic_region"]
        label_2 = "_".join(read_id_2.split("_")[:-1])

        # write to the new tsv file
        with open(outfile, "a") as f:
            f.write(f"{read_id_1}\t{genomic_region_1}\t{label_1}\n")
            f.write(f"{read_id_2}\t{genomic_region_2}\t{label_2}\n")

    


if __name__ == "__main__":
    sys.exit(main())