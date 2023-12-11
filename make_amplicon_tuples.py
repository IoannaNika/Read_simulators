import pandas as pd
import os
import sys
import argparse

def sample_overlapping_positive():
    return

def sample_overlapping_negative():
    return


def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--n', dest = 'num', required=True, type=str, help="number of samples to produce")
    parser.add_argument('--dir', dest = 'data_dir', required=True, type=str, help="data directory to sample from")
    args = parser.parse_args()

    data_directory = args.data_dir
    num_samples = args.num



if __name__ == "__main__":
    sys.exit(main())