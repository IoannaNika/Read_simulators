
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    args = parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())