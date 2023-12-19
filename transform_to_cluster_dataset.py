import pandas as pd
import os
import sys
import argparse
import random
import editdistance
import csv


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--tsv', dest = 'tsv_file', required=True, type=int, help="")
    args = parser.parse_args()