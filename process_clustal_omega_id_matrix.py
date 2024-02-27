

import pandas as pd
import os
import sys
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="process percent identity matrix from clustal omega")
    parser.add_argument('--pim', dest = 'pim', required=True, type=str, help="percentage identity matrix")
    args = parser.parse_args()

    pim = args.pim

    # read in the percent identity matrix with format: 
    # 1: DQ480524.1  100.00   93.29   93.27   93.15   93.70   93.59   93.52   93.43   93.37   93.27   93.47   93.67   93.65   93.60   69.32   69.35   69.01   69.93   69.83   69.98   69.16   70.08   69.81   69.71   70.00   69.95   69.57   69.79   69.79   69.78   69.72   69.87   69.95   69.58   69.45   69.95   69.76   69.36   69.62   69.75   69.99   69.71   69.89   70.03   69.11   69.24   69.56   69.67   69.30   69.91   69.64   69.97   69.79   69.44   69.71   69.68   69.65   69.54   69.71   69.34   69.51   69.34   69.62   69.55   69.67   69.62   69.60   69.72   69.46   69.49   69.79   69.99   70.12   69.90   69.91   69.41   69.91   69.93   70.07   69.78   69.85   69.93   70.06   69.74   69.79   70.07   70.15   70.01   69.95   69.71   69.77   69.89   69.85   69.90   69.58   69.47   69.91   69.73   69.73   69.79   68.77   69.59   69.71   69.87   69.77   69.78   69.85   69.92   69.81   69.79   69.69   69.56   70.02   69.32   69.76   69.54   69.64   69.78   69.59   69.76   69.72   69.85   69.95   69.88   69.85   69.98   70.39   69.87   69.62   69.69   69.83   69.77   69.87   69.86   69.85   70.04   69.67   69.73   69.79   69.91   69.71   69.68   69.94   69.83   69.81   69.61   69.35   69.71   68.90   69.69   69.79   69.58   69.78   69.89   69.73   69.37   69.81   70.08   69.75   70.16   70.20   70.17   69.60   69.79   70.04   69.90   69.83   69.76
    # 2: DQ480523.1   93.29  100.00   95.04   94.51   95.26   95.14   95.09   94.91   94.87   94.91   94.97   95.18   95.14   95.35   69.69   69.73   69.03   69.96   70.16   70.05   69.33   70.25   70.10   69.93   70.43   70.19   70.00   70.21   70.09   70.09   69.98   70.22   69.91   69.78   69.67   70.21   70.12   69.75   69.99   70.15   70.15   70.13   70.11   70.22   69.40   69.51   69.98   70.26   69.61   70.22   70.15   70.01   70.03   69.68   70.10   69.94   69.91   69.92   70.05   69.66   70.04   69.61   70.14   69.95   70.09   69.96   69.78   69.90   69.92   69.95   70.07   70.28   70.26   70.17   70.25   69.85   70.19   70.18   70.40   69.91   70.13   70.30   70.55   70.21   70.17   70.41   70.37   70.32   70.18   70.09   70.38   70.18   70.18   70.39   70.03   69.99   70.27   70.20   70.19   70.19   69.22   70.14   70.14   70.27   69.98   70.19   70.20   70.40   70.05   70.10   69.98   69.92   70.11   69.66   70.14   69.94   70.08   70.22   70.02   70.01   70.10   70.14   70.06   69.91   70.15   70.10   70.56   69.99   70.03   69.85   70.20   70.13   70.10   70.25   70.29   70.41   70.14   70.11   70.20   70.13   70.10   69.92   70.12   70.10   69.97   69.96   69.67   70.03   69.26   70.07   70.26   69.96   70.14   70.10   69.99   69.65   70.04   70.42   70.06   70.23   70.40   70.38   69.79   70.00   70.40   70.22   70.17   70.12

    # read lines from the file
    with open(pim, "r") as f:
        pim_lines = f.readlines()

    seq_dict = {}

    for i in range(len(pim_lines)):
        row = pim_lines[i]
        row = row.split(":")[1].strip()
        seq_name = row.split("  ")[0]
        seq_dict[seq_name] = dict()
        other_seq = row.split("  ")[1:]
        # remove "" from list
        other_seq = [x for x in other_seq if x != ""]
        # to integer 
        other_seq = [float(x.strip()) for x in other_seq]

        for j in range(len(other_seq)):
            # get the sequence name from the row with the same index as the current  column index
            other_seq_name = pim_lines[j].split(":")[1].strip()
            other_seq_name = other_seq_name.split("  ")[0]
            seq_dict[seq_name][other_seq_name] = other_seq[j]

    
    # Get the entry with the highest number of entries with more than 95% identity
    max_count = 0
    max_seq = ""
    for key in seq_dict:
        count = 0
        for key2 in seq_dict[key]:
            if seq_dict[key][key2] > 95:
                count += 1
        if count > max_count:
            max_count = count
            max_seq = key
    
    print(max_seq, max_count)

    # print the entries with which the max_seq has more than 95% identity
    for key in seq_dict[max_seq]:
        if seq_dict[max_seq][key] > 95:
            print(key, seq_dict[max_seq][key])

    



if __name__ == "__main__":
    sys.exit(main())
