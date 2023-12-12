import pandas as pd
import os
import sys
import argparse
import random

def sample_overlapping_positive(genomic_region, sample1):
    
    return

def sample_overlapping_negative(genomic_region, sample1, sample2):
    return


def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--n', dest = 'num', required=True, type=str, help="number of samples to produce")
    parser.add_argument('--dir', dest = 'data_dir', required=True, type=str, help="data directory to sample from")
    parser.add_argument('--out', dest = 'out_dir', required=True, type=str, help="output directory to write samples to")
    args = parser.parse_args()

    data_directory = args.data_dir
    num_samples = args.num
    out_dir = args.out_dir

    # check if output directory exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
       if len(os.listdir(out_dir)) > 0:
           for file in os.listdir(out_dir):
               # if its a directory, skip
                if os.path.isdir(os.path.join(out_dir, file)):
                    continue
                os.remove(os.path.join(out_dir, file))

    # make tsv file to write the samples to
    tsv_file = os.path.join(out_dir, "samples.tsv")

    # make a read directory if it doesn't exist
    read_dir = os.path.join(out_dir, "reads")
    if not os.path.exists(read_dir):
        os.makedirs(read_dir)
    else:
        if len(os.listdir(read_dir)) > 0:
            # make sure it's empty
            for file in os.listdir(read_dir):
                os.remove(os.path.join(read_dir, file))


    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768)]
    

    # get all fasta files in the data directory and subdirectories
    fasta_files = []
    for root, dirs, files in os.walk(data_directory):
        for file in files:
            if file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))

    # remove the sequences.fasta file
    if "sequences.fasta" in fasta_files:
        fasta_files.remove("sequences.fasta")
    info = [(f.split("/")[-1][:-6], "/".join(f.split("/")[:-1])) for f in fasta_files]

    for i in range(int(num_samples/2)):
        # sample a random genomic region with equal probability
        genomic_region = genomic_regions[random.randint(0, len(genomic_regions) - 1)]
        
        # sample two info entries wihtout replacement
        info_sample = random.sample(info, 2)

        sample_overlapping_positive(genomic_region, info_sample[0])
        sample_overlapping_negative(genomic_region, info_sample[0], info_sample[1])

    


    


if __name__ == "__main__":
    sys.exit(main())