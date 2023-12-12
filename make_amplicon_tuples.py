import pandas as pd
import os
import sys
import argparse
import random
import editdistance
import csv


def parse_maf_file(maf_file, read_id):
    maf_dict = {}
    with open(maf_file, "r") as f:
        lines = f.readlines()

    cnt = 0
    while cnt < len(lines):
        if lines[cnt].startswith("a"):
            next_line = lines[cnt+1].strip().split(" ")
            genomic_region = next_line[1].split(":")[1]
            next_line = lines[cnt+2].strip().split(" ")
            sim_read_id ="/".join(next_line[1].split("/")[:-1])
            if genomic_region not in maf_dict.keys():
                maf_dict[genomic_region] = set()
                maf_dict[genomic_region].add(sim_read_id)
            else:
                maf_dict[genomic_region].add(sim_read_id)

            cnt += 4
    return maf_dict

def parse_fastq(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l!= len(lines)):
        read_name = "/".join(lines[l].split(" ")[0].split("/")[-1])[1:]
        l+=1
        read = lines[l]
        l+=2
        # quality_score = lines[l]
        l+=1
        reads[read_name] = read

    return reads

def sample_negative_pair(read_anchor, maf_dict_2, reads2, gr_key):
    # get negative pair of reads from the same genomic region
    read_negative_sid = random.sample(maf_dict_2[gr_key],1)
    reads_negative = reads2[reads_negative]
    ed_neg =  editdistance.eval(read_anchor, reads_negative)
    return (read_negative_sid, reads_negative), ed_neg

def check_if_in_tsv(read_1_sid, read_2_sid, outdir):
    # load tsv file 
    tsv_file = os.path.join(outdir, "samples.tsv")

    out_df = pd.read_csv(tsv_file, sep="\t", header=0)

    if read_1_sid in out_df["read_1"].values and read_2_sid in out_df["read_2"].values:
        print("positive pair is already in tsv file")
        return True
    if read_1_sid in out_df["read_2"].values and read_2_sid in out_df["read_1"].values:
        print("positive pair is already in tsv file")
        return True
    
    return False

def sample_positive_pair(maf_dict_1, reads1, gr_key):
    # sample positive pair of reads from the same genomic region
    reads_positive = random.sample(maf_dict_1[gr_key], 2)
    read_anchor = reads1[reads_positive[0]]
    read_anchor_sid = reads_positive[0]
    read_positive = reads1[reads_positive[1]]
    read_pos_sid = reads_positive[1]
    ed_pos = editdistance.eval(read_anchor, read_positive)
    return (read_anchor_sid, read_anchor), (read_pos_sid, read_positive), ed_pos

def write_to_tsv(gisaid_id_1, gisaid_id_2, read_1_sid, read_2_sid, label, genomic_region, read_length, lineage, outdir):
    complete_read_id_1 = gisaid_id_1 + "_" + read_1_sid 
    complete_read_id_2 = gisaid_id_2 + "_" + read_2_sid
    # Header: read_1 read_2 label genomic_region read_length lineage
    tsv_file = os.path.join(outdir, "samples.tsv")
    # load tsv file
    writer = csv.writer(tsv_file, delimiter="\t")
    writer.writerow([complete_read_id_1, complete_read_id_2, label, genomic_region, read_length, lineage])
    writer.close()
    return 


def write_to_fasta(genome_id, read_sid, read, outdir):
    complete_read_id = genome_id + "_" + read_sid
    read_file = os.path.join(outdir, "reads", complete_read_id + ".fasta")
    with open(read_file, "w") as f:
        f.write(">" + complete_read_id + "\n")
        f.write(read + "\n")
    return

def sample(genomic_region, sample1, sample2, outdir):
    lineage1 = sample1[1].split("/")[-1]
    lineage2 = sample2[1].split("/")[-1]
    maf_file_1 = os.path.join(sample1[1], sample1[0] + ".maf")
    fastq_file_1 = os.path.join(sample1[1], sample1[0] + ".fastq")

    maf_file_2 = os.path.join(sample2[1], sample2[0] + ".maf")
    fastq_file_2 = os.path.join(sample2[1], sample2[0] + ".fastq")

    gr_key  = str(genomic_region[0]) + "_" + str(genomic_region[1])

    maf_dict_1 = parse_maf_file(maf_file_1, sample1[0]) # genomic region -> simulated read id 
    maf_dict_2 = parse_maf_file(maf_file_2, sample2[0])

    reads_1 = parse_fastq(fastq_file_1) # simulated read id -> read
    reads_2 = parse_fastq(fastq_file_2)

    # get positive pair
    (read_anchor_sid, read_anchor), (read_pos_sid, read_positive), ed = sample_positive_pair(maf_dict_1, reads_1, gr_key)
    
    # check if the positive pair is already in the tsv file
    if check_if_in_tsv(read_anchor_sid, read_pos_sid, outdir):
         # resampling positive pair
         (read_anchor_sid, read_anchor), (read_pos_sid, read_positive), ed_pos = sample_positive_pair(maf_dict_1, reads_1, gr_key)

    # get negative pair for the anchor read
    (read_negative_sid, read_negative), ed_neg = sample_negative_pair(read_anchor, maf_dict_2, reads_2, gr_key)

    # check if the negative pair is already in the tsv file
    if check_if_in_tsv(read_anchor_sid, read_negative_sid, outdir):
        # resampling negative pair
        (read_negative_sid, read_negative), ed_neg = sample_negative_pair(read_anchor, maf_dict_2, reads_2, gr_key)
    
    if ed_neg <= ed_pos:
        print("negative sample is too similar to anchor read, writing only positive pair to tsv file")
        write_to_tsv(sample1[0], sample1[0], read_anchor_sid, read_pos_sid, "positive", gr_key, len(read_anchor), len(read_positive), lineage1, lineage1, outdir)
        write_to_fasta(sample1[0], read_anchor_sid, read_anchor, outdir)
        return  1
    
    # write to tsv file both positive and negative pair
    write_to_tsv(sample1[0], sample2[0], read_anchor_sid, read_pos_sid, "positive", gr_key, len(read_anchor), len(read_positive), lineage1, lineage1, outdir)
    write_to_tsv(sample1[0], sample2[0], read_anchor_sid, read_negative_sid, "negative", gr_key, len(read_anchor), len(read_negative), lineage1, lineage2, outdir)
    
    # write to fasta file
    write_to_fasta(sample1[0], read_anchor_sid, read_anchor, outdir)
    write_to_fasta(sample1[0], read_pos_sid, read_positive, outdir)
    write_to_fasta(sample2[0], read_negative_sid, read_negative, outdir)
    return 2
    

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

    if not os.path.exists(tsv_file):
        with open(tsv_file, "w") as f:
            f.write("read_1\tread_2\tlabel\tgenomic_region\tread_length_r1\tread_length_r2\tlineage_r1\tlineage_r2\n")
    else:
        # if it exist remove and create a new one
        os.remove(tsv_file)
        with open(tsv_file, "w") as f:
            f.write("read_1\tread_2\tlabel\tgenomic_region\tread_length_r1\tread_length_r2\\tlineage_r1\tlineage_r2\n")

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

    count = 0
    while count < num_samples:
        # sample a random genomic region with equal probability
        genomic_region = genomic_regions[random.randint(0, len(genomic_regions) - 1)]
        
        # sample two info entries wihtout replacement
        info_sample = random.sample(info, 2)

        # sample the reads
        num_samples = sample(genomic_region, info_sample[0], info_sample[1], out_dir)

        count += num_samples
        print("number of samples: ", count)


if __name__ == "__main__":
    sys.exit(main())