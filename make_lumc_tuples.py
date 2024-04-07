import argparse
import os
import pandas as pd
import editdistance
import numpy as np
from multiprocessing import Process

def rev_complement(read): 
    return  read[::-1].translate(str.maketrans("ACGT", "TGCA"))


def write_fasta_read(reads_path, tsv_file_prefix, read_row):
    if tsv_file_prefix != "":
        tsv_file_prefix = tsv_file_prefix + "_"

    path = reads_path + "/" + tsv_file_prefix + read_row["id"].replace("/", "_") + ".fasta"
    # if it exists, delete it
    if os.path.exists(path):
        os.remove(path)
    
    # create it again
    open(path, 'x').close()

    with open(path, "w") as f:
        f.write(">{}\n".format(tsv_file_prefix + read_row["id"]).replace("/", "_"))
        f.write(read_row["read"])
        f.close()

    return

def per_region_tuples_no_gt(tsv1, outfile, region, reads_path):
    start = region[0]
    end = region[1]

    reads_tsv1 = tsv1[(tsv1["start"] == start) & (tsv1["end"] == end)]
    print(reads_tsv1)

    # sample negative samples
    for i, read1 in reads_tsv1.iterrows():

        if read1["strand"] == "-":
                read1["read"] = rev_complement(read1["read"])
                
        # get all reads that span the region
        for j, read2 in reads_tsv1.iterrows():

            if j <= i: 
                continue

            if read2["strand"] == "-":
                read2["read"] = rev_complement(read2["read"])
                
            ed = editdistance.eval(read1["read"], read2["read"])

            with open(outfile, "a") as f:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("NA", read1["id"].replace("/", "_"), read2["id"].replace("/", "_"), start, end, ed))
                f.close()

            
            # write reads to the /reads folder
            write_fasta_read(reads_path, "", read2)
        write_fasta_read(reads_path, "", read1)

    return 
            
def per_region_tuples_negatives(tsv1, tsv2, outfile, region, reads_path, tsv1_file_prefix, tsv2_file_prefix):
    start = region[0]
    end = region[1]
    
    # get all reads that span the region
    reads_tsv1 = tsv1[(tsv1["start"] == start) & (tsv1["end"] == end)]
    reads_tsv2 = tsv2[(tsv2["start"] == start) & (tsv2["end"] == end)]
    
    # sample negative samples
    for i, read1 in reads_tsv1.iterrows():

        if read1["strand"] == "-":
            read1["read"] = rev_complement(read1["read"])

        # get all reads that span the region
        for j, read2 in reads_tsv2.iterrows():
            
            if read2["strand"] == "-":
                read2["read"] = rev_complement(read2["read"])
                
            ed = editdistance.eval(read1["read"], read2["read"])

            with open(outfile, "a") as f:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("negative", tsv1_file_prefix + "_" + read1["id"].replace("/", "_"), tsv2_file_prefix + "_" + read2["id"].replace("/", "_"), start, end, ed))
                f.close()

            # write reads to the /reads folder
            write_fasta_read(reads_path, tsv2_file_prefix, read2)
        write_fasta_read(reads_path, tsv1_file_prefix, read1)


def make_region_tuples_positives(tsv, outfile, region, reads_path, tsv_file_prefix): 
    start = region[0]
    end = region[1]
    # get all reads that span the region
    reads_tsv = tsv[(tsv["start"] == start) & (tsv["end"] == end)]
    print("make_region_tuples_positives")
    # sample positive samples
    for index, read1 in reads_tsv.iterrows():
        
        if read1["strand"] == "-":
                read1["read"] = read1["read"][::-1].translate(str.maketrans("ACGT", "TGCA"))
            
        #  iterate over the succesive rows of the same dataframe to make all possible pairs
        for jindex, read1_next in reads_tsv.iterrows():
            if jindex <= index:
                continue
            
            if read1_next["strand"] == "-": 
                read1_next["read"] = read1_next["read"][::-1].translate(str.maketrans("ACGT", "TGCA"))
          
            ed = editdistance.eval(read1["read"], read1_next["read"])

            with open(outfile, "a") as f:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("positive", tsv_file_prefix + "_" + read1["id"].replace("/", "_"), tsv_file_prefix + "_" + read1_next["id"].replace("/", "_"), start, end, ed))
                f.close()
            
            # write reads to the /reads folder
            write_fasta_read(reads_path, tsv_file_prefix, read1_next)
        
        write_fasta_read(reads_path, tsv_file_prefix, read1)

    return


def make_all_tuples_no_gt(tsv1, outfile, genomic_regions): 
    # if outfile exists, delete it
    if os.path.exists(outfile):
        os.remove(outfile)
    # create it again
    open(outfile, 'x').close()

    # write header for the tsv file
    with open(outfile, "a") as f:
        f.write("label\tid1\tid2\tstart\tend\tedit_distance\n")
        f.close()
    
    # make a /reads folder 
    reads_path = outfile.split("/")[:-1]
    reads_path = "/".join(reads_path) + "/reads"
    # if it exists delete it and its contents
    if os.path.exists(reads_path):
        os.system("rm -r " + reads_path)

    if not os.path.exists(reads_path):
        os.makedirs(reads_path)

    processes = [Process(target=per_region_tuples_no_gt, args=(tsv1, outfile, region, reads_path)) for region in genomic_regions]

    
    for process in processes:
        process.start()
    for process in processes:
        process.join()

    return

def make_all_tuples(tsv1, tsv2, outfile, genomic_regions, tsv1_file_prefix, tsv2_file_prefix):
    # if outfile exists, delete it
    if os.path.exists(outfile):
        os.remove(outfile)
    # create it again
    open(outfile, 'x').close()

    reads_path = outfile.split("/")[:-1]
    reads_path = "/".join(reads_path) + "/reads"
    # if it exists delete it and its contents
    if os.path.exists(reads_path):
        os.system("rm -r " + reads_path)
        
    if not os.path.exists(reads_path):
        os.makedirs(reads_path)

    # write header for the tsv file
    with open(outfile, "a") as f:
        f.write("label\tid1\tid2\tstart\tend\tedit_distance\n")
        f.close()

    processes = [Process(target=per_region_tuples_negatives, args=(tsv1, tsv2, outfile, region, reads_path, tsv1_file_prefix, tsv2_file_prefix)) for region in genomic_regions]

    for process in processes:
        process.start()
    for process in processes:
        process.join()

    processes = [Process(target=make_region_tuples_positives, args=(tsv1, outfile, region, reads_path, tsv1_file_prefix)) for region in genomic_regions]

    for process in processes:
        process.start()
    for process in processes:
        process.join()
    
    processes = [Process(target=make_region_tuples_positives, args=(tsv2, outfile, region, reads_path, tsv2_file_prefix)) for region in genomic_regions]

    for process in processes:
        process.start()
    for process in processes:
        process.join()

    return


def make_tuples(tsv1, tsv2, n, outfile, singles, genomic_regions):

    # if outfile exists, delete it
    if os.path.exists(outfile):
        os.remove(outfile)
    # create it again
    open(outfile, 'x').close()
    # write header for the tsv file
    with open(outfile, "a") as f:
        if singles == False: 
            f.write("label\tid1\tread1\tid2\tread2\tstart\tend\tedit_distance\n")
        else: 
            f.write("id\tread\tgenomic_regions\tfile\n")

    count = 0
    # strand = ["+", "-"]
    tsvs = [tsv1, tsv2]

    while count < n:        
        # sample uniformly from the genomic regions
        region = genomic_regions[np.random.randint(len(genomic_regions))]
        start = region[0]
        end = region[1]

        # sample uniformly a tsv file
        indx = np.random.randint(len(tsvs))
        tsv = tsvs[indx]
        # sample uniformly a strand
        # s = strand[np.random.randint(len(strand))]

        # find reads that span the region
        reads = tsv[(tsv["start"] == start) & (tsv["end"] == end)]  #& (tsv["strand"] == s))]

        # if there are no reads that span the region, continue
        if reads.shape[0] < 2:
            print("No reads that span the region")
            continue
    
        # sample uniformly two reads
        reads = reads.sample(2)
        read_anchor = reads.iloc[0]
        read_anchor_strand = read_anchor["strand"]
        read_pos = reads.iloc[1]
        read_pos_strand = read_pos["strand"]

        if read_anchor_strand == "-":
            read_anchor["read"] = rev_complement(read_anchor["read"])
        if read_pos_strand == "-":
            read_pos["read"] = rev_complement(read_pos["read"])
            
        ed_pos = editdistance.eval(read_anchor["read"], read_pos["read"])

        # write to tsv file
        with open(outfile, "a") as f:
            if singles == False:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("positive", read_anchor["id"], read_anchor["read"], read_pos["id"], read_pos["read"], read_anchor["start"], read_anchor["end"], ed_pos))
            else: 
                f.write("{}\t{}\t{}\t{}\n".format(read_anchor["id"], read_anchor["read"], str(read_anchor["start"]) + "_" + str(read_anchor["end"]), indx))
        count += 1

        # sample a negative sample for the anchor
        if indx == 0:
            tsv_neg = tsv2
        else:
            tsv_neg = tsv1
        
        # sample uniformly a read from the negative tsv file
        tsv_neg = tsv_neg[(tsv_neg["start"] == start) & (tsv_neg["end"] == end)] # tsv_neg["strand"] == s) &
        
        if tsv_neg.shape[0] == 0:
            print("No reads that span the region")
            continue

        read_neg = tsv_neg.sample(1)
        read_neg = read_neg.iloc[0]
        read_neg_strand = read_neg["strand"]
        if read_neg_strand == "-": 
            read_neg["read"] = rev_complement(read_neg["read"])
        ed_neg = editdistance.eval(read_anchor["read"], read_neg["read"])

        if ed_neg < 1:
            continue
        # write to tsv file
        if ed_neg <= ed_pos:
            with open(outfile, "a") as f:
                if singles == False:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("hard negative", read_anchor["id"], read_anchor["read"],read_neg["id"], read_neg["read"], read_anchor["start"], read_anchor["end"], ed_neg))
                else: 
                    f.write("{}\t{}\t{}\t{}\n".format(read_neg["id"], read_neg["read"], str(read_anchor["start"]) + "_" +  str(read_anchor["end"]), 1-indx))

            count += 1
        else:
            with open(outfile, "a") as f:
                if singles == False:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("negative", read_anchor["id"], read_anchor["read"],read_neg["id"], read_neg["read"], read_anchor["start"], read_anchor["end"], ed_neg))
                else:
                    f.write("{}\t{}\t{}\t{}\n".format(read_neg["id"], read_neg["read"], str(read_anchor["start"]) + "_" + str(read_anchor["end"]), 1-indx))

            count += 1
    return


def main():
    parser = argparse.ArgumentParser(description="Make tuples for LUMC dataset")
    parser.add_argument('--out', dest = 'out', required=True, type=str, help="output tsv file")
    parser.add_argument('--tsv1', dest = 'tsv1', required=True, type=str, help="tsv file with reads from the first dataset")
    parser.add_argument('--tsv2', dest = 'tsv2', required=False, default=None, type=str, help="tsv file with reads from the second dataset")
    parser.add_argument("--n", dest="n", required=False, default=1000, type=int, help="number of tuples to create")
    parser.add_argument("--all", dest="all", required=False, default=False, type=bool, help="True, if all tuples should be created, False, is only the specified number of tuples should be created, when one tsv file is provided, pairs within the same tsv file will be created")
    parser.add_argument('--singles', dest = 'singles', required=False, default=False, type=bool, help="True, if only singles should be created, False, if tuples should be created")
    parser.add_argument('--with_prefix', dest ='with_prefix', required=False, default=False, type=bool, help="include file prefix in read ids e.g. Wuhan/tsv_file_1.tsv, Wuhan will be used as prefix for the read ids to differentiate their origin, implemented only for the --all option")
    args = parser.parse_args()
    
    tsv1_file = args.tsv1
    tsv2_file = args.tsv2
    tsv1_file_prefix = ""
    tsv2_file_prefix = ""

    if args.with_prefix and tsv2_file != None: 
        tsv1_file_prefix = tsv1_file.split("/")[-2]
        tsv2_file_prefix = tsv2_file.split("/")[-2]

    n = args.n
    outfile = args.out
    singles = args.singles

    # Read in the tsv files
    tsv1 = pd.read_csv(tsv1_file, sep="\t")
    tsv1.columns = ["id", "read", "start", "end", "strand"]
    if tsv2_file != None:
        tsv2 = pd.read_csv(tsv2_file, sep="\t")
        tsv2.columns = ["id", "read", "start", "end", "strand"]

    # lumc 
    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    # cerbaresearch nsp12
    # genomic_regions = [(13055,13960), (13803, 14718), (14600, 15599), (15371, 16311)]

    
    if tsv2_file == None and args.all == True: 
        make_all_tuples_no_gt(tsv1, outfile, genomic_regions)

    if args.all == False and tsv2_file != None:
        make_tuples(tsv1, tsv2, n, outfile, singles, genomic_regions)

    if tsv2_file != None and args.all == True:
        make_all_tuples(tsv1, tsv2, outfile, genomic_regions, tsv1_file_prefix, tsv2_file_prefix)



if __name__ == "__main__":
    main()

