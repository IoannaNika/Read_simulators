import argparse
import os
import pandas as pd
import editdistance
import numpy as np



def main():
    parser = argparse.ArgumentParser(description="Map LUMC regions")
    parser.add_argument('--tsv1', dest = 'tsv1', required=True, type=str, help="")
    parser.add_argument('--tsv2', dest = 'tsv2', required=True, type=str, help="")
    parser.add_argument("--n", dest="n", required=True, type=int, help="")
    parser.add_argument('--out', dest = 'out', required=True, type=str, help="")
    args = parser.parse_args()

    tsv1_file = args.tsv1
    tsv2_file = args.tsv2
    n = args.n
    outfile = args.out

    # if outfile exists, delete it
    if os.path.exists(outfile):
        os.remove(outfile)
    # create it again
    open(outfile, 'x').close()
    # write header for the tsv file
    with open(outfile, "a") as f:
        f.write("label\tid1\tread1\tid2\tread2\tstart\tend\tedit_distance\n")

    # Read in the tsv files
    tsv1 = pd.read_csv(tsv1_file, sep="\t")
    tsv1.columns = ["id", "read", "start", "end", "strand"]
    tsv2 = pd.read_csv(tsv2_file, sep="\t")
    tsv2.columns = ["id", "read", "start", "end", "strand"]

    # sample positive and negative tuples from the same regions

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768)]

    count = 0
    strand = ["+", "-"]
    tsvs = [tsv1, tsv2]
    while count < n:


        # from each tsv file find reads that span the same region and are on the same strand
        
        # sample uniformly from the genomic regions
        region = genomic_regions[np.random.randint(len(genomic_regions))]
        start = region[0]
        end = region[1]

        # sample uniformly a tsv file
        indx = np.random.randint(len(tsvs))
        tsv = tsvs[indx]
        # sample uniformly a strand
        s = strand[np.random.randint(len(strand))]

        # find reads that span the region
        reads = tsv[((tsv["start"] == start) & (tsv["end"] == end) & (tsv["strand"] == s))]

        # if there are no reads that span the region, continue
        if reads.shape[0] < 2:
            print("No reads that span the region")
            continue
    
        # sample uniformly two reads
        reads = reads.sample(2)
        read_anchor = reads.iloc[0]
        read_pos = reads.iloc[1]
        ed_pos = editdistance.eval(read_anchor["read"], read_pos["read"])
        # write to tsv file
        with open(outfile, "a") as f:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("positive", read_anchor["id"], read_anchor["read"], read_pos["id"], read_pos["read"], read_anchor["start"], read_anchor["end"], ed_pos))
        count += 1

        # sample a negative sample for the anchor
        if indx == 0:
            tsv_neg = tsv2
        else:
            tsv_neg = tsv1
        
        # sample uniformly a read from the negative tsv file
        tsv_neg = tsv_neg[(tsv_neg["strand"] == s) & (tsv_neg["start"] == start) & (tsv_neg["end"] == end)]
        if tsv_neg.shape[0] == 0:
            print("No reads that span the region")
            continue
        read_neg = tsv_neg.sample(1)
        read_neg = read_neg.iloc[0]
        ed_neg = editdistance.eval(read_anchor["read"], read_neg["read"])
        if ed_neg < 1:
            continue
        # write to tsv file
        if ed_neg <= ed_pos:
            with open(outfile, "a") as f:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("hard negative", read_anchor["id"], read_anchor["read"],read_neg["id"], read_neg["read"], read_anchor["start"], read_anchor["end"], ed_neg))
            count += 1
        else:
            with open(outfile, "a") as f:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("negative", read_anchor["id"], read_anchor["read"],read_neg["id"], read_neg["read"], read_anchor["start"], read_anchor["end"], ed_neg))
            count += 1

if __name__ == "__main__":
    main()


    

    
        



    

    
