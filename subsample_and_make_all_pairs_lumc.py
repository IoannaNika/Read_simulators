import argparse
import sys
import os
import pandas as pd
import editdistance

def write_sequence_to_file(full_read_id, sequence, outdir_read_dir):
    # check if file exists
    if os.path.exists(outdir_read_dir + "/" + full_read_id + ".fasta"):
        os.remove(outdir_read_dir + "/" + full_read_id + ".fasta")

    if not os.path.exists(outdir_read_dir + "/" + full_read_id + ".fasta"):
        with open(outdir_read_dir + "/" + full_read_id + ".fasta", "w") as f:
            f.write(">" + full_read_id + "\n")
            f.write(sequence + "\n")
    return

def main(): 
    parser = argparse.ArgumentParser(description="Subsample reads and make all pairs")
    parser.add_argument('--input_dir', dest = 'input_dir', required=True, type=str, help="path to input file, fasta format")
    parser.add_argument('--n', dest = 'n', required=False, default = 100, type=int, help="coverage to keep")
    parser.add_argument('--seed', dest = 'seed', required=False, default = 1, type=int, help="seed fur shuffling pool of reads to select from")
    parser.add_argument('--gr_start', dest = 'gr_start', required=False, default =-1, type=int, help="datasets_for_genomic_region")
    args = parser.parse_args()

    n = args.n
    reads_tsv_file = args.input_dir + "/reads.tsv"
    reads = pd.read_csv(reads_tsv_file, sep='\t')

    if args.gr_start != -1:
        output_dir = args.input_dir + "/seed_" + str(args.seed) + "/subsampled_reads_" + str(n) + "_" + str(args.gr_start)
    else: 
        output_dir = args.input_dir + "/seed_" + str(args.seed) + "/subsampled_reads_" + str(n)
    # make a directory to store the subsampled reads
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = output_dir + "/reads.tsv"
    # make output file with 'w' mode
    open(output_file, "w").close()

    # header is not present in the reads.tsv file but: 
    # the first column is the read_id and the second column is the sequence the third column is the start position and the fourth column is the end position
    # the fifth column is the strand
    # add header to the reads dataframe

    reads.columns = ["read_id", "sequence", "start", "end", "strand"]
    if args.gr_start != -1:
        print(len(reads))
        print(reads["start"].unique())
        reads = reads[reads["start"] == args.gr_start]
        print(len(reads))

    # subsample the reads and leave only 100 reads per genomic region
    reads = reads.groupby(["start", "end"]).sample(n=args.n, random_state=args.seed).reset_index(drop=True)
    reads.to_csv(output_file, sep='\t', index=False)

    # for all reads in the reads.tsv file, write the read to a file
    for index, row in reads.iterrows():
        read_id = row["read_id"].replace("/", "_")
        sequence = row["sequence"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        full_read_id = read_id
        if strand == "-":
            sequence = sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
        seq_output_dir = output_dir + "/reads"
        # make seq_output_dir if it does not exist
        if not os.path.exists(seq_output_dir):
            os.makedirs(seq_output_dir)
        write_sequence_to_file(full_read_id, sequence, seq_output_dir)

    all_test_pairs = output_dir + "/all_test_pairs.tsv"
    # make output file with 'w' mode
    open(all_test_pairs, "w").close()

    # write header: label	id1	id2	start	end	edit_distance
    with open(all_test_pairs, "w") as f:
        f.write("label\tid1\tid2\tstart\tend\tedit_distance\n")
    
    # for all reads in the reads.tsv file, make all pairs
    for index, row in reads.iterrows():
        read_id = row["read_id"].replace("/", "_")
        sequence = row["sequence"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        for index2, row2 in reads.iterrows():
            read_id2 = row2["read_id"].replace("/", "_")
            sequence2 = row2["sequence"]
            start2 = row2["start"]
            end2 = row2["end"]
            strand2 = row2["strand"]

            if index2 <= index:
                continue

            if start != start2:
                continue

            if read_id == read_id2:
                continue
            
            if strand2 == "-":
                sequence2 = sequence2[::-1].translate(str.maketrans("ATGC", "TACG"))
            if strand == "-":
                sequence = sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
           
            edit_distance = editdistance.eval(sequence, sequence2)
            with open(all_test_pairs, "a") as f:
                f.write("NA\t" + read_id + "\t" + read_id2 + "\t" + str(start) + "\t" + str(end) + "\t" + str(edit_distance) + "\n")

    
if __name__ == "__main__":
    main()

