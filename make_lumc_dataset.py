
import pandas as pd
import os
import sys
import argparse




def map_to_correct_region(start, end):
    # genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
    #                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
    #                     (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
    #                      (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
    #                     (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
    #                      (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
    #                       (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    genomic_regions = [(13055,13960), (13803, 14718), (14600, 15599), (15371, 16311)]
    distances = {}
    
    for i in range(len(genomic_regions)):
        # calculate distance from start points
        dist_start = abs(genomic_regions[i][0] - start)
        distances[i] = dist_start

    # find the key with minimum distance
    min_key = min(distances, key=distances.get)

    return genomic_regions[min_key][0], genomic_regions[min_key][1]


def parse_fastq(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l!= len(lines)):
        read_name = lines[l].split(" ")[0][1:]
        l+=1
        read = lines[l]
        l+=2
        # quality_score = lines[l]
        l+=1
        reads[read_name] = read

    return reads

def main():
    parser = argparse.ArgumentParser(description="Map LUMC regions")
    parser.add_argument('--fastq', dest = 'fastq', required=True, type=str, help="")
    parser.add_argument('--bed', dest = 'bed', required=True, type=str, help="")
    parser.add_argument('--out', dest = 'out', required=True, type=str, help="")
    args = parser.parse_args()

    fastq_file = args.fastq
    bed_file = args.bed
    outfile = args.out

    # if outfile exists, delete it
    if os.path.exists(outfile):
        os.remove(outfile)
    # create it again
    open(outfile, 'a').close()


    # Read in the bed file
    bed = pd.read_csv(bed_file, sep="\t", header=None)
    bed.columns = ["ref_id", "start", "end", "id", "score", "strand"]
    # Read in the fastq file
    reads = parse_fastq(fastq_file)

    for read_name in reads.keys():
        row = bed[bed['id'] == read_name]
        # check if row is empty
        if row.empty:
            continue

        start = int(row['start'].values[0])
        end = int(row['end'].values[0])
        strand = row['strand'].values[0]

            
        # map to correct region
        start, end = map_to_correct_region(start, end)
        # write to tsv file
        with open(outfile, "a") as f:
            f.write(read_name + "\t" + reads[read_name] + "\t" + str(start) + "\t" + str(end) + "\t" + strand +"\n")
    
               
if __name__ == "__main__":
    sys.exit(main())