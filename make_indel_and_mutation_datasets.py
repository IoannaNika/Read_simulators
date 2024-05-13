import argparse
import sys
import os
import random



def make_read_pairs(maf_file, fastq_file, output_directory):
    maf_file_dict = read_maf_file(maf_file)
    fastq_file_dict = read_fastq_file(fastq_file)

    # open the output file
    output_file = open(output_directory + "/all_test_pairs.tsv", "w")
    # write the header: read_1	read_2	label	genomic_region
    output_file.write("read_1\tread_2\tlabel\tgenomic_region\n")

    # make /reads directory within the output directory
    os.system("mkdir {}/reads".format(output_directory))

    # make all read pairs
    for i, s_id in enumerate(fastq_file_dict):
        # get the sequence from the fastq file
        sequence = fastq_file_dict[s_id]
        # get the amplicon from the maf file
        normal_id = maf_file_dict[s_id]

        # write the sequence to a file
        with open("{}/reads/{}.fasta".format(output_directory, normal_id), "w") as read_file:
            read_file.write(">" + normal_id + "\n")
            read_file.write(sequence + "\n")
            read_file.close()

        for j, s_id_2 in enumerate(fastq_file_dict):
            # get the sequence from the fastq file
            sequence_2 = fastq_file_dict[s_id_2]
            # get the amplicon from the maf file
            normal_id_2 = maf_file_dict[s_id_2]

            if s_id == s_id_2:
                continue

            if j <= i:
                continue

            if "modified" in normal_id and "modified" in normal_id_2:
                label = "1"
            else:
                if "modified" not in normal_id and "modified" not in normal_id_2:
                    label = "1"
                else:
                    label = "0"
            if label == "1": 
                continue

            output_file.write(normal_id + "\t" + normal_id_2 + "\t" + label + "\t" + "54_1183" + "\n")

            with open("{}/reads/{}.fasta".format(output_directory, normal_id_2), "w") as read_file:
                read_file.write(">" + normal_id_2 + "\n")
                read_file.write(sequence_2 + "\n")
                read_file.close()
    
    return

def read_maf_file(maf_file):
    maf_file_dict = {}
    # read the maf file
    with open(maf_file, "r") as in_file:
        lines = in_file.readlines()
    
    count = 0
    while count < len(lines):
        if lines[count].startswith("a"):
            # get the amplicon
            normal_id = lines[count+1].split(" ")[1]
            # get the sequence
            s_id =  "/".join(lines[count+2].split(" ")[1].split("/")[:-1])
            maf_file_dict[s_id] = normal_id
        count += 3

    return maf_file_dict


def read_fastq_file(fastq_file):
    # link id to sequence
    fastq_file_dict = {}
    # read the fastq file
    with open(fastq_file, "r") as in_file:
        lines = in_file.readlines()
    
    count = 0

    while count < len(lines):
        # get the id
        s_id = lines[count].replace("\n", "").split(" ")[0][1:-4]
        # get the sequence
        sequence = lines[count+1].replace("\n", "")
        fastq_file_dict[s_id] = sequence
        count += 4
    
    return fastq_file_dict



def simulate_hifi_reads(directory, identifier, cores):
    os.system("pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-SEQUEL.model  --template {}/{}.template --pass-num 10 --prefix {}/{}".format(directory, identifier, directory, identifier))
    os.system("samtools view -bS {}/{}.sam > {}/{}.bam".format(directory, identifier, directory, identifier))
    os.system("ccs {}/{}.bam -j {} {}/{}.fastq".format(directory, identifier, cores, directory, identifier))
    return


def apply_insertion(sequence, number):
    options = ["A", "C", "G", "T"]
    
    for i in range(number):
        # randomly select a position in the sequence
        position = random.randint(0, len(sequence))
        # randomly select a base to insert
        base = random.choice(options)

        # insert the base at the position
        sequence = sequence[:position] + base + sequence[position:]

    return sequence



def apply_deletion(sequence, number):

    for i in range(number):
        # randomly select a position in the sequence
        position = random.randint(0, len(sequence))

        # delete the base at the position
        sequence = sequence[:position] + sequence[position+1:]

    return sequence

def apply_mutation(sequence, number):
    options = ["A", "C", "G", "T"]

    for i in range(number):
        # randomly select a position in the sequence
        position = random.randint(0, len(sequence)-1)
        
        # get the base at the position
        base = sequence[position]

        # randomly select a base to replace the original base
        while True:
            new_base = random.choice(options)
            if new_base != base:
                break

        # mutate the base at the position
        sequence = sequence[:position] + new_base + sequence[position+1:]
    
    return sequence


def main():
    parser = argparse.ArgumentParser(description="Creates amplicon fragments, with the specified number of operations")
    parser.add_argument('--input_amplicon', dest = 'input_amplicon', required=True, type=str, help="path to input amplicon, fasta format")
    parser.add_argument('--operation', dest = 'operation', required=True, type=str, help="operation to perform on amplicon, either 'insertion', deletion or 'mutation'")
    parser.add_argument('--number', dest = 'number', required=True,  type=int, help="number of insertions, deletions or mutations to perform")
    parser.add_argument('--output_file', dest = 'output_file', required=True, type=str, help="path to output file, fasta format")
    parser.add_argument('--coverage', dest = 'coverage', required=False, default=100, type=int, help="coverage  for the dataset")
    parser.add_argument('--cores', dest = 'cores', required=False, default=1, type=int, help="number of cores to use for CCS")
    args = parser.parse_args()

    # read the fasta input file, must contain only one sequence
    with open(args.input_amplicon, "r") as in_file:
        lines = in_file.readlines()
        template_header = lines[0].replace("\n", "")
        template_sequence = "".join(lines[1:]).replace("\n", "")

    # if output file already exists, remove it
    if os.path.exists(args.output_file):
        os.remove(args.output_file)

   
    output_file = open(args.output_file, "w")

    # in the output file, simulate with hald the coverage the original sequence
    half_n = int(args.coverage / 2)
    rest_n = args.coverage - half_n
    for i in range(half_n):
        output_file.write(template_header + "_" + str(i) + "\n")
        output_file.write(template_sequence + "\n")

    # close the file
    output_file.close()

    for i in range(rest_n):
    
        # perform the operation
        if args.operation == "insertion":
            
            modified_sequence = apply_insertion(template_sequence, args.number)

        elif args.operation == "deletion":
            # perform deletion
            modified_sequence = apply_deletion(template_sequence, args.number)
        elif args.operation == "mutation":
            
            modified_sequence = apply_mutation(template_sequence, args.number)
        else:
            print("Invalid operation, must be either 'insertion', 'deletion' or 'mutation'")
            return 1
    
        # write the modified sequence to the output file
        output_file = open(args.output_file, "a")
        output_file.write(template_header + "_modified_" + str(i) + "\n")
        output_file.write(modified_sequence + "\n")
        output_file.close()


    # simulate sequencing errors: pacbio hi-fi reads
    directory =  os.path.dirname(args.output_file)
    identifier = os.path.basename(args.output_file).split(".")[0]
    simulate_hifi_reads(directory, identifier, args.cores)

    # read maf file
    maf_file = directory + "/" + identifier + ".maf"
    # read the fastq file
    fastq_file = directory + "/" + identifier + ".fastq"

    # make all read pairs
    make_read_pairs(maf_file, fastq_file, os.path.join(directory, "dataset"))
   

    
if __name__ == "__main__":
    sys.exit(main())