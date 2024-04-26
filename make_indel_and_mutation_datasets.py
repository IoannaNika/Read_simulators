import argparse
import sys
import os
import random

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
        position = random.randint(0, len(sequence))
        
        # get the base at the position
        base = sequence[position]

        # randomly select a base to replace the original base
        while True:
            new_base = random.choice(options)
            if new_base != base:
                break

        # mutate the base at the position
        sequence = sequence[:position] + base + sequence[position+1:]
    
    return sequence


def main():
    parser = argparse.ArgumentParser(description="Creates amplicon fragments")
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
        template_header = lines[0]
        template_sequence = "".join(lines[1:]).replace("\n", "")

    # if output file already exists, remove it
    if os.path.exists(args.output_file):
        os.remove(args.output_file)

   
    output_file = open(args.output_file, "w")

    # in the output file, simulate with hald the coverage the original sequence
    half_n = args.coverage // 2
    rest_n = args.coverage - half_n
    for i in range(half_n):
        output_file.write(template_header + "_" + str(i) + "\n")
        output_file.write(template_sequence + "\n")

    # close the file
    output_file.close()
    
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
    for i in range(rest_n):
        output_file.write(template_header + "_modified_" + str(i) + "\n")
        output_file.write(modified_sequence + "\n")
    output_file.close()


    # simulate sequencing errors: pacbio hi-fi reads
    directory = os.path.dirname(args.output_file)
    identifier = os.path.basename(args.output_file).split(".")[0]
    simulate_hifi_reads(directory, identifier, args.cores)

    
if __name__ == "__main__":
    sys.exit(main())