import argparse
import sys
import os

def parse_template(template_file):
    with open(template_file, 'r') as f:
        lines = f.readlines()
    
    # >EPI_ISL_408480:54_1183:0
    # TCGTCTCACTTGGTGGCTGGGACGACGCGGCGGCGCGTGCAGCGTGTAGGCTGGGACGTCACTTGGGAGGCTGGTGG
    
    template_dict = {}
    # read two lines at a time
    cnt = 0
    while cnt < len(lines):
        record_info = lines[cnt].strip()
        sequence = lines[cnt+1].strip()
        template_dict[record_info] = sequence
        cnt += 2
       
    return template_dict

def reverse_complement(sequence):
    return sequence[::-1].translate(str.maketrans("ATGC", "TACG"))

def modify_template(template_dict):
    modified_template_dict = {}
    for record_info, sequence in template_dict.items():
        # if the sequence is nubmer between 0 and 4 dont reverse complement
        record_info_digit = int(record_info.split(":")[-1])
        if record_info_digit <= 4:
            record_info = ">+_" + record_info[1:]
            modified_template_dict[record_info] = sequence
        else:
            record_info = ">-" + record_info[1:]
            modified_template_dict[record_info] = reverse_complement(sequence)
    return modified_template_dict


def write_template(modified_template_dict, output_file):
    # if the file exists, delete it
    if os.path.exists(output_file):
        os.remove(output_file)
    
    open(output_file, 'x').close()

    with open(output_file, 'w') as f:
        for record_info, sequence in modified_template_dict.items():
            f.write(record_info + "\n")
            f.write(sequence + "\n")

    return


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="directory with sequences.fasta which contains all sequences to be used for amplicon detection")
    args = parser.parse_args()

    directory = args.dir

    # get all files in directory and subdirectory that end in .template
    template_files = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            for file in os.listdir(os.path.join(root, dir)):
                if file.endswith(".template"):
                    template_files.append(os.path.join(root, dir, file))


    for template_file in template_files:
        template_dict = parse_template(template_file)
        modified_template_dict = modify_template(template_dict)
        write_template(modified_template_dict, template_file)


    
if __name__ == "__main__":
    sys.exit(main())