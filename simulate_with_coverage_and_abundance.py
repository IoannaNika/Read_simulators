

import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
import json

def create_template(directory, s_id, n_reads, bed_file):
    os.system('python make_template_from_fasta-PS.py --dir {} --bed_file {} --seq_id {} --cnt_templates {}'.format(directory, bed_file, s_id, n_reads))
    return

def create_template_with_pf(directory, s_id, n_reads, bed_file, primer_file):
    os.system('python make_template_from_fasta.py --dir {} --bed_file {} --primer_file {} --seq_id {} --cnt_templates {}'.format(directory, bed_file, primer_file, s_id, n_reads))
    return

def simulate_ONT_HQ_reads(directory, identifier, ): 
    os.system("pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-ONT-HQ.model  --template {}/{}.template  --seed 20 --prefix {}/{}".format(directory, identifier, directory, identifier))
    return

def simulate_hifi_reads(directory, identifier, cores):
    os.system("pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-SEQUEL.model  --template {}/{}.template --pass-num 10 --seed 20 --prefix {}/{}".format(directory, identifier, directory, identifier))
    os.system("samtools view -bS {}/{}.sam > {}/{}.bam".format(directory, identifier, directory, identifier))
    os.system("ccs {}/{}.bam -j {} {}/{}.fastq".format(directory, identifier, cores, directory, identifier))
    return

def main():
    parser = argparse.ArgumentParser(description="Simulates pacbio hifi amplicon reads")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    parser.add_argument('--cores', dest = 'cores', default=8, required=False, type=int, help="cpu cores for the simulation")
    parser.add_argument('--strategy', dest = 'strategy', default="pacbio-hifi", required=False, type=str, help="ONT-HQ or pacbio-hifi")   
    parser.add_argument('--ids_and_ab', type=str, required=True, help='A dictionary of ids and their relative abundances (%)')
    parser.add_argument('--coverage', type=int, required=True, help='The coverage to simulate')
    parser.add_argument('--outdir', type=str, required=False,default="", help='The output directory')
    parser.add_argument('--bed_file', type=str, required=True, help='The bed file with primer locations')
    parser.add_argument('--primer_file', type=str, required=False, default=None, help='The primer file, can be ommitted if the bed file contains primer information')
    args = parser.parse_args()

    data_dir = args.dir
    cores = args.cores
    strategy = args.strategy

    with open(args.ids_and_ab) as json_file:
        ids_and_abundances = json.load(json_file)
        
    coverage = args.coverage
    primer_file = args.primer_file
    bed_file = args.bed_file


    # open file to write the ids and their abundances
    out_file = open(args.outdir + "_ids_and_abundances.tsv", "w")
    
    for s_id in ids_and_abundances.keys():
        abundance = ids_and_abundances[s_id]
        ab_fraction = abundance #/ 100
        n_reads = int(ab_fraction * coverage)

        out_file.write(s_id + "\t" + str(n_reads) + "\n")

        if primer_file == None:
            create_template(data_dir, s_id, n_reads, bed_file)
        else:
            create_template_with_pf(data_dir, s_id, n_reads, bed_file, primer_file)
        

        if strategy == "pacbio-hifi":
            simulate_hifi_reads(data_dir, s_id, cores)
            
        if strategy == "ONT-HQ": 
            simulate_ONT_HQ_reads(data_dir, s_id)
        
    

if __name__ == "__main__":
    sys.exit(main())

    

    


