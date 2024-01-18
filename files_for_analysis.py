import os

ref_seq = "data/ON563414.3.fasta"
hifi_reads = "data/test_datasets/pox_B.1.17/IIb_B.1.17/EPI_ISL_15896351"

# for fasta 
os.system("minimap2 -ax map-hifi {} {}.fastq  > {}_analysis.sam".format(ref_seq, hifi_reads, hifi_reads))
os.system("samtools view -@ n -Sb -o {}_analysis.bam {}_analysis.sam".format(hifi_reads, hifi_reads))
os.system("samtools sort -@ n -o {}_analysis.sorted.bam {}_analysis.bam".format(hifi_reads, hifi_reads))
os.system("samtools index {}_analysis.sorted.bam".format(hifi_reads))