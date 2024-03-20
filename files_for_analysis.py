import os

ref_seq = "data/SARS-CoV-2-NC_045513.fa"
hifi_reads = "P0835SEQ00003-HDF6.bc1012--bc1035"
outdir = "data/cerbaresearch/Pacbio-Spike-5-fragment-CDC"

os.system("minimap2 -ax map-hifi {} {}/{}.fastq  > {}/{}_analysis.sam".format(ref_seq, outdir, hifi_reads, outdir, hifi_reads))
os.system("samtools view -@ n -Sb -o {}/{}_analysis.bam {}/{}_analysis.sam".format(outdir, hifi_reads, outdir, hifi_reads))
os.system("samtools sort -@ n -o {}/{}_analysis.sorted.bam {}/{}_analysis.bam".format(outdir, hifi_reads, outdir, hifi_reads))
os.system("samtools index {}/{}_analysis.sorted.bam".format(outdir, hifi_reads))
os.system("bedtools bamtobed -i {}/{}_analysis.sorted.bam > {}/{}_reads.bed".format(outdir,hifi_reads, outdir, hifi_reads))