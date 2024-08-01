# canu-2.2/bin/canu -correct -d Error_corrected_data/hicanu/lumc/03_50/ -p output.canu useGrid=false maxMemory=90g genomeSize=29.8k saveReads=True stopOnLowCoverage=0 -pacbio data/lumc_data/natural_mixtures/03_50/trimmed_03_50.fastq

canu-2.2/bin/canu -d Error_corrected_data/hicanu/lumc/03_50_assembly/ -p output.canu useGrid=false maxMemory=90g genomeSize=29.8k saveReads=True stopOnLowCoverage=0 -pacbio data/lumc_data/natural_mixtures/03_50/trimmed_03_50.fastq

# gzip -d -c Error_corrected_data/hicanu/lumc/03_50/output.canu.correctedReads.fasta.gz > Error_corrected_data/hicanu/lumc/03_50/output.canu.correctedReads.fasta

# bash files_for_analysis.sh data/nCoV-2019.reference.fasta output.canu.correctedReads Error_corrected_data/hicanu/lumc/03_50

