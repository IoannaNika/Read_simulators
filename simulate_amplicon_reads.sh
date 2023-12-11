pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-SEQUEL.model  --template data/hcov_global_2023-11-16_09-28/EG.2.3/EPI_ISL_18127471.template --depth 1000 --pass-num 10
# samtools view -bS sd.sam > sd.bam
# ccs sd.bam -j 8 unzipped.ccs.fastq