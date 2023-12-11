pbsim --strategy templ --method errhmm --errhmm data/pbsim3_models/ERRHMM-SEQUEL.model  --template data/hcov_global_2023-11-16_09-28/B.1/EPI_ISL_445248.template --pass-num 10 --seed 20 --prefix EPI_ISL_445248
samtools view -bS EPI_ISL_445248.sam > EPI_ISL_445248.bam
ccs EPI_ISL_445248.bam -j 8 EPI_ISL_445248.ccs.fastq

