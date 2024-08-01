

ref_seq="$1" 
hifi_reads="$2" 
outdir="$3"

minimap2-2.28_x64-linux/minimap2 -ax map-hifi ${ref_seq} ${outdir}/${hifi_reads}.fastq  > $outdir/${hifi_reads}_analysis.sam
samtools view -@ n -Sb -o ${outdir}/${hifi_reads}_analysis.bam ${outdir}/${hifi_reads}_analysis.sam
samtools sort -@ n -o ${outdir}/${hifi_reads}_analysis.sorted.bam ${outdir}/${hifi_reads}_analysis.bam
bedtools bamtobed -i ${outdir}/${hifi_reads}_analysis.sorted.bam > ${outdir}/${hifi_reads}_reads.bed
samtools index ${outdir}/${hifi_reads}_analysis.sorted.bam