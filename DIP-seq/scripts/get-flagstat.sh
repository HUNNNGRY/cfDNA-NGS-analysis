#!/bin/bash
dst=$1
threads=$2
samples=$3

if [ -d output/$dst/flagstat ];then
	echo "flagstats dir exists""\t"
else
	echo "create flagstats dir... "
	mkdir output/$dst/flagstat
fi

# 1.samtools flagstats
for sample in `cat $samples`; \
do samtools flagstat -@ $threads output/${dst}/unmapped/${sample}_lambda.bam > output/$dst/flagstat/${sample}_lambda.txt ;
samtools flagstat -@ $threads output/${dst}/02bam/${sample}.bam > output/$dst/flagstat/${sample}.txt ;
done 

# 2.counts reads   不需要并行
#对应bowtie2的bam里只包含map上的reads，所以无法根据单一bam的flagstat算mapped ratio
##count hg38
( echo -e "sample""\t""total_reads""\t""mapped_reads""\t""paired_reads""\t""properPaired_reads"; \
for sample in `cat $samples`; \
do awk -v sample=$sample \
' BEGIN {ORS=""} NR==1 {print sample "\t" $1 "\t"; total=$1} NR==5 {print $1 "\t"; map=$1} NR==6 {print $1 "\t"} NR==9 {print $1 "\n"} ' output/${dst}/flagstat/${sample}.txt ;done ) \
>> output/$dst/multiqc/align-hg38-ratio.txt

##count lambda
( echo -e "sample""\t""lambda_reads""\t"; \
for sample in `cat $samples`; \
do awk -v sample=$sample \
' BEGIN {ORS=""} NR==1 {print sample "\t" $1 "\n"} ' output/${dst}/flagstat/${sample}_lambda.txt ;done ) \
>> output/$dst/multiqc/align-lambda-ratio.txt

