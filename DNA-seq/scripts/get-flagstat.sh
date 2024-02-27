#!/bin/bash
dst=$1
samples=$2
##count hg38
( echo -e "sample""\t""total_reads""\t""mapped_reads""\t""paired_reads""\t""properPaired_reads""\t""mapped_ratio"; \
for sample in `cat ${samples}`; \
do awk -v sample=$sample \
' BEGIN {ORS=""} NR==1 {print sample "\t" $1 "\t"; total=$1} NR==5 {print $1 "\t"; map=$1} NR==6 {print $1 "\t"} NR==9 {print $1 "\t" map/total "\n"} ' output/${dst}/flagstat/${sample}.txt ;done ) \
>> output/$dst/multiqc/align-hg38-ratio.txt
