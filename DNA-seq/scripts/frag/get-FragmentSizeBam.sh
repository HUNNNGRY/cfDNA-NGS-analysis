#!/bin/bash

# get short and long fragment bam

id=$1
bam=$2
outdir=$3

mkdir -p ${outdir}/log

#dst=$2
#bam="output/lulab/bam-sorted-deduped/${id}.bam"
min=100
max=150
#outbam="output/lulab/test/${id}-short.bam"
outbam="${outdir}/${id}-short.bam"

(samtools view -H ${bam}; 
samtools view -F 4 -f 2 ${bam} | \
        awk -v min_size=${min} -v max_size=${max} \
        '{{if($9>0){{size=$9}}else{{size=-$9}};if(size>=min_size&&size<=max_size){{print $0}}}}') | \
        samtools view -b > ${outbam} 2> ${outdir}/log/get-short.err

min=151
max=220
#outbam="output/lulab/test/${id}-long.bam"
outbam="${outdir}/${id}-long.bam"
(samtools view -H ${bam}; 
samtools view -F 4 -f 2 ${bam} | \
        awk -v min_size=${min} -v max_size=${max} \
        '{{if($9>0){{size=$9}}else{{size=-$9}};if(size>=min_size&&size<=max_size){{print $0}}}}') | \
        samtools view -b > ${outbam} 2> ${outdir}/log/get-long.err
