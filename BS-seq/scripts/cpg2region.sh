#!/bin/bash
## this scripts convert cpg signal into region signal (methylation beta value)
#workdir=$PWD
#type=$1
dataset=$1
surfix=$2
bed_dir=$3

pro="${bed_dir}/gencodev27_promoter_sort.bed"
gene="${bed_dir}/gencodev27_gene_sort.bed"
cgi="${bed_dir}/CGI_sort.bed"
#bed=$5

for sample in `ls ${dataset}/log`
do 
#input_fname="${sample}${surfix}"
#output_fname="${sample}.${type}"
tmpfile="${dataset}/methylation/${sample}/${sample}.tmp"

echo "start ${sample} at `date`"
zcat ${dataset}/methylation/${sample}/${sample}${surfix} | sort -k1,1 -k2,2n > ${tmpfile}
  
nohup bedtools map -a ${pro} -b ${tmpfile} -c 4 -o mean  > ${dataset}/methylation/${sample}/${sample}.promoter &
nohup bedtools map -a ${gene} -b ${tmpfile} -c 4 -o mean  > ${dataset}/methylation/${sample}/${sample}.gene &
nohup bedtools map -a ${cgi} -b ${tmpfile} -c 4 -o mean  > ${dataset}/methylation/${sample}/${sample}.CGI &

echo "end ${sample} at `date`"
done

rm -rf ${dataset}/methylation/*/*.tmp 
