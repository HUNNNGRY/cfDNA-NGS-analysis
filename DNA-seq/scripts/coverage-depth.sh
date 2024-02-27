#!/bin/bash
dst=$1
threads=$2
samples=$3
if [ -d output/$dst/bam-sorted-deduped ]; then 
	echo "start from dedup bams"
	export bam_dir="output/$dst/bam-sorted-deduped"
else
	echo "start from usual bams"
	export bam_dir="output/$dst/bam"
fi


(echo -e "sample""\t""genome_depth""\t""coverage_depth""\t""coverage""\t""coverage_ratio"; \
cat ${samples} |  parallel -k -I% -j $threads " samtools depth ${bam_dir}/%.bam | awk '{sum+=\$3} END { print \"%\"  \"\t\" sum/3000000000 \"\t\" sum/NR \"\t\" NR \"\t\" NR/3000000000 } ' ")  >> output/$dst/multiqc/coverage-depth.txt


# awk in linux parallel need \
