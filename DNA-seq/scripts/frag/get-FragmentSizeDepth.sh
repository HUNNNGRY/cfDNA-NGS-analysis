#!/bin/bash 
sample_id=$1
outdir=$2
mkdir -p ${outdir}/depth/log

for size in {long,short}
do
bam="${outdir}/${sample_id}-${size}.bam"
gappedbins="genome/hg38.bins.100kb.norepeats.bed"
bins="genome/hg38.bins.100kb.bed"
chromsize="genome/chrom.size"
coverage="${outdir}/depth/${sample_id}_${size}.bed"
log="${outdir}/depth/log/${sample_id}_${size}.log"

bedtools coverage -g ${chromsize} -sorted -a ${gappedbins} -b ${bam} > ${coverage}.tmp ;\
        bedtools map -a ${bins} -b ${coverage}.tmp  -c 7 -o mean > ${coverage} ;\
        rm ${coverage}.tmp > ${log} 2>&1
done
