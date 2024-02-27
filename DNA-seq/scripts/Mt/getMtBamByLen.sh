#!/usr/bin/bash
sample=$1
for min_size in `seq 20 20 200`
do
cd /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq
input="output/lulab/bam-sorted-deduped/${sample}.bam"
#min_size="20"
#declare -i min_size
declare -i max_size=${min_size}+20
echo ${min_size}
echo ${max_size}
echo ${sample}
output="output/lulab/test/MT-bam/${sample}_${min_size}_${max_size}.bam"

(samtools view -H  ${input}; samtools view -f 2 -F 4 ${input} | \
        awk -v min_size=${min_size} -v max_size=${max_size} \
        '{{if($9>0){{size=$9}}else{{size=-$9}};if(size>=min_size&&size<=max_size){{print $0}}}}') | \
        samtools view -b > ${output}
samtools index ${output}
done
