#!/bin/bash
#SBATCH -J multiqc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --output=/data/baopengfei/exOmics/DIP-seq/multiqc.out
#SBATCH --error=/data/baopengfei/exOmics/DIP-seq/multiqc.err


# Get the software
#newPATH=/data/baopengfei/anaconda3/envs/exvariance/bin
#PATH=$newPATH:$PATH

# Get var
workdir=/BioII/lulab_b/baopengfei/projects/exOmics/BS-seq   # GSE*****, change this when processing other seq types
dataset=$1
mkdir output/$dataset/multiqc
outdir="output/$dataset/multiqc"

# RUN
echo "start at `date`"
echo "PATH: /n $PATH"
echo "workdir: /n $workdir"
echo "dataset: /n $dataset"

cd  ${workdir}

## trim_galore qc0 log (optional)
multiqc -f -o ${outdir} -n qc0_multiqc output/$dataset/qc0/*/*

## trim_galore trim log (flter reads ratio)
multiqc -f -o ${outdir} -n trim_multiqc output/$dataset/log/*/*statistics*

## trim_galore qc1 log (remaining reads number)
multiqc -f -o ${outdir} -n qc1_multiqc output/$dataset/qc1/*/*

## bismark align log (mapping ratio)
for i in `cat metadata/$dataset/sample_ids.txt`;
do 
mv -f output/$dataset/log/$i/bismark-report.txt output/$dataset/log/$i/${i}_PE_report.txt 
done 

multiqc -f -o ${outdir} -n align_multiqc output/$dataset/log/*/*_PE_report.txt 

## bismark dedup log 
for i in `cat metadata/$dataset/sample_ids.txt`;
do
mv -f output/$dataset/bam/${i}.deduplication_report.txt output/$dataset/log/$i/${i}.deduplication_report.txt
done

multiqc -f -o ${outdir} -n dedup_multiqc output/$dataset/log/*/*.deduplication_report.txt

## bismark extract log (no useful info)
for i in `cat metadata/$dataset/sample_ids.txt`;
do 
mv -f output/$dataset/methylation/$i/*_splitting_report.txt output/$dataset/log/$i/
done 
multiqc -f -o ${outdir} -n extract_multiqc output/$dataset/log/*/*_splitting_report.txt

## bismark M-bias log (no useful info)
#for i in `cat metadata/$dataset/sample_ids.txt`;
#do 
#mv output/$dataset/methylation/$i/${i}.deduplicated.M-bias.txt output/$dataset/log/$i/
#done 
##multiqc -f -o ${outdir} -n bias_multiqc output/$dataset/log/*/*.deduplicated.M-bias.txt

echo "end at `date`"
