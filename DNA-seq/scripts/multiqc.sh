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
workdir=/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq   # GSE*****, change this when processing other seq types
dataset=$1
outdir="${workdir}/output/${dataset}/multiqc"
mkdir ${outdir}

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

## bowtie2 align log (mapping ratio)
#for i in `cat metadata/$dataset/sample_ids.txt`;
#do 
#mv -f output/$dataset/log/$i/bowtie2.log output/$dataset/log/$i/${i}-bowtie2.log 
#done 
#multiqc -f -o ${outdir} -n align_multiqc output/$dataset/log/*/*bowtie2*

## bwa align log (mapping ratio) 好像没法QC ，multiQC目前也不支持

## gatk dedup log 
multiqc -f -o ${outdir} -n dedup_multiqc output/$dataset/log/*/*dedup*

echo "end at `date`"
