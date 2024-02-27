#!/bin/bash
#SBATCH -J multiqc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --output=/data/baopengfei/exOmics/DIP-seq/multiqc.out
#SBATCH --error=/data/baopengfei/exOmics/DIP-seq/multiqc.err

workdir=/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq   # GSE*****, change this when processing other seq types
dataset=$1
outdir="${workdir}/output/${dataset}/multiqc"

cd  ${workdir}
if [ -d $outdir ];then
	echo "multiqc dir exists"
else
	echo "create multiqc dir..."
	mkdir $outdir	
fi

# RUN
echo "start at `date`"
echo "find multiqc at `which multiqc`"

dataset=$1
for dir in {01fastp,02bam,04bam_dedup}
do
multiqc -f -o ${outdir} -n ${dir}_multiqc ${workdir}/output/${dataset}/$dir/log/*
done

echo "end at `date`"
