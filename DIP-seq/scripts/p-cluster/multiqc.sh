#!/bin/bash
#SBATCH -J multiqc
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --output=/data/baopengfei/exOmics/DIP-seq/multiqc.out
#SBATCH --error=/data/baopengfei/exOmics/DIP-seq/multiqc.err


# Get the software
newPATH=/data/baopengfei/anaconda3/envs/exvariance/bin
PATH=$newPATH:$PATH

# Get var
workdir=/data/baopengfei/exOmics/DIP-seq   # GSE*****, change this when processing other seq types
cd  ${workdir}

# RUN
echo "start at `date`"
echo $PATH

dataset=$1
for dir in {01fastp,02bam,04bam_dedup}
do
cd  ${workdir}/output/${dataset}/${dir}/log
multiqc -f .
cd -
done

#cd  ${workdir}/output/${dataset}/07counts
#multiqc --sample-filters *gene* -n geneqc -o ./ ./
#cd -

# for dataset in `ls /data/baopengfei/exOmics/DIP-seq/output`
# do
# 	for dir in {01fastp,02bam,04bam_dedup,07counts}
# 	do
# 	cd  ${workdir}/output/${dataset}/${dir}/log
# 	multiqc .
# 	cd -
# 	done
# done



echo "end at `date`"
