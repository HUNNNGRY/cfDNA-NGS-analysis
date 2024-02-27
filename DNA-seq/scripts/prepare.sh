#!/bin/bash

workdir=/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq   # GSE*****, change this when processing other seq types
dataset=$1
mode=$2
logdir="${workdir}/log/${dataset}/"
metadir="${workdir}/metadata/${dataset}/"
datadir="${workdir}/data/${dataset}/"
mkdir ${logdir}
mkdir ${metadir}
mkdir ${datadir}

# RUN
echo "start at `date`" 

cd ${workdir}
cp /Share2/home/lulab1/cfDNA-seq/raw/DNA-seq/${dataset}/meta_data/sample_ids.txt ${metadir}/
cp run/PRJNA241385.sh run/${dataset}.sh
chmod 775 run/${dataset}.sh
cp config/PRJNA241385.yaml config/${dataset}.yaml
chmod 775 config/${dataset}.yaml

if [ "$mode" == "se" ] || [ "$mode" == "SE" ]; then 
	for i in `cat ${metadir}/sample_ids.txt`
	do 
	ln -s /Share2/home/lulab1/cfDNA-seq/raw/DNA-seq/${dataset}/fastq/${i}.fastq.gz ${datadir}/${i}.fastq.gz
	done
elif [ "$mode" == "pe" ] || [ "$mode" == "PE" ]; then
	for i in `cat ${metadir}/sample_ids.txt`
	do 
	ln -s /Share2/home/lulab1/cfDNA-seq/raw/DNA-seq/${dataset}/fastq/${i}_1.fastq.gz ${datadir}/${i}_1.fastq.gz
	ln -s /Share2/home/lulab1/cfDNA-seq/raw/DNA-seq/${dataset}/fastq/${i}_2.fastq.gz ${datadir}/${i}_2.fastq.gz
	done	
fi

echo "end at `date`" 