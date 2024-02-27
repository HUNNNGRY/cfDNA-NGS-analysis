#!/bin/bash
#SBATCH -J test
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=/data/baopengfei/exOmics/DIP-seq/log/test/submit_snake.out
#SBATCH --error=/data/baopengfei/exOmics/DIP-seq/log/test/submit_snake.err


# Get the software
newPATH=/data/baopengfei/anaconda3/envs/exvariance/bin
PATH=$newPATH:$PATH

# Get var
workdir=/data/baopengfei/exOmics/DIP-seq   # GSE*****, change this when processing other seq types

dataset="test"   # GSE*****, change this when processing other datasets of DIP-seq

# RUN
cd  ${workdir}

echo "start at `date`"
echo $PATH

snakemake --rerun-incomplete --keep-going --printshellcmds --use-conda --nolock -j 20 --latency-wait 10 --restart-times 3 \
	-s ${workdir}/scripts/snakefile.py \
	--configfile ${workdir}/scripts/config/${dataset}.yaml \
	--cluster-config ${workdir}/scripts/cluster_config/cluster-slurm.json \
	--cluster "sbatch -n {cluster.threads} --job-name {cluster.jobname} --partition {cluster.partition} --output {cluster.output} --error {cluster.error} "  \
	> ${workdir}/log/${dataset}/run_snake.log 2>&1
 
echo "end at `date`"
