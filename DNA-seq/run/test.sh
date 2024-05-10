#!/bin/bash
#BSUB -J lulab
#BSUB -q TEST-1U
#BSUB -o /BioII/lulab_b/baopengfei/projects/%J.log
#BSUB -e /BioII/lulab_b/baopengfei/projects/%J.err
#BSUB -n 2  
#BSUB -R "span[hosts=1]"

# submit to new A cluster
#PATH=/BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin:/BioII/lulab_b/baopengfei/miniconda3/bin:/BioII/lulab_b/baopengfei/biosoft/manta-1.6.0.centos6_x86_64/bin:/BioII/lulab_b/baopengfei/biosoft/gatk-4.2.0.0:/BioII/lulab_b/baopengfei/biosoft/TrimGalore-0.6.6:/BioII/lulab_b/baopengfei/biosoft/ngsplot/bin:/BioII/lulab_b/baopengfei/localperl/bin:/BioII/lulab_b/baopengfei/anaconda3/condabin:/BioII/lulab_b/baopengfei/bin/bin/:/BioII/lulab_b/baopengfei/anaconda3/bin:/BioII/lulab_b/baopengfei/biosoft/TRUST4:/BioII/lulab_b/baopengfei/jdk1.8.0_251/bin:/BioII/lulab_b/baopengfei/.aspera/connect/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.9.6-1-centos_linux64/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.10.5-centos_linux64/bin:/BioI/lulab_b/baopengfei/miniconda3/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/apps/bin:/apps/homer/bin:/apps/ucscKentUtilities:/apps/RSEM/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/texlive/2018/bin/x86_64-linux:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:/BioII/lulab_b/baopengfei/biosoft:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:$PATH
# source ~/.bashrc

# add some commands here
mkdir -p log/$dst
echo -e "start at `date`"

snakemake --rerun-incomplete --keep-going --printshellcmds --reason --use-conda --nolock --latency-wait 200 --restart-times 1 --jobs 200  \
	--snakefile snakemake/DNA-seq-pe.snakemake \
	--configfile config/${dst}.yaml \
	--cluster-config config/cluster-lsf.json \
	--cluster "sbatch -N {cluster.nodesNum} -x {cluster.excludeNodesName} -n {cluster.threads} -J {cluster.jobname} -p {cluster.partition} -o {cluster.output} -e {cluster.error} " \
	> log/${dst}/run-${dst}.log 2>&1
echo "end at `date`"
