#!/bin/bash
#BSUB -J lulab
#BSUB -q TEST-1U

PATH=/BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin:/BioII/lulab_b/baopengfei/miniconda3/bin:/BioII/lulab_b/baopengfei/biosoft/manta-1.6.0.centos6_x86_64/bin:/BioII/lulab_b/baopengfei/biosoft/gatk-4.2.0.0:/BioII/lulab_b/baopengfei/biosoft/TrimGalore-0.6.6:/BioII/lulab_b/baopengfei/biosoft/ngsplot/bin:/BioII/lulab_b/baopengfei/localperl/bin:/BioII/lulab_b/baopengfei/anaconda3/condabin:/BioII/lulab_b/baopengfei/bin/bin/:/BioII/lulab_b/baopengfei/anaconda3/bin:/BioII/lulab_b/baopengfei/biosoft/TRUST4:/BioII/lulab_b/baopengfei/jdk1.8.0_251/bin:/BioII/lulab_b/baopengfei/.aspera/connect/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.9.6-1-centos_linux64/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.10.5-centos_linux64/bin:/BioI/lulab_b/baopengfei/miniconda3/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/apps/bin:/apps/homer/bin:/apps/ucscKentUtilities:/apps/RSEM/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/texlive/2018/bin/x86_64-linux:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:/BioII/lulab_b/baopengfei/biosoft:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:$PATH

source ~/.bashrc


/BioII/lulab_b/baopengfei/biosoft/homer/bin/findMotifsGenome.pl \
	output/lulab/macs2/17CRC-censusPeaks.bed \
	hg38 \
	output/lulab/macs2/homer-peakBed/CRC \
	-p 12 -size 200 -mask  \
> output/lulab/macs2/homer-peakBed/CRC/findMotifsGenome.log 2>&1
