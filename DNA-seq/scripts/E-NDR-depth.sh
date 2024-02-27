#!/bin/bash
#BSUB -J lulab
#BSUB -q TEST-1U
#BSUB -n 18
#BSUB -o /BioII/lulab_b/baopengfei/projects/multi-omics-explore/log/submit-STAR-Fusion-SAMPLE.log 
#BSUB -e /BioII/lulab_b/baopengfei/projects/multi-omics-explore/log/submit-STAR-Fusion-SAMPLE.err

PATH=/BioII/lulab_b/baopengfei/anaconda3/envs/exvariance/bin:/BioII/lulab_b/baopengfei/miniconda3/bin:/BioII/lulab_b/baopengfei/biosoft/manta-1.6.0.centos6_x86_64/bin:/BioII/lulab_b/baopengfei/biosoft/gatk-4.2.0.0:/BioII/lulab_b/baopengfei/biosoft/TrimGalore-0.6.6:/BioII/lulab_b/baopengfei/biosoft/ngsplot/bin:/BioII/lulab_b/baopengfei/localperl/bin:/BioII/lulab_b/baopengfei/anaconda3/condabin:/BioII/lulab_b/baopengfei/bin/bin/:/BioII/lulab_b/baopengfei/anaconda3/bin:/BioII/lulab_b/baopengfei/biosoft/TRUST4:/BioII/lulab_b/baopengfei/jdk1.8.0_251/bin:/BioII/lulab_b/baopengfei/.aspera/connect/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.9.6-1-centos_linux64/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.10.5-centos_linux64/bin:/BioI/lulab_b/baopengfei/miniconda3/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/apps/bin:/apps/homer/bin:/apps/ucscKentUtilities:/apps/RSEM/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/texlive/2018/bin/x86_64-linux:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:/BioII/lulab_b/baopengfei/biosoft:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:$PATH

source ~/.bashrc

cd  /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq
ls -1   /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/bam-sorted-deduped-merged/*-passQC*.bam  >  output/lulab/NDR/single_base_dep-v3/bam.list
for i in {promoter2000TSS2000_bloodTx_unexp,promoter2000TSS2000_bloodTx_fpkm-001-01,promoter2000TSS2000_bloodTx_fpkm-01-5,promoter2000TSS2000_bloodTx_fpkm-5-30,promoter2000TSS2000_bloodTx_fpkm-30,promoter2000exon1end2000_bloodTx_unexp,promoter2000exon1end2000_bloodTx_fpkm-001-01,promoter2000exon1end2000_bloodTx_fpkm-01-5,promoter2000exon1end2000_bloodTx_fpkm-5-30,promoter2000exon1end2000_bloodTx_fpkm-30}
do echo $i ;done | parallel -k -I % -j 10 "samtools depth -a -b ref/gtf/bloodNDR/%.bed -f output/lulab/NDR/single_base_dep-v3/bam.list > output/lulab/NDR/single_base_dep-v3/%.txt"
