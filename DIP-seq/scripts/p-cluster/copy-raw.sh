#!/bin/bash
#SBATCH -J rsync
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --output=/data/baopengfei/lulab_cfmedip_snakemake/%j.out
#SBATCH --error=/data/baopengfei/lulab_cfmedip_snakemake/%j.err

# Get the software
cnodePATH=/BioII/lulab_b/baopengfei/anaconda3/envs/snakemake_medip/bin:/BioII/lulab_b/baopengfei/anaconda3/envs/cfea/bin:/BioII/lulab_b/baopengfei/anaconda3/bin:/BioII/lulab_b/baopengfei/miniconda3/bin:/BioII/lulab_b/baopengfei/anaconda3/condabin:/BioII/lulab_b/baopengfei/bin/bin:/BioII/lulab_b/baopengfei/anaconda3/bin:/BioII/lulab_b/baopengfei/biosoft/TRUST4:/BioII/lulab_b/baopengfei/jdk1.8.0_251/bin:/BioII/lulab_b/baopengfei/.aspera/connect/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.9.6-1-centos_linux64/bin:/BioII/lulab_b/baopengfei/app/sratoolkit.2.10.5-centos_linux64/bin:/BioI/lulab_b/baopengfei/miniconda3/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/apps/bin:/apps/homer/bin:/apps/ucscKentUtilities:/apps/RSEM/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/texlive/2018/bin/x86_64-linux:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin:/BioII/lulab_b/baopengfei/.local/bin:/BioII/lulab_b/baopengfei/bin
export PATH=$cnodePATH:$PATH


echo "start at `date`"
#rsync -av --append-verify --progress /BioII/lulab_b/baopengfei/2020proj/lulab_cfmedip_snakemake/data /data/baopengfei/lulab_cfmedip_snakemake/data
rsync -av --append-verify --progress /BioII/lulab_b/baopengfei/2020proj/lulab_cfmedip_snakemake/ref /data/baopengfei/lulab_cfmedip_snakemake/ &
rsync -av --append-verify --progress /BioII/lulab_b/baopengfei/2020proj/lulab_cfmedip_snakemake/output/lulab_cfmedip_ /data/baopengfei/lulab_cfmedip_snakemake/output/ 

echo "end at `date`"
