cd /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq

# 以下somatic-hg38_其实是下载时的文件夹前缀，比如somatic-hg38_1000g_pon.hg38.vcf.gz指的不是somatic，而是germline的1000g_pon.hg38.vcf.gz
# a panel of normals is simply a vcf of blacklisted sites flagged as recurrent artifacts
# The optional germline resource can be any VCF that contains an AF (population allele frequency) INFO field. The Broad Institute provides a version of gnomAD stripped of all fields except AF. The Broad Institute also 
# provides several panels of normals, but users with a large number (at least 50 or so) may benefit from generating their own panel with CreateSomaticPanelOfNornals
sample="CRC-PKU-10-wgs" 
cd /BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq
gatk Mutect2 \
  -R genome/fasta/hg38.fa \
  -I output/lulab/bam-sorted-deduped-RG/${sample}.bam \
  --germline-resource ref/SNV_ref/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ref/SNV_ref/somatic-hg38_1000g_pon.hg38.vcf.gz \
  --native-pair-hmm-threads 10 \
  -O output/lulab/mutect2/${sample}.vcf.gz

gatk FilterMutectCalls \
   -R genome/fasta/hg38.fa \
   -V output/lulab/mutect2/${sample}.vcf.gz \
   -f-score-beta 1 \
   -O output/lulab/mutect2/${sample}_FilterMutectCalls.vcf.gz 



gatk FilterMutectCalls \
   -R genome/fasta/hg38.fa \
   -V output/lulab/mutect2-vcf/${sample}.vcf.gz \
   -f-score-beta 1 \
   -O output/lulab/vcf-filtered/mutect2/${sample}.vcf.gz 
