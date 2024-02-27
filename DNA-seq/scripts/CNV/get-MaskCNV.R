
# getMaskCNV
# last 210813 by pengfei (not finish yet)

setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
bedtools.path <- "/BioII/lulab_b/baopengfei/anaconda3/bin"

bt.coverage.a <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/hg38.bins.100kb.norepeats.bed"
bt.coverage.b <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/bam-sorted-deduped/CRC-PKU-10-wgs.bam"

bt.map.a <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/genome/hg38.bins.100kb.bed"


options(bedtools.path="/BioII/lulab_b/baopengfei/anaconda3/bin")
suppressPackageStartupMessages(library(bedtoolsr))

# cal mean depth
gap.mean <- bedtoolsr::bt.coverage(a = bt.coverage.a, b = bt.coverage.b, mean = T)

# norm to ratio:  bin mean depth/genome-wide mean depth
gap.mean$V5 <- gap.mean$V4/mean(gap.mean$V4)

write.table(gap.mean,"./CRC-PKU-10-wgs.gap.mean",col.names = F,row.names = F,quote = F,sep = "\t")

tmp <- "./CRC-PKU-10-wgs.gap.mean"
bin.mean <- bedtoolsr::bt.map(a = bt.coverage.a, b = tmp, c = 5, o = mean)

tmp <- read.table("output/lulab/matrix/CNVmaskDepthRatio_matrix_gene.txt",sep = "\t",header = T,row.names = 1)
tmp <- as.matrix(tmp)
quantile(tmp$CRC.PKU.12.wgs, c(.01, .5, .99)) 

mean(tmp$CRC.PKU.10.wgs)


