#plot boxplot for target DNA-seq vcf-filter result


setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

vcf <- read.table("./output/PRJNA687345/vcf-filtered/uniq_vcf.txt",header = F,stringsAsFactors = F)
colnames(vcf) <- c("sample","SNV")
vcf$cov <- as.numeric(lapply(strsplit(vcf$sample,"\\_"),function(x) x[1]))
vcf$cov <- vcf$cov/600

vcf$SNV.1x <- vcf$SNV/vcf$cov
vcf$cov <- as.factor(vcf$cov)

str(vcf)

# plot
library(ggplot2)

ggplot(data=vcf,aes(x=cov,y=SNV,group=cov,fill=cov)) +
  geom_boxplot(show.legend = T,) +
  geom_point(size=1) + 
  scale_y_continuous(trans='log10') + 
  #xlim(c(0,5000000)) +  
  #ylim(c(0,1000000))
  #scale_x_discrete(breaks=c("TMEFF2","NGFR","SEPTIN9","NDRG4","BMP3"),labels=c("TMEFF2","NGFR","SEPTIN9","NDRG4","BMP3")) +
  #scale_y_continuous(oob = scales::squish) +
  #stat_compare_means(comparisons = list(c("CRC","NC"))) +
  #labs(title="in-house cfMeDIP",x="CRC FDA Biomarker",y="transformed counts") +
  #geom_jitter(width = .2) +
  theme_bw()

ggplot(data=vcf,aes(x=cov,y=SNV.1x,group=cov,fill=cov)) +
  geom_boxplot(show.legend = T,) +
  geom_point(size=1) + 
  scale_y_continuous(trans='log10') + 
  #xlim(c(0,5000000)) +  
  #ylim(c(0,1000000))
  #scale_x_discrete(breaks=c("TMEFF2","NGFR","SEPTIN9","NDRG4","BMP3"),labels=c("TMEFF2","NGFR","SEPTIN9","NDRG4","BMP3")) +
  #scale_y_continuous(oob = scales::squish) +
  #stat_compare_means(comparisons = list(c("CRC","NC"))) +
  #labs(title="in-house cfMeDIP",x="CRC FDA Biomarker",y="transformed counts") +
  #geom_jitter(width = .2) +
  theme_bw()
