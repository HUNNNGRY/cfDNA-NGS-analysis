# 
me <- read.table("/BioII/lulab_b/baopengfei/meth.txt",header = T,  sep = "\t")
me$ID <- me$X
me$ID <- sub("Methylation|","",me$ID,fixed = T)
me$ID <- unlist(lapply(strsplit(me$ID,"\\|"),function(x) x[1]))
me$ensg.full <- gsub("promoter_","",me$ID)
ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/gene.bed",header = F, sep = "\t")
head(ref)
table(ref$V4 %in% me$ensg.full)

me$chr <- ref$V1[match(me$ensg.full,ref$V4)]


ref2 <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene_v2.txt",header = T, sep = "\t")
ref2 <- ref2[,c("ensg.full","type","name","ENTREZID")]
library(dplyr)
me <- left_join(me,ref2)
write.table(me,"/BioII/lulab_b/baopengfei/meth_candidate.txt",quote = F,col.names = T,row.names = F,sep = "\t")




# ref2 <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv",header = T, sep = "\t")
# table(ref2$CANCER_TYPE)
# 
# 
# gc <- ref2[ref2$CANCER_TYPE=="ST",]
# summary(gc$X._SAMPLES_COHORT)
# 
# gc <- gc[order(gc$X._SAMPLES_COHORT,decreasing = T),]
# gc <- gc[!duplicated(gc$SYMBOL),]
# gc <- gc[,c("SYMBOL","X._SAMPLES_COHORT")]

help.start() 
library(Xtail)
Xtail