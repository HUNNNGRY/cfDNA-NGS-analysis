## HCC tissue and blood expr
# expr <- read.table("/BioII/lulab_b/baopengfei/projects/tcga/output/diff/GTExwholeblood-TCGA01A_wilcox.txt",sep = "\t",header = T,check.names = F,fill = T)
# expr <- expr[order(-expr$pvalue,expr$log2FoldChange,decreasing = T),]
# 
# table(abs(expr$log2FoldChange)>=5)  # 1340 true
# table(expr$padj<0.01)
# table(expr$padj<0.01 & abs(expr$log2FoldChange)>=5)
# table(expr$padj == 1)
# 
# 
# hcc.tx <- expr[expr$padj<0.01 & expr$log2FoldChange>=5,]
# 
# blood.tx <- expr[expr$padj<0.01 & !is.na(expr$pvalue),]
# table(is.na(blood.tx$pvalue))
# blood.tx <- blood.tx[order(blood.tx$log2FoldChange,blood.tx$pvalue),]
# blood.tx <- blood.tx[!duplicated(rownames(blood.tx)),]
# table(duplicated(rownames(blood.tx)))
# summary(blood.tx$treatMean)
# summary(blood.tx$ctrlMean)
# 
# top <- rownames(blood.tx)[1:1000]
# botom <- rownames(blood.tx)[(nrow(blood.tx)-999):nrow(blood.tx)]
# #tail( rownames(blood.tx))
# 



# TPM ---------------------------------------------------------------------
gtex <- read.table("/BioII/lulab_b/baopengfei/projects/GTEx/wholeblood_TPM.txt",header = T,stringsAsFactors = F,sep = "\t")
colSums(gtex)
gtex$median <- Biobase::rowMedians(as.matrix(gtex))
gtex[1:3,50:52]
gtex <- gtex[order(gtex$median),]
#?rowMedians

tail <- head(rownames(gtex),1000)
top <- tail(rownames(gtex),1000)
#mid30 <- rownames(gtex)[ceiling(0.7*nrow(gtex)):(ceiling(0.7*nrow(gtex))+999)]
#mid60 <- rownames(gtex)[ceiling(0.4*nrow(gtex)):(ceiling(0.4*nrow(gtex))+999)]
mid30 <- rownames(gtex)[(nrow(gtex)-1999):(nrow(gtex)-1000)]
mid60 <- rownames(gtex)[(nrow(gtex)-2999):(nrow(gtex)-1999)]
plot(gtex[mid60,"median"])



# tss <- read.table("./ref/gtf/promoter2000TSS2000.bed",header = F,stringsAsFactors = F)
# tss$ensg <- unlist(lapply(strsplit(tss$V4,"\\."),function(x) x[1]))
# tss$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_","",tss$ensg)
# tss[1:3,]
# tss.top <- tss[tss$ensg %in% top,]
# tss.top <- tss.top[!duplicated(tss.top$ensg),]
# table(duplicated(tss.top$ensg))
# tss.botom <- tss[tss$ensg %in% botom,]
# tss.botom <- tss.botom[!duplicated(tss.botom$ensg),]
# tss.top <- tss.top[1:4]
# tss.botom <- tss.botom[1:4]
# tss.botom[1:3,]
# write.table(tss.top,"./ref/gtf/promoter2000TSS2000_bloodTxTop1000.bed",quote = F,sep = "\t",col.names = F,row.names = F)
# write.table(tss.botom,"./ref/gtf/promoter2000TSS2000_bloodTxTail1000.bed",quote = F,sep = "\t",col.names = F,row.names = F)
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
{
tss <- read.table("./ref/gtf/promoter2000exon1end2000.bed",header = F,stringsAsFactors = F)
tss$ensg <- unlist(lapply(strsplit(tss$V4,"\\."),function(x) x[1]))
tss$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_","",tss$ensg)
tss[1:3,]
tss.top <- tss[tss$ensg %in% top,]
tss.top <- tss.top[!duplicated(tss.top$ensg),]
table(duplicated(tss.top$ensg))
tss.tail <- tss[tss$ensg %in% tail,]
tss.tail <- tss.tail[!duplicated(tss.tail$ensg),]
tss.top <- tss.top[1:4]
tss.tail <- tss.tail[1:4]
tss.top[1:3,]
tss.tail[1:3,]
write.table(tss.top,"./ref/gtf/promoter2000exon1end2000_bloodTxTop1000.bed",quote = F,sep = "\t",col.names = F,row.names = F)
write.table(tss.tail,"./ref/gtf/promoter2000exon1end2000_bloodTxTail1000.bed",quote = F,sep = "\t",col.names = F,row.names = F)
}

{
  tss <- read.table("./ref/gtf/promoter2000TSS2000.bed",header = F,stringsAsFactors = F)
  tss$ensg <- unlist(lapply(strsplit(tss$V4,"\\."),function(x) x[1]))
  tss$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_","",tss$ensg)
  tss[1:3,]
  tss.top <- tss[tss$ensg %in% top,]
  tss.top <- tss.top[!duplicated(tss.top$ensg),]
  table(duplicated(tss.top$ensg))
  tss.tail <- tss[tss$ensg %in% tail,]
  tss.tail <- tss.tail[!duplicated(tss.tail$ensg),]
  tss.top <- tss.top[1:4]
  tss.tail <- tss.tail[1:4]
  tss.top[1:3,]
  tss.tail[1:3,]
  write.table(tss.top,"./ref/gtf/promoter2000TSS2000_bloodTxTop1000.bed",quote = F,sep = "\t",col.names = F,row.names = F)
  write.table(tss.tail,"./ref/gtf/promoter2000TSS2000_bloodTxTail1000.bed",quote = F,sep = "\t",col.names = F,row.names = F)
}

## add 0.3,0.6 percentile
{
  tss <- read.table("./ref/gtf/promoter2000exon1end2000.bed",header = F,stringsAsFactors = F)
  tss$ensg <- unlist(lapply(strsplit(tss$V4,"\\."),function(x) x[1]))
  tss$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_","",tss$ensg)
  tss[1:3,]
  tss.m3 <- tss[tss$ensg %in% mid30,]
  tss.m3 <- tss.m3[!duplicated(tss.m3$ensg),]
  tss.m6 <- tss[tss$ensg %in% mid60,]
  tss.m6 <- tss.m6[!duplicated(tss.m6$ensg),]
  tss.m3 <- tss.m3[1:4]
  tss.m6 <- tss.m6[1:4]
  tss.m3[1:3,]
  tss.m6[1:3,]
  table(duplicated(tss.m6$V4))
  write.table(tss.m3,"./ref/gtf/promoter2000exon1end2000_bloodTxMid30.bed",quote = F,sep = "\t",col.names = F,row.names = F)
  write.table(tss.m6,"./ref/gtf/promoter2000exon1end2000_bloodTxMid60.bed",quote = F,sep = "\t",col.names = F,row.names = F)
}
{
  tss <- read.table("./ref/gtf/promoter2000TSS2000.bed",header = F,stringsAsFactors = F)
  tss$ensg <- unlist(lapply(strsplit(tss$V4,"\\."),function(x) x[1]))
  tss$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_","",tss$ensg)
  tss[1:3,]
  tss.m3 <- tss[tss$ensg %in% mid30,]
  tss.m3 <- tss.m3[!duplicated(tss.m3$ensg),]
  tss.m6 <- tss[tss$ensg %in% mid60,]
  tss.m6 <- tss.m6[!duplicated(tss.m6$ensg),]
  tss.m3 <- tss.m3[1:4]
  tss.m6 <- tss.m6[1:4]
  tss.m3[1:3,]
  tss.m6[1:3,]
  table(duplicated(tss.m6$V4))
  write.table(tss.m3,"./ref/gtf/promoter2000TSS2000_bloodTxMid30.bed",quote = F,sep = "\t",col.names = F,row.names = F)
  write.table(tss.m6,"./ref/gtf/promoter2000TSS2000_bloodTxMid60.bed",quote = F,sep = "\t",col.names = F,row.names = F)
}



# fpkm (newer) --------------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

ref <- read.table("../../GTEx/wholeblood_Tx_fpkm_2021-NC-NDR.txt",sep = "\t",header = T)

## use cut to get desired categories
ref$FPKMblood.f <- cut(as.numeric(ref$FPKMblood), breaks=c(0,0.01,0.1,5,30,100000), labels=c("unexp","fpkm-001-01", "fpkm-01-5","fpkm-5-30", "fpkm-30"),include.lowest=T)
table(ref$FPKMblood.f)
str(ref)
colnames(ref)[3] <- "name"
ref <- ref[,c("name","FPKMblood.f")]

## read ensg-name ref
ref2 <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene_v2.txt",sep = "\t",header = T,stringsAsFactors = F)
ref2 <- ref2[,c("name","ensg.full")]

ref <- dplyr::left_join(ref,ref2)
#ref <- ref[!duplicated(ref$name),]
table(is.na(ref$ensg.full))
ref <- ref[!is.na(ref$ensg.full),]


## read bed ref (tss) 
tss <- read.table("./ref/gtf/promoter2000TSS2000.bed",header = F,stringsAsFactors = F)
n1 <- colnames(tss)
tss$ensg.full <- gsub("promoter300100exon1end_|promoter150TSS50_","",tss$V4)
#dim(tss)
table(duplicated(tss$ensg.full))
tss <- dplyr::left_join(tss,ref)
table(duplicated(tss$V4))
tss <- tss[!duplicated(tss$V4) & !is.na(tss$name),]
table(tss$FPKMblood.f)
#unexp fpkm-001-01   fpkm-01-5   fpkm-5-30     fpkm-30 
#68          38         125          91          59 
#tss <- tss[,n1]
for (i in unique(tss$FPKMblood.f)){
  write.table(tss[tss$FPKMblood.f==i,n1],paste0("./ref/gtf/bloodNDR/promoter2000tss2000_bloodTx_",i,".bed"),quote = F,sep = "\t",col.names = F,row.names = F)
}


## read bed ref (exon1) 
exon1 <- read.table("./ref/gtf/promoter2000exon1end2000.bed",header = F,stringsAsFactors = F)
n1 <- colnames(exon1)
exon1$ensg.full <- gsub("promoter300100exon1end_|promoter150TSS50_","",exon1$V4)
exon1 <- dplyr::left_join(exon1,ref)
exon1 <- exon1[!duplicated(exon1$V4) & !is.na(exon1$name),]
for (i in unique(exon1$FPKMblood.f)){
  write.table(exon1[exon1$FPKMblood.f==i,n1],paste0("./ref/gtf/bloodNDR/promoter2000exon1end2000_bloodTx_",i,".bed"),quote = F,sep = "\t",col.names = F,row.names = F)
}


## get selected candidates for coverage depth plot 
tss <- read.table("./ref/gtf/promoter2000TSS2000.bed",header = F,stringsAsFactors = F)
tss.ensg.full <- gsub("promoter300100exon1end_|promoter150TSS50_","",tss$V4)
exon1 <- read.table("./ref/gtf/promoter2000exon1end2000.bed",header = F,stringsAsFactors = F)
exon1.ensg.full <- gsub("promoter300100exon1end_|promoter150TSS50_","",exon1$V4)

#ENSG00000048140.17
#ENSG00000040487.12
#ENSG00000047936.10
candidates <- c("ENSG00000048140.17","ENSG00000040487.12","ENSG00000047936.10")

dir.create("ref/gtf/candidatesNDR")
write.table(tss[tss.ensg.full %in% candidates,],paste0("./ref/gtf/candidatesNDR/promoter2000tss2000_bloodTx_RFLOO-1202.bed"),quote = F,sep = "\t",col.names = F,row.names = F)
write.table(exon1[exon1.ensg.full %in% candidates,],paste0("./ref/gtf/candidatesNDR/promoter2000exon1end2000_bloodTx_RFLOO-1202.bed"),quote = F,sep = "\t",col.names = F,row.names = F)
