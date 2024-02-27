# inhouse CNV

setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
options(stringsAsFactors = F)

cnv <- read.table("output/lulab/wisecondorx/CNV/sum-score.txt",check.names = F,header = F)
id <- read.table("metadata/lulab/sample_table.txt",check.names = F,header = T, sep = "\t")

cnv <- merge(cnv,id,by.x = "V1", by.y = "data_id")
cnv <- cnv[cnv$overall_QC=="Y" & cnv$group!="STAD",]
cnv <- cnv[,1:3]
colnames(cnv)[1:2] <- c("sample","CNV.score")
cnv$group <- factor(cnv$group,levels = c("CRC","NC"))
mytheme <- theme(plot.title = element_text(size = 16,color="black",hjust = 0.5,face="bold"),
                 axis.title = element_text(size = 16,color ="black",face="bold"), 
                 axis.text = element_text(size= 16,color = "black",face="bold"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 #axis.text.x = element_text( hjust = 1 ), # angle = 45,
                 panel.grid=element_blank(),
                 #legend.position = "top",
                 legend.text = element_text(size= 16),
                 legend.title= element_text(size= 16))


b <- runif(nrow(cnv), -0.15, 0.15)

ggplot(cnv, aes(x=group, y=CNV.score)) +  # , size=cyl
  geom_boxplot(outlier.colour = "white", outlier.size = 0) + 
  geom_point( aes(x = as.numeric(group) + b, y = CNV.score, fill=group), size=6, shape=21, stroke = 1 , color="black" ) +
  scale_fill_manual(values=alpha(c("orange2","steelblue2"), 1))  + 
  stat_compare_means(label.x.npc = 0.5,  bracket.size = .3,tip.length=1,size=10,paired = FALSE,label = "p.format",method = "wilcox.test",method.args = list(alternative = "two.sided")) +
  theme_classic() + mytheme +
  ylab("global CNV score") 

str(cnv)
cnv$CNV.score



#### plot NDR boxplot ####
# {PPP1R16A,RAB25,PRTN3,LSR}
# {PRTN3,RAB25,}
ndr.score <- read.table("output/lulab/NDR/PRTN3.score.summary",check.names = F,header = F, sep = "\t")
colnames(ndr.score) <- c("sample","score")
ndr.score <- merge(ndr.score,id,by.x = "sample", by.y = "data_id")
ndr.score <- ndr.score[ndr.score$overall_QC=="Y" & ndr.score$group!="STAD",]
#ndr.score <- ndr.score[ndr.score$overall_QC=="Y",]
#ndr.score <- ndr.score[ndr.score$overall_QC=="Y" & ndr.score$uniq_reads>=quantile(ndr.score$uniq_reads)[2] & ndr.score$coverage_ratio>=quantile(ndr.score$coverage_ratio)[2] & ndr.score$coverage_depth>=quantile(ndr.score$coverage_depth)[2],]


ndr.score <- ndr.score[,1:3]
ndr.score$score[ndr.score$score=="N.A."] <- 0
ndr.score$score <- as.numeric(ndr.score$score)
str(ndr.score)
#ndr.score$group <- factor(ndr.score$group,levels = c("CRC","STAD","NC"))
ndr.score$group <- factor(ndr.score$group,levels = c("CRC","NC"))
#large.id <- read.table("metadata/lulab/sample_ids_large.txt",check.names = F,header = F, sep = "\t")
#ndr.score <- ndr.score[ndr.score$sample %in% large.id$V1,]
#str(ndr.score)

b <- runif(nrow(ndr.score), -0.15, 0.15)
ggplot(ndr.score, aes(x=group, y=score)) +  # , size=cyl
  geom_boxplot(outlier.colour = "white", outlier.size = 0) + 
  geom_point( aes(x = as.numeric(group) + b, y = score, fill=group), size=6, shape=21, stroke = 1 , color="black") +
  scale_fill_manual(values=alpha(c("orange2","steelblue2"), 1))  +   # c("orange2","firebrick2","steelblue2")
  stat_compare_means(label.x.npc = 0.5, bracket.size = .1, tip.length=1,size=10,
                     paired = FALSE, label = "p.format",
                     #comparisons = list(c("CRC","NC") ) ,  # c("STAD","NC")
                     method = "wilcox.test",method.args = list(alternative = "two.sided")) +  # 
  theme_classic() + mytheme +
  ylab("PRTN3 promoter CNV score") 

# colnames(ndr.score)
# colnames(ndr.score)[2] <- "single.NDR.score"
# colnames(cnv)[2] <- "global.CNV.score"
# 
# t <- merge(cnv,ndr.score,by.x="sample",by.y="sample")
# t <- t[,c(3,2,4)]
# colnames(t)[1] <- c("group")
# write.table(t,"./wgs-roc.txt",sep = "\t",quote = F,col.names = T,row.names = F)
# 
# log2tpm[1:4,1:4]
# tpm <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/dmr/tpm-matrix-promoter-filter-ensg.txt",stringsAsFactors = F, sep = "\t",header = T,row.names = 1,check.names = F)
# log2tpm <- log2(tpm+1)
# tab <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/meta/lulab/sample_table.txt",stringsAsFactors = F, sep = "\t",header = T,row.names = 1,check.names = F)
# tab <- tab[tab$group != "STAD" & tab$overall_QC=="Y",]
# log2tpm <- log2tpm[,match(rownames(tab),colnames(log2tpm))]
# septin9 <- log2tpm[rownames(log2tpm)=="ENSG00000184640",]
# write.table(septin9,"/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/dmr/septin9.txt",sep = "\t",quote = F,row.names = T)
# 
# crc.spec <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/dmr/crc_specific_ensg.txt",stringsAsFactors = F, sep = "\t",check.names = F)
# log2tpm <- log2tpm[match(crc.spec$V1,rownames(log2tpm)),]
# write.table(log2tpm,"/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/dmr/crc.txt",sep = "\t",quote = F,row.names = T)
# write.table(tab$group,"/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/dmr/crc-label.txt",sep = "\t",quote = F,row.names = F,col.names = F)
# ENSG00000184640.17



#### plot NDR coverage ####
cov <- read.table("output/lulab/NDR/PRTN3.relativeDepth.summary",header = T, sep = " ",stringsAsFactors = F,check.names = F)
t <- cov
#t <- read.table("output/lulab/fragment-length/histogram.txt",header = T, sep = "\t",stringsAsFactors = F,check.names = F)  # fill = T,
#colnames(t)
f <- read.table("./metadata/lulab/sample_table.txt",stringsAsFactors = F,header = T, sep = "\t")
colnames(f)
summary(f$coverage_depth)
f <- f[f$overall_QC=="Y" & f$group!="STAD" ,]
table(f$group)
t <- t[,c("TSS",f$data_id)]
all(colnames(t)[2:ncol(t)]==f$data_id)

t$CRC.aver <- apply(t[,grep("CRC",colnames(t))],1,mean) 
t$NC.aver <- apply(t[,grep("NC",colnames(t))],1,mean) 
library(reshape2)
m <- melt(t,id.vars = 1,variable.name = "id")
m <- m[-grep("NC-PKU",m$id),]

library(ggplot2)
m1 <- m[m$id!="CRC.aver" & m$id!="NC.aver",]
large.id <- read.table("metadata/lulab/sample_ids_large.txt",check.names = F,header = F, sep = "\t")
m1 <- m1[m1$id %in% large.id$V1,]
m2 <- m[ m$id=="NC.aver",]
m3 <- m[ m$id=="CRC.aver",]

ggplot() +  # ,fill=id,alpha=0,color=col
  scale_color_manual(values=alpha(c("orange2","steelblue2"), 1))  +  # "1" = "red", "2" = "grey", "3" = "blue"
  geom_line(stat = "summary",data = m1, aes(x=TSS,y=value,group=id),color=alpha("orange2",0.3),size = 0.1) + # 
  geom_line(stat = "summary",data = m2, aes(x=TSS,y=value),color="steelblue2",size = 1) + # 
  geom_line(stat = "summary",data = m3, aes(x=TSS,y=value),color="orange2",size = 1) + # 
  
  geom_vline(xintercept = c(-150,50),color="grey", linetype="dashed", size=0.5) +
  xlim(c(-1000,1000)) +
  ylab("Relative Coverage") +
  xlab("TSS Relative Position (bp)") +
  theme_classic() + mytheme






#### plot RNA expr in tissue and PBMC ####
tissue.expr <- tpm.tissue.m
table(tissue.expr$hgnc_symbol=="PRTN3")
tissue.expr <- tissue.expr[tissue.expr$hgnc_symbol=="PRTN3",]
rownames(tissue.expr) <- tissue.expr$hgnc_symbol
tissue.expr <- tissue.expr[,-1] 

tissue.expr <- as.data.frame(t(tissue.expr))
tissue.expr$group <- "PBMC"
tissue.expr$group[grep("-N",rownames(tissue.expr))] <- "Colon_Normal_Tissue"
tissue.expr$group[grep("-T",rownames(tissue.expr))] <- "CRC_Tissue"

tissue.expr$group2 <- "PBMC"
tissue.expr$group2[grep("-N",rownames(tissue.expr))] <- "Colon_Tissue"
tissue.expr$group2[grep("-T",rownames(tissue.expr))] <- "Colon_Tissue"
tissue.expr$PRTN3 <- log2(tissue.expr$PRTN3 + 1)
colnames(tissue.expr) <- c("log2TPM","group","tissue")

tissue.expr$tissue <- factor(tissue.expr$tissue,levels = c("Colon_Tissue","PBMC"))

mytheme <- theme(plot.title = element_text(size = 16,color="black",hjust = 0.5,face="bold"),
                 axis.title = element_text(size = 16,color ="black",face="bold"), 
                 axis.text = element_text(size= 16,color = "black",face="bold"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 #axis.text.x = element_text( hjust = 1 ), # angle = 45,
                 panel.grid=element_blank(),
                 #legend.position = "top",
                 legend.text = element_text(size= 16),
                 legend.title= element_text(size= 16))
b <- runif(nrow(tissue.expr), -0.15, 0.15)
ggplot(tissue.expr, aes(x=tissue, y=log2TPM)) +  # , size=cyl
  geom_boxplot(outlier.colour = "black", outlier.size = 0) + 
  geom_point( aes(x = as.numeric(tissue) + b, y = log2TPM, fill=tissue), size=6, shape=21, stroke = 1 , color="black") +
  scale_fill_manual(values=alpha(c("orange2","steelblue2"), 1))  +   # c("orange2","firebrick2","steelblue2")
  stat_compare_means(label.x.npc = 0.3, bracket.size = .1, tip.length=1,size=10,
                     paired = FALSE, label = "p.format",
                     #comparisons = list(c("CRC","NC") ) ,  # c("STAD","NC")
                     method = "wilcox.test",method.args = list(alternative = "two.sided")) +  # 
  theme_classic() + mytheme +
  ylab("PRTN3 Expression (log2TPM)") 
