---
title: "5' cfDNA end motif"
author: "bpf"
date: "12/8/2021"
output: html_document
---


# lulab DNA ---------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")

# read func.
a <- function(x){
  # read data and set column names
  #x <- "./output/lulab/motif5/CRC-PKU-10-wgs_motif5.txt" 
  data <- read.table(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  
  data$motif5 <- as.character(lapply(strsplit(data$V1,":", ), function(x) x[1]))
  data$sample <- as.character(lapply(strsplit(data$V1,":", ), function(x) x[2]))

  l <- grep("N",data$motif5,invert = T)
  data <- data[l,c("motif5","sample")]
  sample <- as.character(lapply(strsplit(basename(x),"_" ), function(x) x[1]))
  colnames(data)[2] <- sample
  
  # return the data
  return(data)
}

# get mds by sample
get.mds <- function(x){
  #x <- motif.ratio$`NC-PKU-10-wgs`
  mds <- x
  for(i in 1:length(x)){
    mds[i] <- -x[i]*log(x[i])/log(256)
  }
  return(sum(mds))
}

cmp.boxplot <- function(matrix,point.size=1,point.col="black",point.alpha=0.6,title="",y.lim,x.title="group",y.title="value"){
  library(ggplot2)
  b <- runif(nrow(matrix), -0.2, 0.2)
  #p <- compare
  ggplot(matrix, aes(x = group, y = value))+ 
    labs(y=y.title,x=x.title,title = title)+  
    geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
                 outlier.alpha = 1, outlier.size = 0)+
    geom_point(aes(x=as.numeric(group)+b,y=value),size=point.size,color=point.col, stroke=0.5,alpha=point.alpha) +
    #geom_violin(draw_quantiles = T) +
    geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
    ylim(y.lim) +
    theme_bw() + 
    theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
          axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
          axis.text = element_text(size= 18,color = "black",family="arial"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1 ),
          panel.grid=element_blank(),
          legend.position = "top",
          legend.text = element_text(size= 12),
          legend.title= element_text(size= 12))
}



## sum motif
files <- Sys.glob("./output/lulab/motif5/*_motif5.txt")  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed


# run the anonymous function defined above
motif <- lapply(files, a)
rownames(motif[[1]])
motif <- do.call("cbind", motif)
motif <- motif[,!duplicated(colnames(motif))]
dir.create("./output/lulab/motif5/sum")
write.table(motif,"./output/lulab/motif5/sum/fre.txt",sep = "\t",quote = F,row.names = F,col.names = T)

## get ratio
motif <- read.table("./output/lulab/motif5/sum/fre.txt",sep = "\t",row.names = 1, check.names = F, header = T)
s <- colSums(motif)

motif.ratio <- motif
for (i in 1:nrow(motif)){
  motif.ratio[i,] <- motif[i,]/s
}
# motif.ratio <- motif/s
# head(s)
# motif[1:3,1:3]
write.table(motif.ratio,"./output/lulab/motif5/sum/ratio.txt",sep = "\t",quote = F,row.names = T,col.names = T)

## get motif diversity score (MDS)
#motif.mds <- motif.ratio

motif.mds <- apply(motif.ratio,2,get.mds)
motif.mds.df <- as.data.frame(motif.mds)
motif.mds.df$group <- as.character(lapply(strsplit(rownames(motif.mds.df),"-"),function(x) x[1]))
colnames(motif.mds.df)[1] <- "value"

str(motif.mds.df)

write.table(t(motif.mds),"./output/lulab/motif5/sum/motif_diversity_score.txt",sep = "\t",quote = F,row.names = F,col.names = T)



# filter sample ID
crc.id <- read.delim("../../multi-omics-explore/meta/lulab/paired-wgs-CRC-passQCloose.txt",header = F,stringsAsFactors = F)$V1
stad.id <- read.delim("../../multi-omics-explore/meta/lulab/paired-wgs-STAD-passQCloose.txt",header = F,stringsAsFactors = F)$V1
nc.id <- read.delim("../../multi-omics-explore/meta/lulab/paired-wgs-NC-passQCloose.txt",header = F,stringsAsFactors = F)$V1


# MDS box plot
# motif.mds.df$group <- factor(motif.mds.df$group,levels = c("NC","CRC","STAD"))
# matrix <- motif.mds.df
# point.size <- 1
# point.col <- "black"
# point.alpha <- 0.6
# title <- ""
# y.lim <- c(0.93,0.96)
# x.title <- "group"
# y.title <- "value"

library(ggplot2)
library(ggpubr)
library(ggsci)
# #matrix <- read.table("./output/lulab/motif5/sum/motif_diversity_score.txt",sep = "\t",check.names = F,header = T)
# motif.mds$sample <- rownames(motif.mds)
# motif.mds$group <- as.character(lapply(strsplit(as.character(motif.mds$sample),"\\-"), function(x) x[1]))


matrix <- as.data.frame(motif.mds.df)
matrix <- matrix[rownames(matrix) %in% c(nc.id,crc.id,stad.id),]
matrix$group <- factor(matrix$group,levels = c("NC","CRC","STAD"))
b <- runif(nrow(matrix), -0.2, 0.2)
#p <- compare
ggplot(matrix, aes(x = group, y = value))+  # , fill=group
  labs(y="Motif diversity score")+   # x=x.title,title = title
  geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 0, outlier.size = 0)+
  geom_point(aes(x=as.numeric(group)+b,y=value),size=2,color="black", stroke=0.5,alpha=0.6) +
  #geom_violin(draw_quantiles = T) +
  #geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
  #ylim(c(format(min(matrix$value),digit=2),format(max(matrix$value),digit=2))) +
  #ggsci::scale_fill_nejm()+
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(), # angle = 45, hjust = 1 
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) +
        ggpubr::stat_compare_means(ref.group = "NC",
                                method = "wilcox.test",
                                method.args = list(alternative = "greater"),
                                label = "p.format",size = 8, vjust=0.6 # hide.ns=T,
                                )



# single motif box plot
motif <- read.table("./output/lulab/motif5/sum/ratio.txt",sep = "\t",row.names = 1, check.names = F, header = T)
motif <- motif*100
motif <- motif[,colnames(motif) %in% c(nc.id,crc.id,stad.id)]

# CRC
crc <- read.table("./output/lulab/motif5/sum/fre.diff.STAD.NC.txt",sep = "\t",row.names = 1, check.names = F, header = T)
crc.motif <- rownames(crc)[1:6]
crc.motif <- motif[crc.motif,]
crc.motif$motif <- rownames(crc.motif)
crc.motif <- reshape2::melt(crc.motif,id.var="motif")
crc.motif$group <- as.character(lapply(strsplit(as.character(crc.motif$variable),"\\-"), function(x) x[1]))
#str(crc.motif)
crc.motif$group <- factor(crc.motif$group,levels = c("NC","CRC","STAD"))

wilcox.test(crc.motif[crc.motif$group=="STAD" & crc.motif$motif=="TCCA","value"],crc.motif[crc.motif$group=="NC" & crc.motif$motif=="TCCA","value"],alternative='two.sided')
table(crc.motif$group)

matrix <- crc.motif[crc.motif$group!="CRC",]
library(ggplot2)
library(ggpubr)
library(ggsci)
b <- runif(nrow(matrix), -0.2, 0.2)
#p <- compare
ggplot(matrix, aes(x = motif, y = value, fill=group))+  # , fill=group
  labs(y="Motif frequency(%)")+   # x=x.title,title = title
  geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 0, outlier.size = 0)+
  #geom_point(aes(x=as.numeric(group)+b,y=value),size=2,color="black", stroke=0.5,alpha=0.6) +
  #geom_violin(draw_quantiles = T) +
  #geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
  #ylim(c(format(min(matrix$value),digit=2),format(max(matrix$value),digit=2))) +
  scale_fill_manual(values = c("steelblue4","firebrick2"))+ # "grey50",
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1 ,vjust = 0.5), # 
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) +
        ggpubr::stat_compare_means(#comparisons = list(c("CRC","NC")),# ref.group = "NC",
                                method = "wilcox.test",
                                #method.args = list(alternative = "greater"),
                                label = "p.signif",size = 8, vjust=0.6, hide.ns=F,
                                )




# lulab medip ---------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq")


# read func.
a <- function(x){
  # read data and set column names
  #x <- "./output/lulab/motif5/CRC-PKU-10-wgs_motif5.txt"
  data <- read.table(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  
  data$motif5 <- as.character(lapply(strsplit(data$V1,":", ), function(x) x[1]))
  data$sample <- as.character(lapply(strsplit(data$V1,":", ), function(x) x[2]))
  
  l <- grep("N",data$motif5,invert = T)
  data <- data[l,c("motif5","sample")]
  sample <- as.character(lapply(strsplit(basename(x),"_" ), function(x) x[1]))
  colnames(data)[2] <- sample
  
  # return the data
  return(data)
}

# get mds by sample
get.mds <- function(x){
  #x <- motif.ratio$`NC-PKU-10-wgs`
  mds <- x
  for(i in 1:length(x)){
    mds[i] <- -x[i]*log(x[i])/log(256)
  }
  return(sum(mds))
}

cmp.boxplot <- function(matrix,point.size=1,point.col="black",point.alpha=0.6,title="",y.lim,x.title="group",y.title="value"){
  library(ggplot2)
  b <- runif(nrow(matrix), -0.2, 0.2)
  #p <- compare
  ggplot(matrix, aes(x = group, y = value))+ 
    labs(y=y.title,x=x.title,title = title)+  
    geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
                 outlier.alpha = 1, outlier.size = 0)+
    geom_point(aes(x=as.numeric(group)+b,y=value),size=point.size,color=point.col, stroke=0.5,alpha=point.alpha) +
    #geom_violin(draw_quantiles = T) +
    geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
    ylim(y.lim) +
    theme_bw() + 
    theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
          axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
          axis.text = element_text(size= 18,color = "black",family="arial"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1 ),
          panel.grid=element_blank(),
          legend.position = "top",
          legend.text = element_text(size= 12),
          legend.title= element_text(size= 12))
}



## sum motif
files <- Sys.glob("./output/lulab/motif5/*_motif5.txt")  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed


# run the anonymous function defined above
motif <- lapply(files, a)
rownames(motif[[1]])
motif <- do.call("cbind", motif)
motif <- motif[,!duplicated(colnames(motif))]
dir.create("./output/lulab/motif5/sum")
write.table(motif,"./output/lulab/motif5/sum/fre.txt",sep = "\t",quote = F,row.names = F,col.names = T)


## get ratio
motif <- read.table("./output/lulab/motif5/sum/fre.txt",sep = "\t",row.names = 1, check.names = F, header = T)
s <- colSums(motif)

motif.ratio <- motif
for (i in 1:nrow(motif)){
  motif.ratio[i,] <- motif[i,]/s
}
# motif.ratio <- motif/s
# head(s)
# motif[1:3,1:3]
colSums(motif.ratio)
write.table(motif.ratio,"./output/lulab/motif5/sum/ratio.txt",sep = "\t",quote = F,row.names = T,col.names = T)


## get motif diversity score (MDS)
#motif.mds <- motif.ratio
motif.mds <- apply(motif.ratio,2,get.mds)
motif.mds.df <- as.data.frame(motif.mds)
motif.mds.df$group <- as.character(lapply(strsplit(rownames(motif.mds.df),"-"),function(x) x[1]))
colnames(motif.mds.df)[1] <- "value"

str(motif.mds.df)
write.table(t(motif.mds),"./output/lulab/motif5/sum/motif_diversity_score.txt",sep = "\t",quote = F,row.names = F,col.names = T)



# filter sample ID
crc.id <- read.delim("../../multi-omics-explore/meta/lulab/paired-medip-CRC-passQCloose.txt",header = F,stringsAsFactors = F)$V1
stad.id <- read.delim("../../multi-omics-explore/meta/lulab/paired-medip-STAD-passQCloose.txt",header = F,stringsAsFactors = F)$V1
nc.id <- read.delim("../../multi-omics-explore/meta/lulab/paired-medip-NC-passQCloose.txt",header = F,stringsAsFactors = F)$V1


# MDS box plot
# motif.mds.df$group <- factor(motif.mds.df$group,levels = c("NC","CRC","STAD"))
# matrix <- motif.mds.df
# point.size <- 1
# point.col <- "black"
# point.alpha <- 0.6
# title <- ""
# y.lim <- c(0.93,0.96)
# x.title <- "group"
# y.title <- "value"

library(ggplot2)
library(ggpubr)
library(ggsci)
# #matrix <- read.table("./output/lulab/motif5/sum/motif_diversity_score.txt",sep = "\t",check.names = F,header = T)
# motif.mds$sample <- rownames(motif.mds)
# motif.mds$group <- as.character(lapply(strsplit(as.character(motif.mds$sample),"\\-"), function(x) x[1]))


matrix <- as.data.frame(motif.mds.df)
matrix <- matrix[rownames(matrix) %in% c(nc.id,crc.id,stad.id),]
matrix$group <- factor(matrix$group,levels = c("NC","CRC","STAD"))
matrix$label <- rownames(matrix)
b <- runif(nrow(matrix), -0.2, 0.2)
#p <- compare
ggplot(matrix, aes(x = group, y = value))+  # , fill=group
  labs(y="Motif diversity score")+   # x=x.title,title = title
  geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 0, outlier.size = 0)+
  geom_point(aes(x=as.numeric(group)+b,y=value),size=2,color="black", stroke=0.5,alpha=0.6) +
  ggrepel::geom_text_repel(aes(x=as.numeric(group)+b,y=value,label=label),#nudge_x=.1,nudge_y=.1,
                           #box.padding = 0.5,box.padding = .1,point.padding = .1,
                           size = 3,color="black",max.overlaps = Inf ) + 
  #geom_violin(draw_quantiles = T) +
  #geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
  #ylim(c(format(min(matrix$value),digit=2),format(max(matrix$value),digit=2))) +
  #ggsci::scale_fill_nejm()+
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(), # angle = 45, hjust = 1 
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) +
  ggpubr::stat_compare_means(ref.group = "NC",
                             method = "wilcox.test",
                             method.args = list(alternative = "greater"),
                             label = "p.format",size = 8, vjust=0.6 # hide.ns=T,
  )



# single motif box plot
motif <- read.table("./output/lulab/motif5/sum/ratio.txt",sep = "\t",row.names = 1, check.names = F, header = T)
motif <- motif*100
motif <- motif[,colnames(motif) %in% c(nc.id,crc.id,stad.id)]

# CRC (no diff p<0.05)
# STAD
crc <- read.table("./output/lulab/motif5/sum/fre.diff.STAD.NC.txt",sep = "\t",row.names = 1, check.names = F, header = T)
crc.motif <- rownames(crc)[1:6]
crc.motif <- motif[crc.motif,]
crc.motif$motif <- rownames(crc.motif)
crc.motif <- reshape2::melt(crc.motif,id.var="motif")
crc.motif$group <- as.character(lapply(strsplit(as.character(crc.motif$variable),"\\-"), function(x) x[1]))
#str(crc.motif)
crc.motif$group <- factor(crc.motif$group,levels = c("NC","CRC","STAD"))

#wilcox.test(crc.motif[crc.motif$group=="CRC" & crc.motif$motif=="TCCA","value"],crc.motif[crc.motif$group=="NC" & crc.motif$motif=="TCCA","value"],alternative='two.sided')
table(crc.motif$group)

matrix <- crc.motif[crc.motif$group!="CRC",]
library(ggplot2)
library(ggpubr)
library(ggsci)
b <- runif(nrow(matrix), -0.2, 0.2)
#p <- compare
ggplot(matrix, aes(x = motif, y = value, fill=group))+  # , fill=group
  labs(y="Motif frequency(%)")+   # x=x.title,title = title
  geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 0, outlier.size = 0)+
  #geom_point(aes(x=as.numeric(group)+b,y=value),size=2,color="black", stroke=0.5,alpha=0.6) +
  #geom_violin(draw_quantiles = T) +
  #geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
  #ylim(c(format(min(matrix$value),digit=2),format(max(matrix$value),digit=2))) +
  #scale_fill_manual(values = c("steelblue4","firebrick2"))+ # "grey50",
  ggsci::scale_fill_nejm()+
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1 ,vjust = 0.5), # 
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) +
  ggpubr::stat_compare_means(#comparisons = list(c("CRC","NC")),# ref.group = "NC",
    method = "wilcox.test",
    #method.args = list(alternative = "greater"),
    label = "p.signif",size = 8, vjust=0.6, hide.ns=F,
  )





# public 5hmc (to be continued) ---------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq")


# read func.
a <- function(x){
  # read data and set column names
  #x <- paste0("./output/",dst,"/motif5/SRR6947629_motif5.txt")
  # x <- "./output/GSE112679/motif5/SRR6946743_motif5.txt"
  #print(x)
  data <- read.table(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  
  data$motif5 <- as.character(lapply(strsplit(data$V1,":", ), function(x) x[1]))
  data$sample <- as.character(lapply(strsplit(data$V1,":", ), function(x) x[2]))
  
  l <- grep("N",data$motif5,invert = T)
  data <- data[l,c("motif5","sample")]
  sample <- as.character(lapply(strsplit(basename(x),"_" ), function(x) x[1]))
  colnames(data)[2] <- sample
  
  # return the data
  return(data)
}

# get mds by sample
get.mds <- function(x){
  #x <- motif.ratio$`NC-PKU-10-wgs`
  #print(x)
  mds <- x
  for(i in 1:length(x)){
    mds[i] <- -x[i]*log(x[i])/log(256)
  }
  return(sum(mds))
}

cmp.boxplot <- function(matrix,point.size=1,point.col="black",point.alpha=0.6,title="",y.lim,x.title="group",y.title="value"){
  library(ggplot2)
  b <- runif(nrow(matrix), -0.2, 0.2)
  #p <- compare
  ggplot(matrix, aes(x = group, y = value))+ 
    labs(y=y.title,x=x.title,title = title)+  
    geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
                 outlier.alpha = 1, outlier.size = 0)+
    geom_point(aes(x=as.numeric(group)+b,y=value),size=point.size,color=point.col, stroke=0.5,alpha=point.alpha) +
    #geom_violin(draw_quantiles = T) +
    geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
    ylim(y.lim) +
    theme_bw() + 
    theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
          axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
          axis.text = element_text(size= 18,color = "black",family="arial"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1 ),
          panel.grid=element_blank(),
          legend.position = "top",
          legend.text = element_text(size= 12),
          legend.title= element_text(size= 12))
}


dst <- "GSE112679" # "GSE81314" "GSE89570"  "GSE112679"

## sum motif
files <- Sys.glob(paste0("./output/",dst,"/motif5/*_motif5.txt"))  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed


# run the anonymous function defined above
motif <- lapply(files, a)
motif <- do.call("cbind", motif)
motif <- motif[,!duplicated(colnames(motif))]
motif[1:3,1:3]
dir.create(paste0("./output/",dst,"/motif5/sum"))
#motif$motif <- rownames(motif)
#motif <- motif[,c(ncol(motif),1:(ncol(motif)-1))]
write.table(motif,paste0("./output/",dst,"/motif5/sum/fre.txt"),sep = "\t",quote = F,row.names = F,col.names = T)


## get ratio
motif <- read.table(paste0("./output/",dst,"/motif5/sum/fre.txt"),sep = "\t",row.names = 1, check.names = F, header = T)
s <- colSums(motif)

motif.ratio <- motif
for (i in 1:nrow(motif)){
  motif.ratio[i,] <- motif[i,]/s
}
# motif.ratio <- motif/s
# head(s)
# motif[1:3,1:3]
colSums(motif.ratio)  # should be all 1
motif.ratio.out <- motif.ratio
motif.ratio.out$motif <- gsub("\"","",fixed = T,rownames(motif.ratio.out))
motif.ratio.out <- motif.ratio.out[,c(ncol(motif.ratio.out),1:(ncol(motif.ratio.out)-1))]
write.table(motif.ratio.out,paste0("./output/",dst,"/motif5/sum/ratio.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
motif.ratio.out[1:3,1:3]

## get motif diversity score (MDS)
#motif.mds <- motif.ratio
motif.mds <- apply(motif.ratio,2,get.mds)
motif.mds <- t(motif.mds)
motif.mds <- as.data.frame(motif.mds)
motif.mds.out <- motif.mds
motif.mds.out$motif <- "motif.diversity.score"
motif.mds.out <- motif.mds.out[,c(ncol(motif.mds.out),1:(ncol(motif.mds.out)-1))]
motif.mds.out[,1:3]
write.table(motif.mds.out,paste0("./output/",dst,"/motif5/sum/motif_diversity_score.txt"),sep = "\t",quote = F,row.names = F,col.names = T)







# filter sample ID (to be continued)
crc.id <- read.delim("/Share2/home/lulab1/cfDNA-seq/raw/DIP-seq/GSE81314/meta_data/5hmc_sample_table.txt",header = F,stringsAsFactors = F)$V1


# MDS box plot
library(ggplot2)
library(ggpubr)
library(ggsci)
motif.mds.df <- as.data.frame(motif.mds)
motif.mds.df$group <- as.character(lapply(strsplit(rownames(motif.mds.df),"-"),function(x) x[1]))
colnames(motif.mds.df)[1] <- "value"
str(motif.mds.df)
matrix <- as.data.frame(motif.mds.df)
matrix <- matrix[rownames(matrix) %in% c(nc.id,crc.id,stad.id),]
matrix$group <- factor(matrix$group,levels = c("NC","CRC","STAD"))
matrix$label <- rownames(matrix)
b <- runif(nrow(matrix), -0.2, 0.2)
#p <- compare
ggplot(matrix, aes(x = group, y = value))+  # , fill=group
  labs(y="Motif diversity score")+   # x=x.title,title = title
  geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 0, outlier.size = 0)+
  geom_point(aes(x=as.numeric(group)+b,y=value),size=2,color="black", stroke=0.5,alpha=0.6) +
  ggrepel::geom_text_repel(aes(x=as.numeric(group)+b,y=value,label=label),#nudge_x=.1,nudge_y=.1,
                           #box.padding = 0.5,box.padding = .1,point.padding = .1,
                           size = 3,color="black",max.overlaps = Inf ) + 
  #geom_violin(draw_quantiles = T) +
  #geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
  #ylim(c(format(min(matrix$value),digit=2),format(max(matrix$value),digit=2))) +
  #ggsci::scale_fill_nejm()+
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(), # angle = 45, hjust = 1 
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) +
  ggpubr::stat_compare_means(ref.group = "NC",
                             method = "wilcox.test",
                             method.args = list(alternative = "greater"),
                             label = "p.format",size = 8, vjust=0.6 # hide.ns=T,
  )



# single motif box plot
motif <- read.table(paste0("./output/",dst,"/motif5/sum/ratio.txt"),sep = "\t",row.names = 1, check.names = F, header = T)
motif <- motif*100
motif <- motif[,colnames(motif) %in% c(nc.id,crc.id,stad.id)]

# CRC (no diff p<0.05)
# STAD
crc <- read.table(paste0("./output/",dst,"/motif5/sum/fre.diff.STAD.NC.txt"),sep = "\t",row.names = 1, check.names = F, header = T)
crc.motif <- rownames(crc)[1:6]
crc.motif <- motif[crc.motif,]
crc.motif$motif <- rownames(crc.motif)
crc.motif <- reshape2::melt(crc.motif,id.var="motif")
crc.motif$group <- as.character(lapply(strsplit(as.character(crc.motif$variable),"\\-"), function(x) x[1]))
#str(crc.motif)
crc.motif$group <- factor(crc.motif$group,levels = c("NC","CRC","STAD"))

#wilcox.test(crc.motif[crc.motif$group=="CRC" & crc.motif$motif=="TCCA","value"],crc.motif[crc.motif$group=="NC" & crc.motif$motif=="TCCA","value"],alternative='two.sided')
table(crc.motif$group)

matrix <- crc.motif[crc.motif$group!="CRC",]
library(ggplot2)
library(ggpubr)
library(ggsci)
b <- runif(nrow(matrix), -0.2, 0.2)
#p <- compare
ggplot(matrix, aes(x = motif, y = value, fill=group))+  # , fill=group
  labs(y="Motif frequency(%)")+   # x=x.title,title = title
  geom_boxplot(position=position_dodge(0.5),width=0.5,size=0.4,
               outlier.alpha = 0, outlier.size = 0)+
  #geom_point(aes(x=as.numeric(group)+b,y=value),size=2,color="black", stroke=0.5,alpha=0.6) +
  #geom_violin(draw_quantiles = T) +
  #geom_hline(yintercept = 0, color="grey50",linetype="dashed") + 
  #ylim(c(format(min(matrix$value),digit=2),format(max(matrix$value),digit=2))) +
  #scale_fill_manual(values = c("steelblue4","firebrick2"))+ # "grey50",
  ggsci::scale_fill_nejm()+
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1 ,vjust = 0.5), # 
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) +
  ggpubr::stat_compare_means(#comparisons = list(c("CRC","NC")),# ref.group = "NC",
    method = "wilcox.test",
    #method.args = list(alternative = "greater"),
    label = "p.signif",size = 8, vjust=0.6, hide.ns=F,
  )






# motif logo count mat ----------------------------------------------------

count <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/motif5/sum/fre.txt",header = T,row.names = 1,check.names = F)
count[1:3,1:3]

sample <- "NC-PKU-10-wgs"
count.tbl <- count[,sample,drop=F]
colnames(count.tbl) <- "sample"
count.tbl$motif <- rownames(count.tbl)
#count.tbl$sample <- sample

count.tbl <- as_tibble(count.tbl) %>% 
  # mutate(motif=rownames(count.tbl)) %>% 
  pivot_wider( names_from = "motif",values_from = "sample")

kmer <- 4
seq_df <- data.frame(row.names=1:kmer, pos=0:(kmer-1), "A"=rep(0,kmer), "T"=rep(0,kmer), "C"=rep(0,kmer), "G"=rep(0,kmer) )
seq_fre <- strsplit(colnames(count.tbl),"")
for (i in 1:ncol(count.tbl) ){
  for (j in 1:kmer){
    seq_df[j, seq_fre[[i]][j]] <- count.tbl[,i] 
    #k = k + count.tbl[,i]
  }
}
write.table(seq_df,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/motif5/sum/fre_pwm_NC.txt",quote = F,sep = '\t',row.names = F,col.names = T)
#i <- 1
#j <- 1
#seq_fre[[1]][2]
