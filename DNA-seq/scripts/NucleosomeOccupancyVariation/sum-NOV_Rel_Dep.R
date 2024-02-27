# sum RelativeDepth Of NDR
# last 211125 by bpf

#### sum lulab NDR (newer) ####
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
rm(list = ls())

# install.packages("stringr")
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))

# create function to read in and format data
a <- function(x){
  # read data and set column names
  print(x)
  data <- read.table(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  data <- as_tibble(data)
  sample.name <- as.character(lapply(str_split(string = x,pattern = "_", ), function(x) x[1]))
  sample.name <- as.character(lapply(str_split(string = sample.name,pattern = "/", ), function(x) x[6]))
  
  colnames(data) <- c("chr", "start", "end", "name", "phase","strand",sample.name)
  
  # return the data
  return(data)
}

# get locations of all cn files
{
  files <- Sys.glob("./output/lulab/NDR/depth/*_promoter300100exon1end.bed") 
  # run the anonymous function defined above
  tmp <- lapply(files, a)
  # turn the list of data frames into a single data frame
  tmp <- do.call("cbind", tmp)
  tmp <- tmp[,!duplicated(colnames(tmp))]
  rownames(tmp) <- tmp$name
  tmp <- tmp[,7:ncol(tmp)]
  promoter300100exon1end <- tmp
}
{
  files <- Sys.glob("./output/lulab/NDR/depth/*_promoter150TSS50.bed")  # CRC-PKU-10-wgs
  # run the anonymous function defined above
  tmp <- lapply(files, a)
  
  # turn the list of data frames into a single data frame
  tmp <- do.call("cbind", tmp)
  tmp <- tmp[,!duplicated(colnames(tmp))]
  rownames(tmp) <- tmp$name
  tmp <- tmp[,7:ncol(tmp)]
  promoter150TSS50 <- tmp
}
{
  files <- Sys.glob("./output/lulab/NDR/depth/*_20001000TSS.bed")  # CRC-PKU-10-wgs
  # run the anonymous function defined above
  tmp <- lapply(files, a)
  
  # turn the list of data frames into a single data frame
  tmp <- do.call("cbind", tmp)
  tmp <- tmp[,!duplicated(colnames(tmp))]
  rownames(tmp) <- tmp$name
  tmp <- tmp[,7:ncol(tmp)]
  ref.TSS <- tmp
}
{
  files <- Sys.glob("./output/lulab/NDR/depth/*_exon1end10002000.bed")  # CRC-PKU-10-wgs
  # run the anonymous function defined above
  tmp <- lapply(files, a)
  
  # turn the list of data frames into a single data frame
  tmp <- do.call("cbind", tmp)
  tmp <- tmp[,!duplicated(colnames(tmp))]
  rownames(tmp) <- tmp$name
  tmp <- tmp[,7:ncol(tmp)]
  ref.exon1end <- tmp
}

{ # add in 220115
  files <- Sys.glob("./output/lulab/NDR/depth/*_TSS10002000.bed")  # CRC-PKU-10-wgs
  # run the anonymous function defined above
  tmp <- lapply(files, a)
  
  # turn the list of data frames into a single data frame
  tmp <- do.call("cbind", tmp)
  tmp <- tmp[,!duplicated(colnames(tmp))]
  rownames(tmp) <- tmp$name
  tmp <- tmp[,7:ncol(tmp)]
  ref2.TSS <- tmp
}
{ # add in 220115
  files <- Sys.glob("./output/lulab/NDR/depth/*_20001000exon1end.bed")  # CRC-PKU-10-wgs
  # run the anonymous function defined above
  tmp <- lapply(files, a)
  
  # turn the list of data frames into a single data frame
  tmp <- do.call("cbind", tmp)
  tmp <- tmp[,!duplicated(colnames(tmp))]
  rownames(tmp) <- tmp$name
  tmp <- tmp[,7:ncol(tmp)]
  ref2.exon1end <- tmp
}

# check col names (should be all true !!!)
all(colnames(promoter300100exon1end)==colnames(promoter150TSS50))
all(colnames(promoter150TSS50)==colnames(ref.TSS))
all(colnames(promoter150TSS50)==colnames(ref2.TSS))
ref.TSS <- ref.TSS[,colnames(promoter150TSS50)]
promoter150TSS50 <- promoter150TSS50[,colnames(promoter150TSS50)]
all(colnames(promoter300100exon1end)==colnames(ref.exon1end))
all(colnames(promoter300100exon1end)==colnames(ref2.exon1end))
ref.exon1end <- ref.exon1end[,colnames(promoter300100exon1end)]
ref2.exon1end <- ref2.exon1end[,colnames(promoter300100exon1end)]


# check row names (should be all true !!!)
ref.TSS.1 <- ref.TSS
ref.TSS.1 <- ref.TSS.1[rownames(promoter150TSS50),]
table(is.na(ref.TSS.1))
ref.TSS.2 <- ref2.TSS
ref.TSS.2 <- ref.TSS.2[rownames(promoter150TSS50),]
table(is.na(ref.TSS.2))

rownames(ref.TSS.1) <- rownames(promoter150TSS50)
rownames(ref.TSS.2) <- rownames(promoter150TSS50)
all(rownames(promoter150TSS50)==rownames(ref.TSS.1))
all(rownames(promoter150TSS50)==rownames(ref.TSS.2))


ref.exon1end.1 <- ref.exon1end
rownames(ref.exon1end.1) <- gsub("promoter300exon1end100_","promoter300100exon1end_",rownames(ref.exon1end.1))
ref.exon1end.1 <- ref.exon1end.1[rownames(promoter300100exon1end),]
table(is.na(ref.exon1end.1))
rownames(ref.exon1end.1) <- rownames(promoter300100exon1end)
all(colnames(promoter300100exon1end)==colnames(ref.exon1end.1))

ref.exon1end.2 <- ref2.exon1end
rownames(ref.exon1end.2) <- gsub("promoter300exon1end100_","promoter300100exon1end_",rownames(ref.exon1end.2))
ref.exon1end.2 <- ref.exon1end.2[rownames(promoter300100exon1end),]
table(is.na(ref.exon1end.2))
rownames(ref.exon1end.2) <- rownames(promoter300100exon1end)
all(colnames(promoter300100exon1end)==colnames(ref.exon1end.2))


# cal TSS 
#fill NA by mat col mean
for(i in 1:ncol(ref.TSS.1)) {
  ref.TSS.1[ , i][is.na(ref.TSS.1[ , i])] <- mean(ref.TSS.1[ , i], na.rm=TRUE)
}
for(i in 1:ncol(ref.TSS.2)) {
  ref.TSS.2[ , i][is.na(ref.TSS.2[ , i])] <- mean(ref.TSS.2[ , i], na.rm=TRUE)
}

## cap2 after cal ratio (difference: original paper NDR cor are sing-base resolution)
#table(promoter150TSS50>2)
#hist(promoter150TSS50)
tss <- 0.5*ref.TSS.1 + 0.5*ref.TSS.2
#table(is.na(tss)) # if not fill NA
#FALSE    TRUE 
#5374656   46128
tss <- promoter150TSS50/(tss+0.01)

tss[tss>2] <- 2
hist(as.matrix(tss),xlim = c(0,10),breaks = 10000)
mean(as.matrix(tss),na.rm = T)
median(as.matrix(tss),na.rm = T)
max(as.matrix(tss),na.rm = T)
min(as.matrix(tss),na.rm = T)
#tss.log <- log2(tss+0.1)
#write.table(tss.log,"./output/lulab/matrix/NucOccuVarLog2RelativeRatio_matrix_promoter150TSS50.txt",quote = F,row.names = T,col.names = T,sep = "\t")
write.table(cbind("gene_id"=rownames(tss),tss),"./output/lulab/matrix/NucOccuVarRelDepRatio_matrix_promoter150TSS50.txt",quote = F,row.names = F,col.names = T,sep = "\t")
#table(is.na(tss))
#tss0 <- read.table("./output/lulab/matrix/NucOccuVarRelDepRatio_matrix_promoter150TSS50.txt",sep = "\t",header = T,row.names = 1,check.names = F)
#head(tss0[,1:5])
#cor.test(tss0$`CRC-PKU-10-wgs`,tss$`CRC-PKU-10-wgs`,method = "spearman")

# cal Exon1end  (imputation)
#fill NA by mat col mean
for(i in 1:ncol(ref.exon1end.1)) {
  ref.exon1end.1[ , i][is.na(ref.exon1end.1[ , i])] <- mean(ref.exon1end.1[ , i], na.rm=TRUE)
}
for(i in 1:ncol(ref.exon1end.2)) {
  ref.exon1end.2[ , i][is.na(ref.exon1end.2[ , i])] <- mean(ref.exon1end.2[ , i], na.rm=TRUE)
}
exon1 <- 0.5*ref.exon1end.1 + 0.5*ref.exon1end.2
exon1 <- promoter300100exon1end/(exon1+0.01)

exon1[exon1>2] <- 2
hist(as.matrix(exon1))
mean(as.matrix(exon1),na.rm = T)
median(as.matrix(exon1),na.rm = T)
max(as.matrix(exon1),na.rm = T)
min(as.matrix(exon1),na.rm = T)
#exon1.log <- log2(exon1+0.1)
#write.table(exon1.log,"./output/lulab/matrix/NucOccuVarLog2RelativeRatio_matrix_promoter150TSS50.txt",quote = F,row.names = T,col.names = T,sep = "\t")
write.table(cbind("gene_id"=rownames(exon1),exon1),"./output/lulab/matrix/NucOccuVarRelDepRatio_matrix_promoter300100exon1end.txt",quote = F,row.names = F,col.names = T,sep = "\t")

hist(ref.exon1end.1$`NC-PKU-10-wgs`,breaks = 10000,xlim = c(0,20))
hist(ref.exon1end.1$`CRC-PKU-10-wgs`,breaks = 10000,xlim = c(0,20))


# cal average (ratio)  (deprecated in 220115)
# rownames(tss) <- gsub("promoter150TSS50_|promoter300100exon1end_","",rownames(tss))
# rownames(exon1end) <- gsub("promoter150TSS50_|promoter300100exon1end_","",rownames(exon1end))
# exon1end <- exon1end[rownames(tss),]
# all(rownames(tss)==rownames(exon1end)) # shoud be the same
# aver <- (tss+exon1end)/2
# write.table(cbind("gene_id"=rownames(aver),aver),"./output/lulab/matrix/NucOccuVarRelDepRatio_matrix_aver.txt",quote = F,row.names = F,col.names = T,sep = "\t")



# aver <- data.table::fread("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_aver.txt",header = T,sep = "\t",check.names = F,data.table = F)
# rownames(aver) <- gsub("promoter150TSS50_|promoter300100exon1end_","",perl = T,aver$gene_id)
# aver$gene_id <- NULL
# 
# exon1 <- data.table::fread("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_promoter300100exon1end.txt",header = T,sep = "\t",check.names = F,data.table = F)
# rownames(exon1) <- gsub("promoter150TSS50_|promoter300100exon1end_","",perl = T,exon1$gene_id)
# exon1$gene_id <- NULL
# exon1 <- exon1[rownames(aver),]
#   
# tss <- data.table::fread("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_promoter150TSS50.txt",header = T,sep = "\t",check.names = F,data.table = F)
# rownames(tss) <- gsub("promoter150TSS50_|promoter300100exon1end_","",perl = T,tss$gene_id)
# tss$gene_id <- NULL
# tss <- tss[rownames(aver),]
# #tss <- na.omit(tss)
# 
# head(exon1[,1:4])
# all(colnames(exon1)==colnames(tss))
# all(colnames(aver)==colnames(tss))
# 
# all(rownames(exon1) == rownames(tss))
# all(rownames(aver) == rownames(tss))
# aver2 <- (exon1+tss)/2
# aver[1:4,1:4]
# plot(x = aver$`CRC-PKU-10-wgs`,y = aver2$`CRC-PKU-10-wgs`)
# cor.test(aver$`STAD-PKU-9-wgs`,aver2$`STAD-PKU-9-wgs`)



# cor betweent 2 flanks:  exon1end and tss --------------------------------
exon1_1 <- apply(X = ref.exon1end.1[,grep("-PKU-",colnames(ref.exon1end.1))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
exon1_2 <- apply(X = ref.exon1end.2[,grep("-PKU-",colnames(ref.exon1end.2))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tss_1 <- apply(X = ref.TSS.1[,grep("-PKU-",colnames(ref.TSS.1))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))
tss_2 <- apply(X = ref.TSS.2[,grep("-PKU-",colnames(ref.TSS.2))], MARGIN = 1, FUN = function(x) edgeR::gini(as.numeric(x)))

## read blood high exp ensg
blood <- read.table("ref/gtf/bloodNDR/promoter2000exon1end2000_bloodTxTop1000.bed",header = F,sep = "\t")$V4
blood <- gsub("promoter300100exon1end_","",blood)
blood <- unlist(lapply(strsplit(blood,"|",fixed = T),function(x) x[1]))

df <- rbind(exon1_1=exon1_1,exon1_2=exon1_2,tss_1=tss_1,tss_2=tss_2)
df <- as.data.frame(t(df))
rownames(df) <- gsub("promoter300100exon1end_","",rownames(df))
df <- df[rownames(df) %in% blood,]
df[1:3,]
#head(df)
plot(x=exon1_1,y=exon1_2)
plot(x=tss_1,y=tss_2)

for (i in 1:ncol(df)){
  df[,i][is.na(df[,i])] <- mean(df[,i],na.rm=T)
}

df.m <- reshape2::melt(df)
#df.m <- na.omit(df.m)
summary(df.m$value[df.m$variable=="exon1_2"])

ggplot(df.m,aes(x=variable, y=value, fill=variable)) +
  geom_boxplot(color="black",shape=21,size=1,alpha=0.6)+
  ggsci::scale_fill_nejm()+
  #ylim(c(0.3,0.7))+
  ggpubr::stat_compare_means(paired = T, comparisons=list(c("exon1_1","exon1_2"),c("tss_1","tss_2")) ,
                             #mapping = aes(group=Group), #comparisons = ref.group = "NC",
                             bracket.size = 0.3, size = 6,color="grey30", method = "wilcox.test",hide.ns=T) + # vjust = 2, hjust=1.2,
  theme_minimal() +
  theme(axis.text=element_text(size=24,face="bold"),
        strip.text.x = element_text(size = 28,face="bold", colour = "grey30"),
        strip.text.y = element_text(size = 28,face="bold", colour = "grey30"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=32,face="bold"),
        legend.title=element_text(size=28,face="bold")) +
  ylab("cfDNA Nucleosome Flank Depth Gini") +
  ggtitle("top1000 blood exp") +
  guides(color=guide_legend("Group"))
#看起来只有exon1两个侧翼区域有略微区别，且符合预期



# cor betweent 2 regions:  exon1end and tss -------------------------------
### sample-centric
s <- list()
for(i in colnames(tss)){
  #i <- "ENSG00000278757.1"
  print(i)
  rho <- Hmisc::rcorr(as.numeric(exon1end[,i]),as.numeric(tss[,i]),type ="spearman")
  s[[i]] <- data.frame(sample=i,rho=rho$r["y","x"],p=rho$P["y","x"])
  #print(paste0(i,":",rho))
}
#0.4~0.5
s <- do.call("rbind",s)
hist(s$rho)


### gene-centric
g <- list()
for(i in rownames(tss)){
  #i <- "ENSG00000278757.1"
  print(i)
  rho <- Hmisc::rcorr(as.numeric(exon1end[i,]),as.numeric(tss[i,]),type ="spearman")
  g[[i]] <- data.frame(gene=i,rho=rho$r["y","x"],p=rho$P["y","x"])
  #print(paste0(i,":",rho))
}
#-0.2~1
g <- do.call("rbind",g)
hist(g$rho)
