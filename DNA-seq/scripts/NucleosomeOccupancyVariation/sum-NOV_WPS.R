# sum WPS (and COV) mat

# get wps matrix (tss exon1end average version) ----------------------------------------------------------
options(stringsAsFactors = F)
setwd("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/")

## read func.
a <- function(x){
  # x <- "./results/intermediate/lulab/table/GRCh38/normalize/TSS--STAD-PKU-37-wgs_WPS_normalized.tsv.gz"
  print(x)
  data <- data.table::fread(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  data <- as.data.frame(data)
  
  sample <- gsub("./results/intermediate/lulab/table/GRCh38/normalize/TSS--","",x)
  sample <- gsub("./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--","",sample)
  sample <- gsub("_WPS_normalized.tsv.gz","",sample)
  sample <- gsub("_COV_normalized.tsv.gz","",sample)
  data[[sample]] <- rowMeans(data[,1000:1200])
  data <- data[,c("V1",sample),drop=FALSE]
  print(data[1:2,])
  print(dim(data))
  
  # return the data
  return(data)
}

## read in  mat
for (t in c("WPS","COV")){
  #t <- "WPS"
  files <- Sys.glob(paste0(paste0("./results/intermediate/lulab/table/GRCh38/normalize/TSS--*_",t,"_normalized.tsv.gz")))  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed
  files2 <- Sys.glob(paste0(paste0("./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--*_",t,"_normalized.tsv.gz")))  
  
  # run the anonymous function defined above
  wps.tss <- lapply(files, a)
  wps.tss <- do.call("cbind", wps.tss)
  wps.tss <- wps.tss[,!duplicated(colnames(wps.tss))]
  #wps.tss <- cbind(rownames(wps.tss),wps.tss)
  
  wps.tss$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.tss$V1)
  #wps.tss$V1 <- unlist(lapply(strsplit(wps.tss$V1,".",fixed=T),function(x) x[1]))
  wps.tss <- wps.tss[!duplicated(wps.tss$V1),]
  #rownames(wps.tss) <- wps.tss$V1
  colnames(wps.tss)[1] <- "gene_id"
  
  
  # run the anonymous function defined above
  wps.exon1 <- lapply(files2, a)
  wps.exon1 <- do.call("cbind", wps.exon1)
  wps.exon1 <- wps.exon1[,!duplicated(colnames(wps.exon1))]
  #wps.exon1 <- cbind(rownames(wps.exon1),wps.exon1)
  
  wps.exon1$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.exon1$V1)
  #wps.exon1$V1 <- unlist(lapply(strsplit(wps.exon1$V1,".",fixed=T),function(x) x[1]))
  wps.exon1 <- wps.exon1[!duplicated(wps.exon1$V1),]
  #rownames(wps.exon1) <- wps.exon1$V1
  colnames(wps.exon1)[1] <- "gene_id"
  
  #wps.exon1$gene_id <- gsub("promoter300100exon1end_","",wps.exon1$gene_id)
  
  g <- unique(wps.tss$gene_id,wps.exon1$gene_id)
  rownames(wps.tss) <- wps.tss$gene_id
  rownames(wps.exon1) <- wps.exon1$gene_id
  wps.tss <- wps.tss[g,]
  wps.exon1 <- wps.exon1[g,]
  #head(g,30)
  rownames(wps.tss) <- g
  rownames(wps.exon1) <- g
  wps.tss <- wps.tss[,colnames(wps.tss)!="gene_id"]
  wps.exon1 <- wps.exon1[,colnames(wps.exon1)!="gene_id"]
  
  table(rownames(wps.tss)==rownames(wps.exon1)) # should be all true
  table(colnames(wps.tss)==colnames(wps.exon1)) # should be all true
  
  wps <- cbind(gene_id=rownames(wps.tss),(wps.tss + wps.exon1)/2)
  wps <- na.omit(wps)
  #print(wps[1:3,1:3])
  
  write.table(wps,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/",t,"_matrix_gene-aver.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
}

#/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/WPS_matrix_gene.txt
#/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/COV_matrix_gene.txt


wps <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/WPS_matrix_gene-aver.txt",sep = "\t",header = T,row.names = 1,check.names = F)
wps[1:4,1:4]
cov <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/COV_matrix_gene-aver.txt",sep = "\t",header = T,row.names = 1,check.names = F)
cov[1:4,1:4]
relDep <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_aver.txt",sep = "\t",header = T,row.names = 1,check.names = F)
relDep[1:4,1:4]
relDep <- relDep[rownames(wps),]
CPM <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/CPM_matrix_NOVsum.txt",sep = "\t",header = T,row.names = 1,check.names = F)
CPM[1:4,1:4]
table(duplicated(rownames(CPM)))
rownames(CPM) <- gsub("NOVsum_","",rownames(CPM))
rownames(CPM) <- gsub("|602","",fixed = T,rownames(CPM))
CPM <- CPM[rownames(wps),]



Hmisc::rcorr(wps$`CRC-PKU-10-wgs`,cov$`CRC-PKU-10-wgs`,type ="spearman")  # -0.93  
Hmisc::rcorr(wps$`CRC-PKU-10-wgs`,relDep$`CRC-PKU-10-wgs`,type ="spearman")  # -0.65  
Hmisc::rcorr(cov$`CRC-PKU-10-wgs`,relDep$`CRC-PKU-10-wgs`,type ="spearman")  # 0.7  
Hmisc::rcorr(cov$`CRC-PKU-10-wgs`,CPM$`CRC-PKU-10-wgs`,type ="spearman")  # 0.67  
Hmisc::rcorr(relDep$`CRC-PKU-10-wgs`,CPM$`CRC-PKU-10-wgs`,type ="spearman")  # 0.78  

cor.test(x=wps$`CRC-PKU-10-wgs`,y=cov$`CRC-PKU-10-wgs`,method = "spearman")     # -0.918875
cor.test(x=wps$`CRC-PKU-10-wgs`,y=relDep$`CRC-PKU-10-wgs`)  # -0.5661182
cor.test(x=cov$`CRC-PKU-10-wgs`,y=relDep$`CRC-PKU-10-wgs`)  # 0.630467
cor.test(x=cov$`CRC-PKU-10-wgs`,y=CPM$`CRC-PKU-10-wgs`)  # 0.5262653
cor.test(x=relDep$`CRC-PKU-10-wgs`,y=CPM$`CRC-PKU-10-wgs`)  # 0.4382428

# wps abnormal !!!!




# get wps matrix (tss only version, newer version, 220117) ----------------------------------------------------------
options(stringsAsFactors = F)
setwd("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/")

## read func.
a <- function(x){
  # x <- "./results/intermediate/lulab/table/GRCh38/normalize/TSS--STAD-PKU-37-wgs_WPS_normalized.tsv.gz"
  print(x)
  data <- data.table::fread(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  data <- as.data.frame(data)
  
  sample <- gsub("./results/intermediate/lulab/table/GRCh38/normalize/TSS--","",x)
  #sample <- gsub("./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--","",sample)
  sample <- gsub("_WPS_normalized.tsv.gz","",sample)
  sample <- gsub("_COV_normalized.tsv.gz","",sample)
  data[[sample]] <- rowMeans(data[,1000:1200])
  data <- data[,c("V1",sample),drop=FALSE]
  print(data[1:2,])
  print(dim(data))
  
  # return the data
  return(data)
}

## read in  mat
for (t in c("WPS","COV")){
  #t <- "WPS"
  files <- Sys.glob(paste0(paste0("./results/intermediate/lulab/table/GRCh38/normalize/TSS--*_",t,"_normalized.tsv.gz")))  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed
  #files2 <- Sys.glob(paste0(paste0("./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--*_",t,"_normalized.tsv.gz")))  
  
  # run the anonymous function defined above
  wps.tss <- lapply(files, a)
  wps.tss <- do.call("cbind", wps.tss)
  wps.tss <- wps.tss[,!duplicated(colnames(wps.tss))]
  #wps.tss <- cbind(rownames(wps.tss),wps.tss)
  
  wps.tss$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.tss$V1)
  #wps.tss$V1 <- unlist(lapply(strsplit(wps.tss$V1,".",fixed=T),function(x) x[1]))
  wps.tss <- wps.tss[!duplicated(wps.tss$V1),]
  #rownames(wps.tss) <- wps.tss$V1
  colnames(wps.tss)[1] <- "gene_id"
  
  
  # run the anonymous function defined above
  #wps.exon1 <- lapply(files2, a)
  #wps.exon1 <- do.call("cbind", wps.exon1)
  #wps.exon1 <- wps.exon1[,!duplicated(colnames(wps.exon1))]
  #wps.exon1 <- cbind(rownames(wps.exon1),wps.exon1)
  
  #wps.exon1$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.exon1$V1)
  #wps.exon1$V1 <- unlist(lapply(strsplit(wps.exon1$V1,".",fixed=T),function(x) x[1]))
  #wps.exon1 <- wps.exon1[!duplicated(wps.exon1$V1),]
  #rownames(wps.exon1) <- wps.exon1$V1
  #colnames(wps.exon1)[1] <- "gene_id"
  
  #wps.exon1$gene_id <- gsub("promoter300100exon1end_","",wps.exon1$gene_id)
  
  #g <- unique(wps.tss$gene_id,wps.exon1$gene_id)
  g <- unique(wps.tss$gene_id)
  rownames(wps.tss) <- wps.tss$gene_id
  #rownames(wps.exon1) <- wps.exon1$gene_id
  wps.tss <- wps.tss[g,]
  #wps.exon1 <- wps.exon1[g,]
  #head(g,30)
  rownames(wps.tss) <- g
  #rownames(wps.exon1) <- g
  wps.tss <- wps.tss[,colnames(wps.tss)!="gene_id"]
  #wps.exon1 <- wps.exon1[,colnames(wps.exon1)!="gene_id"]
  
  #table(rownames(wps.tss)==rownames(wps.exon1)) # should be all true
  #table(colnames(wps.tss)==colnames(wps.exon1)) # should be all true
  
  #wps <- cbind(gene_id=rownames(wps.tss),(wps.tss + wps.exon1)/2)
  wps <- cbind(gene_id=rownames(wps.tss),wps.tss)
  wps <- na.omit(wps)
  #print(wps[1:3,1:3])
  
  write.table(wps,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/",t,"_matrix_150TSS50.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
}

#/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/WPS_matrix_gene.txt
#/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/COV_matrix_gene.txt


wps <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/WPS_matrix_gene-aver.txt",sep = "\t",header = T,row.names = 1,check.names = F)
wps[1:4,1:4]
cov <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/COV_matrix_gene-aver.txt",sep = "\t",header = T,row.names = 1,check.names = F)
cov[1:4,1:4]
relDep <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_aver.txt",sep = "\t",header = T,row.names = 1,check.names = F)
relDep[1:4,1:4]
relDep <- relDep[rownames(wps),]
CPM <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/CPM_matrix_NOVsum.txt",sep = "\t",header = T,row.names = 1,check.names = F)
CPM[1:4,1:4]
table(duplicated(rownames(CPM)))
rownames(CPM) <- gsub("NOVsum_","",rownames(CPM))
rownames(CPM) <- gsub("|602","",fixed = T,rownames(CPM))
CPM <- CPM[rownames(wps),]



Hmisc::rcorr(wps$`CRC-PKU-10-wgs`,cov$`CRC-PKU-10-wgs`,type ="spearman")  # -0.93  
Hmisc::rcorr(wps$`CRC-PKU-10-wgs`,relDep$`CRC-PKU-10-wgs`,type ="spearman")  # -0.65  
Hmisc::rcorr(cov$`CRC-PKU-10-wgs`,relDep$`CRC-PKU-10-wgs`,type ="spearman")  # 0.7  
Hmisc::rcorr(cov$`CRC-PKU-10-wgs`,CPM$`CRC-PKU-10-wgs`,type ="spearman")  # 0.67  
Hmisc::rcorr(relDep$`CRC-PKU-10-wgs`,CPM$`CRC-PKU-10-wgs`,type ="spearman")  # 0.78  

cor.test(x=wps$`CRC-PKU-10-wgs`,y=cov$`CRC-PKU-10-wgs`,method = "spearman")     # -0.918875
cor.test(x=wps$`CRC-PKU-10-wgs`,y=relDep$`CRC-PKU-10-wgs`)  # -0.5661182
cor.test(x=cov$`CRC-PKU-10-wgs`,y=relDep$`CRC-PKU-10-wgs`)  # 0.630467
cor.test(x=cov$`CRC-PKU-10-wgs`,y=CPM$`CRC-PKU-10-wgs`)  # 0.5262653
cor.test(x=relDep$`CRC-PKU-10-wgs`,y=CPM$`CRC-PKU-10-wgs`)  # 0.4382428

# wps abnormal !!!!




# get wps matrix (exon1 only version, newer version, 220117) ----------------------------------------------------------
options(stringsAsFactors = F)
setwd("/BioII/lulab_b/baopengfei/gitsoft/cfDNA-2016/")

## read func.
a <- function(x){
  # x <- "./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--STAD-PKU-37-wgs_WPS_normalized.tsv.gz"
  print(x)
  data <- data.table::fread(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  data <- as.data.frame(data)
  
  #sample <- gsub("./results/intermediate/lulab/table/GRCh38/normalize/TSS--","",x)
  sample <- gsub("./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--","",x)
  sample <- gsub("_WPS_normalized.tsv.gz","",sample)
  sample <- gsub("_COV_normalized.tsv.gz","",sample)
  data[[sample]] <- rowMeans(data[,1000:1200])
  data <- data[,c("V1",sample),drop=FALSE]
  print(data[1:2,])
  print(dim(data))
  
  # return the data
  return(data)
}

## read in  mat
for (t in c("WPS","COV")){
  #t <- "WPS"
  #files <- Sys.glob(paste0(paste0("./results/intermediate/lulab/table/GRCh38/normalize/TSS--*_",t,"_normalized.tsv.gz")))  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed
  files <- Sys.glob(paste0(paste0("./results/intermediate/lulab/table/GRCh38/normalize/Exon1end--*_",t,"_normalized.tsv.gz")))  
  
  # run the anonymous function defined above
  wps.tss <- lapply(files, a)
  wps.tss <- do.call("cbind", wps.tss)
  wps.tss <- wps.tss[,!duplicated(colnames(wps.tss))]
  #wps.tss <- cbind(rownames(wps.tss),wps.tss)
  
  wps.tss$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.tss$V1)
  #wps.tss$V1 <- unlist(lapply(strsplit(wps.tss$V1,".",fixed=T),function(x) x[1]))
  wps.tss <- wps.tss[!duplicated(wps.tss$V1),]
  #rownames(wps.tss) <- wps.tss$V1
  colnames(wps.tss)[1] <- "gene_id"
  
  
  # run the anonymous function defined above
  #wps.exon1 <- lapply(files2, a)
  #wps.exon1 <- do.call("cbind", wps.exon1)
  #wps.exon1 <- wps.exon1[,!duplicated(colnames(wps.exon1))]
  #wps.exon1 <- cbind(rownames(wps.exon1),wps.exon1)
  
  #wps.exon1$V1 <- gsub("promoter150TSS50_|promoter300100exon1end_","",wps.exon1$V1)
  #wps.exon1$V1 <- unlist(lapply(strsplit(wps.exon1$V1,".",fixed=T),function(x) x[1]))
  #wps.exon1 <- wps.exon1[!duplicated(wps.exon1$V1),]
  #rownames(wps.exon1) <- wps.exon1$V1
  #colnames(wps.exon1)[1] <- "gene_id"
  
  #wps.exon1$gene_id <- gsub("promoter300100exon1end_","",wps.exon1$gene_id)
  
  #g <- unique(wps.tss$gene_id,wps.exon1$gene_id)
  g <- unique(wps.tss$gene_id)
  rownames(wps.tss) <- wps.tss$gene_id
  #rownames(wps.exon1) <- wps.exon1$gene_id
  wps.tss <- wps.tss[g,]
  #wps.exon1 <- wps.exon1[g,]
  #head(g,30)
  rownames(wps.tss) <- g
  #rownames(wps.exon1) <- g
  wps.tss <- wps.tss[,colnames(wps.tss)!="gene_id"]
  #wps.exon1 <- wps.exon1[,colnames(wps.exon1)!="gene_id"]
  
  #table(rownames(wps.tss)==rownames(wps.exon1)) # should be all true
  #table(colnames(wps.tss)==colnames(wps.exon1)) # should be all true
  
  #wps <- cbind(gene_id=rownames(wps.tss),(wps.tss + wps.exon1)/2)
  wps <- cbind(gene_id=rownames(wps.tss),wps.tss)
  wps <- na.omit(wps)
  #print(wps[1:3,1:3])
  
  write.table(wps,paste0("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/",t,"_matrix_300100exon1end.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
}

#/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/WPS_matrix_gene.txt
#/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/COV_matrix_gene.txt

# wps abnormal !!!!







#
tmp <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/exon1.bed",header = F)
head(tmp)



# test exon1 len distribution -----------------------------------------------------
hist(abs(tmp$V3-tmp$V2),xlim = c(0,1000),breaks = 8000)


#ly <- data.table::fread("/BioII/lulab_b/liyu/projects/multiomics/DNA-seq/output/matrix/promoter.ndr.ave.cap2.withNA.txt", data.table = F,sep = "\t",header = T,check.names = F)
ly <- read.table("/BioII/lulab_b/liyu/projects/multiomics/DNA-seq/output/matrix/promoter.ndr.ave.cap2.withNA.txt",row.names = 1, sep = "\t", header = T, check.names = F)
ly[1:3,1:3]



bpf <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/NucOccuVarRelDepRatio_matrix_promoter150TSS50.txt", row.names = 1, sep = "\t", header = T, check.names = F)
rownames(bpf) <- gsub("promoter150TSS50_","",rownames(bpf))
bpf[1:3,]

li <- list()
for(i in (intersect(colnames(ly),colnames(bpf)))){
print(i)
  ly.tmp <- ly[,i,drop=F]
ly.tmp <- na.omit(ly.tmp)
bpf.tmp <- bpf[,i,drop=F]
bpf.tmp <- bpf.tmp[rownames(ly.tmp),,drop=F]

all(rownames(bpf.tmp) == rownames(ly.tmp))

#hist(bpf.tmp$`CRC-PKU-10-wgs`)
#hist(ly.tmp$`CRC-PKU-10-wgs`)

l <- (bpf.tmp[,1]!=0 & bpf.tmp[,1]!=2) & (ly.tmp[,1]!=0 & ly.tmp[,1]!=2) 
table(l)

#plot(x=bpf.tmp[l,1],y=ly.tmp[l,1])
p <- cor.test(bpf.tmp[,1],ly.tmp[,1]) # ,method = "spearman"
#0.44
#0.5
li[[i]] <- p$estimate
}
li.df <- do.call(rbind,li)
hist(li.df)

