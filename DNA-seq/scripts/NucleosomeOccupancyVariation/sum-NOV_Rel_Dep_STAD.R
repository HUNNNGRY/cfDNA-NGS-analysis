# sum RelativeDepth Of NDR

#### sum in-house NDR and correlation ####
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
rm(list = ls())

# install.packages("stringr")
suppressPackageStartupMessages(library("stringr"))
# get locations of all cn files
files <- Sys.glob("./output/lulab/NDR/depth/*_exon1end10002000.bed")  # CRC-PKU-10-wgs     _promoter150TSS50.bed  _20001000TSS.bed  _promoter300100exon1end.bed  _exon1end10002000.bed

# create function to read in and format data
a <- function(x){
  # read data and set column names
  data <- read.table(x, stringsAsFactors = F, sep = "\t", header=FALSE,check.names = F)
  data <- as_tibble(data)
  sample.name <- as.character(lapply(str_split(string = x,pattern = "_", ), function(x) x[1]))
  sample.name <- as.character(lapply(str_split(string = sample.name,pattern = "/", ), function(x) x[6]))
  
  colnames(data) <- c("chr", "start", "end", "name", "phase","strand",sample.name)
  
  # return the data
  return(data)
}

# run the anonymous function defined above
cnData <- lapply(files, a)

# turn the list of data frames into a single data frame
cnData <- do.call("cbind", cnData)
cnData <- cnData[,!duplicated(colnames(cnData))]
rownames(cnData) <- cnData$name
cnData <- cnData[,7:ncol(cnData)]


# get 4 depth
promoter150TSS50 <- cnData

ref.TSS <- cnData

promoter300100exon1end <- cnData

ref.exon1end <- cnData

# check col names
all(colnames(promoter150TSS50)==colnames(ref.TSS))
promoter150TSS50 <- promoter150TSS50[,colnames(ref.TSS)]

all(colnames(promoter300100exon1end)==colnames(ref.exon1end))
promoter300100exon1end <- promoter300100exon1end[,colnames(ref.TSS)]
ref.exon1end <- ref.exon1end[,colnames(ref.TSS)]

# check rownames
ref.TSS.1 <- ref.TSS
ref.TSS.1 <- ref.TSS.1[rownames(promoter150TSS50),]
table(is.na(ref.TSS.1))
rownames(ref.TSS.1) <- rownames(promoter150TSS50)
ref.TSS.1 <- as.matrix(ref.TSS.1)
all(rownames(promoter150TSS50)==rownames(ref.TSS.1))


ref.exon1end.1 <- ref.exon1end
rownames(ref.exon1end.1) <- sub("promoter300exon1end100","promoter300100exon1end",rownames(ref.exon1end.1))
ref.exon1end.1 <- ref.exon1end.1[rownames(promoter300100exon1end),]
#promoter300100exon1end <- promoter300100exon1end[rownames(promoter150TSS50),]
table(is.na(ref.exon1end.1))
rownames(ref.exon1end.1) <- rownames(promoter300100exon1end)
ref.exon1end.1 <- as.matrix(ref.exon1end.1)
all(colnames(promoter300100exon1end)==colnames(ref.exon1end))


# cal TSS (log2ratio)
tss <- promoter150TSS50/(ref.TSS.1+0.01)
mean(as.matrix(tss),na.rm = T)
median(as.matrix(tss),na.rm = T)
max(as.matrix(tss),na.rm = T)
min(as.matrix(tss),na.rm = T)
hist(as.matrix(tss),breaks = 800000,xlim = c(0,10))
table(as.matrix(tss)==0)
tss.log <- log2(tss+0.01)
mean(as.matrix(tss.log),na.rm = T)
median(as.matrix(tss.log),na.rm = T)
max(as.matrix(tss.log),na.rm = T)
min(as.matrix(tss.log),na.rm = T)
write.table(tss.log,"./output/lulab/matrix/NOVlog2Ratio_matrix_promoter150TSS50.txt",quote = F,row.names = T,col.names = T,sep = "\t")

# cal exon1end (log2ratio)
exon1end <- promoter300100exon1end/(ref.exon1end.1+0.01)
mean(as.matrix(exon1end),na.rm = T)
median(as.matrix(exon1end),na.rm = T)
max(as.matrix(exon1end),na.rm = T)
min(as.matrix(exon1end),na.rm = T)
hist(as.matrix(exon1end),breaks = 800000,xlim = c(0,10))
table(as.matrix(exon1end)==0)
exon1end.log <- log2(exon1end+0.01)
write.table(exon1end.log,"./output/lulab/matrix/NOVlog2Ratio_matrix_promoter300100exon1end.txt",quote = F,row.names = T,col.names = T,sep = "\t")

# cal average (log2ratio)
#cal mean with na removed
#https://stackoverflow.com/questions/9424311/how-to-get-mean-median-and-other-statistics-over-entire-matrix-array-or-dataf
# A <- matrix(c(2,4,3,5), 2)
# B <- matrix(c(6,8,7,9), 2)
# X <- list(tss.log, exon1end.log)
# Y <- do.call(cbind, X)
# Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
# aver <- apply(Y, c(1, 2), mean, na.rm = TRUE)

# aver <- (tss.log+exon1end.log)/2
# write.table(aver,"./output/lulab/matrix/NucOccuVarLog2RelativeRatio_matrix_aver.txt",quote = F,row.names = T,col.names = T,sep = "\t")
# aver <- read.table("./output/lulab/matrix/NucOccuVarLog2RelativeRatio_matrix_aver.txt",header = T,row.names = 1,sep = "\t")
aver <- (tss.log+exon1end.log)/2   # row id not the same order !!!
rownames(tss.log) <- gsub("promoter150TSS50_|promoter300100exon1end_","",rownames(tss.log))
rownames(exon1end.log) <- gsub("promoter150TSS50_|promoter300100exon1end_","",rownames(exon1end.log))
table(rownames(tss.log)==rownames(exon1end.log))
#FALSE  TRUE 
#7787 50501

#write.table(aver,"./output/lulab/matrix/NOVlog2Ratio_matrix_mean.txt",quote = F,row.names = T,col.names = T,sep = "\t")
#aver <- read.table("./output/lulab/matrix/NOVlog2Ratio_matrix_mean.txt",header = T,row.names = 1,sep = "\t")

table(is.na(min)) 
min <- pmin(tss.log,exon1end.log,na.rm = T)
write.table(min,"./output/lulab/matrix/NOVlog2Ratio_matrix_min.txt",quote = F,row.names = T,col.names = T,sep = "\t")
min <- read.table("./output/lulab/matrix/NOVlog2Ratio_matrix_min.txt",header = T,row.names = 1,sep = "\t")
NOVlog2Ratio_matrix_mean.txt

# cor betweent 2 regions:  exon1end.log and tss.log
for(i in colnames(exon1end.log)){
  rho <- cor(exon1end.log[[i]],tss.log[[i]],use = "complete.obs", method ="spearman")
  print(paste0(i,":",rho))
}

# cor betweent 3 regions:  exon1end.log, tss.log and featurecount
exon1end.log <- read.table("./output/lulab/matrix/NOVlog2Ratio_matrix_promoter300100exon1end.txt",header = T,row.names = 1,sep = "\t")
tss.log <- read.table("./output/lulab/matrix/NOVlog2Ratio_matrix_promoter150TSS50.txt",header = T,row.names = 1,sep = "\t")
featurecount <- read.table("./output/lulab/matrix/CPM_matrix_promoter300100exon1end.txt",header = T,row.names = 1,sep = "\t")
#featurecount <- read.table("./output/lulab/matrix/CPM_matrix_promoter150TSS50.txt",header = T,row.names = 1,sep = "\t")
featurecount <- log2(featurecount+1)
rownames(featurecount) <- as.character(lapply(str_split(rownames(featurecount),"\\|"),function(x) x[1]))
rownames(featurecount) <- sub("promoter300100exon1end_|promoter150TSS50_","",rownames(featurecount))
rownames(exon1end.log) <- sub("promoter300100exon1end_|promoter150TSS50_","",rownames(exon1end.log))
rownames(tss.log) <- sub("promoter300100exon1end_|promoter150TSS50_","",rownames(tss.log))

table(rownames(exon1end.log) %in% rownames(tss.log))
all(rownames(exon1end.log)==rownames(tss.log))
table(rownames(exon1end.log)==rownames(featurecount))
all(rownames(exon1end.log)==rownames(featurecount))
exon1end.log <- exon1end.log[rownames(featurecount),]
tss.log <- tss.log[rownames(featurecount),]

all(colnames(exon1end.log)==colnames(tss.log))

for(i in colnames(exon1end.log)){
  rho <- cor(exon1end.log[[i]],featurecount[[i]],use = "complete.obs", method ="spearman")
  print(paste0(i,":",rho))
}

for(i in colnames(exon1end.log)){
  rho <- cor(tss.log[[i]],featurecount[[i]],use = "complete.obs", method ="spearman")
  print(paste0(i,":",rho))
}

#featurecount1 <- featurecount
for(i in colnames(exon1end.log)){
  rho <- cor(featurecount1[[i]],featurecount[[i]],use = "complete.obs", method ="spearman")
  print(paste0(i,":",rho))
}

plot(x=featurecount1$NC.PKU.mix6.wgs,y=featurecount$NC.PKU.mix6.wgs,shape=19,size=0.02,alpha=0.1)
# seems no much difference between featurecount and relative depth





#### TCGA STAD vs GTEx normal blood ####
rm(list = ls())

# GTEx
gtex.meta <- read.table("/BioII/lulab_b/baopengfei/projects/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep = "\t",header = T,check.names = F,fill = T)
table(gtex.meta$SMTS)
gtex.meta <- gtex.meta[gtex.meta$SMTS=="Blood" & gtex.meta$SMTSD=="Whole Blood",]

gtex.data <- data.table::fread("/BioII/lulab_b/baopengfei/projects/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",nThread = 10, stringsAsFactors = F, check.names = F, sep = "\t",data.table = F,header = T, skip = 2)
table(gtex.meta$SAMPID %in% colnames(gtex.data)) # 51 blood samples
gtex.data[1:4,1:4]
rownames(gtex.data) <- gtex.data$Name
gtex.data <- gtex.data[,colnames(gtex.data) %in% gtex.meta$SAMPID]
max(as.matrix(gtex.data),na.rm = T)  # not log value
write.table(colnames(gtex.data),"/BioII/lulab_b/baopengfei/projects/GTEx/wholeblood_id.txt",quote = F,col.names = F,row.names = F,sep = "\t")

# TCGA
TCGA.meta <- read.table("/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-STAD/TCGA-STAD.GDC_phenotype.tsv",sep = "\t",header = T,check.names = F,fill = T)
TCGA.meta[1:4,1:4]
TCGA.meta$type <- as.character(lapply(str_split(TCGA.meta$submitter_id.samples,"-"),function(x) x[4]))
table(TCGA.meta$type)
TCGA.meta <- TCGA.meta[TCGA.meta$type=="01A",]

tcga.data <- data.table::fread("/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-STAD/TCGA-STAD.htseq_fpkm.tsv",nThread = 10, stringsAsFactors = F, check.names = F, sep = "\t",data.table = F,header = T)
tcga.data[1:3,1:3]
rownames(tcga.data) <- tcga.data$Ensembl_ID
table(colnames(tcga.data) %in% TCGA.meta$submitter_id.samples)  # 189 samples
tcga.data <- tcga.data[,colnames(tcga.data) %in% TCGA.meta$submitter_id.samples]
write.table(colnames(tcga.data),"/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-STAD/01A_id.txt",quote = F,col.names = F,row.names = F,sep = "\t")

tcga.data.tpm <- tcga.data/colSums(tcga.data)*10^6
colSums(tcga.data)[1:3]
tcga.data.tpm[1:3,1:3]  # not exact number ?
2.044264*10^6/30392.97
2.849784*10^6/38187.13

# rm dup rowname: ensg
head(as.character(lapply(str_split(rownames(gtex.data),"\\."),function(x) x[1])))
f1 <- duplicated(as.character(lapply(str_split(rownames(gtex.data),"\\."),function(x) x[1])))
table(f1)
gtex.data <- gtex.data[!f1,]
rownames(gtex.data) <- as.character(lapply(str_split(rownames(gtex.data),"\\."),function(x) x[1]))

f2 <- duplicated(as.character(lapply(str_split(rownames(tcga.data.tpm),"\\."),function(x) x[1])))
table(f2)
tcga.data.tpm <- tcga.data.tpm[!f2,]
rownames(tcga.data.tpm) <- as.character(lapply(str_split(rownames(tcga.data.tpm),"\\."),function(x) x[1]))


table(rownames(gtex.data) %in% rownames(tcga.data.tpm))
table( rownames(tcga.data.tpm) %in% rownames(gtex.data))  # 55369 overlap ensg
gtex.data <- gtex.data[rownames(gtex.data) %in% rownames(tcga.data.tpm),]
tcga.data.tpm <- tcga.data.tpm[rownames(gtex.data),]
write.table(gtex.data,"/BioII/lulab_b/baopengfei/projects/GTEx/wholeblood_TPM.txt",quote = F,col.names = T,row.names = T,sep = "\t")
write.table(tcga.data.tpm,"/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-STAD/01A_TPM.txt",quote = F,col.names = T,row.names = T,sep = "\t")
write.table(cbind(tcga.data.tpm,gtex.data),"/BioII/lulab_b/baopengfei/projects/tcga/data/GDC-STAD/wholeblood-01A_TPM.txt",quote = F,col.names = T,row.names = T,sep = "\t")


#### NDR feature selection ####
## STAD tissue and blood expr
expr <- read.table("/BioII/lulab_b/baopengfei/projects/tcga/output/diff/GTExwholeblood-TCGA-STAD-01A_wilcox.txt",sep = "\t",header = T,check.names = F,fill = T)
expr <- expr[order(-expr$pvalue,expr$log2FoldChange,decreasing = T),]
table(abs(expr$log2FoldChange)>=5)  # 1263 true
table(expr$padj<0.01)
table(expr$padj<0.01 & abs(expr$log2FoldChange)>=5)

stad.tx <- expr[expr$padj<0.05 & expr$log2FoldChange>=5,] # 1227
blood.tx <- expr[expr$padj<0.1 & expr$log2FoldChange<=-2,] # 162, 350


## options ---
exon1 <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/promoter300100exon1end/STADvsNC-passQC-edger_exact/STADvsNC-passQC-edger_exact.txt",sep = "\t",header = T,check.names = F,fill = T)
#exon1 <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/promoter300100exon1end/STADvsNC-passQC-deseq2_wald/STADvsNC-passQC-deseq2_wald.txt",sep = "\t",header = T,check.names = F,fill = T)
#exon1 <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/NucOccuVarLog2RelativeRatio-promoter300100exon1end/STADvsNC-passQC-wilcox/STADvsNC-passQC-wilcox.txt",sep = "\t",header = T,check.names = F,fill = T)


f3 <- as.character(lapply(str_split(rownames(exon1),"\\."),function(x) x[1]))
f3 <- as.character(lapply(str_split(f3,"_"),function(x) x[2]))
head(f3)
table(duplicated(f3))
exon1 <- exon1[!duplicated(f3),]
rownames(exon1) <- f3[!duplicated(f3)]


## plasma expr
table(abs(exon1$log2FoldChange)>=1 & exon1$pvalue<=0.05) # 4459
#exon1.candidate <- exon1[abs(exon1$log2FoldChange)>=1 & exon1$pvalue<=0.01,]  # 1146
exon1.candidate <- exon1
exon1.candidate <- exon1.candidate[order(-exon1.candidate$pvalue,exon1.candidate$log2FoldChange,decreasing = T),]
#table(rownames(exon1.candidate)[exon1.candidate$log2FoldChange<0] %in% rownames(stad.tx))  #  depth.ratio: 21, cpm.wilcox: 15, count.deseq2: 2
exon1.candidate.down <- exon1.candidate[exon1.candidate$log2FoldChange<0,]
#exon1.candidate.down <- exon1.candidate.down[rownames(exon1.candidate.down) %in% rownames(stad.tx),]
#exon1.candidate.down$type <- "STAD.plasma.DNA.exon1.down.AND.STAD.tissue.RNA.expr.up"
exon1.candidate.down$type <- "STAD.plasma.DNA.exon1.down"
exon1.candidate.down$name <- rownames(exon1.candidate.down)

#table(rownames(exon1.candidate)[exon1.candidate$log2FoldChange>0] %in% rownames(blood.tx))  # depth.ratio: 1, cpm.wilcox: 7  count.deseq2: 15
exon1.candidate.up <- exon1.candidate[exon1.candidate$log2FoldChange>0,]  # 0 left
#exon1.candidate.up <- exon1.candidate.up[rownames(exon1.candidate.up) %in% rownames(blood.tx),]
#exon1.candidate.up$type <- "STAD.plasma.DNA.exon1.up.AND.STAD.tissue.RNA.expr.down"
exon1.candidate.up$type <- "STAD.plasma.DNA.exon1.up"
exon1.candidate.up$name <- rownames(exon1.candidate.up)

# all.can <- rbind(exon1.candidate.down,tss.candidate.down,exon1.candidate.up,tss.candidate.up)  # watch out dup rownames ! will add 1... at ensg end
all.can <- rbind(exon1.candidate.down,exon1.candidate.up) 
all.can$baseMean <- 2^all.can$baseMean
#all.can$treatMean <- 2^all.can$treatMean
#all.can$ctrlMean <- 2^all.can$ctrlMean
#all.can$STAD.TPM <- 2^expr[all.can$name,"treatMean"]
#all.can$blood.TPM <- 2^expr[all.can$name,"ctrlMean"]
#all.can <- all.can[,c("name","log2FoldChange","pvalue","padj","baseMean","treatMean","ctrlMean","type","STAD.TPM","blood.TPM")]
#colnames(all.can) <- c("ensg","log2FoldChange","pvalue","padj","MeanRatio","STAD.MeanRatio","NC.MeanRatio","trend","STAD.tissue.RNA.expr.TPM","NC.blood.RNA.expr.TPM")
all.can <- all.can[,c("name","log2FoldChange","pvalue","padj","baseMean","type")]
colnames(all.can) <- c("ensg","log2FoldChange","pvalue","padj","MeanRatio","trend")

colnames(all.can)
all.can$region <- unlist(lapply(strsplit(all.can$trend,"\\."), function(x) x[4]))


## add gene meta info
gene <- read.table("/BioII/lulab_b/baopengfei/shared_reference/geneset/encodev27_gene.txt",sep = "\t",header = T,check.names = F,fill = T)
gene <- gene[,c("ensg","type","name","strand")]
gene <- gene[!duplicated(gene$ensg),]
colnames(gene) <- c("ensg","type","hgnc","strand")

exon1.bed <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/promoter300100exon1end.bed",sep = "\t",header = F,check.names = F,fill = T)
exon1.bed$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_|promoter300exon1end100_","",exon1.bed$V4)
exon1.bed$ensg <- unlist(lapply(strsplit(exon1.bed$ensg,"\\."), function(x) x[1]))
exon1.bed <- exon1.bed[!duplicated(exon1.bed$ensg),]
exon1.bed$pos <- paste0(exon1.bed$V1,":",exon1.bed$V2,"-",exon1.bed$V3)
exon1.bed <- exon1.bed[,c("ensg","pos")]
colnames(exon1.bed)[2] <- "exon1.ndr.pos"

exon1.flank.bed <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/NDR/exon1end10002000.bed",sep = "\t",header = F,check.names = F,fill = T)
exon1.flank.bed$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_|promoter300exon1end100_","",exon1.flank.bed$V4)
exon1.flank.bed$ensg <- unlist(lapply(strsplit(exon1.flank.bed$ensg,"\\."), function(x) x[1]))
exon1.flank.bed <- exon1.flank.bed[!duplicated(exon1.flank.bed$ensg),]
exon1.flank.bed$pos <- paste0(exon1.flank.bed$V1,":",exon1.flank.bed$V2,"-",exon1.flank.bed$V3)
exon1.flank.bed <- exon1.flank.bed[,c("ensg","pos")]
colnames(exon1.flank.bed)[2] <- "exon1.flank.pos"
# c(tss.bed,tss.flank.bed,exon1.bed,exon1.flank.bed)
# c("tss.bed","tss.flank.bed","exon1.bed","exon1.flank.bed")
# for(i in c("tss.bed","tss.flank.bed","exon1.bed","exon1.flank.bed")){
#   gene <- dplyr::left_join(gene,i)
# }
#gene <- dplyr::left_join(gene,tss.bed)
#gene <- dplyr::left_join(gene,tss.flank.bed)
gene <- dplyr::left_join(gene,exon1.bed)
gene <- dplyr::left_join(gene,exon1.flank.bed)
all.can <- dplyr::left_join(all.can,gene)


## get expr mat to plot heatmap
# gse.meta <- read.table("./meta/lulab/sample_table.txt",header = T,sep = "\t",stringsAsFactors = F) # row.names = 1,
# table(gse.meta$health_state)
# table(gse.meta$INPUT)
# gse.meta <- gse.meta[gse.meta$INPUT=='input',]
# gse.stad <- gse.meta$acc[gse.meta$health_state=="STAD"]
# gse.hbv <- gse.meta$acc[gse.meta$health_state=="HBV"]
# gse.nc <- gse.meta$acc[gse.meta$health_state=="healthy"]
stad <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-STAD-passQC.txt",header = F,sep = "\t",stringsAsFactors = F)$V1
nc <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-NC-passQC.txt",header = F,sep = "\t",stringsAsFactors = F)$V1

## options ---
gse.exon1 <- read.table("./output/lulab/matrix/CPM_matrix_promoter300100exon1end.txt",header = T,check.names = F,row.names = 1,sep = "\t") # 
#gse.exon1 <- read.table("./output/lulab/matrix/NOVlog2Ratio_matrix_promoter300100exon1end.txt",header = T,check.names = F,row.names = 1,sep = "\t") # 
gse.exon1 <- gse.exon1[,c(stad,nc)]
#gse.exon1[1:3,]
gse.exon1$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_|promoter300exon1end100_","",rownames(gse.exon1))
gse.exon1$ensg <- unlist(lapply(strsplit(gse.exon1$ensg,"\\."), function(x) x[1]))
gse.exon1$region <- "exon1"
#colnames(gse.exon1)[grep("SRR",colnames(gse.exon1))] <- paste0(colnames(gse.exon1)[grep("SRR",colnames(gse.exon1))],".exon1")

# gse.tss <- read.table("./output/lulab/matrix/NucOccuVarLog2RelativeRatio_matrix_promoter150TSS50.txt",header = T,sep = "\t") # row.names = 1,
# gse.tss <- gse.tss[,c(stad,nc)]
# #gse.tss[1:3,]
# gse.tss$ensg <- gsub("promoter150TSS50_|promoter300100exon1end_|promoter300exon1end100_","",rownames(gse.tss))
# gse.tss$ensg <- unlist(lapply(strsplit(gse.tss$ensg,"\\."), function(x) x[1]))
# gse.tss$region <- "tss"
# colnames(gse.tss)[grep("SRR",colnames(gse.tss))] <- paste0(colnames(gse.tss)[grep("SRR",colnames(gse.tss))],".tss")
# head(gse.tss)
colnames(all.can)
str(all.can)
str(gse.exon1)

all.can <- dplyr::left_join(all.can,gse.exon1)
#all.can <- dplyr::left_join(all.can,gse.tss)


## options ---
#write.table(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/STADvsNC_NOV_exon1_candidates_v2.txt",quote = F,row.names = F,col.names = T,sep = "\t")
#write.table(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/deseq2_wald_STADvsNC_NOV_exon1_candidates_v2.txt",quote = F,row.names = F,col.names = T,sep = "\t")
#write.table(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/depthRatio_wilcox_STADvsNC_NOV_exon1_candidates_v2.txt",quote = F,row.names = F,col.names = T,sep = "\t")

write.table(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/edger_exact_STADvsNC_NOV_exon1_all_v2.txt",quote = F,row.names = F,col.names = T,sep = "\t")
rio::export(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/edger_exact_STADvsNC_NOV_exon1_all_v2.xlsx")
#write.table(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/deseq2_wald_STADvsNC_NOV_exon1_all_v2.txt",quote = F,row.names = F,col.names = T,sep = "\t")
#rio::export(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/deseq2_wald_STADvsNC_NOV_exon1_all_v2.xlsx")
#write.table(all.can,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/diff/depthRatio_wilcox_STADvsNC_NOV_exon1_all_v2.txt",quote = F,row.names = F,col.names = T,sep = "\t")




#summary(expr$baseMean)
#summary(2^expr$treatMean)
#summary(2^expr$ctrlMean)
gse.exon1 <- 2^gse.exon1
gse.exon1[1:4,1:4]
g <- "ENSG00000116641"
stad.dep <- as.numeric(gse.exon1[grep(g,rownames(gse.exon1)),stad])
nc.dep <- as.numeric(gse.exon1[grep(g,rownames(gse.exon1)),nc])

hist(stad.dep,breaks = 93,xlim = c(-20,20))
hist(nc.dep,breaks = 93,xlim = c(-20,20))

summary(stad.dep)
summary(nc.dep)

wilcox.test(stad.dep,nc.dep)


# plot candidates boxplot -------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq")
rm(list = ls())

## read meta
stad <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-STAD-passQC.txt",header = F,sep = "\t",stringsAsFactors = F)$V1
nc <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-NC-passQC.txt",header = F,sep = "\t",stringsAsFactors = F)$V1

data <- read.table("output/lulab/matrix/CPM_matrix_promoter300100exon1end.txt", stringsAsFactors = F, sep = "\t", header=T,check.names = F)
data <- data[,c(nc,stad)]
data <- log(data+1)
#hist(as.matrix(data))

## plot 
# ENSG00000116641, ENSG00000019549, ENSG00000243509, ENSG00000075624.13 (actb)
tmp <- data[grep("ENSG00000075624.13",rownames(data)),]

tmp <- reshape2::melt(tmp)
colnames(tmp) <- c("sample","log2CPM")
tmp$group <- unlist(lapply(strsplit(as.character(tmp$sample),"\\-"),function(x) x[1]))
str(tmp)
tmp$group <- factor(tmp$group,levels = c("NC","STAD"))
b <- runif(nrow(tmp), -0.05, 0.05)

ggplot(tmp,aes(x=group,y=log2CPM,fill=group))+
  geom_boxplot()+
  geom_point(aes(x=as.numeric(group)+b,y=log2CPM),size=4, stroke=0.5,alpha=0.5) +
  #geom_point(shape=21,size=2, stroke=0.5,alpha=0.7) +
  theme_classic() + 
  #ggsci::scale_fill_nejm()+
  scale_fill_manual(values = c("steelblue","firebrick"))+
  ggpubr::stat_compare_means(#ref.group = "nc", 
    
    label="p.signif",bracket.size = 0.3, vjust = 1, hjust=-2,size = 16,color="red", method = "wilcox.test",hide.ns=F) +
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black"), #,face="bold"
        legend.title= element_text(size= 26,color = "black",face="bold")) 
