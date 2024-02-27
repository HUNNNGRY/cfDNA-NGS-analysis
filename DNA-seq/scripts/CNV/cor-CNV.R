# compare CNV methods (per sample as a dot, 1v1 methods each box)
# last 211007
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

# methods cor among 4 methods ---------------------------------------------
#methods cor among 4 methods: 
#maskDepth,WisecondorX.CNVkit,WisecondorX.zscore,TPM
cor.m <- "spearman" # pearson


## read samples
crc <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-CRC-passQCstringent.txt",header = F)$V1
stad <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-STAD-passQCstringent.txt",header = F)$V1
nc <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-wgs-NC-passQCstringent.txt",header = F)$V1


## read mat
CNVkit <- read.table("./output/lulab/matrix/CNVkitlog2depth_matrix_gene.txt",header = T,check.names = F,row.names = 1,sep = "\t")
WisecondorX <- read.table("./output/lulab/matrix/CNVzscore_matrix_gene.txt_old",header = T,check.names = F,row.names = 1,sep = "\t")

TPM <- read.table("./output/lulab/matrix/TPM_matrix_gene.txt",header = T,check.names = F,row.names = 1,sep = "\t")
TPM <- log2(TPM+1)
CPM <- read.table("./output/lulab/matrix/CPM_matrix_gene.txt",header = T,check.names = F,row.names = 1,sep = "\t")
CPM <- log2(CPM+1)
maskDepth <- read.table("./output/lulab/matrix/CNVmaskDepthRatio_matrix_gene.txt",header = T,check.names = F,row.names = 1,sep = "\t")

CNVkit[1:3,1:3]
CNVkit <- na.omit(CNVkit)
WisecondorX <- na.omit(WisecondorX)
CNVkit[1:3,1:3]
WisecondorX[1:3,1:3]

rownames(TPM) <- as.character(lapply(str_split(rownames(TPM),"\\|"),function(x) x[1]))
rownames(TPM) <- sub("promoter300100exon1end_|promoter150TSS50_","",rownames(TPM))
rownames(CPM) <- as.character(lapply(str_split(rownames(CPM),"\\|"),function(x) x[1]))
rownames(CPM) <- sub("promoter300100exon1end_|promoter150TSS50_","",rownames(CPM))
CPM[1:3,1:3]

# table(rownames(maskDepth) %in% rownames(CNVkit))
# all(rownames(WisecondorX)==rownames(CNVkit))
# all(colnames(maskDepth)==colnames(CNVkit))
TPM <- TPM[rownames(WisecondorX),c(crc,stad,nc)]
CPM <- CPM[rownames(WisecondorX),c(crc,stad,nc)]
maskDepth <- maskDepth[rownames(WisecondorX),c(crc,stad,nc)]
CNVkit <- CNVkit[rownames(WisecondorX),c(crc,stad,nc)]
WisecondorX <- WisecondorX[,c(crc,stad,nc)]

cor.res <- data.frame(row.names = "i",
                      rho.TPM.CPM="rho.TPM.CPM",
                      rho.TPM.maskDepth="rho.TPM.maskDepth",
                      rho.TPM.CNVkit="rho.TPM.CNVkit",
                      rho.TPM.WisecondorX="rho.TPM.WisecondorX",   
                      rho.CPM.maskDepth="rho.CPM.maskDepth",
                      rho.CPM.CNVkit="rho.CPM.CNVkit",
                      rho.CPM.WisecondorX="rho.CPM.WisecondorX",   
                      rho.maskDepth.CNVkit="rho.maskDepth.CNVkit",
                      rho.maskDepth.WisecondorX="rho.maskDepth.WisecondorX",
                      rho.CNVkit.WisecondorX="rho.CNVkit.WisecondorX")

diy.cor <- function(CPM,TPM,cor.m){
  rho.TPM.CPM <- cor(CPM[[i]],TPM[[i]],use = "complete.obs", method =cor.m)
  return(rho.TPM.CPM)
}
for(i in colnames(maskDepth)){ 
  # rho.maskDepth.TPM <- cor(maskDepth[[i]],TPM[[i]],use = "complete.obs", method =cor.m)  # moderate correalted !
  # rho.maskDepth.WisecondorX <- cor(maskDepth[[i]],WisecondorX[[i]],use = "complete.obs", method =cor.m) # poorly correalted !
  # rho.maskDepth.CNVkit <- cor(maskDepth[[i]],CNVkit[[i]],use = "complete.obs", method =cor.m) # poorly correalted !
  # rho.WisecondorX.CNVkit <- cor(CNVkit[[i]],WisecondorX[[i]],use = "complete.obs", method =cor.m)
  # rho.TPM.CNVkit <- cor(TPM[[i]],CNVkit[[i]],use = "complete.obs", method =cor.m)
  # rho.TPM.WisecondorX <- cor(TPM[[i]],WisecondorX[[i]],use = "complete.obs", method =cor.m)
  
  rho.TPM.CPM=diy.cor(TPM,CPM,cor.m)
  rho.TPM.maskDepth=diy.cor(TPM,maskDepth,cor.m)
  rho.TPM.CNVkit=diy.cor(TPM,CNVkit,cor.m)
  rho.TPM.WisecondorX=diy.cor(TPM,WisecondorX,cor.m)
  rho.CPM.maskDepth=diy.cor(CPM,maskDepth,cor.m)
  rho.CPM.CNVkit=diy.cor(CPM,CNVkit,cor.m)
  rho.CPM.WisecondorX=diy.cor(CPM,WisecondorX,cor.m)
  rho.maskDepth.CNVkit=diy.cor(maskDepth,CNVkit,cor.m)
  rho.maskDepth.WisecondorX=diy.cor(maskDepth,WisecondorX,cor.m)
  rho.CNVkit.WisecondorX=diy.cor(CNVkit,WisecondorX,cor.m)
  
  tmp <- data.frame(row.names = i,
                    rho.TPM.CPM=rho.TPM.CPM,
                    rho.TPM.maskDepth=rho.TPM.maskDepth,
                    rho.TPM.CNVkit=rho.TPM.CNVkit,
                    rho.TPM.WisecondorX=rho.TPM.WisecondorX,   
                    rho.CPM.maskDepth=rho.CPM.maskDepth,
                    rho.CPM.CNVkit=rho.CPM.CNVkit,
                    rho.CPM.WisecondorX=rho.CPM.WisecondorX,   
                    rho.maskDepth.CNVkit=rho.maskDepth.CNVkit,
                    rho.maskDepth.WisecondorX=rho.maskDepth.WisecondorX,
                    rho.CNVkit.WisecondorX=rho.CNVkit.WisecondorX
  )
  cor.res <- rbind(cor.res,tmp)
}


cor.res <- cor.res[c(crc,stad,nc),]
cor.res$sample <- rownames(cor.res)
cor.res.m <- reshape2::melt(cor.res,id.var="sample",value.name = "rho")
cor.res.m$variable <- sub("rho.","",cor.res.m$variable)
colnames(cor.res.m)[2] <- "comparison"
cor.res.m <- as_tibble(cor.res.m)
head(cor.res.m)
cor.res.m$rho <- as.numeric(cor.res.m$rho)
library(ggplot2)
cor.res.m$comparison <- factor(cor.res.m$comparison,levels =c("TPM.CPM","TPM.maskDepth","TPM.CNVkit","TPM.WisecondorX","CPM.maskDepth","CPM.CNVkit","CPM.WisecondorX","maskDepth.CNVkit","maskDepth.WisecondorX","CNVkit.WisecondorX"))


ggplot(cor.res.m, aes(x = comparison, y = rho))+ 
  labs(y="correlation",x= NULL,title = "diff CNV methods cor")+  
  geom_boxplot(position=position_dodge(0.5),color="black",fill="grey",width=0.5,size=0.4,
               outlier.alpha = 1, outlier.size = 0.5)+ 
  theme_bw() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",face="bold",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12))



