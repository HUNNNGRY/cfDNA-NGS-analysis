# compare CNV methods (per sample as a dot, 1v1 methods each box)
# last 211012
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")

# methods cor among 4 methods ---------------------------------------------
#methods cor among 4 methods: 
#maskDepth,WisecondorX.CNVkit,WisecondorX.zscore,GC.CPM
cor.m <- "spearman" # pearson


## read samples
crc <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-medip-CRC-passQCstringent.txt",header = F,stringsAsFactors = F)$V1
stad <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-medip-STAD-passQCstringent.txt",header = F,stringsAsFactors = F)$V1
nc <- read.delim("../../multi-omics-explore/meta/lulab/multiomic-medip-NC-passQCstringent.txt",header = F,stringsAsFactors = F)$V1


## read mat
CPM <- read.table("./output/lulab/matrix/CPM_matrix_promoter.txt",header = T,check.names = F,row.names = 1,sep = "\t")
CPM <- log2(CPM+1)
GC.CPM <- read.table("./output/lulab/matrix/CPMCorrectGC_matrix_promoter.txt",header = T,check.names = F,row.names = 1,sep = "\t")
GC.CPM <- log2(GC.CPM+1)
CPM.5K <- read.table("./output/lulab/matrix/CPM_matrix_promoter5k.txt",header = T,check.names = F,row.names = 1,sep = "\t")
CPM.5K <- log2(CPM.5K+1)
CPM.1K <- read.table("./output/lulab/matrix/CPM_matrix_promoter1kbTSS.txt",header = T,check.names = F,row.names = 1,sep = "\t")
CPM.1K <- log2(CPM.1K+1)
CPM.200b <- read.table("./output/lulab/matrix/CPM_matrix_promoter150TSS50.txt",header = T,check.names = F,row.names = 1,sep = "\t")
CPM.200b <- log2(CPM.200b+1)


CPM.200b[1:3,1:3]

trim.name <- function(GC.CPM){
  rownames(GC.CPM) <- as.character(lapply(strsplit(rownames(GC.CPM),"\\|"),function(x) x[1]))
  rownames(GC.CPM) <- sub("promoter300100exon1end_|promoter150TSS50_|promoter_|promoter1kbTSS_","",rownames(GC.CPM))
  return(GC.CPM)
}

CPM <- trim.name(CPM)
GC.CPM <- trim.name(GC.CPM)
CPM.5K <- trim.name(CPM.5K)
CPM.1K <- trim.name(CPM.1K)
CPM.200b <- trim.name(CPM.200b)


g <- rownames(CPM)
s <- c(crc,stad,nc)

CPM <- CPM[g,s]
GC.CPM <- GC.CPM[g,s]
CPM.5K <- CPM.5K[g,s]
CPM.1K <- CPM.1K[g,s]
CPM.200b <- CPM.200b[g,s]

cor.res <- data.frame(row.names = "i",
                      rho.CPM.GC.CPM=0,
                      rho.CPM.CPM.5K=0,
                      rho.CPM.CPM.1K=0,
                      rho.CPM.CPM.200b=0,

                      rho.GC.CPM.CPM.5K=0,
                      rho.GC.CPM.CPM.1K=0,
                      rho.GC.CPM.CPM.200b=0,
                      
                      rho.CPM.5K.CPM.1K=0,
                      rho.CPM.5K.CPM.200b=0,
                      
                      rho.CPM.1K.CPM.200b=0
                      
)

diy.cor <- function(CPM,GC.CPM,cor.m="spearman"){
  rho.GC.CPM.CPM <- cor(CPM[[i]],GC.CPM[[i]],use = "complete.obs", method =cor.m)
  return(rho.GC.CPM.CPM)
}


for(i in colnames(CPM)){ 
  # rho.CPM.GC.CPM=diy.cor(CPM,GC.CPM)
  # rho.CPM.CPM.5K=diy.cor(CPM,CPM.5K)
  # rho.CPM.CPM.1K=diy.cor(CPM,CPM.1K)
  # rho.CPM.CPM.200b=diy.cor(CPM,CPM.200b)
  # rho.GC.CPM.CPM.5K=diy.cor(GC.CPM,CPM.5K)
  # rho.GC.CPM.CPM.1K=diy.cor(GC.CPM,CPM.1K)
  # rho.GC.CPM.CPM.200b=diy.cor(GC.CPM,CPM.200b)
  # rho.CPM.5K.CPM.1K=diy.cor(CPM.5K,CPM.1K)
  # rho.CPM.5K.CPM.200b=diy.cor(CPM.5K,CPM.200b)
  # rho.CPM.1K.CPM.200b=diy.cor(CPM.1K,CPM.200b)
  tmp <- data.frame(row.names = i,
                    rho.CPM.GC.CPM=diy.cor(CPM,GC.CPM),
                    rho.CPM.CPM.5K=diy.cor(CPM,CPM.5K),
                    rho.CPM.CPM.1K=diy.cor(CPM,CPM.1K),
                    rho.CPM.CPM.200b=diy.cor(CPM,CPM.200b),
                    rho.GC.CPM.CPM.5K=diy.cor(GC.CPM,CPM.5K),
                    rho.GC.CPM.CPM.1K=diy.cor(GC.CPM,CPM.1K),
                    rho.GC.CPM.CPM.200b=diy.cor(GC.CPM,CPM.200b),
                    rho.CPM.5K.CPM.1K=diy.cor(CPM.5K,CPM.1K),
                    rho.CPM.5K.CPM.200b=diy.cor(CPM.5K,CPM.200b),
                    rho.CPM.1K.CPM.200b=diy.cor(CPM.1K,CPM.200b)
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
cor.res.m$comparison <- factor(cor.res.m$comparison,levels =c("CPM.GC.CPM","CPM.CPM.5K","CPM.CPM.1K","CPM.CPM.200b",
                                                              "GC.CPM.CPM.5K","GC.CPM.CPM.1K","GC.CPM.CPM.200b",
                                                              "CPM.5K.CPM.1K","CPM.5K.CPM.200b",
                                                              "CPM.1K.CPM.200b"))



ggplot(cor.res.m, aes(x = comparison, y = rho))+ 
  labs(y="correlation",x= NULL,title = "diff pro meth methods spearman cor")+  
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



