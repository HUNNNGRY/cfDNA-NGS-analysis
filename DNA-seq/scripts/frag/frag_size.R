# fragment size EDA (test)
# last 210721 by pengfei

setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)


# option1: depth ----------------------------------------------------------
s <- read.table("./output/lulab/test/shortFrag_matrix_gene.txt",header = T,row.names = 1, check.names = F) #.not tsv
#s[1:4,1:4]
s.ratio <- s
for(i in 1:ncol(s)){
  s.ratio[i] <- s.ratio[i]/colMeans(s)[i]
}
#colMeans(s)[1:4]
#ratio[1:4,1:4]

l <- read.table("./output/lulab/test/longFrag_matrix_gene.txt",header = T,row.names = 1, check.names = F)
#l[1:4,1:4]

#ratio[1:4,1:4]
l.ratio <- l
for(i in 1:ncol(l)){
  l.ratio[i] <- l.ratio[i]/colMeans(l)[i]
}

# option2: TPM ----------------------------------------------------------
s.ratio <- read.table("./output/lulab/test/TPM_matrix_gene_short.txt",header = T,row.names = 1, check.names = F) #.not tsv
l.ratio <- read.table("./output/lulab/test/TPM_matrix_gene_long.txt",header = T,row.names = 1, check.names = F) #.not tsv
s.ratio <- log2(as.matrix(s.ratio)+0.01)
l.ratio <- log2(as.matrix(l.ratio)+0.01)

# option2: CPM.TMM (final chosen in 1st version multiomics paper)  ----------------------------------------------------------
s.ratio <- read.table("./output/lulab/test/CPM-TMM_matrix_gene_short.txt",header = T,row.names = 1, check.names = F) #.not tsv
l.ratio <- read.table("./output/lulab/test/CPM-TMM_matrix_gene_long.txt",header = T,row.names = 1, check.names = F) #.not tsv
s.ratio <- log2(as.matrix(s.ratio)+0.01)
l.ratio <- log2(as.matrix(l.ratio)+0.01)


# plot ----------------------------------------------------------
ratio <- s.ratio/(l.ratio+0.01)
hist(s.ratio$`CRC-PKU-10-wgs`,breaks = 1000)
hist(l.ratio$`CRC-PKU-10-wgs`,breaks = 1000)
hist(ratio$`CRC-PKU-10-wgs`,breaks = 1000)
table(is.na(ratio))

## write out mat
ratio.out <- cbind(rownames(ratio),ratio)
colnames(ratio.out)[1] <- "gene_id"
write.table(ratio.out,"./output/lulab/matrix/ShortLongFragRatio_matrix_gene.txt",col.names = T,sep = "\t",quote = F)

#summary(ratio)
ratio <- as.matrix(ratio)
std <- apply(ratio,2,sd)

df <- data.frame(stderr=std,group=colnames(ratio))
#all(rownames(df) ==colnames(ratio) )
df$group <- as.character(lapply(strsplit(df$group,"-"), function(x) x[1]))
df$group <- factor(df$group,levels = c("NC","CRC","STAD"))

## FILTTER sample id
sampl <- read.table("./meta/lulab/sample_table.txt",sep = "\t",check.names = F,header = T)
sampl <- sampl$data_id[sampl$wgs_passQCstringent=="Y"]
rownames(df) <- sub("-short|-long","",rownames(df))

df <- df[rownames(df) %in% sampl,]

## plot
library(ggplot2)
library(ggpubr)
library(extrafont)
p <- ggplot(df,aes(x=group,y=stderr,fill=group))+
  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  geom_point(size = 3, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1),color=alpha("black",0.3))+
  #scale_fill_manual(values=alpha(c("steelblue2","orange2","firebrick2"), 0.7)) +
  ggsci::scale_fill_aaas()+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="none",  # "right"
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black", size=24),
    legend.text= element_text(face="bold", color="black", size=24),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))
#?stat_compare_means
#m  <- mean(subset(gene,CancerType==positive[1])$counts)-mean(subset(gene,CancerType==negative[n])$counts)
p <- p+stat_compare_means(ref.group = "NC",
                          method = "wilcox.test",
                          #method.args = list(alternative = "less"),
                          label = "p.format",size = 12, vjust=0.6 # hide.ns=T,
)+labs(x="",y="fragment size stderr", face="bold",fill="Type",title="wilcox-two-tail")
p 



#### lulab Mt boxplot ####
mt.wgs <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_Mt.txt.summary",sep = "\t",row.names = 1,header = T)
mt.wgs <- as.data.frame(t(mt.wgs))
#head(mt.wgs)
mt.wgs$ratio <- mt.wgs$Assigned/(mt.wgs$Assigned+mt.wgs$Unassigned_NoFeatures)*100
#summary(mt.wgs$ratio)
mt.wgs$group <- as.character(lapply(strsplit(rownames(mt.wgs),"\\."), function(x) x[6]))  # 4 For medip, 6 For WGS
mt.wgs$group <- factor(mt.wgs$group,levels = c("NC","CRC","STAD"))

## plot
library(ggplot2)
library(ggpubr)
library(extrafont)
require(scales) # to access break formatting functions

p <- ggplot(mt.wgs,aes(x=group,y=ratio,fill=group))+
  geom_boxplot(alpha = 1, size = 1, position = position_dodge(1.1),outlier.size=0,outlier.alpha = 0)+
  geom_point(size = 3, position = position_jitterdodge(dodge.width=1.1,jitter.width = 1),color=alpha("black",0.3))+
  scale_fill_manual(values=alpha(c("steelblue2","orange2","firebrick2"), 0.7)) +
  #ylim(0,10) + 
  scale_y_log10(breaks = c(0.0000000001,0.01,0.1,1,10), labels = paste0(c(0.00,0.01,0.10,1.00,10.00),"%"))+ 
  #scale_y_continuous(trans = log10_trans(), breaks = c(0.0000000001,0.01,0.1,1,10), labels = paste0(c(0.00,0.01,0.10,1.00,10.00),"%")) + 
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +  # labels = trans_format("log10", math_format(10^.x))
  theme_bw()+
  theme(#legend.position="right",
    legend.position="none",  # "right"
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(size=1, colour = "black"),
    legend.title = element_text(face="bold", color="black", size=24),
    legend.text= element_text(face="bold", color="black", size=24),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=24,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))
#?stat_compare_means
#m  <- mean(subset(gene,CancerType==positive[1])$counts)-mean(subset(gene,CancerType==negative[n])$counts)
p <- p+stat_compare_means(ref.group = "NC",
                          method = "wilcox.test",
                          #method.args = list(alternative = "less"),
                          label = "p.signif",size = 12, vjust=0.6 # hide.ns=T,
)+labs(x="",y="Mt reads ratio (%)", face="bold",fill="Type",title="wilcox-two-tail")
p 





# public 5hmc data --------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")
options(stringsAsFactors = F)
#dst <- "GSE112679" # "GSE81314" "GSE89570"  "GSE112679"
for (dst in c("GSE81314","GSE89570","GSE112679")){
region <- "bin100kb" # "bin100kb", "gene"

s.ratio <- read.table(paste0("./output/",dst,"/matrix/CPM_matrix_",region,"_short.correctGC.txt"),header = T,row.names = 1, check.names = F) #.not tsv
l.ratio <- read.table(paste0("./output/",dst,"/matrix/CPM_matrix_",region,"_long.correctGC.txt"),header = T,row.names = 1, check.names = F) #.not tsv
#s.ratio <- log2(as.matrix(s.ratio)+0.01)
#l.ratio <- log2(as.matrix(l.ratio)+1)

ratio <- s.ratio/(l.ratio+1)
all(gsub("-short|-long","",colnames(s.ratio))==gsub("-short|-long","",colnames(l.ratio)))  # should be true
colnames(ratio) <- gsub("-short|-long","",colnames(s.ratio))
ratio[1:3,1:3]

#hist(s.ratio$`SRR3486314-short`,breaks = 1000)
#hist(l.ratio$`SRR3486314-long`,breaks = 1000)
#hist(ratio$`SRR3486314`,breaks = 1000)
table(is.na(ratio))

## write out mat
ratio.out <- cbind(rownames(ratio),ratio)
colnames(ratio.out)[1] <- "gene_id"
write.table(ratio.out,paste0("./output/",dst,"/matrix/DNA-FragRatio_matrix_",region,".correctGC.txt"),row.names = F,col.names = T,sep = "\t",quote = F)
ratio.out[1:3,1:3]
#hist(as.vector(as.matrix(ratio.out[,2:2000])))
}




# lulab wgs data --------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)
dst <- "lulab" 

## TMM, gc adjusted
for (region in c("bin100kb","gene")){
  #region <- "bin100kb" # "bin100kb", "gene"
  
  #options ---
  #s.ratio <- read.table(paste0("./output/",dst,"/matrix/CPM_matrix_",region,"_short.txt"),header = T,row.names = 1, check.names = F) #.not tsv
  #l.ratio <- read.table(paste0("./output/",dst,"/matrix/CPM_matrix_",region,"_long.txt"),header = T,row.names = 1, check.names = F) #.not tsv
  s.ratio <- read.table(paste0("./output/",dst,"/matrix/CPM-TMM_matrix_",region,"_short.correctGC.txt"),header = T,row.names = 1, check.names = F) #.not tsv
  l.ratio <- read.table(paste0("./output/",dst,"/matrix/CPM-TMM_matrix_",region,"_long.correctGC.txt"),header = T,row.names = 1, check.names = F) #.not tsv
  
  
  
  #s.ratio <- log2(as.matrix(s.ratio)+0.01)
  #l.ratio <- log2(as.matrix(l.ratio)+1)
  
  ratio <- s.ratio/(l.ratio+1)
  all(gsub("-short|-long","",colnames(s.ratio))==gsub("-short|-long","",colnames(l.ratio)))  # should be true
  colnames(ratio) <- gsub("-short|-long","",colnames(s.ratio))
  ratio[1:3,1:3]
  
  #hist(s.ratio$`SRR3486314-short`,breaks = 1000)
  #hist(l.ratio$`SRR3486314-long`,breaks = 1000)
  #hist(ratio$`SRR3486314`,breaks = 1000)
  table(is.na(ratio))
  
  ## write out mat
  ratio.out <- cbind(rownames(ratio),ratio)
  colnames(ratio.out)[1] <- "gene_id"
  #write.table(ratio.out,paste0("./output/",dst,"/matrix/DNA-FragRatio_matrix_",region,".txt"),row.names = F,col.names = T,sep = "\t",quote = F)
  write.table(ratio.out,paste0("./output/",dst,"/matrix/DNA-FragRatio_matrix_",region,".correctGC.txt"),row.names = F,col.names = T,sep = "\t",quote = F)
  
  ratio.out[1:3,1:3]
  #hist(as.vector(as.matrix(ratio.out[,2:2000])))
}
