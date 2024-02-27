# last 210709 by pengfei
getwd()
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

#### seq data and 300 bp saturation ####
getwd()
l <- read.table("../exOmics/DNA-seq/meta/lulab/sample_library.txt",sep = "\t",header = T)
l <- l[l$Yield.Gbases.!="" & !is.na(l$Yield.Gbases.),]
l <- l[,c("data_id","Yield.Gbases.")]

l2 <- read.table("../exOmics/DNA-seq/meta/lulab/sample_table.txt",sep = "\t",header = T)
l2 <- l2[l2$wgs_passQCstringent=="Y",]
l2 <- l2[,c("data_id","estimated.max.saturation")]

l <- dplyr::left_join(l,l2)
#l <- l[l$Yield.Gbases.<20,]
plot(x=l$Yield.Gbases.,y=l$estimated.max.saturation)
cor(l$Yield.Gbases.,l$estimated.max.saturation)

scatter.smooth(l$Yield.Gbases.,l$estimated.max.saturatio)
library(ggplot2)
ggplot(l,aes(x=l$Yield.Gbases.,y=l$estimated.max.saturation))+
  geom_point() +
  geom_smooth()+
  geom_hline(yintercept = c(0.9),linetype=c("dashed")) +
  geom_vline(xintercept = c(10),linetype=c("dashed")) + 
  labs(title = "Saturation-seqdata",x="seq data (G)",y="300bp saturation") +
  # ylim(c(0.5,1))+
  theme_classic() + 
  theme(plot.title = element_text(size = 18,color="black",hjust = 0.5,face="bold",family="arial"),
        axis.title = element_text(size = 18,color ="black",face="bold",family="arial"), 
        axis.text = element_text(size= 18,color = "black",face="bold",family="arial"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12))


#### duplication boxplot ####
d <- read.table("./meta/lulab/sample_table.txt",stringsAsFactors = F,check.names = F, header = T, sep = "\t")
d[1:3,1:3]
colnames(d)
d <- d[,c("data_id","group","nonDup-ratio")]
str(d)
d$group <- factor(d$group,levels = c("NC","CRC","STAD"))
d$duplication.level <- 1-d$`nonDup-ratio`
f.plot <- d
set.seed(1)
b <- runif(nrow(f.plot), -0.2, 0.2)
ggplot(f.plot, aes(x = group, y = duplication.level, fill=group))+   # , fill = Group
  #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
  geom_boxplot(outlier.alpha = 0, outlier.size = 0)  + 
  geom_point(aes(x=as.numeric(group)+b,y=duplication.level),color="black",size=2,stroke=0.5,alpha=0.5) + 
  ggsci::scale_fill_nejm()+
  #scale_color_manual(values = c("passQC"="grey30","faillQC"="firebrick")) + 
  xlab("group") + ylab("Duplication level") +
  geom_vline(xintercept = 0,color="grey50",linetype="dashed") +
  #scale_y_discrete(breaks=QC.factor,labels=QC.lab) +
  guides(color=guide_legend(title="Duplication level")) +
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))


table(d.m$variable)
d.m <- reshape2::melt(d, id.var=c("data_id","group"))
d.m$variable <- factor(d.m$variable,levels = c("duplication.level","nonDup-ratio"))
ggplot(d.m, aes(x = data_id, y = value, fill=variable))+   # , fill = Group
  #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
  geom_bar(stat = "identity",position = "stack",color="black",)  + 
  #geom_point(aes(x=as.numeric(group)+b,y=duplication.level),color="black",size=2,stroke=0.5,alpha=0.5) + 
  #ggsci::scale_fill_nejm()+
  #scale_color_manual(values = c("passQC"="grey30","faillQC"="firebrick")) + 
  xlab("group") + ylab("Duplication level") +
  geom_vline(xintercept = 0,color="grey50",linetype="dashed") +
  #scale_y_discrete(breaks=QC.factor,labels=QC.lab) +
  guides(color=guide_legend(title="Duplication level")) +
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))



#### quality control ####
f <- read.table("./meta/lulab/qc_stats.txt",stringsAsFactors = F,check.names = F, header = T, sep = "\t")
colnames(f) <- gsub("-|_","\\.",colnames(f))
colnames(f)

qc.col <- c("uniq.reads","genome.depth","coverage.ratio","mapped.ratio","estimated.max.saturation.3k","enrichment.score.relH")
qc.col.criteria <- paste0(qc.col,".","criteria")
f.dat <- f[,qc.col]
f.criteria <- f[,qc.col.criteria]
f.plot <- log10(f.dat/f.criteria)

f.plot <- tidyr::gather(f.plot,key = "QC", value = "log10ratio", 1:ncol(f.plot))
f.plot$id <- f$data.id
#unique(f.plot$QC)
QC.factor <- qc.col
QC.lab <- c("Criteria 6: Uniq reads",
            "Criteria 5: Genome.depth",
            "Criteria 4: Coverage ratio",
            "Criteria 3: Mapped ratio",
            "Criteria 2: Saturation",
            "Criteria 1: Enrich relH score")
f.plot$QC <- factor(f.plot$QC,levels = QC.factor)
f.plot$point.col <- "passQC"
f.plot$point.col[f.plot$log10ratio<(0)] <- "faillQC" 
f.plot$point.col[f.plot$QC=="enrichment.score.relH"] <- "faillQC"
f.plot$point.col[f.plot$log10ratio<(0) & f.plot$QC=="enrichment.score.relH"] <- "passQC"

#table(f.plot$point.col)

set.seed(1)
b <- runif(nrow(f.plot), -0.2, 0.2)
ggplot(f.plot, aes(x = log10ratio, y = QC))+   # , fill = Group
  #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
  geom_boxplot(outlier.alpha = 1, outlier.size = 0)  + 
  geom_point(aes(x=log10ratio,y=as.numeric(QC)+b,color=point.col),size=2,stroke=0.5,alpha=0.5) + 
  scale_color_manual(values = c("passQC"="grey30","faillQC"="firebrick")) + 
  xlab("Log10ratio((Read pairs or Ratio)/Criteria))") + ylab("") +
  geom_vline(xintercept = 0,color="grey50",linetype="dashed") +
  scale_y_discrete(breaks=QC.factor,labels=QC.lab) +
  guides(color=guide_legend(title="Quality control")) +
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))


#### plot cfDNA fragment size ####
t <- read.table("output/lulab/fragment-length/histogram.txt",header = T, sep = "\t",stringsAsFactors = F,check.names = F)  # fill = T,
#colnames(t)
f <- read.table("./meta/lulab/sample_table.txt",stringsAsFactors = F,header = T, sep = "\t")
#f <- f[f$overall_QC=="Y" & f$group!="STAD",]
f <- f[f$wgs_passQC=="Y",]


table(f$group)
t <- t[,c("insertion-size",f$data_id)]
all(colnames(t)[2:ncol(t)]==f$data_id)

t.freq <- t[,2:ncol(t)]
s <- colSums(t.freq,na.rm = T)

#TEST: t.freq <- matrix(c(1:20),nrow = 4,ncol = 5)
t.freq <- sweep(t.freq, 2, colSums(t.freq,na.rm = T), FUN="/")   # scale(t.freq, center=FALSE, scale=colSums(t.freq))

t <- cbind(t[1:401,1],t.freq[1:401,])  # colnames(t) %in% f$data_id
colnames(t)[1] <- "length"
t$CRC.aver <- apply(t[,grep("CRC",colnames(t))],1,mean) 
t$STAD.aver <- apply(t[,grep("STAD",colnames(t))],1,mean) 

t$NC.aver <- apply(t[,grep("NC",colnames(t))],1,mean) 
library(reshape2)
m <- melt(t,id.vars = 1,variable.name = "id")
#m$col <- "grey"
#m$col[grep("CRC.aver",m$id)] <- "firebrick"
#m$col[grep("NC.aver",m$id)] <- "steelblue"

#rownames(t) <- t$insertion.size
#unique(m$id)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
unique(m$id)
#str(m)
m$Group <- unlist(lapply(strsplit(as.character(m$id),"-",fixed = T),function(x) x[1]))
#table(m$Group)
m1 <- m[m$id!="CRC.aver" & m$id!="NC.aver" & m$id!="STAD.aver",]
m2 <- m[m$id=="CRC.aver" | m$id=="NC.aver" | m$id=="STAD.aver",]
m2$Group <- factor(m2$Group,levels = c("NC.aver","CRC.aver","STAD.aver"))

p0 <- ggplot() +  # ,fill=id,alpha=0,color=col
  scale_color_manual(values=alpha(c("steelblue2","#6B4E04","orange2"), 1),labels=c("HD","CRC","STAD"))  +  # "steelblue2","firebrick","orange2"
  geom_line(stat = "summary",data = m1, aes(x=length,y=value,group=id),color=alpha("grey",0.5),size = 0.1) + # 
  geom_line(stat = "summary",data = m2, aes(x=length,y=value,group=id,color=id),size = 1) + # 

  geom_vline(xintercept = c(167,334),color="grey30", linetype="dashed", size=0.8) +
  #geom_text(x=0.150, y=0.28, label="167 bp",size=10) +
  annotate(geom="text", x=120, y=0.03, label="167 bp",color="grey30",size=8)+
  xlim(c(0,400)) +
  ylim(c(0,0.04)) + 
  ylab("Frequency") +
  xlab("Fragment size (bp)") +
  theme_classic() +
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_text(size=18,face="bold")) +
  guides(color=guide_legend("Group"))
pdf("./DNA-fragLen.pdf",width = 11,height = 8)
p0
dev.off()

## CRCvsNC
m1 <- m[m$id!="CRC.aver" & m$id!="NC.aver" & m$id!="STAD.aver" & m$Group!="STAD",]
m2 <- m[(m$id=="CRC.aver" | m$id=="NC.aver") & m$Group!="STAD",]
m2$Group <- factor(m2$Group,levels = c("NC.aver","CRC.aver"))

p0 <- ggplot() +  # ,fill=id,alpha=0,color=col
  scale_color_manual(values=alpha(c("steelblue2","#6B4E04"), 1),labels=c("HD","CRC"))  +  # "steelblue2","firebrick","orange2"
  geom_line(stat = "summary",data = m1, aes(x=length,y=value,group=id),color=alpha("grey",0.5),size = 0.1) + # 
  geom_line(stat = "summary",data = m2, aes(x=length,y=value,group=id,color=id),size = 1) + # 
  
  geom_vline(xintercept = c(167,334),color="grey30", linetype="dashed", size=0.8) +
  #geom_text(x=0.150, y=0.28, label="167 bp",size=10) +
  annotate(geom="text", x=120, y=0.03, label="167 bp",color="grey30",size=8)+
  xlim(c(0,400)) +
  ylim(c(0,0.04)) + 
  ylab("Frequency") +
  xlab("Fragment size (bp)") +
  theme_classic() +
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_text(size=18,face="bold")) +
  guides(color=guide_legend("Group"))
pdf("./DNA-fragLen-CRCvsNC.pdf",width = 11,height = 8)
p0
dev.off()



#### plot cfDNA conc. ####
conc <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/meta/lulab/sample_table.txt",stringsAsFactors = F,check.names = F,sep = "\t",header = T) 
table(conc$plasma=="-")
table(conc$DNA_concentration)

conc <- conc[!conc$plasma=="-" & !conc$plasma=="" & !conc$DNA_concentration=="-" & !conc$DNA_concentration=="",]
conc$plasma <- as.numeric(conc$plasma)
conc$DNA_concentration <- as.numeric(conc$DNA_concentration)

conc$conc <- conc$DNA_concentration/conc$plasma
#conc <- conc[conc$group!="STAD",]
colnames(conc)
age <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/meta/lulab/age.txt",stringsAsFactors = F,check.names = F,sep = "\t",header = T) 
conc <- merge(conc,age,by.x="data_id",by.y="data_ID")
conc$group.x[conc$group.x=="CRC"] <- "CRC" 
conc$group.x[conc$group.x=="STAD"] <- "STAD" 
conc$group.x[conc$group.x=="NC"] <- "HD" 
conc$group.x <- factor(conc$group.x,levels = c("HD","CRC","STAD"))
library(ggplot2)
library(ggExtra)
#conc <- conc[conc$group.x!="STAD",]

p <- ggplot(conc, aes(x=age, y=conc, color=factor(group.x))) +  # , size=cyl
  geom_point(size=4) +
  #scale_fill_discrete(RColorBrewer::brewer.pal(3, "Set2")[1:3]) +  
  #scale_fill_discrete(c("orange2","steelblue2")) +  
  #scale_color_manual(values=alpha(c("orange2","steelblue2","firebrick2"), 1))  + 
  scale_color_manual(values=alpha(c("steelblue2","#6B4E04","orange2"), 1),labels=c("HD","CRC","STAD"))  +  # "1" = "red", "2" = "grey", "3" = "blue"
  
  #theme(legend.position="none")
  theme_classic() +
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_text(size=18,face="bold")) +
  ylab("cfDNA concentrition (ng/ml)") +
  guides(color=guide_legend("Group"))

p1 <- ggMarginal(p, type = "density",color="black",size=4, groupColour=F,groupFill=T,alpha=0.3)  #

pdf("./DNA-concentration.pdf",width = 10,height = 8)
p1
dev.off()

# ggplot(conc, aes(x=age, y=conc,fill=factor(group.x))) +  # , size=cyl
#      geom_point( size=6, shape=21, stroke = 1) +
#      #scale_fill_discrete(RColorBrewer::brewer.pal(3, "Set2")[1:3]) +  
#      #scale_fill_discrete(c("orange2","steelblue2")) +  
#      scale_fill_manual(values=alpha(c("steelblue2","#6B4E04","orange2"), 1))  +  #"firebrick2"
#      #scale_color_manual(values = c("black")) +
#      #theme(legend.position="none")
#      theme_classic() +
#      theme(axis.text=element_text(size=16,face="bold"),
#                      axis.title=element_text(size=18,face="bold"),
#                      legend.text=element_text(size=20,face="bold"),
#                      legend.title=element_text(size=18,face="bold")) +
#      ylab("cfDNA concentrition (ng/ml)") 




#### plot coverage hist ####
dataset <- "lulab"   # 
region <- "gene"  # gene or promoter
datatype <- "wgs"  # medip or wgs
#outFile <- paste0("output/",datatype,"/",dataset,"/ssGSEA-MtPathway-",region,".csv")
TPM <- read.table(paste0("data/",datatype,"/",dataset,"/matrix/TPM_matrix_",region,".txt"),stringsAsFactors = F,header = T,row.names = 1, check.names = F)
r <- as.character(lapply(strsplit(rownames(TPM),"\\."), function(x) x[1]))
#table(duplicated(r))
TPM <- TPM[!duplicated(r),]
rownames(TPM) <- as.character(lapply(strsplit(rownames(TPM),"\\."), function(x) x[1]))

mt <- read.table("/BioII/lulab_b/baopengfei/projects/methylation/ref/Mt.txt",sep = "\t",stringsAsFactors = F,header = T,row.names=NULL)
mt <- unique(mt$ensembl_gene_id)
rRNA <- read.table("/BioII/lulab_b/baopengfei/projects/methylation/ref/rRNA.txt",sep = "\t",header = T)
rRNA <- unique(rRNA$ensembl_gene_id)
ribo <- read.table("/BioII/lulab_b/baopengfei/projects/methylation/ref/ensembl_entrez_pathway.txt",sep = "\t",stringsAsFactors = F,header = T,row.names=NULL)
Ribosome <- ribo[ribo$DESCRPTION=="Ribosome",]
Ribosome.biogenesis <-  ribo[ribo$DESCRPTION=="Ribosome biogenesis in eukaryotes",]
ribo <- unique(c(Ribosome$ensembl_gene_id,Ribosome.biogenesis$ensembl_gene_id))
  
TPM.1 <- TPM[rownames(TPM) %in% ribo,]
TPM.1 <- TPM.1/1000000
hist(as.matrix(TPM.1),breaks = 200, xlim = c(0,0.004))


#### plot regions assign ratio ####
regions <- c("promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE",'retroposon')
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
# setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")

tmp <- read.table("./output/lulab/matrix/count_matrix_gene.txt.summary",sep = "\t",check.names = F,header = T,stringsAsFactors = F,row.names = 1)
tmp <- as.data.frame(t(tmp))
tmp$sample <- unlist(lapply(strsplit(rownames(tmp),"/"),function(x) x[4]))
tmp$sample <- gsub(".bam","",tmp$sample,fixed = T)
tmp$region <- "gene"

#for(i in c("promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE",'retroposon')){
for(i in regions){
  tmp1 <- read.table(paste0("./output/lulab/matrix/count_matrix_",i,".txt.summary"),sep = "\t",check.names = F,header = T,stringsAsFactors = F,row.names = 1)
  tmp1 <- as.data.frame(t(tmp1))
  tmp1$sample <- unlist(lapply(strsplit(rownames(tmp1),"/"),function(x) x[4]))
  tmp1$sample <- gsub(".bam","",tmp1$sample,fixed = T)
  tmp1$region <- i
  tmp <- rbind(tmp,tmp1)
}
tmp$assign.ratio <- tmp$Assigned/(tmp$Assigned+tmp$Unassigned_NoFeatures)
tmp$group <- unlist(lapply(strsplit(tmp$sample,"-|_|\\."),function(x) x[1]))
tmp <- tmp[,c("assign.ratio","region","sample","group")]
read.ratio <- tmp
read.ratio$region <- factor(read.ratio$region,levels=c("gene",regions))

library(ggpubr)
library(ggplot2)
library(RColorBrewer)
ggbarplot(read.ratio, x = "region", y = "assign.ratio", fill = "region",
          add = "mean_sd", palette = colorRampPalette( brewer.pal(8, "Dark2"))(12),
          position = position_dodge(0.7)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  ylab("reads count ratio") + 
  #scale_x_discrete(breaks=regions,labels=regions.simplify) +
  facet_wrap(nrow = 1,ncol = 3,facets = group~.,strip.position = "top") +
  #scale_fill_brewer(palette = colorRampPalette(brewer.pal(9, "Blues"))(10)) +
  #ggsci::scale_fill_jama(palette=color) +
  coord_flip() +
  #theme_classic( ) +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) 

write.table(read.ratio,"./output/lulab/matrix/count_matrix_allregions.sum",sep = '\t',row.names = F,col.names = T,quote = F)

tmp <- read.table("./output/lulab/matrix/count_matrix_promoter300100exon1end.txt.summary",sep = "\t",check.names = F,header = T,stringsAsFactors = F,row.names = 1)
tmp <- as.data.frame(t(tmp))
tmp$sample <- unlist(lapply(strsplit(rownames(tmp),"/"),function(x) x[4]))
tmp$sample <- gsub(".bam","",tmp$sample,fixed = T)
tmp$region <- "promoter300100exon1end"
tmp$sum <- (tmp$Assigned+tmp$Unassigned_NoFeatures)
#hist(tmp$Assigned[read.ratio$region=="promoter300100exon1end"])
ggplot(tmp,aes(x=Assigned,y=sum))+
  geom_point()
cor.test(tmp$Assigned,tmp$sum)
#cor 0.98
