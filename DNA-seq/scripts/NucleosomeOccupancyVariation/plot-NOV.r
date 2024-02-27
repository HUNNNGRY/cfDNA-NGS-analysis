# plot NDR rel depth region
# last 210901 by pengfei
getwd()
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)



# exon1end and TSS coverage (new: 202112) ---------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)
## read sample id
#id <- Sys.glob("output/lulab/bam-sorted-deduped-merged/CRC-passQCloose.sort.bam")
#id <- unlist(lapply(strsplit( basename(id), "\\_"),function(x) x[1]))
bams <- read.delim2("./output/lulab/NDR/single_base_dep-v2/bam.list",header = F)$V1
id <- basename(bams)
id <- gsub(".bam","",id,fixed = T)


## read depth file
### exon1 top expressed rank
#region <- "tss"  # "exon1end" ,"tss"
#group <- "NC"  # "NC","STAD" "CRC"

get.plotCoverage.mat <- function(x){
# x <- "fpkm-001-01"
read <- function(x){
  #x <- "output/lulab/NDR/single_base_dep-v2/promoter2000exon1end2000_bloodTx_fpkm-001-01.txt"
  e1.top <- data.table::fread(x, header=F, sep="\t")
  #e1.top[1:3,]
  colnames(e1.top) <- c("chr","pos",id)
  l <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(x))
  l <- gsub(".txt","",l)
  e1.top$Group <- l
  return(e1.top)
}

dat0 <- lapply(files,read)
dat0 <- do.call("rbind",dat0)


e1.top = as.data.frame(dat0[dat0$Group==x,])

#filter pass qc sample
e1.top <- e1.top[,c(1,2,grep("passQCstringent",colnames(e1.top)))]
#id <- colnames(e1.top)[3:5]
head(e1.top)
e1.top$sum.depth <- e1.top[[paste0(group,"-passQCstringent")]]
e1.top = e1.top[,c("chr","pos","sum.depth")]


## read bed file
{ # exon1 top
  bed <- read.table(paste0("ref/gtf/bloodNDR/promoter2000",region,"2000_bloodTx_",x,".bed"),header = F,sep = "\t")
  head(bed)
}

#table(duplicated(bed$V4))


#head(bed,3)
bed$len <- abs(bed$V2-bed$V3)
#table(bed$len)
colnames(bed) <- c("chr","left","right","ensg","score","str","len")
bed <- bed[, c("chr","left","right","ensg","len")]
bed$ensg <- gsub("promoter300100exon1end_|promoter150TSS50_","",bed$ensg) 

### get tss and exon1end pos
tss.exon1 <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
#head(tss.exon1,3)
tss.exon1 <- tss.exon1[,c("ensg","strand","tss","exon1end")]
bed <- dplyr::left_join(bed,tss.exon1)
bed$ensg <- unlist(lapply(strsplit(bed$ensg,"\\."), function(x) x[1]))


## expand bed record to sing base 
n <- 1
bed.expan <- {}
for(n in 1:nrow(bed)){
  gene <- bed$ensg[n]
  ch <- bed$chr[n]
  num <- bed$len[n]
  l <- bed$left[n] 
  tss.tmp <- bed$tss[n]
  exon1end.tmp <- bed$exon1end[n]
  strand.tmp <- bed$strand[n]
  tmp <- data.frame(ensg=rep(gene,num),chr=ch,pos=seq(l+1,l+num),tss=tss.tmp,exon1end=exon1end.tmp,strand=strand.tmp)
  bed.expan <- rbind(bed.expan,tmp)
}


## add ensg info to depth table
e1.top.m <- dplyr::left_join(bed.expan,e1.top)
head(e1.top.m,2)
e1.top.m$rel.exon1end <- e1.top.m$pos-e1.top.m[[region]]
e1.top.m$rel.exon1end[e1.top.m$strand=="-"] <- e1.top.m$exon1end[e1.top.m$strand=="-"]-e1.top.m$pos[e1.top.m$strand=="-"]
#summary(e1.top.m$rel.exon1end)
e1.top.m$flank <- "Y"
e1.top.m$flank[e1.top.m$rel.exon1end>=-1000 & e1.top.m$rel.exon1end<=1000] <- "N"
#table(e1.top.m$flank)
library(tidyverse)
flank.df <- e1.top.m %>% dplyr::group_by(ensg,flank) %>% dplyr::summarise(flank.dep=median(sum.depth))
flank.df <- flank.df[flank.df$flank=="Y",]
flank.df <- as.data.frame(flank.df)
flank.df <- flank.df[,c("ensg","flank.dep")]
#head(flank.df,3)

e1.top.m <- dplyr::left_join(e1.top.m[e1.top.m$flank=="N",],flank.df)
e1.top.m$rel.dep <- e1.top.m$sum.depth/e1.top.m$flank.dep
#summary(e1.top.m$rel.dep)
e1.top.m$rel.dep[e1.top.m$rel.dep>2] <- 2  # trim to 2
#dim(e1.top.m)
#e1.top.m[1:3,]
#hist(e1.top.m[e1.top.m$rel.exon1end=="-1000","rel.dep"])

## get median from 1000 tx
tx.df <- e1.top.m %>% dplyr::group_by(rel.exon1end) %>% dplyr::summarise(median.tx=mean(rel.dep))
#head(tx.df,3)
#dim(tx.df)
#hist(tx.df$median.tx)

e1.top.m <- dplyr::left_join(e1.top.m,tx.df)
e1.top.m$Group <- x
return(e1.top.m)
#head(e1.top.m,3)
#summary(e1.top.m$flank)
}

library(ggplot2)

## plot exon1end NC
group <- "NC"
region <- "exon1end"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat1 <- lapply(l0,get.plotCoverage.mat)
#head(dat1[[1]])
dat1 <- do.call("rbind",dat1)
#table(dat1$Group)
dat1$Group <- factor(dat1$Group,levels = c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) 
a <- ggplot(dat1, aes(x = rel.exon1end, y = median.tx, color = Group)) + # group=ensg, 
  geom_line(alpha=.9,size=1.5) +
  #scale_color_brewer(palette = "Blues") +
  #scale_color_manual(name = "Blood Expr", 
  #                   values=c("unexp"="steelblue4","fpkm-001-01"="steelblue3","fpkm-01-5"="steelblue2","fpkm-5-30"="steelblue1","fpkm-30"="steelblue"),
  #                   labels=c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) +
  ggsci::scale_color_jama()+
  ylim(c(0.5,1.2)) + 
  geom_vline(xintercept = c(-300,-100),color="grey50",linetype="dashed") +
  labs(title = "Exon1End",x="Relative to Exon1End",y="Relative depth") +
  
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black"), #,face="bold
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))
a
cowplot::save_plot(filename = "./exon1end-ndr-v4.pdf",base_height = 10,base_width = 15,plot = a)


## plot tss NC
group <- "NC"
region <- "tss"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat2 <- lapply(l0,get.plotCoverage.mat)
dat2 <- do.call("rbind",dat2)

dat2$Group <- factor(dat2$Group,levels = c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) 
b <- ggplot(dat2[dat2$Group=="unexp" | dat2$Group=="fpkm-5-30",],aes(x = rel.exon1end, y = median.tx, color = Group)) + # group=ensg, 
  geom_line(alpha=.9,size=2) +
  #scale_color_brewer(palette = "Blues") +
  # scale_color_manual(name = "Blood Expr", 
  #                    values=c("unexp"="steelblue4","fpkm-001-01"="steelblue3","fpkm-01-5"="steelblue2","fpkm-5-30"="steelblue1","fpkm-30"="steelblue"),
  #                    labels=c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) +
  ggsci::scale_color_nejm()+
  ylim(c(0.5,1.2)) + 
  geom_vline(xintercept = c(-150,50),color="grey50",linetype="dashed") +
  labs(title = "TSS",x="Relative to TSS",y="Relative depth") +
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black"), #,face="bold
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))
b
cowplot::save_plot(filename = "./tss-ndr-v4.pdf",base_height = 10,base_width = 15,plot = b)



## def func.
get.disease.NDR.geneWise.aver.cov <- function(exon1,tss,disease){
dat1 <- exon1
dat2 <- tss
dat1$region <- "Exon1"
dat2$region <- "TSS"
dat11 <- dat1[dat1$rel.exon1end>=-300 & dat1$rel.exon1end<=-100,c("ensg","rel.exon1end","rel.dep","Group","region")]
dat22 <- dat2[dat2$rel.exon1end>=-150 & dat2$rel.exon1end<=50,c("ensg","rel.exon1end","rel.dep","Group","region")]
g <- intersect(dat11$ensg,dat22$ensg)
dat11 <- dat11[dat11$ensg %in% g,]
dat22 <- dat22[dat22$ensg %in% g,]
dat11 <- dat11[!duplicated(dat11$ensg),]
dat <- rbind(dat11,dat22)
library(dplyr)
dat <- dat %>% dplyr::group_by(ensg,Group,region) %>% summarise(rel.dep.mean.region=mean(rel.dep,na.rm=T))
#hist(dat$rel.dep.mean.region)
#table(is.na(rbind(dat11,dat22)$median.tx))
#table(is.na(dat$median.tx))
#table(dat$region)

dat.c <- reshape2::dcast(data = dat, formula = ensg + Group ~ region, mean ,value.var = "rel.dep.mean.region")
#(data = dat[,c("Description","type","p.adjust")],Description ~ type,mean,value.var=c("p.adjust"))
dat.c$Group <- as.character(dat.c$Group)
#table(duplicated(dat.c$ensg))
#str(dat.c)
#hist(dat.c$TSS,breaks = 100)
#table(duplicated(dat.c$Exon1))
dat.c <- dat.c[dat.c$Group=="unexp" | dat.c$Group=="fpkm-5-30",]
dat.c$disease <- disease
return(dat.c)
}

## plot exon1end and TSS (NC)
dat.c.nc <- get.disease.NDR.geneWise.aver.cov(dat1,dat2,"NC")

p0 <- ggplot(dat.c.nc,aes(x=Exon1, y=TSS)) +
  geom_density_2d(color="grey50",size=0.3,alpha=0.5)+
  #stat_summary(fun.data=) + 
  geom_point(aes(fill=Group,color=Group),shape=21,size=2)+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  geom_smooth(method='lm', formula= y~x,color="black") +
  ggpubr::stat_cor(label.x = 0.15,size=10,label.y = 1.6)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  ggpubr::stat_regline_equation(label.x = 0.15,size=10,label.y = 1.5)+ #this means at 30th unit regresion line equation will be shown
  # xlim(c(0.30,0.4))+
  # ylim(c(0.35,0.48))+
  theme_classic() +
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_text(size=18,face="bold")) +
  ylab("cfDNA Nucleosome Occupancy") +
  guides(color=guide_legend("Group"))
p0
p <- ggExtra::ggMarginal(p0, type = "density",color="black",size=4, groupColour=F,groupFill=T,alpha=0.6)
p

pdf("./cfDNA-TSS-exon1end.pdf",width = 10,height = 8)
p
dev.off()


# add group compare boxplot for blood highly expressed gene
## plot exon1end and TSS (CRC)
group <- "CRC"
region <- "exon1end"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat3 <- lapply(l0,get.plotCoverage.mat)
dat3 <- do.call("rbind",dat3)
dat3$Group <- factor(dat3$Group,levels = c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) 

region <- "tss"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat4 <- lapply(l0,get.plotCoverage.mat)
dat4 <- do.call("rbind",dat4)
dat4$Group <- factor(dat4$Group,levels = c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) 

dat.c.crc <- get.disease.NDR.geneWise.aver.cov(dat3,dat4,"CRC")

## plot exon1end and TSS (STAD)
group <- "STAD"
region <- "exon1end"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat5 <- lapply(l0,get.plotCoverage.mat)
dat5 <- do.call("rbind",dat5)
dat5$Group <- factor(dat5$Group,levels = c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) 

region <- "tss"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat6 <- lapply(l0,get.plotCoverage.mat)
dat6 <- do.call("rbind",dat6)
dat6$Group <- factor(dat6$Group,levels = c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")) 

dat.c.stad <- get.disease.NDR.geneWise.aver.cov(dat5,dat6,"STAD")



## combine and plot box
dat <- rbind(dat.c.nc,dat.c.crc,dat.c.stad)
colnames(dat)
dat.m <- reshape2::melt(dat,id.var=c("ensg","Group","disease"))
p <- ggplot(dat.m,aes(x=variable, y=value, fill=Group)) +
  geom_boxplot(color="black",shape=21,size=1,alpha=0.6)+
  ggsci::scale_fill_nejm()+
  ggpubr::stat_compare_means(mapping = aes(group=Group), #comparisons = ref.group = "NC", 
                             bracket.size = 0.3, size = 6,color="grey30", method = "wilcox.test",hide.ns=T) + # vjust = 2, hjust=1.2,
  theme_minimal() +
  theme(axis.text=element_text(size=24,face="bold"),
        strip.text.x = element_text(size = 28,face="bold", colour = "grey30"),
        strip.text.y = element_text(size = 28,face="bold", colour = "grey30"),
        axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=32,face="bold"),
        legend.title=element_text(size=28,face="bold")) +
  ylab("cfDNA Nucleosome Occupancy Difference") +
  guides(color=guide_legend("Group"))+
  facet_grid(disease~.)
pdf("./cfDNA-TSS-exon1end-diff-boxplot.pdf",width = 10,height = 12)
p
dev.off()

head(dat,3)

## plot all coverage depth
plot.tss.cov.depth <- function(dat){
  #dat <- dat3
  #colnames(dat)
  li <- c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")
  #regions <- c("tss","exon1")
  
  dat <- dat[,c("rel.exon1end","median.tx","Group")] # ,"region"
  #table(duplicated(dat))
  dat <- dat[!duplicated(dat),]
  
  for (i in li){
    l <- dat$Group==i
    tmp <- lowess(dat[l,"rel.exon1end"],dat[l,"median.tx"],f = 0.05)
    dat[l,"rel.exon1end"] <-  tmp$x
    dat[l,"median.tx"] <-  tmp$y
  }

  dat$Group <- factor(dat$Group,levels = rev(li)) 
  p <- ggplot(dat[dat$Group=="unexp" | dat$Group=="fpkm-5-30",],aes(x = rel.exon1end, y = median.tx, color = Group)) + # group=ensg, 
    geom_line(alpha=.9,size=2) +
    geom_rect(xmin=-150,ymin=0,xmax=50,ymax=1.5,fill = "grey",color="white",alpha=0.0002)+
    geom_vline(xintercept = c(-150,50),color="grey30", linetype="dashed", size=0.8) +
    #geom_text(x=0.150, y=0.28, label="167 bp",size=10) +
    annotate(geom="text", x=-240, y=0.56, label="-150",color="grey30",size=8)+
    annotate(geom="text", x=120, y=0.56, label="+50",color="grey30",size=8)+
    annotate("rect", xmin = -149, xmax = 49, ymin = 0.5, ymax = 1.2,
             alpha = .4,fill = "grey90")+
    labs(title = "TSS",x="Relative to TSS",y="Relative depth") +
    ylim(c(0.5,1.2)) + 
    ggsci::scale_color_nejm()+
    theme_classic() + 
    theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
          axis.title = element_text(size = 28,color ="black",face="bold"), 
          axis.text = element_text(size= 24,color = "black"), #,face="bold
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
          axis.text.y = element_text( hjust = 0 ), # angle = 45,
          panel.grid=element_blank(),
          #legend.position = "top",
          legend.text = element_text(size= 26,color = "black",face="bold"),
          legend.title= element_text(size= 26,color = "black",face="bold"))
  print(p)
}
plot.exon1.cov.depth <- function(dat){
  #dat <- dat1
  #str(dat)
  li <- c("unexp","fpkm-001-01","fpkm-01-5","fpkm-5-30","fpkm-30")
  #regions <- c("tss","exon1")
  
  dat <- dat[,c("rel.exon1end","median.tx","Group")] # ,"region"
  #table(duplicated(dat))
  dat <- dat[!duplicated(dat),]
  
  for (i in li){
    l <- dat$Group==i
    tmp <- lowess(dat[l,"rel.exon1end"],dat[l,"median.tx"],f = 0.05)
    dat[l,"rel.exon1end"] <-  tmp$x
    dat[l,"median.tx"] <-  tmp$y
  }
  
  dat$Group <- factor(dat$Group,levels = rev(li))
  p <- ggplot(dat[dat$Group=="unexp" | dat$Group=="fpkm-5-30",],aes(x = rel.exon1end, y = median.tx, color = Group)) + # group=ensg, 
    geom_line(alpha=.9,size=2) +
    geom_rect(xmin=-300,ymin=0,xmax=-100,ymax=1.5,fill = "grey",color="white",alpha=0.0002)+
    geom_vline(xintercept = c(-300,-100),color="grey30", linetype="dashed", size=0.8) +
    #geom_text(x=0.150, y=0.28, label="167 bp",size=10) +
    annotate(geom="text", x=-240, y=0.56, label="-300",color="grey30",size=8)+
    annotate(geom="text", x=120, y=0.56, label="-100",color="grey30",size=8)+
    annotate("rect", xmin = -299, xmax = -101, ymin = 0.5, ymax = 1.2,
             alpha = .4,fill = "grey90")+
    labs(title = "Exon1",x="Relative to Exon1",y="Relative depth") +
    ylim(c(0.5,1.2)) + 
    ggsci::scale_color_nejm()+
    theme_classic() + 
    theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
          axis.title = element_text(size = 28,color ="black",face="bold"), 
          axis.text = element_text(size= 24,color = "black"), #,face="bold
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
          axis.text.y = element_text( hjust = 0 ), # angle = 45,
          panel.grid=element_blank(),
          #legend.position = "top",
          legend.text = element_text(size= 26,color = "black",face="bold"),
          legend.title= element_text(size= 26,color = "black",face="bold"))
  print(p)
}

pdf("./tss-NC-cov.pdf",width = 11,height = 6)
plot.tss.cov.depth(dat2)
dev.off()
pdf("./exon1-NC-cov.pdf",width = 11,height = 6)
plot.exon1.cov.depth(dat1)
dev.off()
pdf("./tss-CRC-cov.pdf",width = 11,height = 6)
plot.tss.cov.depth(dat4)
dev.off()
pdf("./exon1-CRC-cov.pdf",width = 11,height = 6)
plot.exon1.cov.depth(dat3)
dev.off()
pdf("./tss-STAD-cov.pdf",width = 11,height = 6)
plot.tss.cov.depth(dat6)
dev.off()
pdf("./exon1-STAD-cov.pdf",width = 11,height = 6)
plot.exon1.cov.depth(dat5)
dev.off()




# plot&cmp candidates coverage (GI ML candidates example, newest version) --------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)
## read sample id
#id <- Sys.glob("output/lulab/bam-sorted-deduped-merged/CRC-passQCloose.sort.bam")
#id <- unlist(lapply(strsplit( basename(id), "\\_"),function(x) x[1]))
bams <- read.delim2("./output/lulab/NDR/single_base_dep-v2-candidates/bam.list",header = F)$V1
id <- basename(bams)
id <- gsub(".sort.bam","",id,fixed = T)


## define mat func.
get.plotCoverage.mat.v2 <- function(x=("fpkm-001-01"),disease=c("NC"|"CRC"|"STAD"),region=c("tss"|"exon1end"),chr=c("chr1"),pos=c(100,1000)){
  #x <- "RFLOO-1202"
  #region <- "exon1end"
  #disease <- "CRC"
  read <- function(x){
    #x <- "output/lulab/NDR/single_base_dep-v2/promoter2000exon1end2000_bloodTx_fpkm-001-01.txt"
    e1.top <- data.table::fread(x, header=F, sep="\t")
    #e1.top[1:3,]
    colnames(e1.top) <- c("chr","pos",id)
    l <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(x))
    l <- gsub(".txt","",l)
    e1.top$Group <- l
    return(e1.top)
  }
  
  dat0 <- lapply(files,read)
  dat0 <- do.call("rbind",dat0)

  e1.top = as.data.frame(dat0[dat0$Group==x,])
  
  #filter pass qc sample
  e1.top <- e1.top[,c(1,2,grep(paste0(disease,"-passQCstringent"),colnames(e1.top)))]
  #id <- colnames(e1.top)[3:5]
  head(e1.top)
  e1.top$sum.depth <- e1.top[[paste0(disease,"-passQCstringent")]]
  e1.top = e1.top[,c("chr","pos","sum.depth")]
  e1.top <- e1.top[e1.top$chr==chr & e1.top$pos>=pos[1] & e1.top$pos<=pos[2],]
  
  ## read bed file
  { # exon1 top
    bed <- read.table(paste0("ref/gtf/candidatesNDR/promoter2000",region,"2000_bloodTx_",x,".bed"),header = F,sep = "\t")
    head(bed)
  }
  
  #table(duplicated(bed$V4))
  
  
  #head(bed,3)
  bed$len <- abs(bed$V2-bed$V3)
  #table(bed$len)
  colnames(bed) <- c("chr","left","right","ensg","score","str","len")
  bed <- bed[, c("chr","left","right","ensg","len")]
  bed$ensg <- gsub("promoter300100exon1end_|promoter150TSS50_","",bed$ensg) 
  
  ### get tss and exon1end pos
  tss.exon1 <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
  #head(tss.exon1,3)
  tss.exon1 <- tss.exon1[,c("ensg","strand","tss","exon1end")]
  bed <- dplyr::left_join(bed,tss.exon1)
  bed$ensg <- unlist(lapply(strsplit(bed$ensg,"\\."), function(x) x[1]))
  
  
  ## expand bed record to sing base 
  n <- 1
  bed.expan <- {}
  for(n in 1:nrow(bed)){
    gene <- bed$ensg[n]
    ch <- bed$chr[n]
    num <- bed$len[n]
    l <- bed$left[n] 
    tss.tmp <- bed$tss[n]
    exon1end.tmp <- bed$exon1end[n]
    strand.tmp <- bed$strand[n]
    tmp <- data.frame(ensg=rep(gene,num),chr=ch,pos=seq(l+1,l+num),tss=tss.tmp,exon1end=exon1end.tmp,strand=strand.tmp)
    bed.expan <- rbind(bed.expan,tmp)
  }
  bed.expan <- bed.expan[bed.expan$chr==chr & bed.expan$pos>=pos[1] & bed.expan$pos<=pos[2],]
  
  
  ## add ensg info to depth table
  e1.top.m <- dplyr::left_join(bed.expan,e1.top)
  #unique(e1.top.m$ensg)
  e1.top.m$rel.exon1end <- e1.top.m$pos-e1.top.m[[region]]
  e1.top.m$rel.exon1end[e1.top.m$strand=="-"] <- e1.top.m$exon1end[e1.top.m$strand=="-"]-e1.top.m$pos[e1.top.m$strand=="-"]
  #summary(e1.top.m$rel.exon1end)
  e1.top.m$flank <- "Y"
  e1.top.m$flank[e1.top.m$rel.exon1end>=-1000 & e1.top.m$rel.exon1end<=1000] <- "N"
  #table(e1.top.m$flank)
  library(tidyverse)
  flank.df <- e1.top.m %>% dplyr::group_by(ensg,flank) %>% dplyr::summarise(flank.dep=median(sum.depth))
  flank.df <- flank.df[flank.df$flank=="Y",]
  flank.df <- as.data.frame(flank.df)
  flank.df <- flank.df[,c("ensg","flank.dep")]
  #head(flank.df,3)
  
  e1.top.m <- dplyr::left_join(e1.top.m[e1.top.m$flank=="N",],flank.df)
  e1.top.m$rel.dep <- e1.top.m$sum.depth/e1.top.m$flank.dep
  #summary(e1.top.m$rel.dep)
  e1.top.m$rel.dep[e1.top.m$rel.dep>2] <- 2  # trim to 2
  #dim(e1.top.m)
  #e1.top.m[1:3,]
  #hist(e1.top.m[e1.top.m$rel.exon1end=="-1000","rel.dep"])
  
  ## get median from 1000 tx
  tx.df <- e1.top.m %>% dplyr::group_by(rel.exon1end) %>% dplyr::summarise(median.tx=mean(rel.dep))
  #head(tx.df,3)
  #dim(tx.df)
  #hist(tx.df$median.tx)
  
  e1.top.m <- dplyr::left_join(e1.top.m,tx.df)
  e1.top.m$Group <- x
  e1.top.m$Disease <- disease
  
  return(e1.top.m)
  #head(e1.top.m,3)
  #summary(e1.top.m$flank)
}


## plot exon1end NC
#unique(dat1$ensg)
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2-candidates/promoter2000*2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
l0 <- unique(l0)

#TSPAN17	ENSG00000048140.17 chr5	gene	176647387	176659054
#SLC66A1	ENSG00000040487.12 chr1	gene	19312326	19329300    
#ROS1	ENSG00000047936.10 chr6	gene	117288300	117425855
ensg <- "ENSG00000048140"  # "ENSG00000040487" "ENSG00000048140" "ENSG00000047936"
chr <- "chr5"
pos <- c(0,1766590540)

region <- "tss"  # tss  exon1end


disease <- "NC"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2-candidates/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
dat1 <- get.plotCoverage.mat.v2(x = l0, disease = disease, region = region, chr = chr, pos = pos)

disease <- "CRC"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2-candidates/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
dat2 <- get.plotCoverage.mat.v2(x = l0, disease = disease, region = region, chr = chr, pos = pos)

disease <- "STAD"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2-candidates/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
dat3 <- get.plotCoverage.mat.v2(x = l0, disease = disease, region = region, chr = chr, pos = pos)


dat <- rbind(dat1,dat2,dat3)
colnames(dat)


## lowess smooth
dat <- dat[,c("ensg","rel.exon1end","median.tx","Group","Disease")] # ,"region"
#table(duplicated(dat))
dat <- dat[!duplicated(dat),]
unique(dat$ensg)

for (i in l0){
  for(j in c("NC","CRC","STAD")){
    for(k in c(ensg)){
      l <- (dat$Group==i & dat$Disease==j & dat$ensg==k)
      tmp <- lowess(dat[l,"rel.exon1end"],dat[l,"median.tx"],f = 0.05)
      dat[l,"rel.exon1end"] <-  tmp$x
      dat[l,"median.tx"] <-  tmp$y 
    }
  }
}


## plot exon1end
a <- ggplot(dat, aes(x = rel.exon1end, y = median.tx, color = Disease)) + # group=ensg, 
  geom_line(alpha=.9,size=1.5) +
  ggsci::scale_color_jama()+
  ylim(c(0,2)) + 
  geom_vline(xintercept = c(-300,-100),color="grey50",linetype="dashed") +
  labs(title = paste0(ensg,": Exon1end"),x="",y="Relative depth") +
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black"), #,face="bold
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))
a

## plot tss
b <- ggplot(dat, aes(x = rel.exon1end, y = median.tx, color = Disease)) + # group=ensg, 
  geom_line(alpha=.9,size=1.5) +
  ggsci::scale_color_jama()+
  ylim(c(0,2)) + 
  geom_vline(xintercept = c(-150,50),color="grey50",linetype="dashed") +
  labs(title = paste0(ensg,": TSS"),x="",y="Relative depth") +
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black"), #,face="bold
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))
b


#cowplot::save_plot(filename = "./exon1end-ndr-v4-candidates.pdf",base_height = 10,base_width = 15,plot = a)



# plot&cmp candidates coverage (CRC example, older version) --------------------------------------------
#PRTN3 Blood Promoter
options(stringsAsFactors = F)
getwd()
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

files <- Sys.glob("./output/lulab/NDR/archive/*-PRTN3.relativeDepth")  # CRC-PKU-10-wgs
samples <- unlist(lapply(strsplit(files,"/",fixed = T),function(x) x[6]))
samples <- gsub("-PRTN3.relativeDepth","",samples)

read <- function(x){
    #x <- "STAD-PKU-8-wgs"
    tmp <- read.table(paste0("./output/lulab/NDR/archive/",x,"-PRTN3.relativeDepth"),sep = "\t",header = T,check.names = F)
    colnames(tmp)[2] <- "Relative.Depth"
    tmp$Sample <- x
  return(tmp)
}
dat <- lapply(samples,read)
dat <- do.call("rbind",dat)
hist(dat$Relative.Depth)
dat$Group <- unlist(lapply(strsplit(dat$Sample,"-",fixed = T),function(x) x[1]))
dat$Group <- factor(dat$Group,levels = c("NC","CRC","STAD"))
str(dat)

library(dplyr)
dat1 <- dat %>% group_by(TSS,Group) %>% summarise(Mean.Depth=mean(Relative.Depth,na.rm=T),Median.Depth=median(Relative.Depth,na.rm=T))
#median(dat[dat$TSS=="-1000" & dat$Group=="CRC","Relative.Depth"],na.rm=T)


a <- ggplot(dat1[dat1$Group!="STAD",],aes(TSS, Mean.Depth, color = Group)) + # group=ensg, 
  geom_line(alpha=.9,size=1.5) +
  # scale_color_manual(name = "Blood Expr", 
  #                    values=c("blood.top1000"="steelblue4","blood.mid30"="steelblue3","blood.mid60"="steelblue2","blood.tail1000"="steelblue1"),
  #                    labels=c("top","mid30%","mid60%","tail")) +
  
  #ylim(c(0.5,1.2)) + 
  geom_rect(xmin=-150,ymin=0,xmax=50,ymax=1.5,fill = "grey",color="white",alpha=0.0002)+
  geom_vline(xintercept = c(-150,50),color="grey30", linetype="dashed", size=0.8) +
  #geom_text(x=0.150, y=0.28, label="167 bp",size=10) +
  annotate(geom="text", x=-240, y=0.26, label="-150",color="grey30",size=8)+
  annotate(geom="text", x=120, y=0.26, label="+50",color="grey30",size=8)+
  annotate("rect", xmin = -149, xmax = +49, ymin = 0, ymax = 2,
           alpha = .4,fill = "grey90")+
  labs(title = "CRC Biomarker: PRTN3 TSS (Blood Tx)",x="Relative Distance to TSS",y="Relative Depth") +
  scale_color_manual(values=alpha(c("steelblue2","#6B4E04","orange2"), 1),labels=c("HD","CRC","STAD"))  +  # "1" = "red", "2" = "grey", "3" = "blue"
  theme_classic() + 
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 28,color ="black",face="bold"), 
        axis.text = element_text(size= 24,color = "black"), #,face="bold
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
        axis.text.y = element_text( hjust = 0 ), # angle = 45,
        panel.grid=element_blank(),
        #legend.position = "top",
        legend.text = element_text(size= 26,color = "black",face="bold"),
        legend.title= element_text(size= 26,color = "black",face="bold"))
a
#cowplot::save_plot(filename = "./tss-ndr-v2.pdf",base_height = 12,base_width = 17,plot = a)
pdf("./NDR-CRC-eg.pdf",width = 11,height = 6)
a
dev.off()





# test flank correlation (new: 20220115) ---------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)
bams <- read.delim2("./output/lulab/NDR/single_base_dep-v2/bam.list",header = F)$V1
id <- basename(bams)
id <- gsub(".bam","",id,fixed = T)


## read depth file
### exon1 top expressed rank
#region <- "tss"  # "exon1end" ,"tss"
#group <- "NC"  # "NC","STAD" "CRC"

get.plotCoverage.mat <- function(x){
  # x <- "fpkm-001-01"
  read <- function(x){
    #x <- "output/lulab/NDR/single_base_dep-v2/promoter2000exon1end2000_bloodTx_fpkm-001-01.txt"
    e1.top <- data.table::fread(x, header=F, sep="\t")
    #e1.top[1:3,]
    colnames(e1.top) <- c("chr","pos",id)
    l <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(x))
    l <- gsub(".txt","",l)
    e1.top$Group <- l
    return(e1.top)
  }
  
  dat0 <- lapply(files,read)
  dat0 <- do.call("rbind",dat0)
  
  
  e1.top = as.data.frame(dat0[dat0$Group==x,])
  
  #filter pass qc sample
  e1.top <- e1.top[,c(1,2,grep("passQCstringent",colnames(e1.top)))]
  #id <- colnames(e1.top)[3:5]
  head(e1.top)
  e1.top$sum.depth <- e1.top[[paste0(group,"-passQCstringent")]]
  e1.top = e1.top[,c("chr","pos","sum.depth")]
  
  
  ## read bed file
  { # exon1 top
    bed <- read.table(paste0("ref/gtf/bloodNDR/promoter2000",region,"2000_bloodTx_",x,".bed"),header = F,sep = "\t")
    head(bed)
  }
  
  #table(duplicated(bed$V4))
  
  
  #head(bed,3)
  bed$len <- abs(bed$V2-bed$V3)
  #table(bed$len)
  colnames(bed) <- c("chr","left","right","ensg","score","str","len")
  bed <- bed[, c("chr","left","right","ensg","len")]
  bed$ensg <- gsub("promoter300100exon1end_|promoter150TSS50_","",bed$ensg) 
  
  ### get tss and exon1end pos
  tss.exon1 <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/tss-exon1.gtf",sep = "\t",header = T)
  #head(tss.exon1,3)
  tss.exon1 <- tss.exon1[,c("ensg","strand","tss","exon1end")]
  bed <- dplyr::left_join(bed,tss.exon1)
  bed$ensg <- unlist(lapply(strsplit(bed$ensg,"\\."), function(x) x[1]))
  
  
  ## expand bed record to sing base 
  n <- 1
  bed.expan <- {}
  for(n in 1:nrow(bed)){
    gene <- bed$ensg[n]
    ch <- bed$chr[n]
    num <- bed$len[n]
    l <- bed$left[n] 
    tss.tmp <- bed$tss[n]
    exon1end.tmp <- bed$exon1end[n]
    strand.tmp <- bed$strand[n]
    tmp <- data.frame(ensg=rep(gene,num),chr=ch,pos=seq(l+1,l+num),tss=tss.tmp,exon1end=exon1end.tmp,strand=strand.tmp)
    bed.expan <- rbind(bed.expan,tmp)
  }
  
  
  ## add ensg info to depth table
  e1.top.m <- dplyr::left_join(bed.expan,e1.top)
  head(e1.top.m,2)
  e1.top.m$rel.exon1end <- e1.top.m$pos-e1.top.m[[region]]
  e1.top.m$rel.exon1end[e1.top.m$strand=="-"] <- e1.top.m$exon1end[e1.top.m$strand=="-"]-e1.top.m$pos[e1.top.m$strand=="-"]
  #summary(e1.top.m$rel.exon1end)
  e1.top.m$flank <- "Y"
  #e1.top.m$flank[e1.top.m$rel.exon1end>=-1000 & e1.top.m$rel.exon1end<=1000] <- "N"
  e1.top.m$flank[e1.top.m$rel.exon1end>=-2000 & e1.top.m$rel.exon1end<=2000] <- "N"
  #table(e1.top.m$flank)
  library(tidyverse)
  flank.df <- e1.top.m %>% dplyr::group_by(ensg,flank) %>% dplyr::summarise(flank.dep=median(sum.depth))
  flank.df <- flank.df[flank.df$flank=="Y",]
  flank.df <- as.data.frame(flank.df)
  flank.df <- flank.df[,c("ensg","flank.dep")]
  #head(flank.df,3)
  
  e1.top.m <- dplyr::left_join(e1.top.m[e1.top.m$flank=="N",],flank.df)
  e1.top.m$rel.dep <- e1.top.m$sum.depth/e1.top.m$flank.dep
  #summary(e1.top.m$rel.dep)
  e1.top.m$rel.dep[e1.top.m$rel.dep>2] <- 2  # trim to 2
  #dim(e1.top.m)
  #e1.top.m[1:3,]
  #hist(e1.top.m[e1.top.m$rel.exon1end=="-1000","rel.dep"])
  
  ## get median from 1000 tx
  tx.df <- e1.top.m %>% dplyr::group_by(rel.exon1end) %>% dplyr::summarise(median.tx=mean(rel.dep))
  #head(tx.df,3)
  #dim(tx.df)
  #hist(tx.df$median.tx)
  
  e1.top.m <- dplyr::left_join(e1.top.m,tx.df)
  e1.top.m$Group <- x
  return(e1.top.m)
  #head(e1.top.m,3)
  #summary(e1.top.m$flank)
}

library(ggplot2)

## plot exon1end NC
group <- "NC"
region <- "exon1end"
files <- Sys.glob(paste0("output/lulab/NDR/single_base_dep-v2/promoter2000",region,"2000_bloodTx_*.txt"))  # 001-01
l0 <- gsub("promoter2000exon1end2000_bloodTx_|promoter2000tss2000_bloodTx_","",basename(files))
l0 <- gsub(".txt","",l0)
dat1 <- lapply(l0,get.plotCoverage.mat)

dat1 <- do.call("rbind",dat1)

library(dplyr)
dat1$test <- "N"
dat1$test[dat1$rel.exon1end<=-1000] <- "-2k_-1k"
dat1$test[dat1$rel.exon1end>=1000] <- "+1k_+2k"
table(dat1$Group)
#dat1[dat1$Group=="N",]

dat1[1:4,]
dat.test <- dat1 %>% dplyr::group_by(ensg,test,Group) %>% dplyr::summarise(mean.dep=mean(rel.dep))  # fpkm-30,  unexp
#dat.test <- dat.test[dat.test$Group=="fpkm-30",]  # unexp

#dat.test.acast <- reshape2::acast(dat.test,formula = Group+test~ensg)
#dat.test$gini <- edgeR::gini(as.numeric(dat.test$mean.dep))


dat.test <- dat.test[dat.test$Group=="fpkm-30",]  # unexp
dat.test.cast <- reshape2::dcast(data = dat.test, formula = ensg~test)
plot(dat.test.cast$`-2k_-1k`,dat.test.cast$`+1k_+2k`,xlab="-2k_-1k",ylab="+1k_+2k",main="bloox fpkm-30 tx")
text(x=0.9,y=1.2,labels="pearson: 0.6397",cex=2)
cor.test(dat.test.cast$`-2k_-1k`,dat.test.cast$`+1k_+2k`,method = "pearson")  # 0.5621
