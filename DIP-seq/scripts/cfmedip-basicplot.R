# plot basic cfmedip 
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")

#### duplication boxplot ####
d <- read.table("./meta/lulab/sample_table.txt",stringsAsFactors = F,check.names = F, header = T, sep = "\t")
d[1:3,1:3]
colnames(d)
d <- d[,c("data_id","group","nonDup-ratio","medip_passQC")]
str(d)
d$group <- factor(d$group,levels = c("NC","CRC","STAD"))
colnames(d)[3] <- "non-duplication"
d$duplication <- 1-d$`non-duplication`
d <- d[d$medip_passQC=="Y",]
d$data_id <- factor(d$data_id,levels = unique(d[order(d$duplication),"data_id"]))

f.plot <- d
set.seed(1)
b <- runif(nrow(f.plot), -0.2, 0.2)
ggplot(f.plot, aes(x = group, y = duplication, fill=group))+   # , fill = Group
  #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
  geom_boxplot(outlier.alpha = 0, outlier.size = 0)  + 
  geom_point(aes(x=as.numeric(group)+b,y=duplication),color="black",size=2,stroke=0.5,alpha=0.5) + 
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

#str(d.m)
#table(d.m$variable)
d.m <- reshape2::melt(d, id.var=c("data_id","group","medip_passQC"))
d.m$variable <- factor(d.m$variable,levels = c("duplication","non-duplication"))
ggplot(d.m, aes(x = data_id, y = value, fill=variable))+   # , fill = Group
  #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
  geom_bar(stat = "identity",position = "stack",color="black")  + 
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


#### plot saturation and relH boxplot ####
f <- read.table("./meta/lulab/archive/sample_table_20210601.txt",stringsAsFactors = F,header = T, sep = "\t")
#f <- f[f$overall_QC=="Y" & f$group!="STAD",]
f <- f[f$overall_QC=="Y" & f$antibody!=1,]
table(f$antibody)
mytheme <- theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
                 axis.title = element_text(size = 28,color ="black",face="bold"), 
                 axis.text = element_text(size= 24,color = "black",face="bold"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text( hjust = 1 ), # angle = 45,
                 panel.grid=element_blank(),
                 #legend.position = "top",
                 legend.text = element_text(size= 26,color = "black",face="bold"),
                 legend.title= element_text(size= 26,color = "black",face="bold"))
colnames(f)
f <- f[,c("data_id","group","enrichment.score.relH","enrichment.score.GoGe")]  # estimated.max.saturation
colnames(f) <- c("sample","group","relH","GoGe")  # saturation
f <- reshape2::melt(f,id.var = c("sample",'group') )
colnames(f)[3] <- "QC"

b <- runif(nrow(f), -0.15, 0.15)
ggplot(f,aes(QC,value,fill=QC)) + 
  geom_boxplot(outlier.colour = "white", outlier.size = 0.1) + 
  geom_point( aes(x = as.numeric(QC) + b, y = value, fill=QC) ,size=2, shape=21,alpha=0.2, stroke = 1,color="black") +
  ggsci::scale_fill_jama() +
  #scale_fill_manual(values=alpha(c("steelblue2","orange2","firebrick2"), 1))  + 
  theme_classic() + mytheme +
  geom_hline(yintercept = c(2.2,1.4),color=c("steelblue2","orange2"),linetype='dotted',cex=0.5,pch=1) +
  ylim(c(0,4)) + 
  geom_text(aes(0,2.2,label = "2.2", vjust = 0, hjust=2))


#### plot antibody boxplot ####
f <- read.table("./meta/lulab/sample_table.txt",stringsAsFactors = F,header = T, sep = "\t")
#f <- f[f$overall_QC=="Y" & f$group!="STAD",]
f <- f[f$medip_passQC=="Y",]
#f$log10libyield <- log10(f$lib_yield)
#colnames(f)
# nonDup.ratio , coverage_ratio
#str(f)
f$antibody <- paste0(f$antibody,"x")
f$antibody <- factor(f$antibody,levels = c("1x","1.5x","2x"))

#table(f$antibody)
mytheme <- theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
                 axis.title = element_text(size = 28,color ="black",face="bold"), 
                 axis.text = element_text(size= 24,color = "black",face="bold"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text( hjust = 1 ), # angle = 45,
                 panel.grid=element_blank(),
                 #legend.position = "top",
                 legend.text = element_text(size= 26,color = "black",face="bold"),
                 legend.title= element_text(size= 26,color = "black",face="bold"))

ggplot(f, aes(x = antibody, y = coverage_ratio, fill=antibody))+   # , fill = Group
  #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
  geom_boxplot( # position=position_dodge(0.5),width=0.5,size=0.4,
    outlier.alpha = 1, outlier.size = 0.5)  + xlab("Antibody") + ylab("Genome Coverage Ratio") +
  theme_classic() + mytheme +
  scale_fill_manual(values = c("grey","steelblue","orange2")) +
  ggpubr::stat_compare_means(ref.group = "1x",
                     label = "p.format",
                     method = "wilcox.test",
                     size = 6,
                     hide.ns = T) + 
  guides(fill=guide_legend(title="Antibody"))



#### plot passQC boxplot ####
f <- read.table("./meta/lulab/qc_stats.txt",stringsAsFactors = F,check.names = F, header = T, sep = "\t")

f <- f[f$medip_passQCloose!="excluded",]
colnames(f) <- gsub("-|_","\\.",colnames(f))

colnames(f)
qc.col <- c("mapped.ratio","coverage.ratio","uniq.reads","enrichment.score.relH","enrichment.score.GoGe","estimated.max.saturation")
qc.col.criteria <- paste0(qc.col,".","criteria")
f.dat <- f[,qc.col]
f.criteria <- f[,qc.col.criteria]
f.plot <- log10(f.dat/f.criteria)
  
# f.dat[1:3,1:3]
# f.criteria[1:3,1:3]
# f.plot[1:3,1:3]
f.plot <- tidyr::gather(f.plot,key = "QC", value = "log10ratio", 1:ncol(f.plot))
f.plot$id <- f$data.id
#unique(f.plot$QC)
QC.factor <- c("uniq.reads","coverage.ratio","mapped.ratio",
               "enrichment.score.relH","enrichment.score.GoGe","estimated.max.saturation")
QC.lab <- c("Criteria 6: Uniq reads",
                "Criteria 5: Coverage ratio",
                "Criteria 4: Mapped ratio",
                "Criteria 3: Enrich relH score",
                "Criteria 2: Enrich GoGe score",
                "Criteria 1: Saturation")
f.plot$QC <- factor(f.plot$QC,levels = QC.factor)
f.plot$point.col <- "passQC"
f.plot$point.col[f.plot$log10ratio<(0)] <- "faillQC" 

table(f.plot$point.col)

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







#### plot genome cov barplot ####
# e <- read.table("./shaozhen.txt",stringsAsFactors = F,header = T, sep = "\t")
# colnames(e)
# e$Value <- as.numeric(e$Value)
# mytheme <- theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
#                  axis.title = element_text(size = 28,color ="black",face="bold"), 
#                  axis.text = element_text(size= 20,color = "black",face="bold"),
#                  panel.grid.minor.y = element_blank(),
#                  panel.grid.minor.x = element_blank(),
#                  axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
#                  panel.grid=element_blank(),
#                  #legend.position = "top",
#                  legend.text = element_text(size= 26,color = "black",face="bold"),
#                  legend.title= element_text(size= 26,color = "black",face="bold"),
#                   strip.text = element_text(size=32,face="bold"))
# 
# ggplot(e, aes(x = Result, y = Value , fill=Kit))+   #
#   #labs(y="Score",x= NULL,title = "inhouse-15pairCRCtissue-ESTIMATE")+  
#   geom_bar(stat = "summary", position = "dodge",color="black")  + 
#   #xlab("Antibody") + ylab("Genome Coverage Ratio") +
#   theme_classic() + mytheme  + 
#   scale_fill_manual(values = c("grey","steelblue","orange2")) +
#   scale_y_continuous(breaks=c(0,1,2),labels = c(0,1,2)) + 
#   facet_wrap(~ Sample,scales = "free",ncol = 1) 
  






##############################################
#barplot or boxplot of 22 regions
#dmr distribution on genomic regions
#last 210823
##############################################
getwd()
#setwd("/BioII/lulab_b/baopengfei/projects/methylation")

#### meta barplot of regions ####
regions <- c("gene","promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE",'retroposon')
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")

tmp <- read.table("./output/lulab/matrix/count_matrix_gene.txt.summary",sep = "\t",check.names = F,header = T,stringsAsFactors = F,row.names = 1)
tmp <- as.data.frame(t(tmp))
tmp$sample <- unlist(lapply(strsplit(rownames(tmp),"/"),function(x) x[4]))
tmp$sample <- gsub(".bam","",tmp$sample,fixed = T)
tmp$region <- "gene"
for(i in c("promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE",'retroposon')){
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
read.ratio$region <- factor(read.ratio$region,levels=regions)

region.meta <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/summary-regions.txt",header = T,stringsAsFactors = F)
region.meta <- region.meta[region.meta$region %in% regions,]
region.meta$region <- factor(region.meta$region,levels=regions) # 
region.meta$coverage <- region.meta$coverage/1000000
table(region.meta$region)


#### 1.plot region total number barplot #### 
library(ggplot2)
library(tidyverse)
library(ggtext)
library(glue)
ggplot(data=region.meta,aes(x=region,y=number)) +
  geom_bar(show.legend = T,stat = "summary") +
  #geom_point(size=1) + 
  scale_x_discrete(breaks=regions,labels=regions) +
  scale_y_continuous(trans='log10') + 
  #xlim(c(0,5000000)) +  
  #ylim(c(0,1000000))
  #xlab( ) + 
  ylab("region total number") +
  #scale_y_continuous(oob = scales::squish) +
  #stat_compare_means(comparisons = list(c("CRC","NC"))) +
  #labs(title="in-house cfMeDIP",x="CRC FDA Biomarker",y="transformed counts") +
  #geom_jitter(width = .2) +
  #scale_x_discrete(labels= function(x) highlight(x, c("promoter","CpG_island"), "red")) + 
  #theme(axis.text.x=element_markdown()) 
  theme_bw() + 
  theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.2, hjust = 1)) 


#### 2.relative coverage barplot #### 
ggplot(data=region.meta,aes(x=region,y=coverage)) +
  geom_bar(show.legend = T,stat = "summary") +
  #geom_point(size=1) + 
  scale_x_discrete(breaks=regions,labels=regions) +
  scale_y_continuous(trans='log10') + 
  #xlim(c(0,5000000)) +  
  #ylim(c(0,1000000))
  #xlab( ) + 
  ylab("total coverage (Mb)") +
  #scale_y_continuous(oob = scales::squish) +
  #stat_compare_means(comparisons = list(c("CRC","NC"))) +
  #labs(title="in-house cfMeDIP",x="CRC FDA Biomarker",y="transformed counts") +
  #geom_jitter(width = .2) +
  #scale_x_discrete(labels= function(x) highlight(x, c("promoter","CpG_island"), "red")) + 
  #theme(axis.text.x=element_markdown()) 
  theme_bw() + 
  theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) +
  theme(axis.text.x = element_text(angle = 90,vjust=0.2, hjust = 1)) 


#### 3.reads distribution boxplot #### 
# mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5,face="bold"),
#                  axis.title = element_text(size = 12,color ="black",face="bold"), 
#                  axis.text = element_text(size= 12,color = "black",face="bold"),
#                  panel.grid.minor.y = element_blank(),
#                  panel.grid.minor.x = element_blank(),
#                  axis.text.x = element_text(angle = 45, hjust = 1 ),
#                  panel.grid=element_blank(),
#                  legend.position = "top",
#                  legend.text = element_text(size= 12),
#                  legend.title= element_text(size= 12))
# ggplot(data=read.ratio,aes(x=region,y=assign.ratio,fill=group)) +
#   geom_boxplot(show.legend = T) + 
#   #scale_fill_manual(values=alpha(c("orange2","steelblue2"), 1))  + 
#   ggsci::scale_fill_jama() +
# 
#   #scale_x_discrete(breaks=regions,labels=regions) +
#   #scale_y_continuous(trans='log10') + 
#   #xlim(c(0,5000000)) +  
#   #ylim(c(0,1000000))
#   #xlab( ) + 
#   ylab("reads coverage ratio") +  # reads assigned ratio
#   #scale_y_continuous(oob = scales::squish) +
#   #stat_compare_means(comparisons = list(c("CRC","NC"))) +
#   #labs(title="in-house cfMeDIP",x="CRC FDA Biomarker",y="transformed counts") +
#   #geom_jitter(width = .2) +
#   #scale_x_discrete(labels= function(x) highlight(x, c("promoter","CpG_island"), "red")) + 
#   #theme(axis.text.x=element_markdown()) 
#   theme_classic() + mytheme +
#   #theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) +
#   #theme(axis.text.x = element_text(angle = 90,vjust=0.2, hjust = 1)) + 
#   ggpubr::stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test")

## version1
bp <- ggbarplot(read.ratio, x = "region", y = "assign.ratio", fill = "group",
                palette = c("steelblue","chocolate","firebrick"), add = "mean_sd",
                position = position_dodge(0.7))
stat.test <- read.ratio %>%
  group_by(region) %>%
  rstatix::wilcox_test(assign.ratio ~ group)
stat.test <- stat.test %>%
  add_xy_position(x = "region", fun = "mean_sd", dodge = 1)
regions.simplify <- c("gene","promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")
bp +  stat_pvalue_manual(stat.test, label = "p",hide.ns=T, size=8, tip.length = 0.03,bracket.nudge.y = 0.4) +   # label = "p.adj.signif"
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  ylab("reads count ratio") + 
  scale_x_discrete(breaks=regions,labels=regions.simplify) +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold")) 

## version2
regions.simplify <- c("gene","promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")

ggbarplot(read.ratio, x = "region", y = "assign.ratio", fill = "region",
                add = "mean_sd", palette = colorRampPalette(brewer.pal(8, "Dark2"))(10),
                position = position_dodge(0.7)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  ylab("reads count ratio") + 
  scale_x_discrete(breaks=regions,labels=regions.simplify) +
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


#### 4. final: dmr number barplot #### 
# old version
ggplot(data=dmr,aes(x=region,y=number,fill=trend)) +
  geom_bar(show.legend = T,stat = "summary",position = "dodge") + #stack
  #geom_point(size=1) + 
  scale_x_discrete(breaks=regions,labels=regions) +
  ylab("DMR number") +  # dmr edgeR p < 0.01
  theme_bw() + 
  theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) +
  theme(axis.text.x = element_text(angle = 25,vjust=0.2, hjust = 1)) 

# new
regions <- c("gene","promoter","promoter5k","promoter150TSS50","promoter300exon1end100","gene_exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")
dmr <-  read.table("/BioII/lulab_b/baopengfei/projects/methylation/data/lulab/diff/summary-dmr.txt",header = T,stringsAsFactors = F)
str(dmr)
dmr$region <- factor(dmr$region,levels=regions)
#dmr.ratio <- dmr
#dmr.ratio$hyper.ratio <- dmr.ratio$hyper/region.meta$number
#dmr.ratio$hypo.ratio <- dmr.ratio$hypo/region.meta$number
#dmr.ratio <- dmr.ratio[,c(1,4,5)]
dmr <- reshape2::melt(dmr)
colnames(dmr) <- c("region","disease","method","trend","number")
#dmr.ratio <- melt(dmr.ratio)
#colnames(dmr.ratio) <- c("region","trend","number")
dmr.crc <- dmr[dmr$disease=="CRC",]
dmr.stad <- dmr[dmr$disease=="STAD",]


library(ggpubr)
library(rstatix)
bp <- ggbarplot(dmr.crc, x = "region", y = "number", fill = "trend",
                palette = c("steelblue","chocolate","firebrick"), #add = "mean_sd",   # dmr只能barplot无法检验分组显著性
                position = position_dodge(0.8))

regions.simplify <- c("gene","promoter","promoter5k","-150TSS+50","-300exon1end+100","exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")

bp +  #stat_pvalue_manual(stat.test, label = "p.adj.signif",hide.ns=T, size=8, tip.length = 0.03,bracket.nudge.y = 0.4) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  ylab("DMR number") + 
  scale_x_discrete(breaks=regions,labels=regions.simplify) +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold")) 


#### 6.enrich boxplot ####
#regions <- c("gene","promoter","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")
#regions <- c("gene","promoter","promoter5k","promoter150TSS50","promoter300exon1end100","gene_exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")
peak <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/macs2/summary-peak.txt")
colnames(peak) <- c("region","sample","a","b","c","d")
table(peak$region)
peak$region <- factor(peak$region,levels=regions)

qc <- read.table("meta/lulab/sample_table.txt",sep = "\t",stringsAsFactors = F,header = T)
qc <- qc$data_id[qc$medip_passQCstringent == "Y"]
peak <- peak[peak$sample %in% qc,]
peak$enrich <- peak$a*peak$d/((peak$c-peak$a)*peak$b)
peak$log2enrich <- log2(peak$enrich)

peak$group <- "CRC"
peak$group[grep(pattern = "NC",x = peak$sample)] <- "NC"
peak$group[grep(pattern = "STAD",x = peak$sample)] <- "STAD"
#peak <- peak[peak$group != "STAD",]
peak$group <- factor(peak$group,levels=c("NC","CRC","STAD"))

peak <- peak[peak$region %in% regions,]
peak$group[grep(pattern = "STAD",x = peak$sample)] <- "STAD"

library(ggplot2)
ggplot(data=peak, aes(x=region,y=log2enrich,fill=group)) +
  geom_boxplot(show.legend = T) +
  scale_fill_manual(values=alpha(c("orange2","steelblue2","firebrick2"), 1))  + 
  scale_x_discrete(breaks=regions,labels=regions) +
  ylab("peak log2(odds ratio)") +  # reads assigned ratio
  theme_classic() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold")) + 
  stat_compare_means(aes(group = group,label = ..p.signif..) ,bracket.size = 0.3, vjust = 2,size = 5,color="red", method = "wilcox.test",hide.ns=F)


#### 7.enrich boxplot (CRC + STAD + NC) #### 
regions <- c("gene","promoter","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")

peak <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/macs2/summary-peak.txt")
colnames(peak) <- c("region","sample","a","b","c","d")
qc <- read.table("meta/lulab/archive/sample_QC.txt",sep = "\t",stringsAsFactors = F,header = T)
qc <- qc$medip_data_id[qc$medip_overall_QC == "Y"]
peak <- peak[peak$sample %in% qc,]
peak$enrich <- peak$a*peak$d/((peak$c-peak$a)*peak$b)
peak$log2enrich <- log2(peak$enrich)

peak$group <- "CRC"
peak$group[grep(pattern = "NC",x = peak$sample)] <- "NC"
peak$group[grep(pattern = "STAD",x = peak$sample)] <- "STAD"
#peak <- peak[peak$group != "STAD",]
peak$group <- factor(peak$group,levels = c("NC","CRC","STAD") )
peak <- peak[peak$region %in% regions,]
peak$region <- factor(peak$region,levels=regions)

mytheme <-  theme(plot.title = element_text(size = 16,color="black",hjust = 0.5,face="bold"),
                  axis.title = element_text(size = 16,color ="black",face="bold"), 
                  axis.text = element_text(size= 16,color = "black",face="bold"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1 ),
                  panel.grid=element_blank(),
                  legend.position = "top",
                  legend.text = element_text(size= 16),
                  legend.title= element_text(size= 16))

box <- ggplot(data=peak, aes(x=group,y=log2enrich,fill=group)) +
  facet_wrap(~region)+
  geom_boxplot(show.legend = T) +
  scale_fill_manual(values=alpha(c("steelblue2","orange2","firebrick1"), 1))  + 
  #geom_point(size=1) + 
  #scale_x_discrete(breaks=group,labels=group) +
  #scale_y_continuous(trans='log10') + 
  #xlim(c(0,5000000)) +  
  #ylim(c(0,1000000))
  #xlab( ) + 
  ylab("peak log2(odds ratio)") +  # reads assigned ratio
  #scale_y_continuous(oob = scales::squish) +
  #stat_compare_means(comparisons = list(c("CRC","NC"))) +
  #labs(title="in-house cfMeDIP",x="CRC FDA Biomarker",y="transformed counts") +
  #geom_jitter(width = .2) +
  #scale_x_discrete(labels= function(x) highlight(x, c("promoter","CpG_island"), "red")) + 
  #theme(axis.text.x=element_markdown()) 
  theme_classic() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold")) + 
  #theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) +
  #theme(axis.text.x = element_text(angle = 90,vjust=0.2, hjust = 1)) + 
  stat_compare_means(aes(label = ..p.signif..),size = 5,color="red", ref.group = "NC", method = "wilcox.test",hide.ns=T)  # group = group, ref.group = "NC",

ggsave(filename = "./enrich.pdf",width = 10, height = 12, plot = box,path ="./",limitsize = F)


#### 8. final in-house data: enrich barplot: ggpubr ####
# in-house data
#https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
#regions <- c("gene","promoter","promoter5k","promoter150TSS50","promoter300exon1end100","gene_exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")
regions <- c("gene","promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE",'retroposon')
peak <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/lulab/macs2/summary-peak.txt")
colnames(peak) <- c("region","sample","a","b","c","d")
table(peak$region)
peak$region <- factor(peak$region,levels=regions)

qc <- read.table("meta/lulab/sample_table.txt",sep = "\t",stringsAsFactors = F,header = T)
qc <- qc$data_id[qc$medip_passQCstringent == "Y"]
peak <- peak[peak$sample %in% qc,]
peak$enrich <- peak$a*peak$d/((peak$c-peak$a)*peak$b)
peak$log2enrich <- log2(peak$enrich)

peak$group <- "CRC"
peak$group[grep(pattern = "NC",x = peak$sample)] <- "NC"
peak$group[grep(pattern = "STAD",x = peak$sample)] <- "STAD"
#peak <- peak[peak$group != "STAD",]
peak$group <- factor(peak$group,levels=c("NC","CRC","STAD"))

peak <- peak[peak$region %in% regions,]
#peak <- peak[peak$group!="CRC",]

library(ggpubr)
library(rstatix)
bp <- ggbarplot(peak, x = "region", y = "log2enrich", fill = "group",
                palette = c("steelblue","chocolate","firebrick"), add = "mean_sd",
                position = position_dodge(0.7))

stat.test <- peak %>%
  group_by(region) %>%
  rstatix::wilcox_test(log2enrich ~ group)

stat.test <- stat.test %>%
  add_xy_position(x = "region", fun = "mean_sd", dodge = 1)

regions.simplify <- c("gene","promoter","CDS","gene_exon1","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")

bp +  stat_pvalue_manual(stat.test, label = "p",hide.ns=T, size=8, tip.length = 0.03,bracket.nudge.y = 0.4) +   # label = "p.adj.signif"
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  ylab("peak log2(odds ratio)") + 
  scale_x_discrete(breaks=regions,labels=regions.simplify) +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold")) 


#### 8. final public data: enrich barplot: ggpubr ####
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/")
regions <- c("gene","promoter","promoter5k","promoter150TSS50","promoter300exon1end100","gene_exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")
#GSE87096,GSE109202
peak <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DIP-seq/output/GSE109202/macs2/summary-peak.txt")
colnames(peak) <- c("region","sample","a","b","c","d")
table(peak$region)
peak$region <- factor(peak$region,levels=regions)

crc.qc <- read.table("meta/GSE109202/sample_ids_medip_CRC.txt",sep = "\t",stringsAsFactors = F,header = F)$V1
nc.qc <- read.table("meta/GSE109202/sample_ids_medip_NC.txt",sep = "\t",stringsAsFactors = F,header = F)$V1
peak <- peak[peak$sample %in% c(crc.qc,nc.qc),]
peak$enrich <- peak$a*peak$d/((peak$c-peak$a)*peak$b)
peak$log2enrich <- log2(peak$enrich)

peak$group <- "CRC"
peak$group[peak$sample %in% nc.qc] <- "NC"

peak$group <- factor(peak$group,levels=c("NC","CRC"))

peak <- peak[peak$region %in% regions,]

library(ggpubr)
library(rstatix)
bp <- ggbarplot(peak, x = "region", y = "log2enrich", fill = "group",
                palette = c("steelblue","chocolate"), add = "mean_sd",
                position = position_dodge(0.8))

stat.test <- peak %>%
  group_by(region) %>%
  rstatix::wilcox_test(log2enrich ~ group)

stat.test <- stat.test %>%
  add_xy_position(x = "region", fun = "mean_sd", dodge = 1)

regions.simplify <- c("gene","promoter","promoter5k","-150TSS+50","-300exon1end+100","exon1","CDS","UTR3","UTR5","CpG_island","LINE","SINE","retroposon")

bp +  stat_pvalue_manual(stat.test, label = "p",hide.ns=T, size=8, tip.length = 0.03,bracket.nudge.y = 0.4) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  ylab("peak log2(odds ratio)") + 
  scale_x_discrete(breaks=regions,labels=regions.simplify) +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black",face="bold"), 
        axis.text = element_text(size= 16,color = "black",face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "right",
        legend.text = element_text(size= 24,face="bold"),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold")) 


