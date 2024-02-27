# correct WGS GC (test)
# last 220111 by bpf
#correct from count mat not CPM !!! 


 
#library(PopSV)


# get GC corrected for wgs bins short&long-----------------------------------------------
ref <- read.table("/BioII/lulab_b/baopengfei/gitsoft/delfi_scripts-master/hg38.masked.bed",header = F)
colnames(ref) <- c("chr","start","end")
str(ref)
hg38 <- PopSV::getGC(bins.df = ref, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

## get short&long count GC adjusted
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
short <- read.table("./output/lulab/matrix/count_matrix_bin100kb_short.txt",header = T, row.names = 1, check.names = F)
long <- read.table("./output/lulab/matrix/count_matrix_bin100kb_long.txt",header = T, row.names = 1, check.names = F)

long.short <- cbind(long,short)

correct <- function(sample.id){
# sample.id <- "NC-PKU-10-wgs-long"
# "NC-PKU-10-wgs-short"
print(sample.id)
test <- cbind(rownames(long.short),long.short[[sample.id]])
test <- as.data.frame(test)
colnames(test) <- c("bin","bc")  # bin100kb_chr1:29555235-34534399|4979164
test$bin <- gsub("bin100kb_","",test$bin)
test$chr <- unlist(lapply(strsplit(test$bin,":",fixed = T), function(x) x[1]))
test$start <- as.numeric(unlist(lapply(strsplit(test$bin,":|-",perl = T), function(x) x[2])))
test$end <- as.numeric(unlist(lapply(strsplit(test$bin,":|-|\\|",perl = T), function(x) x[3])))
test <- test[,c("chr","start","end","bc")]
test$bc <- as.numeric(test$bc)
#str(test)

#correct.gc <- PopSV::correct.GC(bc.f = test, gc.df = hg38)  # , outfile.prefix = "./pop"
gc.df <- hg38
bc.f <- test

gc.class = cut(gc.df$GCcontent, breaks = seq(0, 1, 0.02), include.lowest = TRUE)
samp.ii = unlist(tapply(1:length(gc.class), gc.class, function(e) e[sample(1:length(e),
                                                                           min(c(length(e), 500)))]))
if(!is.data.frame(bc.f) & length(bc.f)==1 & is.character(bc.f)){
  bc.df = utils::read.table(bc.f, as.is = TRUE, header = TRUE)
} else {
  bc.df = bc.f
}

if(any(bc.df$chr != gc.df$chr) | any(bc.df$start != gc.df$start)){
  bc.df = merge(bc.df, gc.df[,c('chr','start','end','GCcontent')])
  bc.df = bc.df[order(as.character(bc.df$chr), bc.df$start),]
} else {
  bc.df$GCcontent = gc.df$GCcontent
}

lo = stats::loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/stats::predict(lo, newdata = bc.df)
bc.df$bc = round(bc.df$bc, digits = 2)
if (any(bc.df$bc < 0, na.rm = TRUE))
  bc.df$bc[bc.df$bc < 0] = 0
bc.df$GCcontent = NULL
colnames(bc.df)[4] <- sample.id
return(bc.df)
}

dat <- lapply(colnames(long.short),correct)
dat <- do.call("cbind",dat)
dat <- dat[,!duplicated(colnames(dat))]
dat <- dat[,!(colnames(dat) %in% c("chr","start","end"))]
dat <- round(dat)   # keep interger

dat.l <- cbind(rownames(long.short),dat[,grep("long",colnames(dat))])
colnames(dat.l)[1] <- "gene_id"
dat.s <- cbind(rownames(long.short),dat[,grep("short",colnames(dat))])
colnames(dat.s)[1] <- "gene_id"

#options ---
write.table(dat.l,"./output/lulab/matrix/count_matrix_bin100kb_long.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(dat.s,"./output/lulab/matrix/count_matrix_bin100kb_short.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)

# plot(x=bc.df$bc,y=test$bc)
# cor.test(bc.df$bc,test$bc)
# 
# plot(x=gc.df$GCcontent, y=test$bc)
# plot(x=gc.df$GCcontent, y=bc.df$bc)







# get GC corrected for wgs genes&promoter ----------------------------------------------

#options ---
#ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/gene.bed",header = F)
ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/promoter.bed",header = F)

ref <- ref[,1:4]
colnames(ref) <- c("chr","start","end","ENSG")
str(ref)
hg38 <- PopSV::getGC(bins.df = ref, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
hg38$pos <- paste0(hg38$chr,":",hg38$start,"-",hg38$end)
table(duplicated(hg38$pos))  # 13 TURE
hg38$ID <- paste0(hg38$pos,":",hg38$ENSG)
table(duplicated(hg38$ID))  # 0 TURE


## get short&long count GC adjusted
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

#options ---
#gene <- read.table("./output/lulab/matrix/count_matrix_gene.txt",header = T, row.names = 1, check.names = F)
gene <- read.table("./output/lulab/matrix/count_matrix_promoter.txt",header = T, row.names = 1, check.names = F)

g <- unlist(lapply(strsplit(rownames(gene),"|",fixed = T), function(x) x[1]))

correct <- function(sample.id){
  # sample.id <- "CRC-PKU-10-wgs"
  # "NC-PKU-10-wgs-short"
  print(sample.id)
  pos <- hg38$pos[match(g,hg38$ENSG)]
  # chr <- hg38$chr[match(g,hg38$ENSG)]
  # start <- hg38$start[match(g,hg38$ENSG)]
  # end <- hg38$end[match(g,hg38$ENSG)]
  
  
  test <- cbind(pos,gene[[sample.id]])
  #table(pos==rownames(g))
  test <- as.data.frame(test)
  colnames(test) <- c("gene","bc")  # bin100kb_chr1:29555235-34534399|4979164
  #test$gene <- gsub("bin100kb_","",test$gene)
  test$chr <- unlist(lapply(strsplit(test$gene,":",fixed = T), function(x) x[1]))
  test$start <- as.numeric(unlist(lapply(strsplit(test$gene,":|-",perl = T), function(x) x[2])))
  test$end <- as.numeric(unlist(lapply(strsplit(test$gene,":|-|\\|",perl = T), function(x) x[3])))
  test <- test[,c("chr","start","end","bc")]
  test$bc <- as.numeric(test$bc)
  #str(test)
  
  #correct.gc <- PopSV::correct.GC(bc.f = test, gc.df = hg38)  # , outfile.prefix = "./pop"
  gc.df <- hg38
  bc.f <- test
  
  gc.class = cut(gc.df$GCcontent, breaks = seq(0, 1, 0.02), include.lowest = TRUE)
  samp.ii = unlist(tapply(1:length(gc.class), gc.class, function(e) e[sample(1:length(e),
                                                                             min(c(length(e), 500)))]))
  if(!is.data.frame(bc.f) & length(bc.f)==1 & is.character(bc.f)){
    bc.df = utils::read.table(bc.f, as.is = TRUE, header = TRUE)
  } else {
    bc.df = bc.f
  }
  
  if(any(bc.df$chr != gc.df$chr) | any(bc.df$start != gc.df$start)){
    bc.df = merge(bc.df, gc.df[,c('chr','start','end','GCcontent')])
    bc.df = bc.df[order(as.character(bc.df$chr), bc.df$start),]
  } else {
    bc.df$GCcontent = gc.df$GCcontent
  }
  
  lo = stats::loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
  bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/stats::predict(lo, newdata = bc.df)
  bc.df$bc = round(bc.df$bc, digits = 2)
  if (any(bc.df$bc < 0, na.rm = TRUE))
    bc.df$bc[bc.df$bc < 0] = 0
  bc.df$GCcontent = NULL
  colnames(bc.df)[4] <- sample.id
  return(bc.df)
}

dat <- lapply(colnames(gene),correct)
dat1 <- do.call("cbind",dat)
dat1 <- dat1[,!duplicated(colnames(dat1))]
# head(ref)
# head(gene[,1:4])
# head(dat1[,1:4])
# tail(ref)
# tail(gene[,1:4])
# tail(dat1[,1:4])  # seems order not changed
dat1 <- dat1[,!(colnames(dat1) %in% c("chr","start","end"))]
dat1 <- round(dat1)   # keep interger
row.names(dat1) <- g
dat.l <- cbind(rownames(gene),dat1)
colnames(dat.l)[1] <- "gene_id"

#options ---
#write.table(dat.l,"./output/lulab/matrix/count_matrix_gene.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(dat.l,"./output/lulab/matrix/count_matrix_promoter.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)

##check cor and scatter plot
cor.test(dat.l$`CRC-PKU-12-wgs`,gene$`CRC-PKU-12-wgs`)
plot(x=hg38$GCcontent,y=gene$`CRC-PKU-10-wgs`)
plot(x=hg38$GCcontent,y=dat.l$`CRC-PKU-10-wgs`)



# get GC corrected for wgs genes short&long----------------------------------------------
#ref <- read.table("/BioII/lulab_b/baopengfei/gitsoft/delfi_scripts-master/hg38.masked.bed",header = F)
ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/gene.bed",header = F)
ref <- ref[,1:4]
colnames(ref) <- c("chr","start","end","ENSG")
str(ref)
hg38 <- PopSV::getGC(bins.df = ref, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
hg38$pos <- paste0(hg38$chr,":",hg38$start,"-",hg38$end)
table(duplicated(hg38$pos))  # 13 TURE
hg38$ID <- paste0(hg38$pos,":",hg38$ENSG)
table(duplicated(hg38$ID))  # 0 TURE


## get short&long count GC adjusted
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")

#options ---
gene <- read.table("./output/lulab/matrix/count_matrix_gene_long.txt",header = T, row.names = 1, check.names = F)
#gene <- read.table("./output/lulab/matrix/count_matrix_gene_short.txt",header = T, row.names = 1, check.names = F)

g <- unlist(lapply(strsplit(rownames(gene),"|",fixed = T), function(x) x[1]))

correct <- function(sample.id){
  # sample.id <- "CRC-PKU-10-wgs"
  # "NC-PKU-10-wgs-short"
  print(sample.id)
  pos <- hg38$pos[match(g,hg38$ENSG)]
  # chr <- hg38$chr[match(g,hg38$ENSG)]
  # start <- hg38$start[match(g,hg38$ENSG)]
  # end <- hg38$end[match(g,hg38$ENSG)]
  
  
  test <- cbind(pos,gene[[sample.id]])
  #table(pos==rownames(g))
  test <- as.data.frame(test)
  colnames(test) <- c("gene","bc")  # bin100kb_chr1:29555235-34534399|4979164
  #test$gene <- gsub("bin100kb_","",test$gene)
  test$chr <- unlist(lapply(strsplit(test$gene,":",fixed = T), function(x) x[1]))
  test$start <- as.numeric(unlist(lapply(strsplit(test$gene,":|-",perl = T), function(x) x[2])))
  test$end <- as.numeric(unlist(lapply(strsplit(test$gene,":|-|\\|",perl = T), function(x) x[3])))
  test <- test[,c("chr","start","end","bc")]
  test$bc <- as.numeric(test$bc)
  #str(test)
  
  #correct.gc <- PopSV::correct.GC(bc.f = test, gc.df = hg38)  # , outfile.prefix = "./pop"
  gc.df <- hg38
  bc.f <- test
  
  gc.class = cut(gc.df$GCcontent, breaks = seq(0, 1, 0.02), include.lowest = TRUE)
  samp.ii = unlist(tapply(1:length(gc.class), gc.class, function(e) e[sample(1:length(e),
                                                                             min(c(length(e), 500)))]))
  if(!is.data.frame(bc.f) & length(bc.f)==1 & is.character(bc.f)){
    bc.df = utils::read.table(bc.f, as.is = TRUE, header = TRUE)
  } else {
    bc.df = bc.f
  }
  
  if(any(bc.df$chr != gc.df$chr) | any(bc.df$start != gc.df$start)){
    bc.df = merge(bc.df, gc.df[,c('chr','start','end','GCcontent')])
    bc.df = bc.df[order(as.character(bc.df$chr), bc.df$start),]
  } else {
    bc.df$GCcontent = gc.df$GCcontent
  }
  
  lo = stats::loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
  bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/stats::predict(lo, newdata = bc.df)
  bc.df$bc = round(bc.df$bc, digits = 2)
  if (any(bc.df$bc < 0, na.rm = TRUE))
    bc.df$bc[bc.df$bc < 0] = 0
  bc.df$GCcontent = NULL
  colnames(bc.df)[4] <- sample.id
  return(bc.df)
}

dat <- lapply(colnames(gene),correct)
dat1 <- do.call("cbind",dat)
dat1 <- dat1[,!duplicated(colnames(dat1))]
# head(ref)
# head(gene[,1:4])
# head(dat1[,1:4])
# tail(ref)
# tail(gene[,1:4])
# tail(dat1[,1:4])  # seems order not changed
dat1 <- dat1[,!(colnames(dat1) %in% c("chr","start","end"))]
dat1 <- round(dat1)   # keep interger
row.names(dat1) <- g
dat.l <- cbind(rownames(gene),dat1)
colnames(dat.l)[1] <- "gene_id"

#options ---
write.table(dat.l,"./output/lulab/matrix/count_matrix_gene_long.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)
#write.table(dat.l,"./output/lulab/matrix/count_matrix_gene_short.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)

##check cor and scatter plot
cor.test(dat.l$`CRC-PKU-12-wgs`,gene$`CRC-PKU-12-wgs`)
plot(x=hg38$GCcontent,y=gene$`CRC-PKU-10-wgs`)
plot(x=hg38$GCcontent,y=dat.l$`CRC-PKU-10-wgs`)


# get GC corrected for wgs&medip promoter&gene (newer) ----------------------------------------------
#get gene&promoter count GC adjusted

for(d in c("DIP","DNA")){
print(d)
#d <- "DNA"
setwd(paste0("/BioII/lulab_b/baopengfei/projects/exOmics/",d,"-seq/"))

for (region in c("gene","promoter")){
#region <- "promoter"
print(region)
  
#ref <- read.table("/BioII/lulab_b/baopengfei/shared_reference/hg38/promoter.bed",header = F)
ref <- read.table(paste0("/BioII/lulab_b/baopengfei/shared_reference/hg38/",region,".bed"),header = F)
ref <- ref[,1:4]
colnames(ref) <- c("chr","start","end","ENSG")
str(ref)

hg38 <- PopSV::getGC(bins.df = ref, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
hg38$pos <- paste0(hg38$chr,":",hg38$start,"-",hg38$end)
table(duplicated(hg38$pos))  # 13 TURE
hg38$ID <- paste0(hg38$pos,":",hg38$ENSG)
table(duplicated(hg38$ID))  # 0 TURE


#gene <- read.table("./output/lulab/matrix/count_matrix_promoter.txt",header = T, row.names = 1, check.names = F)
gene <- read.table(paste0("./output/lulab/matrix/count_matrix_",region,".txt"),header = T, row.names = 1, check.names = F)
g <- unlist(lapply(strsplit(rownames(gene),"|",fixed = T), function(x) x[1]))
#g <- gsub("promoter_","",g)

# if(region=="promoter"){
#   hg38$ENSG <- paste0("promoter_",hg38$ENSG)
# }else{
#   hg38$ENSG <- gsub("promoter_","",hg38$ENSG)
# }
correct <- function(sample.id){
  # sample.id <- "CRC-PKU-10-me"
  # "NC-PKU-10-wgs-short"
  print(sample.id)
  
  pos <- hg38$pos[match(g,hg38$ENSG)]
  # chr <- hg38$chr[match(g,hg38$ENSG)]
  # start <- hg38$start[match(g,hg38$ENSG)]
  # end <- hg38$end[match(g,hg38$ENSG)]
  
  
  test <- cbind(pos,gene[[sample.id]])
  #table(pos==rownames(g))
  test <- as.data.frame(test)
  colnames(test) <- c("gene","bc")  # bin100kb_chr1:29555235-34534399|4979164
  #test$gene <- gsub("bin100kb_","",test$gene)
  test$chr <- unlist(lapply(strsplit(test$gene,":",fixed = T), function(x) x[1]))
  test$start <- as.numeric(unlist(lapply(strsplit(test$gene,":|-",perl = T), function(x) x[2])))
  test$end <- as.numeric(unlist(lapply(strsplit(test$gene,":|-|\\|",perl = T), function(x) x[3])))
  test <- test[,c("chr","start","end","bc")]
  test$bc <- as.numeric(test$bc)
  #str(test)
  
  #correct.gc <- PopSV::correct.GC(bc.f = test, gc.df = hg38)  # , outfile.prefix = "./pop"
  gc.df <- hg38
  bc.f <- test
  
  gc.class = cut(gc.df$GCcontent, breaks = seq(0, 1, 0.02), include.lowest = TRUE)
  samp.ii = unlist(tapply(1:length(gc.class), gc.class, function(e) e[sample(1:length(e),
                                                                             min(c(length(e), 500)))]))
  if(!is.data.frame(bc.f) & length(bc.f)==1 & is.character(bc.f)){
    bc.df = utils::read.table(bc.f, as.is = TRUE, header = TRUE)
  } else {
    bc.df = bc.f
  }
  
  if(any(bc.df$chr != gc.df$chr) | any(bc.df$start != gc.df$start)){
    bc.df = merge(bc.df, gc.df[,c('chr','start','end','GCcontent')])
    bc.df = bc.df[order(as.character(bc.df$chr), bc.df$start),]
  } else {
    bc.df$GCcontent = gc.df$GCcontent
  }
  
  lo = stats::loess(bc ~ GCcontent, data = bc.df[samp.ii, ])
  bc.df$bc = mean(bc.df$bc, na.rm = TRUE) * bc.df$bc/stats::predict(lo, newdata = bc.df)
  bc.df$bc = round(bc.df$bc, digits = 2)
  if (any(bc.df$bc < 0, na.rm = TRUE))
    bc.df$bc[bc.df$bc < 0] = 0
  bc.df$GCcontent = NULL
  colnames(bc.df)[4] <- sample.id
  return(bc.df)
}

dat <- lapply(colnames(gene),correct)
dat1 <- do.call("cbind",dat)
dat1 <- dat1[,!duplicated(colnames(dat1))]
# head(ref)
# head(gene[,1:4])
# head(dat1[,1:4])
# tail(ref)
# tail(gene[,1:4])
# tail(dat1[,1:4])  # seems order not changed
dat1 <- dat1[,!(colnames(dat1) %in% c("chr","start","end"))]
dat1 <- round(dat1)   # keep interger
row.names(dat1) <- g
dat.l <- cbind(rownames(gene),dat1)
colnames(dat.l)[1] <- "gene_id"

#options ---
#write.table(dat.l,"./output/lulab/matrix/count_matrix_promoter.correctGC.txt",sep = "\t",quote = F,row.names = F,col.names = T)
write.table(dat.l,paste0("./output/lulab/matrix/count_matrix_",region,".correctGC.txt"),sep = "\t",quote = F,row.names = F,col.names = T)

}
}


##check cor and scatter plot
print(cor.test(dat.l$`CRC-PKU-12-me`,gene$`CRC-PKU-12-me`))
print(cor.test(dat.l$`CRC-PKU-12-wgs`,gene$`CRC-PKU-12-wgs`))
plot(x=hg38$GCcontent,y=gene$`CRC-PKU-10-me`)
plot(x=hg38$GCcontent,y=dat.l$`CRC-PKU-10-me`)
