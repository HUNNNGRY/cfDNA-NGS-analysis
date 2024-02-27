# microbe
# last 220114 by bpf
#NA should be replaced with 0

# sum mat -----------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/")
options(stringsAsFactors = F)

files <- Sys.glob("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/microbe/report/*.txt")
read.mb <- function(x){
  #x <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/microbe/report/STAD-PKU-10-wgs.txt"
  
  tmp <- read.table(x,sep = "\t")
  
  ## only keep genus
  #table(tmp$V4)
  tmp <- tmp[tmp$V4=="G",]
  
  tmp <- tmp[,c(1,6)]
  sample <- gsub(".txt","",basename(x))
  print(sample)
  colnames(tmp) <- c("percent","feature")
  tmp[["sample"]] <- sample
  #tmp$V6
  return(tmp)
}
dat <- lapply(files,read.mb)
dat <- do.call("rbind",dat)

f <- unique(dat$feature)
read.mb2 <- function(x){
  #x <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/microbe/report/STAD-PKU-10-wgs.txt"
  
  tmp <- read.table(x,sep = "\t")
  
  ## only keep genus
  #table(tmp$V4)
  tmp <- tmp[tmp$V4=="G",]
  tmp <- tmp[,c(6,1)]
  sample <- gsub(".txt","",basename(x))
  print(sample)
  colnames(tmp) <- c("feature",sample)
  #tmp[["sample"]] <- sample
  rownames(tmp) <- tmp$feature
  tmp <- tmp[f,]
  rownames(tmp) <- f
  tmp$feature <- f
  
  #table(is.na(tmp$`STAD-PKU-10-wgs`))
  #tmp$V6
  return(tmp)
}
dat <- lapply(files,read.mb2)
dat <- do.call("cbind",dat)
dat <- dat[,!duplicated(colnames(dat))] 
dat[is.na(dat)] <- 0  ## for Mb Mat, NA should be replaced with 0

dat$feature <- gsub(" ","",dat$feature)
#dat.a <- reshape2::acast(data = dat,formula =  percent ~ sample, value.var="feature")


write.table(dat,"./output/lulab/matrix/DNA-Mb_matrix_genus.txt",sep = "\t",quote = F,row.names = F,col.names = T)

