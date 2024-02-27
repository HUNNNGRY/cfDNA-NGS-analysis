# Mt 



# wgs Mt.assign.ratio ------------------------------------------------------------
setwd("/BioII/lulab_b/baopengfei/projects/exOmics/")

crc <- read.delim("../multi-omics-explore/meta/lulab/paired-wgs-CRC-passQCloose.txt",header = F,stringsAsFactors = F)$V1
nc <- read.delim("../multi-omics-explore/meta/lulab/paired-wgs-NC-passQCloose.txt",header = F,stringsAsFactors = F)$V1
stad <- read.delim("../multi-omics-explore/meta/lulab/paired-wgs-STAD-passQCloose.txt",header = F,stringsAsFactors = F)$V1
wgs <- c(crc,nc,stad)


readin <- function(x){
print(x)
#x <- "/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/100_120_count_matrix_Mt.txt.summary"
#mt.wgs <- read.table("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/count_matrix_Mt.txt.summary",sep = "\t",row.names = 1,header = T)
mt.wgs <- read.table(x,sep = "\t",row.names = 1,header = T)

mt.wgs <- as.data.frame(t(mt.wgs))
#head(mt.wgs)
mt.wgs$ratio <- mt.wgs$Assigned/(mt.wgs$Assigned+mt.wgs$Unassigned_NoFeatures)*100
#summary(mt.wgs$ratio)
mt.wgs$group <- unlist(lapply(strsplit(rownames(mt.wgs),"\\."), function(x) x[6]))  # 6 For WGS
mt.wgs$group <- factor(mt.wgs$group,levels = c("NC","CRC","STAD"))
mt.wgs[1:3,1:3]
rownames(mt.wgs) <- gsub("output.lulab.bam.sorted.deduped.","",rownames(mt.wgs))
rownames(mt.wgs) <- gsub("output.lulab.test.MT.bam.","",rownames(mt.wgs))
rownames(mt.wgs) <- gsub(".bam","",rownames(mt.wgs))
rownames(mt.wgs) <- gsub(".","-",rownames(mt.wgs),fixed = T)
rownames(mt.wgs) <- unlist(lapply(strsplit(rownames(mt.wgs),"_",fixed=T),function(x) x[1]))

table(wgs %in% rownames(mt.wgs))

mt.wgs <- as.data.frame(t(mt.wgs))
rownames(mt.wgs)

mt.wgs <- mt.wgs[rownames(mt.wgs)=="ratio",]

f <- unlist(lapply(strsplit(x,"/",fixed=T),function(x) x[11]))
f <- gsub("_count_matrix_Mt.txt.summary","",f)
f <- gsub("count_matrix_Mt.txt.summary","",f)

rownames(mt.wgs) <- paste0("Mt.assign.ratio","_",f)
mt.wgs <- mt.wgs[,wgs]
#dir.create("./data/DNA-Mt")
#write.table(mt.wgs,"./data/DNA-Mt/ML.txt",sep = "\t",quote = F,row.names = T,col.names = T)
mt.wgs <- cbind(rownames(mt.wgs),mt.wgs)
colnames(mt.wgs)[1] <- "feature"
return(mt.wgs)
}

files <- Sys.glob("/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/*count_matrix_Mt.txt.summary")
dat <- lapply(files,readin)
dat <- do.call("rbind",dat)


write.table(dat,"/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/output/lulab/matrix/MtAssignRatio_matrix_Mt.txt",sep = "\t",quote = F,row.names = F,col.names = T)



# Mt.Size.Ratio -----------------------------------------------------------
#already in gene-centric fragsize matrix (chrM part)



# Mt.Size.Assign.Ratio ---------------------------------------------------------------------


